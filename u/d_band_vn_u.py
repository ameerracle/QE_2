from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import numpy as np


DEFAULT_BASE_DIR = Path("/scratch/anizami/QE_2/u/PDOS/VN_slab")
DEFAULT_OUTPUT_NAME = "d_band_VN_plus_U.csv"


def get_fermi_energy(nscf_file: Path) -> float | None:
    if not nscf_file.exists():
        return None

    float_pattern = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][-+]?\d+)?")

    with nscf_file.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            lowered = line.lower()
            if "fermi" not in lowered and "highest occupied" not in lowered:
                continue

            numbers = float_pattern.findall(line)
            if not numbers:
                continue

            if "spin" in lowered and len(numbers) >= 2:
                try:
                    return float(np.mean([float(numbers[0]), float(numbers[1])]))
                except ValueError:
                    continue

            try:
                return float(numbers[-1])
            except ValueError:
                continue
    return None


def calculate_dband_center(pdos_file: Path, e_fermi: float) -> float | None:
    try:
        data = np.loadtxt(pdos_file)
    except Exception as exc:
        print(f"Error reading {pdos_file}: {exc}")
        return None

    if data.ndim != 2 or data.shape[1] < 2:
        print(f"Skipping {pdos_file}: unexpected PDOS shape {getattr(data, 'shape', None)}")
        return None

    energy = data[:, 0] - e_fermi

    ncols = data.shape[1]
    if ncols == 7:
        # Non-spin: E, ldos, pdos(m=-2..2)
        dos = data[:, 2:].sum(axis=1)
    elif ncols == 13:
        # Spin-polarized: E, ldosup, ldosdw, pdosup/dw for the five d components.
        dos = data[:, 3::2].sum(axis=1) + data[:, 4::2].sum(axis=1)
    elif ncols > 3:
        # Conservative fallback for other projwfc layouts: keep only projected columns.
        dos = data[:, 2:].sum(axis=1)
    else:
        dos = data[:, 1]

    numerator = np.trapz(dos * energy, energy)
    denominator = np.trapz(dos, energy)

    if denominator == 0:
        return 0.0

    return float(numerator / denominator)


def find_nscf_file(base_dir: Path) -> Path | None:
    candidates = sorted(base_dir.glob("*nscf*.out"))
    if candidates:
        return candidates[0]

    fallback = base_dir / "VN_slab_nscf.out"
    if fallback.exists():
        return fallback

    return None


def find_pdos_files(base_dir: Path, element: str) -> list[Path]:
    pattern = f"slab_pdos.pdos_atm#*({element})_wfc#*(d)"
    return sorted(base_dir.glob(pattern), key=lambda path: int(path.name.split("#")[1].split("(")[0]))


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compute d-band center values for the VN +U slab PDOS."
    )
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=DEFAULT_BASE_DIR,
        help=f"Directory containing the VN slab NSCF output and PDOS files (default: {DEFAULT_BASE_DIR})",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Optional output CSV path. Defaults to d_band_VN_plus_U.csv inside --base-dir.",
    )
    args = parser.parse_args()

    base_dir = args.base_dir.expanduser().resolve()
    if not base_dir.exists():
        print(f"Base directory not found: {base_dir}")
        return 1

    nscf_path = find_nscf_file(base_dir)
    if nscf_path is None:
        print(f"Could not find an NSCF output in {base_dir}")
        return 1

    e_fermi = get_fermi_energy(nscf_path)
    if e_fermi is None:
        print(f"Could not parse Fermi energy from {nscf_path}")
        return 1

    pdos_files = find_pdos_files(base_dir, element="V")
    if not pdos_files:
        print(f"No V d-PDOS files found in {base_dir}")
        return 1

    print("--- VN (+U) d-Band Center Analysis ---")
    print(f"Source: {nscf_path}")
    print(f"Fermi Energy (Ef): {e_fermi:.4f} eV")
    print("-" * 40)

    rows: list[dict[str, float | int]] = []
    for v_index, pdos_file in enumerate(pdos_files, start=1):
        try:
            atom_id = int(pdos_file.name.split("#")[1].split("(")[0])
        except (IndexError, ValueError):
            print(f"Skipping unrecognized file name: {pdos_file.name}")
            continue

        db_center = calculate_dband_center(pdos_file, e_fermi)
        if db_center is None:
            continue

        print(f"  V {v_index:3} (slab atom {atom_id:3}): d-band center = {db_center:8.4f} eV")
        rows.append(
            {
                "v_atom_index": v_index,
                "slab_atom_index": atom_id,
                "d_band_center_eV": round(db_center, 6),
            }
        )

    if not rows:
        print("No d-band centers were computed.")
        return 1

    output_path = args.output if args.output is not None else base_dir / DEFAULT_OUTPUT_NAME
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["v_atom_index", "slab_atom_index", "d_band_center_eV"],
        )
        writer.writeheader()
        writer.writerows(rows)

    average = np.mean([row["d_band_center_eV"] for row in rows])
    print(f"  Average: {average:.4f} eV -> saved to {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())