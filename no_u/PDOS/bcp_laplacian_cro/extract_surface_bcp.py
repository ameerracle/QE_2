#!/usr/bin/env python3
"""Extract surface BCPs and Laplacian values from Critic2 .cro outputs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

BOHR_TO_ANG = 0.52917721092
DEFAULT_KEYWORDS = ("(3,-1)",)
DEFAULT_KEYWORDS_STR = ";".join(DEFAULT_KEYWORDS)
TF_PREFACTOR = (3.0 / 10.0) * (3.0 * 3.141592653589793 ** 2) ** (2.0 / 3.0)


def parse_keywords(raw: str) -> list[str]:
    parts = [part.strip().lower() for part in raw.split(";")]
    return [part for part in parts if part]


def cro_files_from_input(path: Path) -> list[Path]:
    if path.is_dir():
        return sorted(path.rglob("*.cro"))
    return [path]


def parse_cro(path: Path, keywords: list[str], min_z_ang: float, max_z_ang: float) -> list[dict[str, float | int | str]]:
    results: list[dict[str, float | int | str]] = []
    in_section = False
    current: dict[str, object] | None = None

    def flush_current() -> None:
        nonlocal current
        if current is None:
            return
        cp_type = current.get("type")
        coords = current.get("coords")
        lap = current.get("lap")
        rho = current.get("rho")
        if not isinstance(cp_type, str) or not isinstance(coords, tuple) or lap is None or rho is None:
            current = None
            return
        cp_type_lower = cp_type.lower()
        if cp_type_lower not in keywords:
            current = None
            return
        x_bohr, y_bohr, z_bohr = coords
        z_ang = z_bohr * BOHR_TO_ANG
        if z_ang < min_z_ang or z_ang > max_z_ang:
            current = None
            return
        rho_f = float(rho)
        lap_f = float(lap)
        g_ked = TF_PREFACTOR * rho_f ** (5.0 / 3.0) + lap_f / 6.0
        v_ped = -lap_f / 4.0 - g_ked
        h_ted = g_ked + v_ped
        abs_v_over_g = abs(v_ped) / g_ked if g_ked != 0 else float("inf")
        results.append(
            {
                "cp_no": int(current["cp_no"]),
                "cp_type": cp_type,
                "x_bohr": x_bohr,
                "y_bohr": y_bohr,
                "z_bohr": z_bohr,
                "x_ang": x_bohr * BOHR_TO_ANG,
                "y_ang": y_bohr * BOHR_TO_ANG,
                "z_ang": z_ang,
                "rho": rho_f,
                "laplacian": lap_f,
                "G_ked": g_ked,
                "V_ped": v_ped,
                "H_ted": h_ted,
                "abs_V_over_G": abs_v_over_g,
            }
        )
        current = None

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not in_section:
                if line.startswith("* Additional properties at the critical points"):
                    in_section = True
                continue

            if line.startswith("+ Critical point no."):
                flush_current()
                parts = line.strip().split()
                cp_no = int(parts[-1])
                current = {"cp_no": cp_no, "type": None, "coords": None, "lap": None, "rho": None}
                continue

            if current is None:
                continue

            if "Type :" in line:
                current["type"] = line.split(":", 1)[1].strip()
                continue

            if "Cartesian coordinates (bohr):" in line:
                coords_text = line.split(":", 1)[1].strip().split()
                if len(coords_text) >= 3:
                    coords = tuple(float(value) for value in coords_text[:3])
                    current["coords"] = coords
                continue

            if "Laplacian (del2 f):" in line:
                value = float(line.split(":", 1)[1].strip())
                current["lap"] = value
                continue

            if line.strip().startswith("Field value (f):"):
                value = float(line.split(":", 1)[1].strip())
                current["rho"] = value
                continue

    flush_current()
    return results


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract surface BCP coordinates and Laplacian values from Critic2 .cro files."
    )
    parser.add_argument(
        "input",
        nargs="?",
        default="no_u/PDOS/bcp_laplacian_cro",
        help=".cro file or directory to scan (default: no_u/PDOS/bcp_laplacian_cro).",
    )
    parser.add_argument(
        "--min-z-ang",
        type=float,
        default=15.5,
        help="Minimum z (Angstrom) to keep a critical point (default: 15.5).",
    )
    parser.add_argument(
        "--max-z-ang",
        type=float,
        default=20.5,
        help="Maximum z (Angstrom) to keep a critical point (default: 20.5).",
    )
    parser.add_argument(
        "--type-keywords",
        default=DEFAULT_KEYWORDS_STR,
        help="Semicolon-separated substrings that mark a BCP (default: (3,-1)).",
    )
    parser.add_argument(
        "--output",
        default="surface_bcp_laplacian.csv",
        help="Output CSV file (default: surface_bcp_laplacian.csv).",
    )
    args = parser.parse_args()

    input_path = Path(args.input).expanduser()
    keywords = parse_keywords(args.type_keywords)
    output_path = Path(args.output).expanduser()

    fieldnames = [
        "cp_no",
        "cp_type",
        "x_bohr",
        "y_bohr",
        "z_bohr",
        "x_ang",
        "y_ang",
        "z_ang",
        "rho",
        "laplacian",
        "G_ked",
        "V_ped",
        "H_ted",
        "abs_V_over_G",
    ]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    total = 0
    for cro_path in cro_files_from_input(input_path):
        entries = parse_cro(cro_path, keywords, args.min_z_ang, args.max_z_ang)
        if not entries:
            continue
        entries.sort(key=lambda row: int(row["cp_no"]))
        stem = cro_path.stem
        csv_path = cro_path.parent / f"{stem}_bcp_laplacian.csv"
        with csv_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(entries)
        print(f"  {stem}: {len(entries)} BCP(s) -> {csv_path}")
        total += len(entries)

    print(f"Wrote {total} surface BCP(s) across {output_path.parent}/")


if __name__ == "__main__":
    main()
