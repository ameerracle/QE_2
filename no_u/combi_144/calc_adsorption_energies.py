#!/usr/bin/env python3
"""Compute adsorption energies from QE output files.

Formula:
E_adsorption = E_combi - E_adsorbate - E_slab

Defaults are set for this repository layout but can be overridden by CLI args.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from ase.io import read

RY_TO_EV = 13.605693009


def energy_from_qe(path: Path) -> tuple[float, str]:
    """Read final QE energy in eV using ASE, with text fallback for robustness."""
    for kwargs in ({"results_required": False}, {}):
        try:
            atoms = read(str(path), index=-1, format="espresso-out", **kwargs)
            return float(atoms.get_potential_energy()), "ase"
        except Exception:
            pass

    text = path.read_text(errors="ignore")
    matches = re.findall(r"!\s+total energy\s+=\s+([-+0-9.eE]+)\s+Ry", text)
    if matches:
        return float(matches[-1]) * RY_TO_EV, "text-fallback"

    raise RuntimeError(f"No energy found in file: {path}")


def pick_combi_energy(candidates: list[Path]) -> tuple[float, str, Path]:
    """Pick first candidate that provides a valid energy (preferred order supplied)."""
    errors = []
    for path in candidates:
        try:
            energy, src = energy_from_qe(path)
            return energy, src, path
        except Exception as exc:
            errors.append(f"{path}: {type(exc).__name__}")
    raise RuntimeError(" ; ".join(errors))


def parse_args() -> argparse.Namespace:
    here = Path(__file__).resolve()
    repo_root = here.parents[2]
    default_combi = repo_root / "no_u/combi_144"
    default_ads = repo_root / "adsorbates"
    default_slab = repo_root / "no_u/slab/output_slab_144"
    default_output = default_combi / "adsorption_energies_from_outs.csv"

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--combi-dir", type=Path, default=default_combi)
    parser.add_argument("--ads-dir", type=Path, default=default_ads)
    parser.add_argument("--slab-dir", type=Path, default=default_slab)
    parser.add_argument("--output", type=Path, default=default_output)
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    combi_dir: Path = args.combi_dir
    ads_dir: Path = args.ads_dir
    slab_dir: Path = args.slab_dir
    output: Path = args.output

    ads_energy: dict[str, tuple[float, str, Path]] = {}
    for p in sorted(ads_dir.glob("*/espresso.pwo")):
        adsorbate = p.parent.name
        e, src = energy_from_qe(p)
        ads_energy[adsorbate] = (e, src, p)

    slab_energy: dict[str, tuple[float, str, Path]] = {}
    for p in sorted(slab_dir.glob("*_slab_relax/*_slab.out")):
        slab = p.name.replace("_slab.out", "")
        e, src = energy_from_qe(p)
        slab_energy[slab] = (e, src, p)

    combi_candidates: dict[str, list[Path]] = {}
    for p in sorted(combi_dir.glob("*_combi_relax/*_combi_restart.out")):
        key = p.parent.name.replace("_combi_relax", "")
        combi_candidates.setdefault(key, []).append(p)
    for p in sorted(combi_dir.glob("*_combi_relax/*_combi.out")):
        key = p.parent.name.replace("_combi_relax", "")
        combi_candidates.setdefault(key, []).append(p)
    for p in sorted(combi_dir.glob("*_combi.out")):
        key = p.name.replace("_combi.out", "")
        combi_candidates.setdefault(key, []).append(p)

    rows: list[dict[str, str | float]] = []
    skipped: list[tuple[str, str]] = []

    for key, candidates in sorted(combi_candidates.items()):
        if "_" not in key:
            skipped.append((key, "bad_key"))
            continue

        slab, adsorbate = key.split("_", 1)
        if slab not in slab_energy or adsorbate not in ads_energy:
            reasons = []
            if slab not in slab_energy:
                reasons.append("missing_slab")
            if adsorbate not in ads_energy:
                reasons.append("missing_ads")
            skipped.append((key, "+".join(reasons)))
            continue

        try:
            e_combi, src_combi, combi_file = pick_combi_energy(candidates)
        except Exception as exc:
            skipped.append((key, f"no_combi_energy ({exc})"))
            continue

        e_ads, src_ads, _ = ads_energy[adsorbate]
        e_slab, src_slab, _ = slab_energy[slab]
        e_adsorption = e_combi - e_ads - e_slab

        rows.append(
            {
                "slab": slab,
                "adsorbate": adsorbate,
                "E_adsorption_eV": e_adsorption,
                "E_combi_eV": e_combi,
                "E_ads_eV": e_ads,
                "E_slab_eV": e_slab,
                "combi_file": str(combi_file),
                "src_combi": src_combi,
                "src_ads": src_ads,
                "src_slab": src_slab,
            }
        )

    print("RESULTS")
    print("slab adsorbate E_adsorption_eV E_combi_eV E_ads_eV E_slab_eV source file")
    for row in rows:
        print(
            f"{row['slab']} {row['adsorbate']} "
            f"{row['E_adsorption_eV']:.6f} {row['E_combi_eV']:.6f} "
            f"{row['E_ads_eV']:.6f} {row['E_slab_eV']:.6f} "
            f"{row['src_combi']} {row['combi_file']}"
        )

    output.write_text(
        "slab,adsorbate,E_adsorption_eV,E_combi_eV,E_ads_eV,E_slab_eV,"
        "combi_file,src_combi,src_ads,src_slab\n"
        + "\n".join(
            f"{r['slab']},{r['adsorbate']},{r['E_adsorption_eV']:.10f},"
            f"{r['E_combi_eV']:.10f},{r['E_ads_eV']:.10f},{r['E_slab_eV']:.10f},"
            f"{r['combi_file']},{r['src_combi']},{r['src_ads']},{r['src_slab']}"
            for r in rows
        )
        + "\n"
    )
    print(f"\nWROTE {output}")

    if rows:
        best = min(rows, key=lambda r: float(r["E_adsorption_eV"]))
        print(
            f"BEST_BINDING {best['slab']} {best['adsorbate']} "
            f"{best['E_adsorption_eV']:.6f} eV"
        )

    if skipped:
        print("\nSKIPPED")
        for key, reason in skipped:
            print(key, reason)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
