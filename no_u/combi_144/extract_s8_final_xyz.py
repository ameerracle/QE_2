#!/usr/bin/env python3
"""
Extract final relaxed combination structures from QE outputs and write XYZ files.

Selection priority per system:
1) {tag}_combi_re.out
2) Adsorbate-specific relaxed outputs (for example: {tag}_combi_re_Li2S8.out, re_Li2S8.out)
3) Newest relaxed variant among {tag}_combi_re_*.out and re_*.out
4) {tag}_combi_restartN.out (highest N first: restart4 > restart3 > ...)
5) {tag}_combi_restart.out
6) {tag}_combi.out

This script uses ASE's QE output reader (ase.io.read, format='espresso-out')
to take the last structure from the selected output file.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from ase.io import read, write


DEFAULT_ADS = ["S8", "Li2S4", "Li2S8"]


def pick_priority_output(run_dir: Path, tag: str) -> Path | None:
    """Return the highest-priority QE output that exists for a given tag."""
    exact_re = run_dir / f"{tag}_combi_re.out"
    if exact_re.exists():
        return exact_re

    ads_hint = tag.split("_")[-1] if "_" in tag else tag
    preferred_re = [
        run_dir / f"{tag}_combi_re_{ads_hint}.out",
        run_dir / f"{tag}_combi_re_{ads_hint.lower()}.out",
        run_dir / f"re_{ads_hint}.out",
        run_dir / f"re_{ads_hint.lower()}.out",
    ]
    for candidate in preferred_re:
        if candidate.exists():
            return candidate

    # Handles variants such as *_combi_re_Li2S8.out and picks the newest one.
    re_candidates = list(run_dir.glob(f"{tag}_combi_re_*.out")) + list(run_dir.glob("re_*.out"))
    if re_candidates:
        re_candidates.sort(key=lambda p: p.stat().st_mtime, reverse=True)
        return re_candidates[0]

    restart_n = []
    for candidate in run_dir.glob(f"{tag}_combi_restart*.out"):
        m = re.fullmatch(rf"{re.escape(tag)}_combi_restart(\d+)\.out", candidate.name)
        if m:
            restart_n.append((int(m.group(1)), candidate))

    if restart_n:
        restart_n.sort(key=lambda x: x[0], reverse=True)
        return restart_n[0][1]

    restart = run_dir / f"{tag}_combi_restart.out"
    if restart.exists():
        return restart

    base = run_dir / f"{tag}_combi.out"
    if base.exists():
        return base

    return None


def read_final_atoms(qe_out: Path):
    """Read the final structure from a QE pw.x output using ASE."""
    try:
        return read(str(qe_out), format="espresso-out", index=-1)
    except Exception:
        # Fallback to ASE auto-detection if explicit QE parser fails.
        return read(str(qe_out), index=-1)


def extract_structures(
    base_dir: Path, out_dir: Path, ads_list: list[str], extract_all: bool = False
) -> tuple[int, int, int]:
    """Extract final structures in *_combi_relax folders."""
    if extract_all:
        run_dirs = sorted(d for d in base_dir.glob("*_combi_relax") if d.is_dir())
    else:
        run_dirs = []
        for ads in ads_list:
            run_dirs.extend(d for d in base_dir.glob(f"*_{ads}_combi_relax") if d.is_dir())
        run_dirs = sorted(set(run_dirs))

    if not run_dirs:
        print(f"[INFO] No matching '*_ADS_combi_relax' folders found in: {base_dir}")
        return 0, 0, 0

    out_dir.mkdir(parents=True, exist_ok=True)

    ok = 0
    skipped = 0
    failed = 0

    for run_dir in run_dirs:
        tag = run_dir.name.removesuffix("_combi_relax")
        chosen_out = pick_priority_output(run_dir, tag)

        if chosen_out is None:
            print(f"[SKIP] {tag}: none of priority outputs found")
            skipped += 1
            continue

        try:
            atoms = read_final_atoms(chosen_out)
            xyz_path = out_dir / f"{tag}_combi_final.xyz"
            write(str(xyz_path), atoms)
            print(f"[OK] {tag}: {chosen_out.name} -> {xyz_path.name}")
            ok += 1
        except Exception as exc:
            print(f"[FAIL] {tag}: could not parse {chosen_out.name} ({exc})")
            failed += 1

    return ok, skipped, failed


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Extract final combination structures from QE outputs and write XYZ files"
    )
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Directory containing *_<ads>_combi_relax folders",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "PDOS" / "final_combi",
        help="Output directory for generated XYZ files",
    )
    parser.add_argument(
        "--ads",
        nargs="+",
        default=DEFAULT_ADS,
        help="Adsorbate names to extract (default: S8 Li2S4 Li2S8)",
    )
    parser.add_argument(
        "-all",
        "--all",
        action="store_true",
        help="Extract from all *_combi_relax folders (ignores --ads filter)",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()

    ok, skipped, failed = extract_structures(
        args.base_dir, args.out_dir, args.ads, extract_all=args.all
    )
    print("\nSummary")
    print(f"  extracted: {ok}")
    print(f"  skipped:   {skipped}")
    print(f"  failed:    {failed}")


if __name__ == "__main__":
    main()
