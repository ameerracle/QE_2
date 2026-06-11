#!/usr/bin/env python3
"""
bader_slab.py -- Bader charge analysis for TMN slab systems only
Run from: /lustre10/scratch/anizami/QE_2/

Workflow:
  1. Runs `bader <cube>` in each slab directory
  2. Parses ACF.dat for each slab
  3. Outputs one CSV per metal slab

Atom layout (slab, 144 atoms):
  1-72:   metal
  73-144: N
"""

import argparse
import subprocess
import numpy as np
import pandas as pd
from pathlib import Path

# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------
BADER_BIN = "/lustre10/scratch/anizami/QE_2/bader"
BASE_DIR = Path("/lustre10/scratch/anizami/QE_2")
PDOS_BASE = BASE_DIR / "no_u/PDOS"

METALS = ["ScN", "TiN", "VN", "NbN", "ZrN"]
N_SLAB = 144
N_METAL = 72
N_N = 72


def find_cube_file(directory, preferred_name=None):
    """Return the cube file in a directory, preferring an exact filename when given."""
    directory = Path(directory)
    if not directory.exists():
        return None

    if preferred_name is not None:
        preferred = directory / preferred_name
        if preferred.exists():
            return preferred

    cubes = sorted(directory.glob("*.cube"))
    if len(cubes) == 1:
        return cubes[0]
    if cubes:
        return cubes[0]
    return None


def run_bader(cube_path, acf_label=None, force=False):
    """Run bader on a cube file. Returns path to (renamed) ACF file."""
    cube_path = Path(cube_path)
    label = acf_label or cube_path.stem
    acf_final = cube_path.parent / f"ACF_{label}.dat"

    if acf_final.exists():
        if force:
            acf_final.unlink()
            print(f"  ACF_{label}.dat exists, forcing rerun: {cube_path.name}")
        else:
            print(f"  ACF_{label}.dat exists, skipping: {cube_path.name}")
            return acf_final

    stale = cube_path.parent / "ACF.dat"
    if stale.exists():
        stale.unlink()

    print(f"  Running bader on {cube_path.name} ...")
    result = subprocess.run(
        [BADER_BIN, str(cube_path)],
        cwd=str(cube_path.parent),
        capture_output=True,
        text=True
    )
    if result.returncode != 0:
        print(f"  ERROR: bader failed on {cube_path}")
        print(result.stderr)
        return None

    raw_acf = cube_path.parent / "ACF.dat"
    if not raw_acf.exists():
        print(f"  ERROR: ACF.dat not produced in {cube_path.parent}")
        return None

    raw_acf.rename(acf_final)
    return acf_final


def parse_acf(acf_path):
    """
    Parse ACF.dat, return numpy array of charges (1-indexed rows -> 0-indexed array).
    """
    charges = []
    with open(acf_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("-"):
                continue
            parts = line.split()
            if not parts[0].isdigit():
                continue
            charges.append(float(parts[4]))
    return np.array(charges)


def analyze_slab(metal, force=False):
    """
    For a given metal slab (e.g. 'ScN'), run bader and output charges.
    """
    print(f"\n[{metal}]")

    slab_dir = PDOS_BASE / f"{metal}_slab"
    slab_cube = find_cube_file(slab_dir, preferred_name=f"{metal}_charge.cube")

    if slab_cube is None or not slab_cube.exists():
        print(f"  ERROR: slab cube not found in {slab_dir}")
        return None

    acf_slab = run_bader(slab_cube, acf_label=metal, force=force)

    if acf_slab is None:
        return None

    q_slab = parse_acf(acf_slab)

    if len(q_slab) != N_SLAB:
        print(f"  WARNING: slab ACF has {len(q_slab)} atoms, expected {N_SLAB}")

    rows = []

    for i in range(N_METAL):
        rows.append({
            "atom_index": i + 1,
            "element": metal.replace("N", ""),
            "region": "metal",
            "q_slab": q_slab[i],
        })

    for i in range(N_METAL, N_SLAB):
        rows.append({
            "atom_index": i + 1,
            "element": "N",
            "region": "N",
            "q_slab": q_slab[i],
        })

    df = pd.DataFrame(rows)

    out_path = slab_dir / f"{metal}_slab_bader.csv"
    df.to_csv(out_path, index=False, float_format="%.4f")
    print(f"  Saved: {out_path}")
    print(df.groupby("region")[["q_slab"]].mean().round(4))

    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Bader analysis for slab systems")
    parser.add_argument("--force", action="store_true", help="Delete existing ACF_*.dat files and rerun Bader")
    parser.add_argument(
        "--metals",
        nargs="+",
        choices=METALS,
        default=METALS,
        help="Metals to process. Default: all metals.",
    )
    args = parser.parse_args()

    for metal in args.metals:
        analyze_slab(metal, force=args.force)
    print("\nAll done.")
