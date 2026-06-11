"""
bader_combi.py -- Bader charge analysis for combi systems
Run from: /lustre10/scratch/anizami/QE_2/

Workflow:
  1. Runs `bader <cube>` in combi directory
  2. Parses ACF.dat for combi system
  3. Outputs one CSV per metal+adsorbate combo

Atom layout (combi, 150 atoms):
  1-72:   metal
  73-144: N
  145-148: S  (4 atoms for Li2S4; varies for Li2S8/S8)
  149-150: Li (2 atoms for Li2S4; varies for Li2S8/S8)

Adsorbate atom counts:
  Li2S4: 4S + 2Li = 6
  Li2S8: 8S + 2Li = 10
  S8:    8S + 0Li = 8
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

Z_SURFACE_CUTOFF = 2.0   # angstrom -- atoms within this of max_z = surface layer

N_SLAB = 144
N_METAL = 72
N_N = 72

ADSORBATE_ATOMS = {
    "Li2S4": {"S": 4, "Li": 2},
    "Li2S8": {"S": 8, "Li": 2},
    "S8":    {"S": 8, "Li": 0},
}

METALS = ["ScN", "TiN", "VN", "NbN", "ZrN"]
METAL_SYMBOL = {
    "ScN": "Sc", "TiN": "Ti", "VN": "V", "NbN": "Nb", "ZrN": "Zr"
}

COMBI_BASE  = PDOS_BASE


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

# --------------------------------------------------------------------------
# Bader runner
# --------------------------------------------------------------------------

def run_bader(cube_path, acf_label=None, force=False):
    """Run bader on a cube file. Returns path to (renamed) ACF file."""
    cube_path = Path(cube_path)
    label     = acf_label or cube_path.stem   # e.g. "Li2S8" or "ScN_Li2S4"
    acf_final = cube_path.parent / f"ACF_{label}.dat"

    if acf_final.exists():
        if force:
            acf_final.unlink()
            print(f"  ACF_{label}.dat exists, forcing rerun: {cube_path.name}")
        else:
            print(f"  ACF_{label}.dat exists, skipping: {cube_path.name}")
            return acf_final

    # remove any stale generic ACF.dat before running
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

# --------------------------------------------------------------------------
# Parsers
# --------------------------------------------------------------------------

def parse_acf(acf_path):
    """
    Parse ACF.dat, return numpy array of charges (1-indexed rows -> 0-indexed array).
    ACF.dat format:
      # X Y Z CHARGE MIN_DIST ATOMIC_VOL
      1  x  y  z  charge  ...
      ...
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


def parse_cube_atoms(cube_path):
    """
    Parse atom symbols and z-coordinates from cube file header.
    Returns (symbols list, z_coords array) both length N_atoms.

    Cube format:
      line 1: comment
      line 2: comment
      line 3: natoms  ox  oy  oz
      line 4-6: grid lines
      line 7+: Z_atomic  charge  x  y  z  (one per atom, in Bohr)
    """
    BOHR_TO_ANG = 0.529177

    # Atomic number -> symbol map (only elements we care about)
    Z_TO_SYM = {
        3: "Li", 7: "N", 16: "S",
        21: "Sc", 22: "Ti", 23: "V",
        40: "Zr", 41: "Nb"
    }

    symbols = []
    z_coords = []

    with open(cube_path) as f:
        lines = f.readlines()

    natoms = abs(int(lines[2].split()[0]))

    for i in range(6, 6 + natoms):
        parts = lines[i].split()
        z_atomic = int(float(parts[0]))
        sym = Z_TO_SYM.get(z_atomic, f"Z{z_atomic}")
        z_bohr = float(parts[4])
        z_coords.append(z_bohr * BOHR_TO_ANG)
        symbols.append(sym)

    return symbols, np.array(z_coords)


def get_surface_mask(z_coords, cutoff=Z_SURFACE_CUTOFF):
    """Return boolean mask for atoms within cutoff of max z."""
    return z_coords >= (np.max(z_coords) - cutoff)


# --------------------------------------------------------------------------
# Main analysis
# --------------------------------------------------------------------------

def analyze_system(metal_slab, adsorbate, force=False):
    """
    For a given metal slab (e.g. 'ScN') and adsorbate (e.g. 'Li2S4'),
    run bader on combi system and output charges.
    """
    system_key = f"{metal_slab}_{adsorbate}"
    print(f"\n[{system_key}]")

    ads_cfg = ADSORBATE_ATOMS[adsorbate]
    n_s  = ads_cfg["S"]
    n_li = ads_cfg["Li"]
    n_ads = n_s + n_li

    # Paths
    combi_dir  = COMBI_BASE / system_key
    combi_cube = find_cube_file(combi_dir, preferred_name=f"{system_key}_charge.cube")

    if combi_cube is None or not combi_cube.exists():
        print(f"  ERROR: combi cube not found in expected location")
        return None

    # Run bader
    acf_combi = run_bader(combi_cube, acf_label=system_key, force=force)

    if acf_combi is None:
        return None

    # Parse charges
    q_combi = parse_acf(acf_combi)

    # Verify atom counts
    if len(q_combi) != N_SLAB + n_ads:
        print(f"  WARNING: combi ACF has {len(q_combi)} atoms, expected {N_SLAB + n_ads}")

    # Parse atom positions from combi cube for z-filter
    symbols_combi, z_combi = parse_cube_atoms(combi_cube)

    # Surface layer mask (slab atoms only, indices 0..143)
    z_slab_only = z_combi[:N_SLAB]
    surf_mask   = get_surface_mask(z_slab_only)

    # Build rows
    rows = []

    # --- Surface metal atoms ---
    metal_sym = METAL_SYMBOL[metal_slab]
    for i in range(N_METAL):
        if not surf_mask[i]:
            continue
        rows.append({
            "atom_index": i + 1,
            "element":    metal_sym,
            "region":     "surface_metal",
            "q_combi":    q_combi[i],
        })

    # --- Surface N atoms ---
    for i in range(N_METAL, N_SLAB):
        if not surf_mask[i]:
            continue
        rows.append({
            "atom_index": i + 1,
            "element":    "N",
            "region":     "surface_N",
            "q_combi":    q_combi[i],
        })

    # --- Adsorbate S atoms ---
    for j in range(n_s):
        i_combi = N_SLAB + j
        rows.append({
            "atom_index": i_combi + 1,
            "element":    "S",
            "region":     "adsorbate_S",
            "q_combi":    q_combi[i_combi],
        })

    # --- Adsorbate Li atoms ---
    for j in range(n_li):
        i_combi = N_SLAB + n_s + j
        rows.append({
            "atom_index": i_combi + 1,
            "element":    "Li",
            "region":     "adsorbate_Li",
            "q_combi":    q_combi[i_combi],
        })

    df = pd.DataFrame(rows)

    # Output CSV
    out_path = combi_dir / f"{system_key}_bader.csv"
    df.to_csv(out_path, index=False, float_format="%.4f")
    print(f"  Saved: {out_path}")
    print(df.groupby("region")[["q_combi"]].mean().round(4))

    return df


# --------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Bader analysis for slab + adsorbate systems")
    parser.add_argument("--force", action="store_true", help="Delete existing ACF_*.dat files and rerun Bader")
    parser.add_argument(
        "--metals",
        nargs="+",
        choices=METALS,
        default=METALS,
        help="Metals to process. Default: all metals.",
    )
    parser.add_argument(
        "--adsorbates",
        nargs="+",
        choices=list(ADSORBATE_ATOMS.keys()),
        default=list(ADSORBATE_ATOMS.keys()),
        help="Adsorbates to process. Default: all adsorbates.",
    )
    args = parser.parse_args()

    print(f"Processing combi systems only (slab and adsorbate files not available)")
    for metal in args.metals:
        for ads in args.adsorbates:
            analyze_system(metal, ads, force=args.force)
    print("\nAll done.")
