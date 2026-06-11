#!/usr/bin/env python3
"""
bader_analysis_VN_U.py -- Bader charge analysis for VN+U slab + adsorbate systems
Run from: /lustre10/scratch/anizami/QE_2/

Uses plot_num=0 cube files (valence charge with PAW augmentation).
delta_q = q_combi - q_slab - q_adsorbate (nspin2).
No pandas dependency.
"""

import argparse
import subprocess
import numpy as np
import csv
from pathlib import Path

BASE_DIR = Path("/lustre10/scratch/anizami/QE_2")
BADER_BIN = BASE_DIR / "bader"
SLAB_DIR = BASE_DIR / "u/PDOS/VN_slab"
ADS_BASE = BASE_DIR / "adsorbates/PDOS/nspin2"
COMBI_BASE = BASE_DIR / "u/PDOS"
OUTPUT_DIR = BASE_DIR / "bader_output"
ACF_DIR = OUTPUT_DIR / "ACF"

N_SLAB = 144
N_METAL = 72
Z_SURFACE_CUTOFF = 2.0

ADSORBATE_ATOMS = {
    "Li2S4": {"S": 4, "Li": 2},
    "Li2S8": {"S": 8, "Li": 2},
    "S8":    {"S": 8, "Li": 0},
}

COMBI_DIR_MAP = {
    "Li2S4": "VN_Li2S4_adsorbate",
    "Li2S8": "VN_Li2S8_adsorbate",
    "S8":    "VN_S8_adsorbate",
}

Z_TO_SYM = {
    3: "Li", 7: "N", 16: "S",
    21: "Sc", 22: "Ti", 23: "V",
    40: "Zr", 41: "Nb"
}

BOHR_TO_ANG = 0.529177


def find_cube_file(directory, preferred_name=None):
    directory = Path(directory)
    if not directory.exists():
        return None
    if preferred_name is not None:
        preferred = directory / preferred_name
        if preferred.exists():
            return preferred
    cubes = sorted(directory.glob("*_charge_plot0.cube"))
    if cubes:
        return cubes[0]
    cubes = sorted(directory.glob("*_charge.cube"))
    if cubes:
        return cubes[0]
    return None


def run_bader(cube_path, acf_label, output_dir, force=False):
    cube_path = Path(cube_path)
    acf_final = output_dir / f"ACF_{acf_label}.dat"
    if acf_final.exists():
        if force:
            acf_final.unlink()
            print(f"  ACF_{acf_label}.dat exists, forcing rerun: {cube_path.name}")
        else:
            print(f"  ACF_{acf_label}.dat exists, skipping: {cube_path.name}")
            return acf_final
    stale = cube_path.parent / "ACF.dat"
    if stale.exists():
        stale.unlink()
    print(f"  Running bader on {cube_path.name} ...")
    result = subprocess.run(
        [str(BADER_BIN), str(cube_path)],
        cwd=str(cube_path.parent),
        capture_output=True, text=True
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
    return z_coords >= (np.max(z_coords) - cutoff)


def get_acf(cube_path, acf_label, force=False):
    acf_final = ACF_DIR / f"ACF_{acf_label}.dat"
    if acf_final.exists() and not force:
        print(f"  ACF_{acf_label}.dat exists, skipping bader")
        return acf_final
    if cube_path is not None and cube_path.exists():
        return run_bader(cube_path, acf_label, ACF_DIR, force=force)
    if acf_final.exists():
        return acf_final
    print(f"  ERROR: no cube and no ACF for {acf_label}")
    return None


def analyze_VN_U_system(adsorbate, force=False):
    system_key = f"VN_U_{adsorbate}"
    print(f"\n[{system_key}]")

    ads_cfg = ADSORBATE_ATOMS[adsorbate]
    n_s = ads_cfg["S"]
    n_li = ads_cfg["Li"]
    n_ads = n_s + n_li

    combi_dirname = COMBI_DIR_MAP[adsorbate]
    combi_dir = COMBI_BASE / combi_dirname
    combi_cube = find_cube_file(combi_dir, preferred_name=f"{combi_dirname}_charge_plot0.cube")

    slab_cube = find_cube_file(SLAB_DIR, preferred_name="VN_slab_charge_plot0.cube")

    ads_dir = ADS_BASE / adsorbate
    ads_cube = find_cube_file(ads_dir, preferred_name=f"{adsorbate}_charge_plot0.cube")

    acf_combi = get_acf(combi_cube, f"VN_{adsorbate}_u_combi", force=force)
    acf_slab = get_acf(slab_cube, "VN_slab_u", force=force)
    acf_ads = get_acf(ads_cube, f"{adsorbate}_nspin2", force=force)

    if any(x is None for x in [acf_combi, acf_slab, acf_ads]):
        return None

    q_combi = parse_acf(acf_combi)
    q_slab = parse_acf(acf_slab)
    q_ads = parse_acf(acf_ads)

    if len(q_combi) != N_SLAB + n_ads:
        print(f"  WARNING: combi ACF has {len(q_combi)} atoms, expected {N_SLAB + n_ads}")
    if len(q_slab) != N_SLAB:
        print(f"  WARNING: slab ACF has {len(q_slab)} atoms, expected {N_SLAB}")
    if len(q_ads) != n_ads:
        print(f"  WARNING: ads ACF has {len(q_ads)} atoms, expected {n_ads}")

    if slab_cube is not None and slab_cube.exists():
        _, z_slab = parse_cube_atoms(slab_cube)
        surf_mask = get_surface_mask(z_slab)
    else:
        print(f"  WARNING: no slab cube for z-coords, using all atoms as surface")
        surf_mask = np.ones(N_SLAB, dtype=bool)

    rows = []

    for i in range(N_METAL):
        if not surf_mask[i]:
            continue
        rows.append({
            "atom_index": i + 1,
            "element": "V",
            "region": "surface_metal",
            "q_combi": q_combi[i],
            "q_slab": q_slab[i],
            "q_ads": 0.0,
            "delta_q": q_combi[i] - q_slab[i],
        })

    for i in range(N_METAL, N_SLAB):
        if not surf_mask[i]:
            continue
        rows.append({
            "atom_index": i + 1,
            "element": "N",
            "region": "surface_N",
            "q_combi": q_combi[i],
            "q_slab": q_slab[i],
            "q_ads": 0.0,
            "delta_q": q_combi[i] - q_slab[i],
        })

    for j in range(n_s):
        i_combi = N_SLAB + j
        rows.append({
            "atom_index": i_combi + 1,
            "element": "S",
            "region": "adsorbate_S",
            "q_combi": q_combi[i_combi],
            "q_slab": 0.0,
            "q_ads": q_ads[j],
            "delta_q": q_combi[i_combi] - q_ads[j],
        })

    for j in range(n_li):
        i_combi = N_SLAB + n_s + j
        rows.append({
            "atom_index": i_combi + 1,
            "element": "Li",
            "region": "adsorbate_Li",
            "q_combi": q_combi[i_combi],
            "q_slab": 0.0,
            "q_ads": q_ads[n_s + j],
            "delta_q": q_combi[i_combi] - q_ads[n_s + j],
        })

    out_path = OUTPUT_DIR / f"VN_U_{adsorbate}_bader.csv"
    fieldnames = ["atom_index", "element", "region", "q_combi", "q_slab", "q_ads", "delta_q"]
    float_cols = {"q_combi", "q_slab", "q_ads", "delta_q"}
    with open(out_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            formatted = {}
            for k, v in row.items():
                if k in float_cols:
                    formatted[k] = f"{v:.4f}"
                else:
                    formatted[k] = v
            writer.writerow(formatted)
    print(f"  Saved: {out_path}")

    # Print summary
    regions = {}
    for row in rows:
        r = row["region"]
        if r not in regions:
            regions[r] = {"q_combi": [], "q_slab": [], "q_ads": [], "delta_q": []}
        regions[r]["q_combi"].append(row["q_combi"])
        regions[r]["q_slab"].append(row["q_slab"])
        regions[r]["q_ads"].append(row["q_ads"])
        regions[r]["delta_q"].append(row["delta_q"])
    for r, vals in regions.items():
        print(f"  {r}: q_combi={np.mean(vals['q_combi']):.4f}, q_slab={np.mean(vals['q_slab']):.4f}, q_ads={np.mean(vals['q_ads']):.4f}, delta_q={np.mean(vals['delta_q']):.4f}")
    return rows


def main():
    parser = argparse.ArgumentParser(description="Run Bader analysis for VN+U slab + adsorbate systems")
    parser.add_argument("--adsorbates", nargs="+", choices=list(ADSORBATE_ATOMS.keys()),
                        default=list(ADSORBATE_ATOMS.keys()))
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    ACF_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Processing VN+U with adsorbates: {args.adsorbates}")
    for ads in args.adsorbates:
        analyze_VN_U_system(ads, force=args.force)
    print("\nDone.")


if __name__ == "__main__":
    main()
