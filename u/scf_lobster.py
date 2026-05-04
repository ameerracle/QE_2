#!/usr/bin/env python3
"""
Generate and optionally run QE SCF + NSCF (LOBSTER-prep) for VN combined structures
with DFT+U enabled.

Uses structures from:
  /lustre10/scratch/anizami/QE_2/u/combi/final_xyz

Writes to:
  /lustre10/scratch/anizami/QE_2/u/lobster
"""

from __future__ import annotations

import argparse
import os
import subprocess
from pathlib import Path

from ase.io import read

XYZ_DIR = Path("/lustre10/scratch/anizami/QE_2/u/combi/final_xyz")
RUN_ROOT = Path("/lustre10/scratch/anizami/QE_2/u/lobster")
PSEUDO_DIR = "/lustre10/scratch/anizami/QE_2/PAW_pslib"
ECUTWFC = 45.0
ECUTRHO = 450.0
HUBBARD_U_V_3D = 2.50

PSEUDOS = {
    "Ti": "Ti.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "V": "V.pbe-spnl-kjpaw_psl.1.0.0.UPF",
    "Sc": "Sc.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Nb": "Nb.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Zr": "Zr.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "N": "N.pbe-n-kjpaw_psl.1.0.0.UPF",
    "Li": "Li.pbe-s-kjpaw_psl.1.0.0.UPF",
    "S": "S.pbe-n-kjpaw_psl.1.0.0.UPF",
    "C": "C.pbe-n-kjpaw_psl.1.0.0.UPF",
}

TARGET_SLAB = "VN"


def ensure_directory_exists(path: Path, label: str) -> None:
    if not path.exists():
        path.mkdir(parents=True, exist_ok=True)
        print(f"Created missing {label}: {path}")


def calculate_min_nbnd(atoms) -> int:
    total_orbitals = 0
    for atom in atoms:
        if atom.symbol == "V":
            total_orbitals += 9
        elif atom.symbol in {"N", "Li", "S", "C"}:
            total_orbitals += 4
        else:
            total_orbitals += 4
    return int(total_orbitals * 1.1)


def ordered_species(tag: str, atoms) -> list[str]:
    slab = tag.split("_", 1)[0]
    metal = slab.rstrip("N")
    present = set(atoms.get_chemical_symbols())
    species: list[str] = []

    for sym in [metal, "N", "Li", "S", "C"]:
        if sym in present:
            species.append(sym)

    for sym in atoms.get_chemical_symbols():
        if sym not in species:
            species.append(sym)

    return species


def write_qe_input(file_path: Path, atoms, tag: str, degauss: float, k_pts_str: str, calc_type: str) -> None:
    nat = len(atoms)
    species = ordered_species(tag, atoms)
    ntyp = len(species)
    k_pts = k_pts_str.split()

    if len(k_pts) != 6:
        raise ValueError("K grid must contain 6 integers: kx ky kz sx sy sz")

    if "V" not in species:
        raise ValueError(f"This script is VN-only; V not found in species for {tag}")
    v_index = species.index("V") + 1

    with open(file_path, "w", encoding="utf-8") as f:
        f.write("&CONTROL\n")
        f.write(f"  calculation = '{calc_type}'\n")
        f.write(f"  prefix = '{tag}'\n")
        f.write(f"  pseudo_dir = '{PSEUDO_DIR}'\n")
        f.write("  outdir = './tmp/'\n")
        f.write("  disk_io = 'high'\n")
        f.write("  verbosity = 'high'\n")
        f.write("  wf_collect = .true.\n")
        f.write("/\n")

        f.write("&SYSTEM\n")
        f.write(f"  ibrav = 0, nat = {nat}, ntyp = {ntyp}\n")
        f.write(f"  ecutwfc = {ECUTWFC:.1f}, ecutrho = {ECUTRHO:.1f}\n")

        if calc_type == "nscf":
            f.write(f"  nbnd = {calculate_min_nbnd(atoms)}\n")
            f.write("  nosym = .true.\n")
            f.write("  noinv = .true.\n")

        f.write(f"  occupations = 'smearing', smearing = 'cold', degauss = {degauss:.5f}\n")
        f.write("  nspin = 2\n")
        f.write(f"  starting_magnetization({v_index}) = 1.0\n")
        f.write("/\n")

        f.write("&ELECTRONS\n")
        f.write("  conv_thr = 1.0d-7\n")
        f.write("  mixing_beta = 0.25\n")
        if calc_type == "nscf":
            f.write("  diago_david_ndim = 6\n")
            f.write("  diago_thr_init = 1.0d-4\n")
            f.write("  electron_maxstep = 300\n")
        f.write("/\n")

        f.write("\nCELL_PARAMETERS (angstrom)\n")
        for vec in atoms.get_cell():
            f.write(f"  {vec[0]:15.10f} {vec[1]:15.10f} {vec[2]:15.10f}\n")

        f.write("\nATOMIC_SPECIES\n")
        for sym in species:
            if sym not in PSEUDOS:
                raise KeyError(f"Missing pseudopotential for '{sym}'")
            f.write(f"  {sym:<3} 1.0  {PSEUDOS[sym]}\n")

        f.write("\nATOMIC_POSITIONS (angstrom)\n")
        for atom in atoms:
            f.write(f"  {atom.symbol:<3} {atom.position[0]:14.9f} {atom.position[1]:14.9f} {atom.position[2]:14.9f}\n")

        f.write("\nK_POINTS (automatic)\n")
        f.write(f"  {k_pts[0]} {k_pts[1]} {k_pts[2]} {k_pts[3]} {k_pts[4]} {k_pts[5]}\n")

        f.write("\nHUBBARD {atomic}\n")
        f.write(f"U V-3d {HUBBARD_U_V_3D}\n")


def write_pp_ae_density(file_path: Path, tag: str) -> None:
    with open(file_path, "w", encoding="utf-8") as f:
        f.write("&INPUTPP\n")
        f.write(f"  prefix = '{tag}'\n")
        f.write("  outdir = './tmp/'\n")
        f.write("  plot_num = 21\n")
        f.write(f"  filplot = '{tag}_ae_rho'\n")
        f.write("/\n")
        f.write("&PLOT\n")
        f.write("  iflag = 3\n")
        f.write("  output_format = 6\n")
        f.write(f"  fileout = '{tag}_AE_density.cube'\n")
        f.write("/\n")


def run_chain(run_dir: Path, tag: str, np: int) -> None:
    is_slurm = "SLURM_JOB_ID" in os.environ
    pw_cmd = f"srun pw.x -nk 2" if is_slurm else f"mpirun -np {np} pw.x -nk 2"
    pp_cmd = f"srun pp.x" if is_slurm else f"mpirun -np {np} pp.x"

    steps = [
        ("SCF", f"{pw_cmd} < {tag}_scf.in > {tag}_scf.out"),
        ("NSCF_LOBSTER", f"{pw_cmd} < {tag}_nscf.in > {tag}_nscf.out"),
        ("PP_AE", f"{pp_cmd} < {tag}_pp21.in > {tag}_pp21.out"),
    ]
    for label, cmd in steps:
        print(f"  [{label}] {cmd}")
        subprocess.run(cmd, shell=True, cwd=str(run_dir), check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate/run VN (+U) LOBSTER preprocessing")
    parser.add_argument("--slab", choices=["VN"], default="VN", help="Only VN is supported in this +U workflow")
    parser.add_argument("--ads", nargs="+", default=["Li2S4", "S8", "Li2S8"], help="Adsorbates to include")
    parser.add_argument("--degauss", type=float, default=0.015)
    parser.add_argument("--k-scf", default="4 4 1 0 0 0")
    parser.add_argument("--np", type=int, default=64)
    parser.add_argument("--run", action="store_true")
    args = parser.parse_args()

    ensure_directory_exists(XYZ_DIR, "xyz directory")
    ensure_directory_exists(RUN_ROOT, "run root")
    ensure_directory_exists(Path(PSEUDO_DIR), "pseudo directory")

    matched = 0
    for xyz in sorted(XYZ_DIR.glob("*_final.xyz")):
        tag = xyz.name.removesuffix("_final.xyz")

        if tag.endswith("_slab"):
            continue

        slab, ads = tag.split("_", 1)
        if slab != TARGET_SLAB or args.slab != TARGET_SLAB:
            continue
        if ads not in args.ads:
            continue

        atoms = read(str(xyz))
        run_dir = RUN_ROOT / f"{tag}_combi"
        run_dir.mkdir(parents=True, exist_ok=True)
        (run_dir / "tmp").mkdir(exist_ok=True)

        write_qe_input(run_dir / f"{tag}_scf.in", atoms, tag, args.degauss, args.k_scf, "scf")
        write_qe_input(run_dir / f"{tag}_nscf.in", atoms, tag, args.degauss, args.k_scf, "nscf")
        write_pp_ae_density(run_dir / f"{tag}_pp21.in", tag)

        matched += 1
        print(f"Generated 2-step VN (+U) inputs for {tag}")
        if args.run:
            run_chain(run_dir, tag, args.np)

    if matched == 0:
        print("No matching VN combined xyz files found.")


if __name__ == "__main__":
    main()
