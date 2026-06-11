#!/usr/bin/env python3
"""
Generate and optionally run QE SCF/NSCF for LOBSTER workflow with DFT+U (VN only).

Modified for Nibi/Non-Lustre: Uses local ./tmp for I/O.
"""

from __future__ import annotations

import argparse
import os
import subprocess
from pathlib import Path

from ase.io import read

TMP_DIR = "./tmp"

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

HUBBARD_U_V_3D = 2.50
ECUTWFC_DEFAULT = 60.0
ECUTRHO_DEFAULT = 480.0

TARGET_SLAB = "VN"

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

def parse_kgrid(text: str) -> list[int]:
    vals = [int(x) for x in text.split()]
    if len(vals) != 6:
        raise ValueError("K grid must contain 6 integers")
    return vals

def ensure_directory_exists(path: Path, label: str) -> None:
    if not path.exists():
        path.mkdir(parents=True, exist_ok=True)
        print(f"Created {label}: {path}")

def write_scf_input(
    file_path: Path, atoms, tag: str, pseudo_dir: str, ecutwfc: float, ecutrho: float, degauss: float, k_scf: list[int],
) -> None:
    nat = len(atoms)
    species = ordered_species(tag, atoms)
    ntyp = len(species)

    with open(file_path, "w", encoding="utf-8") as f:
        f.write("&CONTROL\n")
        f.write("  calculation = 'scf'\n")
        f.write(f"  prefix = '{tag}'\n")
        f.write("  restart_mode = 'restart'\n")
        f.write("  max_seconds = 148000\n")
        f.write(f"  pseudo_dir = '{pseudo_dir}'\n")
        f.write(f"  outdir = '{TMP_DIR}'\n")
        f.write("  disk_io = 'medium'\n")
        f.write("  verbosity = 'low'\n")
        f.write("  wf_collect = .true.\n")
        f.write("/\n")

        f.write("&SYSTEM\n")
        f.write(f"  ibrav = 0, nat = {nat}, ntyp = {ntyp}\n")
        f.write(f"  ecutwfc = {ecutwfc:.1f}, ecutrho = {ecutrho:.1f}\n")
        f.write(f"  occupations = 'smearing', smearing = 'cold', degauss = {degauss:.5f}\n")
        f.write("  nspin = 2\n")
        f.write("  starting_magnetization(1) = 1.5\n")
        f.write("  starting_magnetization(2) = -0.15\n")
        f.write("/\n")

        f.write("&ELECTRONS\n")
        f.write("  conv_thr = 1.0d-6\n")
        f.write("  mixing_beta = 0.20\n")
        f.write("  mixing_mode = 'local-TF'\n")
        f.write("  electron_maxstep = 250\n")
        f.write("  diago_david_ndim = 12\n")
        f.write("/\n")

        f.write("\nCELL_PARAMETERS (angstrom)\n")
        for vec in atoms.get_cell():
            f.write(f"  {vec[0]:15.10f} {vec[1]:15.10f} {vec[2]:15.10f}\n")

        f.write("\nATOMIC_SPECIES\n")
        for sym in species:
            f.write(f"  {sym:<3} 1.0  {PSEUDOS[sym]}\n")

        f.write("\nATOMIC_POSITIONS (angstrom)\n")
        for atom in atoms:
            f.write(f"  {atom.symbol:<3} {atom.position[0]:14.9f} {atom.position[1]:14.9f} {atom.position[2]:14.9f}\n")

        f.write("\nK_POINTS (automatic)\n")
        f.write(f"  {k_scf[0]} {k_scf[1]} {k_scf[2]} {k_scf[3]} {k_scf[4]} {k_scf[5]}\n")

        f.write("\nHUBBARD {atomic}\n")
        f.write(f"U V-3d {HUBBARD_U_V_3D}\n")

def write_nscf_input(
    file_path: Path, atoms, tag: str, pseudo_dir: str, ecutwfc: float, ecutrho: float, k_nscf: list[int],
) -> None:
    nat = len(atoms)
    species = ordered_species(tag, atoms)
    ntyp = len(species)
    nbnd = calculate_min_nbnd(atoms)

    with open(file_path, "w", encoding="utf-8") as f:
        f.write("&CONTROL\n")
        f.write("  calculation = 'nscf'\n")
        f.write(f"  prefix = '{tag}'\n")
        f.write(f"  pseudo_dir = '{pseudo_dir}'\n")
        f.write(f"  outdir = '{TMP_DIR}'\n")
        f.write("  disk_io = 'medium'\n")
        f.write("  verbosity = 'low'\n")
        f.write("  wf_collect = .true.\n")
        f.write("/\n")

        f.write("&SYSTEM\n")
        f.write(f"  ibrav = 0, nat = {nat}, ntyp = {ntyp}\n")
        f.write(f"  ecutwfc = {ecutwfc:.1f}, ecutrho = {ecutrho:.1f}\n")
        f.write(f"  nbnd = {nbnd}\n")
        f.write("  nosym = .true.\n")
        f.write("  noinv = .true.\n")
        f.write("  occupations = 'smearing'\n")
        f.write("  smearing = 'cold', degauss = 0.02\n")
        f.write("  nspin = 2\n")
        f.write("  starting_magnetization(1) = 1.5\n")
        f.write("  starting_magnetization(2) = -0.15\n")
        f.write("/\n")

        f.write("&ELECTRONS\n")
        f.write("  conv_thr = 1.0d-6\n")
        f.write("  diago_david_ndim = 6\n")
        f.write("  diago_thr_init = 1.0d-4\n")
        f.write("  electron_maxstep = 300\n")
        f.write("/\n")

        f.write("\nCELL_PARAMETERS (angstrom)\n")
        for vec in atoms.get_cell():
            f.write(f"  {vec[0]:15.10f} {vec[1]:15.10f} {vec[2]:15.10f}\n")

        f.write("\nATOMIC_SPECIES\n")
        for sym in species:
            f.write(f"  {sym:<3} 1.0  {PSEUDOS[sym]}\n")

        f.write("\nATOMIC_POSITIONS (angstrom)\n")
        for atom in atoms:
            f.write(f"  {atom.symbol:<3} {atom.position[0]:14.9f} {atom.position[1]:14.9f} {atom.position[2]:14.9f}\n")

        f.write("\nK_POINTS (automatic)\n")
        f.write(f"  {k_nscf[0]} {k_nscf[1]} {k_nscf[2]} {k_nscf[3]} {k_nscf[4]} {k_nscf[5]}\n")

        f.write("\nHUBBARD {atomic}\n")
        f.write(f"U V-3d {HUBBARD_U_V_3D}\n")

def write_pp_ae_density(file_path: Path, tag: str) -> None:
    with open(file_path, "w", encoding="utf-8") as f:
        f.write("&INPUTPP\n")
        f.write(f"  prefix = '{tag}'\n")
        f.write(f"  outdir = '{TMP_DIR}'\n")
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
    parser = argparse.ArgumentParser(description="VN (+U) LOBSTER workflow")
    
    parser.add_argument("--ads", nargs="+", default=["Li2S4", "S8", "Li2S8"])
    
    script_path = Path(__file__).resolve().parent
    parser.add_argument("--xyz-dir", type=Path, default=script_path / "combi" / "final_xyz")
    parser.add_argument("--run-root", type=Path, default=script_path / "lobster")
    parser.add_argument("--pseudo-dir", default=str(script_path.parent / "PAW_pslib"))
    
    parser.add_argument("--ecutwfc", type=float, default=ECUTWFC_DEFAULT)
    parser.add_argument("--ecutrho", type=float, default=ECUTRHO_DEFAULT)
    parser.add_argument("--degauss", type=float, default=0.02)
    parser.add_argument("--k-scf", default="4 4 1 0 0 0")
    parser.add_argument("--k-nscf", default="8 8 1 0 0 0")
    parser.add_argument("--np", type=int, default=64)
    parser.add_argument("--run", action="store_true")
    args = parser.parse_args()

    ensure_directory_exists(args.xyz_dir, "xyz directory")
    ensure_directory_exists(args.run_root, "run root")
    ensure_directory_exists(Path(args.pseudo_dir), "pseudo directory")

    k_scf = parse_kgrid(args.k_scf)
    k_nscf = parse_kgrid(args.k_nscf)

    matched = 0
    for xyz in sorted(args.xyz_dir.glob("*_final.xyz")):
        stem = xyz.stem
        if not stem.endswith("_final"):
            continue
        tag = stem.removesuffix("_final")
        
        if tag.endswith("_slab"):
            continue
        slab, ads = tag.split("_", 1)
        if slab != TARGET_SLAB or ads not in args.ads:
            continue

        atoms = read(str(xyz))
        run_dir = args.run_root / f"{tag}_combi"
        run_dir.mkdir(parents=True, exist_ok=True)
        (run_dir / "tmp").mkdir(exist_ok=True)

        write_scf_input(
            run_dir / f"{tag}_scf.in",
            atoms,
            tag,
            args.pseudo_dir,
            args.ecutwfc,
            args.ecutrho,
            args.degauss,
            k_scf,
        )
        write_nscf_input(
            run_dir / f"{tag}_nscf.in",
            atoms,
            tag,
            args.pseudo_dir,
            args.ecutwfc,
            args.ecutrho,
            k_nscf,
        )
        write_pp_ae_density(run_dir / f"{tag}_pp21.in", tag)

        matched += 1
        if args.run:
            run_chain(run_dir, tag, args.np)

    if matched == 0:
        print("No matching VN combined xyz files found.")

if __name__ == "__main__":
    main()
