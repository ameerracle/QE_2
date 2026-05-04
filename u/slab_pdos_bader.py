#!/usr/bin/env python3
"""
Generate and optionally run QE PDOS/Bader workflow for relaxed VN/TiN slabs.
Updated to read relaxed structures from .xyz files.
"""

from __future__ import annotations

import argparse
import os
import subprocess
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import read

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

SLABS = ["VN", "TiN"]
SCF_CONV_THR = "1.0d-7"
NSCF_CONV_THR = "1.0d-7"
# Hubbard U values
HUBBARD_U = {
    "V": 2.50,
    "Ti": 2.50 # Added default for Ti if needed
}

def parse_kgrid(text: str) -> list[int]:
    vals = [int(x) for x in text.split()]
    if len(vals) != 6:
        raise ValueError("K grid must contain 6 integers: kx ky kz sx sy sz")
    return vals

def ensure_directory_exists(path: Path, label: str) -> None:
    if not path.exists():
        path.mkdir(parents=True, exist_ok=True)
        print(f"Created missing {label}: {path}")

def write_scf_input(file_path: Path, atoms: Atoms, slab_name: str, pseudo_dir: str, ecutwfc: float, ecutrho: float, degauss: float, k_scf: list[int]) -> None:
    metal = "".join([c for c in slab_name if not c.isdigit() and c != 'N'])
    u_val = HUBBARD_U.get(metal, 2.50)
    
    with open(file_path, "w", encoding="utf-8") as f:
        f.write("&CONTROL\n")
        f.write("  calculation = 'scf'\n")
        f.write(f"  prefix = '{slab_name}_slab'\n")
        f.write(f"  pseudo_dir = '{pseudo_dir}'\n")
        f.write("  outdir = './tmp/'\n")
        f.write("  disk_io = 'medium'\n")
        f.write("  verbosity = 'high'\n")
        f.write("/\n")

        f.write("&SYSTEM\n")
        f.write(f" ibrav = 0, nat = {len(atoms)}, ntyp = 2\n")
        f.write(f" ecutwfc = {ecutwfc:.1f}, ecutrho = {ecutrho:.1f}\n")
        f.write(f" occupations = 'smearing', smearing = 'cold', degauss = {degauss:.5f}\n")
        f.write("  nspin = 2\n")
        f.write("  starting_magnetization(1) = 1.0\n")
        f.write("/\n")

        f.write("&ELECTRONS\n")
        f.write(f"  conv_thr = {SCF_CONV_THR}\n")
        f.write("  mixing_beta = 0.2\n")
        f.write(f" mixing_mode = 'local-TF'\n") # Essential for slab-vacuum interfaces
        f.write(f" mixing_fixed_ns = 10\n")
        f.write("  electron_maxstep = 250\n")
        f.write("  diago_david_ndim = 4\n")
        f.write("/\n")

        f.write("\nCELL_PARAMETERS (angstrom)\n")
        for vec in atoms.get_cell():
            f.write(f"  {vec[0]:15.10f} {vec[1]:15.10f} {vec[2]:15.10f}\n")

        f.write("\nATOMIC_SPECIES\n")
        f.write(f"  {metal:<3} 1.0  {PSEUDOS[metal]}\n")
        f.write(f"  N   1.0  {PSEUDOS['N']}\n")

        f.write("\nATOMIC_POSITIONS (angstrom)\n")
        for atom in atoms:
            f.write(f"  {atom.symbol:<3} {atom.position[0]:14.9f} {atom.position[1]:14.9f} {atom.position[2]:14.9f}\n")

        f.write("\nK_POINTS (automatic)\n")
        f.write(f"  {k_scf[0]} {k_scf[1]} {k_scf[2]} {k_scf[3]} {k_scf[4]} {k_scf[5]}\n")

        f.write("\nHUBBARD {atomic}\n")
        f.write(f"U {metal}-3d {u_val}\n")

def write_nscf_input(file_path: Path, atoms: Atoms, slab_name: str, pseudo_dir: str, ecutwfc: float, ecutrho: float, k_nscf: list[int]) -> None:
    metal = "".join([c for c in slab_name if not c.isdigit() and c != 'N'])
    u_val = HUBBARD_U.get(metal, 2.50)

    with open(file_path, "w", encoding="utf-8") as f:
        f.write("&CONTROL\n")
        f.write("  calculation = 'nscf'\n")
        f.write(f"  prefix = '{slab_name}_slab'\n")
        f.write(f"  pseudo_dir = '{pseudo_dir}'\n")
        f.write("  outdir = './tmp/'\n")
        f.write("  disk_io = 'low'\n")
        f.write("  verbosity = 'high'\n")
        f.write("/\n")

        f.write("&SYSTEM\n")
        f.write(f"  ibrav = 0, nat = {len(atoms)}, ntyp = 2\n")
        f.write(f"  ecutwfc = {ecutwfc:.1f}, ecutrho = {ecutrho:.1f}\n")
        f.write("  occupations = 'tetrahedra_opt'\n")
        f.write("  nspin = 2\n")
        f.write("  starting_magnetization(1) = 1.0\n")
        f.write("/\n")

        f.write("&ELECTRONS\n")
        f.write(f"  conv_thr = {NSCF_CONV_THR}\n")
        f.write("  diago_david_ndim = 6\n")
        f.write("  diago_thr_init = 1.0d-4\n")
        f.write("  electron_maxstep = 300\n")
        f.write("/\n")

        f.write("\nCELL_PARAMETERS (angstrom)\n")
        for vec in atoms.get_cell():
            f.write(f"  {vec[0]:15.10f} {vec[1]:15.10f} {vec[2]:15.10f}\n")

        f.write("\nATOMIC_SPECIES\n")
        f.write(f"  {metal:<3} 1.0  {PSEUDOS[metal]}\n")
        f.write(f"  N   1.0  {PSEUDOS['N']}\n")

        f.write("\nATOMIC_POSITIONS (angstrom)\n")
        for atom in atoms:
            f.write(f"  {atom.symbol:<3} {atom.position[0]:14.9f} {atom.position[1]:14.9f} {atom.position[2]:14.9f}\n")

        f.write("\nK_POINTS (automatic)\n")
        f.write(f"  {k_nscf[0]} {k_nscf[1]} {k_nscf[2]} 0 0 0\n")

        f.write("\nHUBBARD {atomic}\n")
        f.write(f"U {metal}-3d {u_val}\n")

def write_projwfc_input(file_path: Path, slab_name: str, filpdos: str, delta_e: float, e_min: float, e_max: float) -> None:
    with open(file_path, "w", encoding="utf-8") as f:
        f.write("&PROJWFC\n")
        f.write(f"  prefix = '{slab_name}_slab'\n")
        f.write("  outdir = './tmp/'\n")
        f.write("  ngauss = -1\n")
        f.write(f"  DeltaE = {delta_e:.4f}\n")
        f.write(f"  Emin = {e_min:.2f}\n")
        f.write(f"  Emax = {e_max:.2f}\n")
        f.write(f"  filpdos = '{filpdos}'\n")
        f.write("/\n")

def write_pp_input(file_path: Path, slab_name: str, filplot: str, fileout: str) -> None:
    with open(file_path, "w", encoding="utf-8") as f:
        f.write("&INPUTPP\n")
        f.write(f"  prefix = '{slab_name}_slab'\n")
        f.write("  outdir = './tmp/'\n")
        f.write("  plot_num = 21\n")
        f.write(f"  filplot = '{filplot}'\n")
        f.write("/\n")
        f.write("&PLOT\n")
        f.write("  iflag = 3\n")
        f.write("  output_format = 6\n")
        f.write(f"  fileout = '{fileout}'\n")
        f.write("/\n")

def run_chain(run_dir: Path, slab_name: str, np: int, nk: int) -> None:
    pref = f"{slab_name}_slab"
    cmds = [
        ("SCF", f"pw.x -nk {nk} < {pref}_scf.in > {pref}_scf.out"),
        ("PP", f"pp.x < {pref}_pp.in > {pref}_pp.out"),
        ("NSCF", f"pw.x -nk {nk} < {pref}_nscf.in > {pref}_nscf.out"),
        ("PROJWFC", f"projwfc.x < {pref}_projwfc.in > {pref}_projwfc.out"),
    ]

    for label, cmd_base in cmds:
        cmd = f"srun {cmd_base}" if "SLURM_JOB_ID" in os.environ else f"mpirun -np {np} {cmd_base}"
        print(f"  [{label}] {cmd}")
        subprocess.run(cmd, shell=True, cwd=str(run_dir), check=True)

def prepare_single_slab(
    slab_name: str,
    input_root: Path,
    run_root: Path,
    pseudo_dir: str,
    **kwargs,
) -> None:
    # UPDATED: Search for the .xyz file instead of .out
    xyz_path = input_root / f"{slab_name}_slab_final.xyz"
    
    if not xyz_path.exists():
        print(f"[SKIP] Missing XYZ structure: {xyz_path}")
        return

    # Use ASE to read the XYZ file
    # Note: Ensure your XYZ files have the unit cell information (Extended XYZ)
    try:
        atoms = read(str(xyz_path))
    except Exception as e:
        print(f"[ERROR] Could not read {xyz_path}: {e}")
        return

    run_dir = run_root / f"{slab_name}_slab"
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "tmp").mkdir(exist_ok=True)

    pref = f"{slab_name}_slab"
    write_scf_input(run_dir / f"{pref}_scf.in", atoms, slab_name, pseudo_dir, kwargs["ecutwfc"], kwargs["ecutrho"], kwargs["degauss"], kwargs["k_scf"])
    write_nscf_input(run_dir / f"{pref}_nscf.in", atoms, slab_name, pseudo_dir, kwargs["ecutwfc"], kwargs["ecutrho"], kwargs["k_nscf"])
    write_projwfc_input(run_dir / f"{pref}_projwfc.in", slab_name, kwargs["filpdos"], kwargs["deltae"], kwargs["emin"], kwargs["emax"])
    write_pp_input(run_dir / f"{pref}_pp.in", slab_name, f"{pref}_rho", f"{pref}_charge.cube")

    print(f"Generated {slab_name} slab (+U) inputs in {run_dir}")
    if kwargs["run"]:
        run_chain(run_dir, slab_name, kwargs["np"], kwargs["nk"])

def main() -> None:
    parser = argparse.ArgumentParser(description="Generate/run VN/TiN slab (+U) PDOS workflow from XYZ")
    parser.add_argument("--structure", choices=SLABS, default="VN")
    # UPDATED: New default path for your XYZ files
    parser.add_argument("--input-root", type=Path, default=Path("/lustre10/scratch/anizami/QE_2/u/combi/final_xyz"))
    parser.add_argument("--run-root", type=Path, default=Path("/lustre10/scratch/anizami/QE_2/u/PDOS"))
    parser.add_argument("--pseudo-dir", default="/lustre10/scratch/anizami/QE_2/PAW_pslib")
    parser.add_argument("--ecutwfc", type=float, default=60.0)
    parser.add_argument("--ecutrho", type=float, default=480.0)
    parser.add_argument("--degauss", type=float, default=0.015)
    parser.add_argument("--k-scf", default="4 4 1 0 0 0")
    parser.add_argument("--k-nscf", default="8 8 1 0 0 0")
    parser.add_argument("--filpdos", default="slab_pdos")
    parser.add_argument("--deltae", type=float, default=0.05)
    parser.add_argument("--emin", type=float, default=-12.0)
    parser.add_argument("--emax", type=float, default=8.0)
    parser.add_argument("--np", type=int, default=64)
    parser.add_argument("--nk", type=int, default=2)
    parser.add_argument("--run", action="store_true")
    args = parser.parse_args()

    ensure_directory_exists(args.input_root, "input root")
    ensure_directory_exists(args.run_root, "run root")
    ensure_directory_exists(Path(args.pseudo_dir), "pseudo directory")

    k_scf = parse_kgrid(args.k_scf)
    k_nscf = parse_kgrid(args.k_nscf)

    prepare_single_slab(
        args.structure,
        args.input_root,
        args.run_root,
        args.pseudo_dir,
        ecutwfc=args.ecutwfc,
        ecutrho=args.ecutrho,
        degauss=args.degauss,
        k_scf=k_scf,
        k_nscf=k_nscf,
        filpdos=args.filpdos,
        deltae=args.deltae,
        emin=args.emin,
        emax=args.emax,
        run=args.run,
        np=args.np,
        nk=args.nk,
    )

if __name__ == "__main__":
    main()