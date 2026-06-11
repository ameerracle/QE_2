#!/usr/bin/env python3
"""
Generate and optionally run QE PDOS workflow for adsorbates only.
Uses final structures from adsorbates/final_xyz/
Output directory: /lustre10/scratch/anizami/QE_2/adsorbates/PDOS
"""

from __future__ import annotations
import argparse
import os
import shutil
import subprocess
from pathlib import Path
from ase.io import read

PSEUDOS = {
    "Li": "Li.pbe-s-kjpaw_psl.1.0.0.UPF",
    "S": "S.pbe-n-kjpaw_psl.1.0.0.UPF",
    "C": "C.pbe-n-kjpaw_psl.1.0.0.UPF"
}

ADSORBATES = ["Li2S", "Li2S2", "Li2S4", "Li2S6", "Li2S8", "S8"]
FULL_PDOS_ADS = {"Li2S4", "Li2S8", "S8"}
FIXED_NK = 1
SCF_CONV_THR = 1.0e-6
NSCF_CONV_THR = 1.0e-6

def parse_kgrid(text: str) -> list[int]:
    return [int(x) for x in text.split()]

def ordered_species(atoms) -> list[str]:
    present = set(atoms.get_chemical_symbols())
    species = []
    for sym in ["Li", "S", "C"]:
        if sym in present:
            species.append(sym)
    for sym in atoms.get_chemical_symbols():
        if sym not in species:
            species.append(sym)
    return species

def format_conv_thr(value: float) -> str:
    return f"{value:.1e}".replace("e", "d")

def write_scf_input(file_path: Path, atoms, tag: str, pseudo_dir: str, ecutwfc: float, ecutrho: float, degauss: float, k_scf: list[int], nspin: int = 1):
    nat = len(atoms)
    species = ordered_species(atoms)
    ntyp = len(species)

    with open(file_path, "w", encoding="utf-8") as f:
        f.write("&CONTROL\n")
        f.write("  calculation = 'scf'\n")
        f.write(f"  prefix = '{tag}'\n")
        f.write(f"  pseudo_dir = '{pseudo_dir}'\n")
        f.write("  outdir = './tmp/'\n")
        f.write("  disk_io = 'medium'\n")
        f.write("  verbosity = 'high'\n")
        f.write("/\n")
        f.write("&SYSTEM\n")
        f.write(f"  ibrav = 0, nat = {nat}, ntyp = {ntyp}\n")
        f.write(f"  ecutwfc = {ecutwfc:.1f}, ecutrho = {ecutrho:.1f}\n")
        f.write(f"  occupations = 'smearing', smearing = 'cold', degauss = {degauss:.5f}\n")
        f.write(f"  nspin = {nspin}\n")
        if nspin == 2:
            for i, sym in enumerate(species, 1):
                f.write(f"  starting_magnetization({i}) = 0.01\n")
        f.write("/\n")
        f.write("&ELECTRONS\n")
        f.write(f"  conv_thr = {format_conv_thr(SCF_CONV_THR)}\n")
        if nspin == 2:
            f.write("  mixing_beta = 0.15\n")
            f.write("  mixing_mode = 'local-TF'\n")
        else:
            f.write("  mixing_beta = 0.25\n")
        f.write("  electron_maxstep = 250\n")
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

def write_nscf_input(file_path, atoms, tag, pseudo_dir, ecutwfc, ecutrho, k_nscf, nspin=1):
    nat = len(atoms)
    species = ordered_species(atoms)
    ntyp = len(species)

    with open(file_path, "w", encoding="utf-8") as f:
        f.write("&CONTROL\n")
        f.write("  calculation = 'nscf'\n")
        f.write(f"  prefix = '{tag}'\n")
        f.write(f"  pseudo_dir = '{pseudo_dir}'\n")
        f.write("  outdir = './tmp/', disk_io = 'low', verbosity = 'high'\n/\n")
        f.write("&SYSTEM\n")
        f.write(f"  ibrav = 0, nat = {nat}, ntyp = {ntyp}\n")
        f.write(f"  ecutwfc = {ecutwfc:.1f}, ecutrho = {ecutrho:.1f}\n")
        f.write("  occupations = 'smearing'\n")
        f.write("  smearing = 'cold', degauss = 0.015\n")
        f.write(f"  nspin = {nspin}\n")
        if nspin == 2:
            for i, sym in enumerate(species, 1):
                f.write(f"  starting_magnetization({i}) = 0.01\n")
        f.write("/\n")
        f.write("&ELECTRONS\n")
        f.write(f"  conv_thr = {format_conv_thr(NSCF_CONV_THR)}\n")
        f.write("  diago_david_ndim = 6\n")
        f.write("  diago_thr_init = 1.0d-4\n")
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

def write_projwfc_input(file_path, tag, filpdos, delta_e, e_min, e_max):
    with open(file_path, "w", encoding="utf-8") as f:
        f.write("&PROJWFC\n")
        f.write(f"  prefix = '{tag}'\n")
        f.write("  outdir = './tmp/'\n")
        f.write("  ngauss = -1\n")
        f.write(f"  DeltaE = {delta_e:.4f}\n")
        f.write(f"  Emin = {e_min:.2f}\n")
        f.write(f"  Emax = {e_max:.2f}\n")
        f.write(f"  filpdos = '{filpdos}'\n")
        f.write("/\n")

def write_pp_input(file_path, tag, filplot, fileout, plot_num=21):
    with open(file_path, "w", encoding="utf-8") as f:
        f.write("&INPUTPP\n")
        f.write(f"  prefix = '{tag}'\n")
        f.write("  outdir = './tmp/'\n")
        f.write(f"  plot_num = {plot_num}\n")
        f.write(f"  filplot = '{filplot}'\n")
        f.write("/\n")
        f.write("&PLOT\n")
        f.write("  iflag = 3\n")
        f.write("  output_format = 6\n")
        f.write(f"  fileout = '{fileout}'\n")
        f.write("/\n")

def run_chain(run_dir, tag, np, step, do_pdos, nspin=1):
    steps = []
    
    if step in ["scf", "full"]:
        steps.append(("SCF", f"srun pw.x -nk {FIXED_NK} < {tag}_scf.in > {tag}_scf.out" if "SLURM_JOB_ID" in os.environ else f"mpirun -np {np} pw.x -nk {FIXED_NK} < {tag}_scf.in > {tag}_scf.out"))
        if nspin == 2:
            steps.append(("PP_plot0", f"srun pp.x < {tag}_pp_plot0.in > {tag}_pp_plot0.out" if "SLURM_JOB_ID" in os.environ else f"mpirun -np {np} pp.x < {tag}_pp_plot0.in > {tag}_pp_plot0.out"))
        steps.append(("PP", f"srun pp.x < {tag}_pp.in > {tag}_pp.out" if "SLURM_JOB_ID" in os.environ else f"mpirun -np {np} pp.x < {tag}_pp.in > {tag}_pp.out"))
    
    if do_pdos and step in ["nscf", "full"]:
        steps.append(("NSCF", f"srun pw.x -nk {FIXED_NK} < {tag}_nscf.in > {tag}_nscf.out" if "SLURM_JOB_ID" in os.environ else f"mpirun -np {np} pw.x -nk {FIXED_NK} < {tag}_nscf.in > {tag}_nscf.out"))
        steps.append(("PROJWFC", f"srun projwfc.x < {tag}_projwfc.in > {tag}_projwfc.out" if "SLURM_JOB_ID" in os.environ else f"mpirun -np {np} projwfc.x < {tag}_projwfc.in > {tag}_projwfc.out"))

    if not steps:
        print(f"  No steps to run for {tag} (step={step}, do_pdos={do_pdos})")
        return

    for label, cmd in steps:
        print(f"  [{label}] {cmd}")
        subprocess.run(cmd, shell=True, cwd=str(run_dir), check=True)

def main():
    parser = argparse.ArgumentParser(description="Generate/run PAW PDOS workflow for adsorbates")
    parser.add_argument("--ads", nargs="+", default=ADSORBATES, help="Adsorbates to process")
    parser.add_argument("--xyz-dir", type=Path, default=Path("/lustre10/scratch/anizami/QE_2/adsorbates/final_xyz"))
    parser.add_argument("--run-root", type=Path, default=Path("/lustre10/scratch/anizami/QE_2/adsorbates/PDOS"))
    parser.add_argument("--local-tf", action="store_true", help="Use local test folder (for debugging)")
    parser.add_argument("--pseudo-dir", default="/lustre10/scratch/anizami/QE_2/PAW_pslib")
    parser.add_argument("--ecutwfc", type=float, default=60.0)
    parser.add_argument("--ecutrho", type=float, default=480.0)
    parser.add_argument("--degauss", type=float, default=0.015)
    parser.add_argument("--k-scf", default="1 1 1 0 0 0")
    parser.add_argument("--k-nscf", default="1 1 1 0 0 0")
    parser.add_argument("--filpdos", default="ads_pdos")
    parser.add_argument("--deltae", type=float, default=0.05)
    parser.add_argument("--emin", type=float, default=-12.0)
    parser.add_argument("--emax", type=float, default=8.0)
    parser.add_argument("--np", type=int, default=64)
    parser.add_argument("--nspin", type=int, choices=[1, 2], default=1, help="Spin polarization (1=nspin1, 2=nspin2)")
    parser.add_argument("--mode", choices=["auto", "full", "bader-only"], default="auto",
                        help="auto: PDOS only for Li2S4/Li2S8/S8; full: PDOS for all; bader-only: no PDOS")
    parser.add_argument("--step", choices=["scf", "nscf", "full"], default="full",
                        help="Which step to run: scf (SCF+PP only), nscf (NSCF+PROJWFC only, assumes SCF done), full (both)")
    parser.add_argument("--run", action="store_true")
    args = parser.parse_args()

    k_scf, k_nscf = parse_kgrid(args.k_scf), parse_kgrid(args.k_nscf)
    
    if args.local_tf:
        args.xyz_dir = Path("/tmp/qe_test/adsorbates/final_xyz")
        args.run_root = Path("/tmp/qe_test/adsorbates/PDOS")
        args.pseudo_dir = "/tmp/qe_test/pseudos"
    
    if args.nspin == 2:
        args.run_root = args.run_root / "nspin2"

    tags = []
    for xyz in sorted(args.xyz_dir.glob("*.xyz")):
        tag = xyz.stem
        if tag not in args.ads:
            continue
        tags.append(tag)

    if not tags:
        return print("No matching xyz files found.")

    for tag in tags:
        xyz_file = args.xyz_dir / f"{tag}.xyz"
        atoms = read(str(xyz_file))
        run_dir = args.run_root / tag
        run_dir.mkdir(parents=True, exist_ok=True)

        tmp_dir = run_dir / "tmp"
        if tmp_dir.exists():
            shutil.rmtree(tmp_dir)
        for save_dir in run_dir.glob("*.save"):
            if save_dir.is_dir():
                shutil.rmtree(save_dir)
        tmp_dir.mkdir(exist_ok=True)

        write_scf_input(run_dir / f"{tag}_scf.in", atoms, tag, args.pseudo_dir, args.ecutwfc, args.ecutrho, args.degauss, k_scf, nspin=args.nspin)
        
        ads_name = tag
        if args.mode == "full":
            do_pdos = True
        elif args.mode == "bader-only":
            do_pdos = False
        else:
            do_pdos = ads_name in FULL_PDOS_ADS
        
        if args.step in ["scf", "full"]:
            if args.nspin == 2:
                write_pp_input(run_dir / f"{tag}_pp_plot0.in", tag, f"{tag}_rho_plot0", f"{tag}_charge_plot0.cube", plot_num=0)
            write_pp_input(run_dir / f"{tag}_pp.in", tag, f"{tag}_rho", f"{tag}_charge.cube", plot_num=21)
        
        if do_pdos and args.step in ["nscf", "full"]:
            write_nscf_input(run_dir / f"{tag}_nscf.in", atoms, tag, args.pseudo_dir, args.ecutwfc, args.ecutrho, k_nscf, nspin=args.nspin)
            write_projwfc_input(run_dir / f"{tag}_projwfc.in", tag, args.filpdos, args.deltae, args.emin, args.emax)
        
        print(f"Generated inputs for {tag} (mode={args.mode}, step={args.step}, do_pdos={do_pdos}, nspin={args.nspin}) in {run_dir}")
        if args.run:
            run_chain(run_dir, tag, args.np, args.step, do_pdos, nspin=args.nspin)

if __name__ == "__main__":
    main()
