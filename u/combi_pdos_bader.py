#!/usr/bin/env python3
"""
Generate and optionally run QE SCF/NSCF/PDOS/Bader workflow for combined slab+adsorbate
structures with DFT+U enabled (VN only).

Uses final structures from:
  /lustre10/scratch/anizami/QE_2/u/combi/final_xyz

Writes run folders under:
  /lustre10/scratch/anizami/QE_2/u/PDOS
"""

from __future__ import annotations

import argparse
import os
import subprocess
from collections import defaultdict
from pathlib import Path

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

TARGET_SLAB = "VN"
FULL_PDOS_ADS = {"Li2S4", "Li2S8", "S8"}
FIXED_NK = 2
SCF_CONV_THR = 1.0e-7
NSCF_CONV_THR = 1.0e-7
HUBBARD_U_V_3D = 2.50


def parse_kgrid(text: str) -> list[int]:
    vals = [int(x) for x in text.split()]
    if len(vals) != 6:
        raise ValueError("K grid must contain 6 integers: kx ky kz sx sy sz")
    return vals


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


def format_conv_thr(value: float) -> str:
    return f"{value:.1e}".replace("e", "d")


def ensure_directory_exists(path: Path, label: str) -> None:
    if not path.exists():
        path.mkdir(parents=True, exist_ok=True)
        print(f"Created missing {label}: {path}")


def write_scf_input(
    file_path: Path,
    atoms,
    tag: str,
    pseudo_dir: str,
    ecutwfc: float,
    ecutrho: float,
    degauss: float,
    k_scf: list[int],
) -> None:
    nat = len(atoms)
    species = ordered_species(tag, atoms)
    ntyp = len(species)

    metal = tag.split("_", 1)[0].rstrip("N")
    if metal != "V":
        raise ValueError(f"This script is VN-only; got metal={metal}")

    metal_index = species.index("V") + 1

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
        f.write("  nspin = 2\n")
        f.write(f"  starting_magnetization({metal_index}) = 1.0\n")
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
        for sym in species:
            if sym not in PSEUDOS:
                raise KeyError(f"Missing pseudopotential for '{sym}'")
            f.write(f"  {sym:<3} 1.0  {PSEUDOS[sym]}\n")

        f.write("\nATOMIC_POSITIONS (angstrom)\n")
        for atom in atoms:
            f.write(f"  {atom.symbol:<3} {atom.position[0]:14.9f} {atom.position[1]:14.9f} {atom.position[2]:14.9f}\n")

        f.write("\nK_POINTS (automatic)\n")
        f.write(f"  {k_scf[0]} {k_scf[1]} {k_scf[2]} {k_scf[3]} {k_scf[4]} {k_scf[5]}\n")

        f.write("\nHUBBARD {atomic}\n")
        f.write(f"U V-3d {HUBBARD_U_V_3D}\n")


def write_nscf_input(
    file_path: Path,
    atoms,
    tag: str,
    pseudo_dir: str,
    ecutwfc: float,
    ecutrho: float,
    k_nscf: list[int],
) -> None:
    nat = len(atoms)
    species = ordered_species(tag, atoms)
    ntyp = len(species)

    metal = tag.split("_", 1)[0].rstrip("N")
    if metal != "V":
        raise ValueError(f"This script is VN-only; got metal={metal}")

    metal_index = species.index("V") + 1

    with open(file_path, "w", encoding="utf-8") as f:
        f.write("&CONTROL\n")
        f.write("  calculation = 'nscf'\n")
        f.write(f"  prefix = '{tag}'\n")
        f.write(f"  pseudo_dir = '{pseudo_dir}'\n")
        f.write("  outdir = './tmp/'\n")
        f.write("  disk_io = 'low'\n")
        f.write("  verbosity = 'high'\n")
        f.write("/\n")

        f.write("&SYSTEM\n")
        f.write(f"  ibrav = 0, nat = {nat}, ntyp = {ntyp}\n")
        f.write(f"  ecutwfc = {ecutwfc:.1f}, ecutrho = {ecutrho:.1f}\n")
        f.write("  occupations = 'tetrahedra_opt'\n")
        f.write("  nspin = 2\n")
        f.write(f"  starting_magnetization({metal_index}) = 1.0\n")
        f.write("/\n")

        f.write("&ELECTRONS\n")
        f.write(f"  conv_thr = {format_conv_thr(NSCF_CONV_THR)}\n")
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
        f.write(f"  {k_nscf[0]} {k_nscf[1]} {k_nscf[2]} {k_nscf[3]} {k_nscf[4]} {k_nscf[5]}\n")

        f.write("\nHUBBARD {atomic}\n")
        f.write(f"U V-3d {HUBBARD_U_V_3D}\n")


def write_projwfc_input(file_path: Path, tag: str, filpdos: str, delta_e: float, e_min: float, e_max: float) -> None:
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


def write_pp_input(file_path: Path, tag: str, filplot: str, fileout: str) -> None:
    with open(file_path, "w", encoding="utf-8") as f:
        f.write("&INPUTPP\n")
        f.write(f"  prefix = '{tag}'\n")
        f.write("  outdir = './tmp/'\n")
        f.write("  plot_num = 21\n")
        f.write(f"  filplot = '{filplot}'\n")
        f.write("/\n")
        f.write("&PLOT\n")
        f.write("  iflag = 3\n")
        f.write("  output_format = 6\n")
        f.write(f"  fileout = '{fileout}'\n")
        f.write("/\n")


def run_chain(run_dir: Path, tag: str, np: int, step: str, do_pdos: bool) -> None:
    steps = []

    if step in ["scf", "full"]:
        pw_cmd = (
            f"srun pw.x -nk {FIXED_NK} < {tag}_scf.in > {tag}_scf.out"
            if "SLURM_JOB_ID" in os.environ
            else f"mpirun -np {np} pw.x -nk {FIXED_NK} < {tag}_scf.in > {tag}_scf.out"
        )
        pp_cmd = (
            f"srun pp.x < {tag}_pp.in > {tag}_pp.out"
            if "SLURM_JOB_ID" in os.environ
            else f"mpirun -np {np} pp.x < {tag}_pp.in > {tag}_pp.out"
        )
        steps.append(("SCF", pw_cmd))
        steps.append(("PP", pp_cmd))

    if do_pdos and step in ["nscf", "full"]:
        nscf_cmd = (
            f"srun pw.x -nk {FIXED_NK} < {tag}_nscf.in > {tag}_nscf.out"
            if "SLURM_JOB_ID" in os.environ
            else f"mpirun -np {np} pw.x -nk {FIXED_NK} < {tag}_nscf.in > {tag}_nscf.out"
        )
        proj_cmd = (
            f"srun projwfc.x < {tag}_projwfc.in > {tag}_projwfc.out"
            if "SLURM_JOB_ID" in os.environ
            else f"mpirun -np {np} projwfc.x < {tag}_projwfc.in > {tag}_projwfc.out"
        )
        steps.append(("NSCF", nscf_cmd))
        steps.append(("PROJWFC", proj_cmd))

    if not steps:
        print(f"  No steps to run for {tag} (step={step}, do_pdos={do_pdos})")
        return

    for label, cmd in steps:
        print(f"  [{label}] {cmd}")
        subprocess.run(cmd, shell=True, cwd=str(run_dir), check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate/run VN (+U) PAW PDOS workflow")
    parser.add_argument("--ads", nargs="+", default=["Li2S4", "S8", "Li2S8"])
    parser.add_argument("--xyz-dir", type=Path, default=Path("/lustre10/scratch/anizami/QE_2/u/combi/final_xyz"))
    parser.add_argument("--run-root", type=Path, default=Path("/lustre10/scratch/anizami/QE_2/u/PDOS"))
    parser.add_argument("--pseudo-dir", default="/lustre10/scratch/anizami/QE_2/PAW_pslib")
    parser.add_argument("--ecutwfc", type=float, default=60.0)
    parser.add_argument("--ecutrho", type=float, default=480.0)
    parser.add_argument("--degauss", type=float, default=0.015)
    parser.add_argument("--k-scf", default="4 4 1 0 0 0")
    parser.add_argument("--k-nscf", default="8 8 1 0 0 0")
    parser.add_argument("--filpdos", default="combi_pdos")
    parser.add_argument("--deltae", type=float, default=0.05)
    parser.add_argument("--emin", type=float, default=-12.0)
    parser.add_argument("--emax", type=float, default=8.0)
    parser.add_argument("--np", type=int, default=64)
    parser.add_argument(
        "--mode",
        choices=["auto", "full", "bader-only"],
        default="auto",
        help="auto: PDOS only for Li2S4/Li2S8/S8; full: PDOS for all; bader-only: no PDOS",
    )
    parser.add_argument(
        "--step",
        choices=["scf", "nscf", "full"],
        default="full",
        help="scf: SCF+PP, nscf: NSCF+PROJWFC, full: both",
    )
    parser.add_argument("--run", action="store_true")
    args = parser.parse_args()

    ensure_directory_exists(args.xyz_dir, "xyz directory")
    ensure_directory_exists(args.run_root, "run root")
    ensure_directory_exists(Path(args.pseudo_dir), "pseudo directory")

    k_scf = parse_kgrid(args.k_scf)
    k_nscf = parse_kgrid(args.k_nscf)

    groups: dict[str, list[Path]] = defaultdict(list)
    for xyz in sorted(args.xyz_dir.glob("*_final.xyz")):
        stem = xyz.stem
        if not stem.endswith("_final"):
            continue
        tag = stem.removesuffix("_final")

        if tag.endswith("_slab"):
            continue
        slab, ads = tag.split("_", 1)
        if slab != TARGET_SLAB:
            continue
        if ads not in args.ads:
            continue
        groups[tag].append(xyz)

    if not groups:
        print("No matching VN combined xyz files found.")
        return

    for tag in sorted(groups):
        xyz_file = groups[tag][0]
        atoms = read(str(xyz_file))

        run_dir = args.run_root / tag
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

        ads_name = tag.split("_", 1)[1]
        if args.mode == "full":
            do_pdos = True
        elif args.mode == "bader-only":
            do_pdos = False
        else:
            do_pdos = ads_name in FULL_PDOS_ADS

        if args.step in ["scf", "full"]:
            write_pp_input(run_dir / f"{tag}_pp.in", tag, f"{tag}_rho", f"{tag}_charge.cube")

        if do_pdos and args.step in ["nscf", "full"]:
            write_nscf_input(
                run_dir / f"{tag}_nscf.in",
                atoms,
                tag,
                args.pseudo_dir,
                args.ecutwfc,
                args.ecutrho,
                k_nscf,
            )
            write_projwfc_input(
                run_dir / f"{tag}_projwfc.in",
                tag,
                args.filpdos,
                args.deltae,
                args.emin,
                args.emax,
            )

        print(f"Generated inputs for {tag} (VN-only, +U, mode={args.mode}, step={args.step}, do_pdos={do_pdos})")
        if args.run:
            run_chain(run_dir, tag, args.np, args.step, do_pdos)


if __name__ == "__main__":
    main()
