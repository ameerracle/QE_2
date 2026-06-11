#!/usr/bin/env python3
"""
Generate QE PDOS workflow for adsorbate-only (Li/S) extracted from relaxed VN+U combi structures.
Outputs to /scratch/anizami/QE_2/u/PDOS/adsorbates_combi/
"""

import argparse
import os
import subprocess
from pathlib import Path
from ase.io import read
from ase import Atoms

PSEUDOS = {
    "Li": "Li.pbe-s-kjpaw_psl.1.0.0.UPF",
    "S": "S.pbe-n-kjpaw_psl.1.0.0.UPF"
}

SYSTEMS = ["VN_Li2S4", "VN_Li2S8", "VN_S8"]

SCF_CONV_THR = "1.0d-7"
NSCF_CONV_THR = "1.0d-7"


def parse_kgrid(text: str) -> list[int]:
    vals = [int(x) for x in text.split()]
    if len(vals) != 6:
        raise ValueError("K grid must contain 6 integers: kx ky kz sx sy sz")
    return vals


def read_combi_structure(xyz_path: Path):
    return read(str(xyz_path), format="xyz", index=-1)


VACUUM_BOX = 10.0

def filter_adsorbate_atoms(atoms):
    adsorbate_indices = [i for i, atom in enumerate(atoms) if atom.symbol in ["Li", "S"]]
    if not adsorbate_indices:
        return None
    adsorbate_atoms = Atoms(
        symbols=[atoms[i].symbol for i in adsorbate_indices],
        positions=[atoms[i].position for i in adsorbate_indices],
        cell=[VACUUM_BOX, VACUUM_BOX, VACUUM_BOX],
        pbc=True
    )
    adsorbate_atoms.center()
    return adsorbate_atoms


def write_scf_input(file_path, atoms, system_name, pseudo_dir, ecutwfc, ecutrho, degauss, k_scf):
    elements = sorted(set(atoms.get_chemical_symbols()))
    ntyp = len(elements)

    with open(file_path, "w", encoding="utf-8") as f:
        f.write(f"&CONTROL\n  calculation = 'scf'\n  prefix = '{system_name}_adsorbate'\n")
        f.write(f"  pseudo_dir = '{pseudo_dir}'\n  outdir = './tmp/', disk_io = 'medium', verbosity = 'high'\n/\n")
        f.write(f"&SYSTEM\n  ibrav = 0, nat = {len(atoms)}, ntyp = {ntyp}\n")
        f.write(f"  ecutwfc = {ecutwfc:.1f}, ecutrho = {ecutrho:.1f}\n")
        f.write(f"  occupations = 'smearing', smearing = 'gaussian', degauss = {degauss:.5f}, nspin = 2\n")
        f.write("  starting_magnetization(1) = 0.01\n")
        f.write("/\n")
        f.write(f"&ELECTRONS\n  conv_thr = {SCF_CONV_THR}\n  mixing_beta = 0.15\n")
        f.write("  mixing_mode = 'local-TF'\n")
        f.write("  electron_maxstep = 250\n/\n")
        f.write("\nCELL_PARAMETERS (angstrom)\n")
        for vec in atoms.get_cell():
            f.write(f"  {vec[0]:15.10f} {vec[1]:15.10f} {vec[2]:15.10f}\n")
        f.write("\nATOMIC_SPECIES\n")
        for elem in elements:
            f.write(f"  {elem:<3} 1.0  {PSEUDOS[elem]}\n")
        f.write("\nATOMIC_POSITIONS (angstrom)\n")
        for atom in atoms:
            f.write(f"  {atom.symbol:<3} {atom.position[0]:14.9f} {atom.position[1]:14.9f} {atom.position[2]:14.9f}\n")
        f.write(f"\nK_POINTS (automatic)\n  {k_scf[0]} {k_scf[1]} {k_scf[2]} {k_scf[3]} {k_scf[4]} {k_scf[5]}\n")


def write_nscf_input(file_path, atoms, system_name, pseudo_dir, ecutwfc, ecutrho, k_nscf):
    elements = sorted(set(atoms.get_chemical_symbols()))
    ntyp = len(elements)

    with open(file_path, "w", encoding="utf-8") as f:
        f.write(f"&CONTROL\n  calculation = 'nscf'\n  prefix = '{system_name}_adsorbate'\n")
        f.write(f"  pseudo_dir = '{pseudo_dir}'\n  outdir = './tmp/', disk_io = 'low', verbosity = 'high'\n/\n")
        f.write(f"&SYSTEM\n  ibrav = 0, nat = {len(atoms)}, ntyp = {ntyp}\n")
        f.write(f"  ecutwfc = {ecutwfc:.1f}, ecutrho = {ecutrho:.1f}\n")
        f.write("  occupations = 'smearing', smearing = 'gaussian', degauss = 0.005, nspin = 2\n")
        f.write("  starting_magnetization(1) = 0.01\n")
        f.write("/\n")
        f.write(f"&ELECTRONS\n  conv_thr = {NSCF_CONV_THR}\n")
        f.write("  diago_david_ndim = 6\n  diago_thr_init = 1.0d-4\n  electron_maxstep = 300\n/\n")
        f.write("\nCELL_PARAMETERS (angstrom)\n")
        for vec in atoms.get_cell():
            f.write(f"  {vec[0]:15.10f} {vec[1]:15.10f} {vec[2]:15.10f}\n")
        f.write("\nATOMIC_SPECIES\n")
        for elem in elements:
            f.write(f"  {elem:<3} 1.0  {PSEUDOS[elem]}\n")
        f.write("\nATOMIC_POSITIONS (angstrom)\n")
        for atom in atoms:
            f.write(f"  {atom.symbol:<3} {atom.position[0]:14.9f} {atom.position[1]:14.9f} {atom.position[2]:14.9f}\n")
        f.write(f"\nK_POINTS (automatic)\n  {k_nscf[0]} {k_nscf[1]} {k_nscf[2]} 0 0 0\n")


def write_projwfc_input(file_path: Path, system_name: str, filpdos: str, delta_e: float, e_min: float, e_max: float):
    with open(file_path, "w", encoding="utf-8") as f:
        f.write("&PROJWFC\n")
        f.write(f"  prefix = '{system_name}_adsorbate'\n")
        f.write("  outdir = './tmp/'\n")
        f.write("  ngauss = -1\n")
        f.write(f"  DeltaE = {delta_e:.4f}\n")
        f.write(f"  Emin = {e_min:.2f}\n")
        f.write(f"  Emax = {e_max:.2f}\n")
        f.write(f"  filpdos = '{filpdos}'\n")
        f.write("/\n")


def run_chain(run_dir: Path, system_name: str, np: int, nk: int):
    pref = f"{system_name}_adsorbate"
    cmds = [
        ("SCF", f"pw.x -nk {nk} < {pref}_scf.in > {pref}_scf.out"),
        ("NSCF", f"pw.x -nk {nk} < {pref}_nscf.in > {pref}_nscf.out"),
        ("PROJWFC", f"projwfc.x < {pref}_projwfc.in > {pref}_projwfc.out")
    ]

    for label, cmd_base in cmds:
        cmd = f"srun {cmd_base}" if "SLURM_JOB_ID" in os.environ else f"mpirun -np {np} {cmd_base}"
        print(f"  [{label}] {cmd}")
        subprocess.run(cmd, shell=True, cwd=str(run_dir), check=True)


def prepare_single_system(system_name: str, input_root: Path, run_root: Path, pseudo_dir: str, **kwargs):
    xyz_file = input_root / f"{system_name}_final.xyz"
    if not xyz_file.exists():
        return print(f"[SKIP] Missing structure: {xyz_file}")

    atoms = read_combi_structure(xyz_file)
    adsorbate = filter_adsorbate_atoms(atoms)

    if adsorbate is None or len(adsorbate) == 0:
        return print(f"[SKIP] No Li/S atoms found in {system_name}")

    run_dir = run_root / f"{system_name}_adsorbate"
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "tmp").mkdir(exist_ok=True)

    pref = f"{system_name}_adsorbate"
    write_scf_input(run_dir / f"{pref}_scf.in", adsorbate, system_name, pseudo_dir, kwargs['ecutwfc'], kwargs['ecutrho'], kwargs['degauss'], kwargs['k_scf'])
    write_nscf_input(run_dir / f"{pref}_nscf.in", adsorbate, system_name, pseudo_dir, kwargs['ecutwfc'], kwargs['ecutrho'], kwargs['k_nscf'])
    write_projwfc_input(run_dir / f"{pref}_projwfc.in", system_name, kwargs['filpdos'], kwargs['deltae'], kwargs['emin'], kwargs['emax'])

    print(f"Generated: {system_name}_adsorbate ({len(adsorbate)} atoms)")

    if kwargs['run']:
        run_chain(run_dir, system_name, kwargs['np'], kwargs['nk'])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--system", choices=SYSTEMS)
    parser.add_argument("--input-root", type=Path, default=Path("/scratch/anizami/QE_2/u/combi/final_xyz"))
    parser.add_argument("--run-root", type=Path, default=Path("/scratch/anizami/QE_2/u/PDOS/adsorbates_combi"))
    parser.add_argument("--pseudo-dir", default="/scratch/anizami/QE_2/PAW_pslib")
    parser.add_argument("--ecutwfc", type=float, default=60.0)
    parser.add_argument("--ecutrho", type=float, default=480.0)
    parser.add_argument("--degauss", type=float, default=0.005)
    parser.add_argument("--k-scf", default="1 1 1 0 0 0")
    parser.add_argument("--k-nscf", default="1 1 1 0 0 0")
    parser.add_argument("--filpdos", default="adsorbate_pdos")
    parser.add_argument("--deltae", type=float, default=0.05)
    parser.add_argument("--emin", type=float, default=-12.0)
    parser.add_argument("--emax", type=float, default=8.0)
    parser.add_argument("--np", type=int, default=8)
    parser.add_argument("--nk", type=int, default=2)
    parser.add_argument("--run", action="store_true")
    args = parser.parse_args()

    targets = [args.system] if args.system else SYSTEMS
    k_scf = parse_kgrid(args.k_scf)
    k_nscf = parse_kgrid(args.k_nscf)

    for system in targets:
        prepare_single_system(system, args.input_root, args.run_root, args.pseudo_dir,
                            ecutwfc=args.ecutwfc, ecutrho=args.ecutrho, degauss=args.degauss,
                            k_scf=k_scf, k_nscf=k_nscf, filpdos=args.filpdos,
                            deltae=args.deltae, emin=args.emin, emax=args.emax,
                            run=args.run, np=args.np, nk=args.nk)


if __name__ == "__main__":
    main()
