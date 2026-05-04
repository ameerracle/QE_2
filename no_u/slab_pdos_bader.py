#!/usr/bin/env python3
"""
Generate and optionally run QE PDOS workflow for relaxed slabs.
Outputs to the centralized /no_u/PDOS directory.
"""

import argparse
import os
import subprocess
from pathlib import Path
from ase.io import read

# Strictly PAW kjpaw_psl files
PSEUDOS = {
    "Ti": "Ti.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "V":  "V.pbe-spnl-kjpaw_psl.1.0.0.UPF",
    "Sc": "Sc.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Nb": "Nb.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Zr": "Zr.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "N":  "N.pbe-n-kjpaw_psl.1.0.0.UPF"
}

SLABS = ["TiN", "VN", "ScN", "NbN", "ZrN"]
SCF_CONV_THR = "1.0d-7"
NSCF_CONV_THR = "1.0d-7"

def parse_kgrid(text: str) -> list[int]:
    vals = [int(x) for x in text.split()]
    if len(vals) != 6:
        raise ValueError("K grid must contain 6 integers: kx ky kz sx sy sz")
    return vals

def read_relaxed_atoms(qe_out: Path):
    try:
        return read(str(qe_out), format="espresso-out", index=-1)
    except Exception:
        return read(str(qe_out), index=-1)

def write_scf_input(file_path, atoms, slab_name, pseudo_dir, ecutwfc, ecutrho, degauss, k_scf):
    metal = slab_name.rstrip("N")
    with open(file_path, "w", encoding="utf-8") as f:
        f.write(f"&CONTROL\n  calculation = 'scf'\n  prefix = '{slab_name}_slab'\n")
        f.write(f"  pseudo_dir = '{pseudo_dir}'\n  outdir = './tmp/', disk_io = 'medium', verbosity = 'high'\n/\n")
        f.write(f"&SYSTEM\n  ibrav = 0, nat = {len(atoms)}, ntyp = 2\n")
        f.write(f"  ecutwfc = {ecutwfc:.1f}, ecutrho = {ecutrho:.1f}\n")
        f.write(f"  occupations = 'smearing', smearing = 'cold', degauss = {degauss:.5f}, nspin = 1\n/\n")
        f.write(f"&ELECTRONS\n  conv_thr = {SCF_CONV_THR}\n  mixing_beta = 0.25\n")
        f.write("  electron_maxstep = 250\n  diago_david_ndim = 4\n/\n")
        f.write("\nCELL_PARAMETERS (angstrom)\n")
        for vec in atoms.get_cell(): f.write(f"  {vec[0]:15.10f} {vec[1]:15.10f} {vec[2]:15.10f}\n")
        f.write(f"\nATOMIC_SPECIES\n  {metal:<3} 1.0  {PSEUDOS[metal]}\n  N    1.0  {PSEUDOS['N']}\n")
        f.write("\nATOMIC_POSITIONS (angstrom)\n")
        for atom in atoms: f.write(f"  {atom.symbol:<3} {atom.position[0]:14.9f} {atom.position[1]:14.9f} {atom.position[2]:14.9f}\n")
        f.write(f"\nK_POINTS (automatic)\n  {k_scf[0]} {k_scf[1]} {k_scf[2]} {k_scf[3]} {k_scf[4]} {k_scf[5]}\n")

def write_nscf_input(file_path, atoms, slab_name, pseudo_dir, ecutwfc, ecutrho, k_nscf):
    metal = slab_name.rstrip("N")
    with open(file_path, "w", encoding="utf-8") as f:
        f.write(f"&CONTROL\n  calculation = 'nscf'\n  prefix = '{slab_name}_slab'\n")
        f.write(f"  pseudo_dir = '{pseudo_dir}'\n  outdir = './tmp/', disk_io = 'low', verbosity = 'high'\n/\n")
        f.write(f"&SYSTEM\n  ibrav = 0, nat = {len(atoms)}, ntyp = 2\n")
        f.write(f"  ecutwfc = {ecutwfc:.1f}, ecutrho = {ecutrho:.1f}\n")
        f.write("  occupations = 'tetrahedra_opt', nspin = 1\n/\n")
        f.write(f"&ELECTRONS\n  conv_thr = {NSCF_CONV_THR}\n")
        f.write("  diago_david_ndim = 6\n  diago_thr_init = 1.0d-4\n  electron_maxstep = 300\n/\n")
        f.write("\nCELL_PARAMETERS (angstrom)\n")
        for vec in atoms.get_cell(): f.write(f"  {vec[0]:15.10f} {vec[1]:15.10f} {vec[2]:15.10f}\n")
        f.write(f"\nATOMIC_SPECIES\n  {metal:<3} 1.0  {PSEUDOS[metal]}\n  N    1.0  {PSEUDOS['N']}\n")
        f.write("\nATOMIC_POSITIONS (angstrom)\n")
        for atom in atoms: f.write(f"  {atom.symbol:<3} {atom.position[0]:14.9f} {atom.position[1]:14.9f} {atom.position[2]:14.9f}\n")
        f.write(f"\nK_POINTS (automatic)\n  {k_nscf[0]} {k_nscf[1]} {k_nscf[2]} 0 0 0\n")

def write_projwfc_input(file_path: Path, slab_name: str, filpdos: str, delta_e: float, e_min: float, e_max: float):
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

def write_pp_input(file_path: Path, slab_name: str, filplot: str, fileout: str):
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

def run_chain(run_dir: Path, slab_name: str, np: int, nk: int):
    pref = f"{slab_name}_slab"
    cmds = [
        ("SCF", f"pw.x -nk {nk} < {pref}_scf.in > {pref}_scf.out"),
        ("PP", f"pp.x < {pref}_pp.in > {pref}_pp.out"),
        ("NSCF", f"pw.x -nk {nk} < {pref}_nscf.in > {pref}_nscf.out"),
        ("PROJWFC", f"projwfc.x < {pref}_projwfc.in > {pref}_projwfc.out")
    ]
    
    for label, cmd_base in cmds:
        cmd = f"srun {cmd_base}" if "SLURM_JOB_ID" in os.environ else f"mpirun -np {np} {cmd_base}"
        print(f"  [{label}] {cmd}")
        subprocess.run(cmd, shell=True, cwd=str(run_dir), check=True)

def prepare_single_slab(slab_name: str, input_root: Path, run_root: Path, pseudo_dir: str, **kwargs):
    # Old structure: output_slab_144/TiN_slab_relax/TiN_slab.out
    qe_out = input_root / f"{slab_name}_slab_relax" / f"{slab_name}_slab.out"
    if not qe_out.exists():
        return print(f"[SKIP] Missing relax output: {qe_out}")

    atoms = read_relaxed_atoms(qe_out)
    # New centralized structure: no_u/PDOS/TiN_slab/
    run_dir = run_root / f"{slab_name}_slab"
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "tmp").mkdir(exist_ok=True)

    pref = f"{slab_name}_slab"
    write_scf_input(run_dir / f"{pref}_scf.in", atoms, slab_name, pseudo_dir, kwargs['ecutwfc'], kwargs['ecutrho'], kwargs['degauss'], kwargs['k_scf'])
    write_nscf_input(run_dir / f"{pref}_nscf.in", atoms, slab_name, pseudo_dir, kwargs['ecutwfc'], kwargs['ecutrho'], kwargs['k_nscf'])
    write_projwfc_input(run_dir / f"{pref}_projwfc.in", slab_name, kwargs['filpdos'], kwargs['deltae'], kwargs['emin'], kwargs['emax'])
    write_pp_input(run_dir / f"{pref}_pp.in", slab_name, f"{pref}_rho", f"{pref}_charge.cube")

    if kwargs['run']:
        run_chain(run_dir, slab_name, kwargs['np'], kwargs['nk'])

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--structure", choices=SLABS)
    parser.add_argument("--input-root", type=Path, default=Path("/lustre10/scratch/anizami/QE_2/no_u/slab/output_slab_144"))
    parser.add_argument("--run-root", type=Path, default=Path("/lustre10/scratch/anizami/QE_2/no_u/PDOS"))
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

    targets = [args.structure] if args.structure else SLABS
    k_scf = parse_kgrid(args.k_scf)
    k_nscf = parse_kgrid(args.k_nscf)

    for slab in targets:
        prepare_single_slab(slab, args.input_root, args.run_root, args.pseudo_dir, 
                            ecutwfc=args.ecutwfc, ecutrho=args.ecutrho, degauss=args.degauss, 
                            k_scf=k_scf, k_nscf=k_nscf, filpdos=args.filpdos, 
                            deltae=args.deltae, emin=args.emin, emax=args.emax, 
                            run=args.run, np=args.np, nk=args.nk)

if __name__ == "__main__":
    main()
