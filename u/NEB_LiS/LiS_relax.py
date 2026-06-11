#!/usr/bin/env python3
"""
Relax all .xyz structures in this folder using Quantum ESPRESSO (pw.x).
Updated with stabilized convergence parameters for magnetic energy materials.

Run:
    python Li_relax.py --xyz NbN_Li_combi1_144.xyz
    python Li_relax.py --pattern "*combi1*.xyz"
    python Li_relax.py --all
"""

import argparse
import subprocess
import os
from pathlib import Path
from ase.io import read

PSEUDOS = {
    'Ti': 'ti_pbe_v1.4.uspp.F.UPF',
    'V':  'v_pbe_v1.4.uspp.F.UPF',
    'Sc': 'sc_pbe_v1.uspp.F.UPF',
    'Nb': 'nb_pbe_v1.uspp.F.UPF',
    'Zr': 'zr_pbe_v1.uspp.F.UPF',
    'N':  'n_pbe_v1.2.uspp.F.UPF',
    'Li': 'li_pbe_v1.4.uspp.F.UPF',
    'S':  's_pbe_v1.4.uspp.F.UPF'
}

WORK_DIR = Path(__file__).resolve().parent
PSEUDO_DIR = (WORK_DIR.parent.parent / "USPP").resolve()

ECUTWFC = 45.0
ECUTRHO = 450.0
HUBBARD_U = {
    "V": 2.50,
}
Z_FIX = 13.0
FIX_INDICES = {144, 145}

def get_unique_elements(atoms):
    """Return element symbols preserving their first occurrence order."""
    seen = set()
    unique = []
    for sym in atoms.get_chemical_symbols():
        if sym not in seen:
            seen.add(sym)
            unique.append(sym)
    return unique

def generate_and_run(xyz_file: Path, dry_run: bool = False):
    """Generate QE input and optionally run a relax job for one .xyz file."""
    base = xyz_file.stem
    atoms = read(str(xyz_file))
    cell = atoms.get_cell()
    unique_elements = get_unique_elements(atoms)
    ntyp = len(unique_elements)

    run_dir = WORK_DIR / f"{base}_relax"
    run_dir.mkdir(exist_ok=True)

    pwi_file = run_dir / f"{base}.in"
    pwo_file = run_dir / f"{base}.out"

    print(f"\n--- Preparing relax for {xyz_file.name} ---")
    print(f"nat = {len(atoms)}, ntyp = {ntyp}, elements = {unique_elements}")
    print(f"run directory: {run_dir}")

    with open(pwi_file, "w", encoding="utf-8") as f:
        f.write("&CONTROL\n")
        f.write("  calculation = 'relax'\n")
        f.write("  restart_mode = 'from_scratch'\n")
        f.write(f"  prefix = '{base}'\n")
        f.write(f"  pseudo_dir = '{PSEUDO_DIR}'\n")
        f.write("  outdir = './tmp/'\n")
        f.write("  disk_io = 'low'\n")
        f.write("  verbosity = 'low'\n")
        f.write("  forc_conv_thr = 0.00156\n")
        f.write("  nstep = 160\n")
        f.write("/\n")

        f.write("&SYSTEM\n")
        f.write(f"  ibrav = 0, nat = {len(atoms)}, ntyp = {ntyp}\n")
        f.write(f"  ecutwfc = {ECUTWFC}, ecutrho = {ECUTRHO}\n")
        f.write("  occupations = 'smearing', smearing = 'cold', degauss = 0.02\n")
        f.write("  nspin = 2\n")
        f.write("  starting_magnetization(1) = 1.45\n")
        f.write("  starting_magnetization(2) = -0.2\n")
        f.write("  vdw_corr = 'dft-d3', dftd3_version = 4\n")
        f.write("/\n")

        f.write("&ELECTRONS\n")
        f.write("  conv_thr = 1.0d-6\n")
        f.write("  mixing_beta = 0.20\n")
        f.write("  mixing_ndim = 12\n")
        f.write("  mixing_mode = 'local-TF'\n")
        f.write("  electron_maxstep = 250\n")
        f.write("  diagonalization = 'rmm-davidson'\n")
        f.write("/\n")

        f.write("&IONS\n")
        f.write("  ion_dynamics = 'bfgs'\n")
        f.write("  upscale = 2\n")
        f.write("  pot_extrapolation = 'second_order'\n")
        f.write("  wfc_extrapolation = 'second_order'\n")
        f.write("/\n")

        f.write("\nCELL_PARAMETERS (angstrom)\n")
        for vec in cell:
            f.write(f"  {vec[0]:15.10f} {vec[1]:15.10f} {vec[2]:15.10f}\n")

        f.write("\nATOMIC_SPECIES\n")
        for sym in unique_elements:
            if sym not in PSEUDOS:
                raise KeyError(f"No pseudopotential defined for element '{sym}'")
            f.write(f"  {sym:<4} 1.0  {PSEUDOS[sym]}\n")

        f.write("\nATOMIC_POSITIONS (angstrom)\n")
        n_fixed = 0
        for atom in atoms:
            if atom.z < Z_FIX or atom.index in FIX_INDICES:
                fix = "0 0 0"
                n_fixed += 1
            else:
                fix = "1 1 1"
            f.write(
                f"  {atom.symbol:<4} "
                f"{atom.position[0]:14.9f} {atom.position[1]:14.9f} {atom.position[2]:14.9f}"
                f"  {fix}\n"
            )

        f.write("\nK_POINTS (automatic)\n")
        f.write("  4 4 1 0 0 0\n")

        if "V" in unique_elements:
            f.write("\nHUBBARD {atomic}\n")
            f.write(f"U V-3d {HUBBARD_U['V']}\n")

    print(f"Fixed {n_fixed} atoms (z < {Z_FIX:.1f} Angstrom), {len(atoms) - n_fixed} atoms free.")

    cmd = f"srun pw.x -nk 2 < {pwi_file.name} > {pwo_file.name}"

    if dry_run:
        print(f"[dry-run] would run: {cmd}")
        return

    print(f"Launching: {cmd}")
    try:
        subprocess.run(cmd, shell=True, cwd=str(run_dir), check=True)
        print(f"OK: {base} finished successfully.")
    except subprocess.CalledProcessError as exc:
        print(f"ERROR: {base} failed (exit code {exc.returncode}).")

def main():
    parser = argparse.ArgumentParser(description="Relax one or more .xyz structures.")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--xyz", default=None)
    parser.add_argument("--pattern", default=None)
    parser.add_argument("--all", action="store_true")
    args = parser.parse_args()

    selected_modes = sum([bool(args.xyz), bool(args.pattern), bool(args.all)])
    if selected_modes != 1:
        parser.error("Choose exactly one mode: --xyz, --pattern, or --all")

    if args.xyz:
        xyz_files = [WORK_DIR / args.xyz] if not Path(args.xyz).is_absolute() else [Path(args.xyz)]
    elif args.pattern:
        xyz_files = sorted(WORK_DIR.glob(args.pattern))
    else:
        xyz_files = sorted(WORK_DIR.glob("*.xyz"))

    xyz_files = [p for p in xyz_files if p.suffix.lower() == ".xyz"]

    if not xyz_files:
        print("No matching .xyz files found.")
        return

    for xyz_file in xyz_files:
        if xyz_file.exists():
            generate_and_run(xyz_file, dry_run=args.dry_run)

if __name__ == "__main__":
    main()
