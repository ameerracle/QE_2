#!/usr/bin/env python3
"""
Relax UMA-relaxed slabs using Quantum ESPRESSO (pw.x) with DFT-D3.
Fixes bottom layer, relaxes top layers.

Run: python slab_144_relax.py [--structure TiN|VN|ScN|NbN|ZrN]
"""
import os
import subprocess
import argparse
import numpy as np
from pathlib import Path
from ase.io import read

# --- Settings & Constants ---
PSEUDOS = {
    "Ti": "ti_pbe_v1.4.uspp.F.UPF",
    "V":  "v_pbe_v1.4.uspp.F.UPF",
    "Sc": "sc_pbe_v1.uspp.F.UPF",
    "Nb": "nb_pbe_v1.uspp.F.UPF",
    "Zr": "zr_pbe_v1.uspp.F.UPF",
    "N":  "n_pbe_v1.2.uspp.F.UPF",
}

SLAB_DIR = Path(".")
RUN_ROOT = Path(".")
PSEUDO_DIR = "/scratch/anizami/QE_2/USPP/"
SLAB_NAMES = ["TiN", "VN", "ScN", "NbN", "ZrN"]
ECUTWFC = 45.0
ECUTRHO = 450.0


def generate_and_run(name: str):
    """Generate QE input and run relaxation for a single slab."""
    metal = name.rstrip("N")
    
    # Input: UMA-relaxed structure
    xyz_file = SLAB_DIR / f"{name}_144_slab_uma.xyz"
    if not xyz_file.exists():
        print(f"Skipping {name}: {xyz_file} not found.")
        return

    atoms = read(str(xyz_file))
    cell = atoms.get_cell()
    
    # Output directory
    run_dir = RUN_ROOT / f"{name}_slab_relax"
    run_dir.mkdir(exist_ok=True)
    
    pwi_file = run_dir / f"{name}_slab.in"
    pwo_file = run_dir / f"{name}_slab.out"
    
    print(f"\n--- Preparing {name} Slab Relaxation (UMA pre-relaxed) ---")
    print(f"Input: {xyz_file}")
    print(f"Output dir: {run_dir}")

    with open(pwi_file, "w") as f:
        # Control
        f.write(f"&CONTROL\n")
        f.write(f"  calculation = 'relax'\n")
        f.write(f"  prefix = '{name}_slab'\n")
        f.write(f"  pseudo_dir = '{PSEUDO_DIR}'\n")
        f.write(f"  outdir = './tmp/'\n")
        f.write(f"  disk_io = 'low'\n")
        f.write(f"  verbosity = 'low'\n")
        f.write(f"  forc_conv_thr = 0.000778\n")
        f.write(f"/\n")

        # System
        f.write(f"&SYSTEM\n")
        f.write(f"  ibrav = 0, nat = {len(atoms)}, ntyp = 2\n")
        f.write(f"  ecutwfc = {ECUTWFC}, ecutrho = {ECUTRHO}\n")
        f.write(f"  occupations = 'smearing', smearing = 'cold', degauss = 0.015\n")
        f.write(f"  vdw_corr = 'dft-d3', dftd3_version = 4\n")
        # no spin polarization
        f.write(f"/\n")

        # Electrons
        f.write(f"&ELECTRONS\n")
        f.write(f"  conv_thr = 1.0d-6\n")
        f.write(f"  mixing_beta = 0.25\n")
        f.write(f"  electron_maxstep = 285\n")
        f.write(f"/\n")

        # Ions
        f.write(f"&IONS\n  ion_dynamics = 'bfgs'\n/\n")

        # Cell parameters
        f.write("\nCELL_PARAMETERS (angstrom)\n")
        for vec in cell:
            f.write(f"  {vec[0]:15.10f} {vec[1]:15.10f} {vec[2]:15.10f}\n")

        # Atomic species
        f.write(f"\nATOMIC_SPECIES\n")
        f.write(f"  {metal}  1.0  {PSEUDOS[metal]}\n")
        f.write(f"  N      0.0  {PSEUDOS['N']}\n")

        # Atomic positions with constraints
        # Fix bottom 2 layers (z < 13.0 Å), relax top 2 layers
        z_threshold = 13.0
        
        f.write("\nATOMIC_POSITIONS (angstrom)\n")
        n_fixed = 0
        for atom in atoms:
            if atom.z < z_threshold:
                fix = "0 0 0"
                n_fixed += 1
            else:
                fix = "1 1 1"
            
            f.write(f"  {atom.symbol:2} {atom.position[0]:14.9f} {atom.position[1]:14.9f} {atom.position[2]:14.9f} {fix}\n")
        
        print(f"Fixed {n_fixed} atoms (z < {z_threshold:.2f} Å)")

        # K-points
        f.write(f"\nK_POINTS (automatic)\n  4 4 1 0 0 0\n")

    print(f"Launching pw.x for {name}...")
    
    cmd = f" pw.x -nk 6 < {pwi_file.name} > {pwo_file.name}"
    print(f"Command: {cmd}\n")
    
    try:
        subprocess.run(cmd, shell=True, cwd=str(run_dir), check=True)
        print(f"✓ Successfully finished {name}.")
    except Exception as e:
        print(f"✗ Error during {name} run: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Relax UMA-pre-relaxed slabs using QE with DFT-D3."
    )
    parser.add_argument(
        "--structure",
        choices=SLAB_NAMES,
        help="Which metal nitride to run (if omitted, runs all)."
    )
    args = parser.parse_args()

    if args.structure:
        generate_and_run(args.structure)
    else:
        for target in SLAB_NAMES:
            generate_and_run(target)

    print("\nDone!")
