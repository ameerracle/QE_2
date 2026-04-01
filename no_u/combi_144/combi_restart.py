#!/usr/bin/env python3
"""
Restart script for QE relaxations. 
Uses restart_mode='restart' to continue from the .save folder in the tmp directory.
"""
import subprocess
import argparse
from pathlib import Path
from ase.io import read

# --- Settings & Constants (Identical to original) ---
PSEUDOS = {
    "Ti": "ti_pbe_v1.4.uspp.F.UPF",
    "V":  "v_pbe_v1.4.uspp.F.UPF",
    "Sc": "sc_pbe_v1.uspp.F.UPF",
    "Nb": "nb_pbe_v1.uspp.F.UPF",
    "Zr": "zr_pbe_v1.uspp.F.UPF",
    "N":  "n_pbe_v1.2.uspp.F.UPF",
    "Li": "li_pbe_v1.4.uspp.F.UPF",
    "S":  "s_pbe_v1.4.uspp.F.UPF",
}

COMBI_DIR  = Path(".")
RUN_ROOT   = Path(".")
PSEUDO_DIR = "/scratch/anizami/QE_2/USPP/"

SLAB_NAMES = ["TiN", "VN", "ScN", "NbN", "ZrN"]
ADS_NAMES  = ["Li2S", "Li2S2", "Li2S4", "Li2S6", "Li2S8", "S8"]

ECUTWFC = 45.0
ECUTRHO = 450.0
Z_FIX   = 13.0 

def generate_and_restart(slab: str, ads: str, dry_run: bool = False):
    """Generate QE restart input and launch."""
    metal = slab.rstrip("N")
    tag   = f"{slab}_{ads}"

    # Use the original XYZ to get the structure (QE will actually update 
    # positions from the .save file, but the input needs the atom list)
    xyz_file = COMBI_DIR / f"{tag}_combi_144.xyz"
    if not xyz_file.exists():
        print(f"Skipping {tag}: {xyz_file} not found.")
        return

    atoms = read(str(xyz_file))
    cell  = atoms.get_cell()
    seen, unique_elements = set(), []
    for sym in atoms.get_chemical_symbols():
        if sym not in seen:
            seen.add(sym)
            unique_elements.append(sym)

    ntyp = len(unique_elements)

    # Output directory and files (Points to the EXISTING run directory)
    run_dir  = RUN_ROOT / f"{tag}_combi_relax"
    pwi_file = run_dir / f"{tag}_combi_restart.in"
    pwo_file = run_dir / f"{tag}_combi_restart.out"

    if not run_dir.exists():
        print(f"Warning: {run_dir} does not exist. Restart might fail if no .save found.")

    print(f"\n--- Restarting {tag} Relaxation ---")
    print(f"Using existing data in: {run_dir}/tmp/")

    with open(pwi_file, "w") as f:
        # &CONTROL
        f.write("&CONTROL\n")
        f.write(f"  calculation = 'relax'\n")
        f.write(f"  restart_mode = 'restart'\n")
        f.write(f"  prefix = '{tag}_combi'\n")
        f.write(f"  pseudo_dir = '{PSEUDO_DIR}'\n")
        f.write(f"  outdir = './tmp/'\n")
        f.write(f"  disk_io = 'low'\n")
        f.write(f"  verbosity = 'low'\n")
        f.write(f"  forc_conv_thr = 0.00156\n")
        f.write(f"  nstep = 160\n")
        f.write("/\n")

        # &SYSTEM
        f.write("&SYSTEM\n")
        f.write(f"  ibrav = 0, nat = {len(atoms)}, ntyp = {ntyp}\n")
        f.write(f"  ecutwfc = {ECUTWFC}, ecutrho = {ECUTRHO}\n")
        f.write(f"  occupations = 'smearing', smearing = 'cold', degauss = 0.025\n")
        f.write(f"  vdw_corr = 'dft-d3', dftd3_version = 4\n")
        f.write("/\n")

        # &ELECTRONS
        f.write("&ELECTRONS\n")
        f.write(f"  conv_thr = 7.5d-6\n")
        f.write(f"  mixing_beta = 0.25\n")
        f.write(f"  electron_maxstep = 250\n")
        f.write("/\n")

        # &IONS
        f.write("&IONS\n")
        f.write(f"  ion_dynamics = 'bfgs'\n")
        f.write(f"  upscale = 750\n")
        f.write("/\n")

        # CELL_PARAMETERS
        f.write("\nCELL_PARAMETERS (angstrom)\n")
        for vec in cell:
            f.write(f"  {vec[0]:15.10f} {vec[1]:15.10f} {vec[2]:15.10f}\n")

        # ATOMIC_SPECIES  (dynamic — only elements actually present)
        f.write("\nATOMIC_SPECIES\n")
        for sym in unique_elements:
            f.write(f"  {sym:<4} 1.0  {PSEUDOS[sym]}\n")

        # ATOMIC_POSITIONS
        f.write("\nATOMIC_POSITIONS (angstrom)\n")
        for atom in atoms:
            fix = "0 0 0" if atom.z < Z_FIX else "1 1 1"
            f.write(f"  {atom.symbol:<4} {atom.position[0]:14.9f} {atom.position[1]:14.9f} {atom.position[2]:14.9f}  {fix}\n")

        # K_POINTS
        f.write("\nK_POINTS (automatic)\n")
        f.write("  4 4 1 0 0 0\n")

    cmd = f"srun pw.x -nk 4 < {pwi_file.name} > {pwo_file.name}"

    if dry_run:
        print(f"[dry-run] Input written to {pwi_file}")
        print(f"[dry-run] Would run: {cmd}")
        return

    print(f"Launching srun pw.x for {tag} (restart) ...")
    print(f"Command: {cmd}\n")

    try:
        subprocess.run(cmd, shell=True, cwd=str(run_dir), check=True)
        print(f"✓ {tag} restart finished successfully.")
    except subprocess.CalledProcessError as e:
        print(f"✗ {tag} restart failed (exit code {e.returncode}).")
    except Exception as e:
        print(f"✗ {tag} restart error: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Restart QE relaxations.")
    parser.add_argument("--structure", choices=SLAB_NAMES, default=None)
    parser.add_argument("--adsorbate", choices=ADS_NAMES, default=None)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    slabs = [args.structure] if args.structure else SLAB_NAMES
    adses = [args.adsorbate] if args.adsorbate else ADS_NAMES

    for slab in slabs:
        for ads in adses:
            generate_and_restart(slab, ads, dry_run=args.dry_run)