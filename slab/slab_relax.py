import os
import subprocess
import argparse
import numpy as np
from pathlib import Path
from ase.io import read

# --- Settings & Constants ---
HUBBARD_MAP = {
    "Ti": (5.70, 5.60),  # (Surface U, Bulk U)
    "V":  (5.15, 5.09),
    "Sc": (3.54, 3.56),
    "Nb": (3.42, 3.51),
    "Zr": (2.95, 3.03),
}

PSEUDOS = {
    "Ti": "ti_pbe_v1.4.uspp.F.UPF",
    "V":  "v_pbe_v1.4.uspp.F.UPF",
    "Sc": "sc_pbe_v1.uspp.F.UPF",
    "Nb": "nb_pbe_v1.uspp.F.UPF",
    "Zr": "zr_pbe_v1.uspp.F.UPF",
    "N":  "n_pbe_v1.2.uspp.F.UPF",
}

SLAB_DIR = Path("./input_slab")
RUN_ROOT = Path(".")
PSEUDO_DIR = "/scratch/anizami/QE_2/USPP/"
SLAB_NAMES = ["TiN", "VN", "ScN", "NbN", "ZrN"]

def generate_and_run(name: str):
    metal = name.rstrip("N")
    u_surf, u_bulk = HUBBARD_MAP[metal]
    manifold = "4d" if metal in ["Nb", "Zr"] else "3d"
    
    xyz_file = SLAB_DIR / f"{name}_144_slab.xyz"
    if not xyz_file.exists():
        print(f"Skipping {name}: {xyz_file} not found.")
        return

    # Load atoms - this file has Tag 1 and Tag 2 in the 5th column
    atoms = read(str(xyz_file))
    cell = atoms.get_cell()
    
    # Create directory for this specific metal
    run_dir = RUN_ROOT / f"{name}_relax_production"
    run_dir.mkdir(exist_ok=True)
    
    # Define file names
    pwi_file = run_dir / f"{name}_slab.in"
    pwo_file = run_dir / f"{name}_slab.out"
    
    print(f"--- Preparing {name} Slab (144 atoms) ---")
    print(f"Surface U: {u_surf} | Bulk U: {u_bulk} | Manifold: {manifold}")

    with open(pwi_file, "w") as f:
        # 1. CONTROL
        f.write(f"&CONTROL\n")
        f.write(f"  calculation = 'relax'\n")
        f.write(f"  prefix = '{name}_slab'\n")
        f.write(f"  pseudo_dir = '{PSEUDO_DIR}'\n")
        f.write(f"  outdir = './tmp/'\n")
        f.write(f"  disk_io = 'low'\n")
        f.write(f"  verbosity = 'low'\n")
        f.write(f"/\n")

        # 2. SYSTEM
        f.write(f"&SYSTEM\n")
        f.write(f"  ibrav = 0, nat = {len(atoms)}, ntyp = 3\n")
        f.write(f"  ecutwfc = 45.0, ecutrho = 450.0\n")
        f.write(f"  occupations = 'smearing', smearing = 'cold', degauss = 0.01\n")
        # VN is special: must be magnetic
        if name == "VN":
            f.write(f"  nspin = 2\n")
            f.write(f"  starting_magnetization(1) = 1.0\n") # Surface V
            f.write(f"  starting_magnetization(2) = 1.0\n") # Bulk V
        f.write(f"/\n")

        # 3. ELECTRONS
        f.write(f"&ELECTRONS\n")
        f.write(f"  conv_thr = 1.0d-6\n")
        f.write(f"  mixing_beta = 0.3\n")
        f.write(f"  electron_maxstep = 125\n")
        f.write(f"/\n")

        # 4. IONS
        f.write(f"&IONS\n  ion_dynamics = 'bfgs'\n/\n")

        # 5. CELL_PARAMETERS
        f.write("\nCELL_PARAMETERS (angstrom)\n")
        for vec in cell:
            f.write(f"  {vec[0]:15.10f} {vec[1]:15.10f} {vec[2]:15.10f}\n")

        # 6. ATOMIC_SPECIES
        f.write(f"\nATOMIC_SPECIES\n")
        f.write(f"  {metal}1  1.0  {PSEUDOS[metal]}\n") # Surface type
        f.write(f"  {metal}2  1.0  {PSEUDOS[metal]}\n") # Bulk type
        f.write(f"  N   0.0  {PSEUDOS['N']}\n")

        # 7. ATOMIC_POSITIONS (Logic using Tags)
        f.write("\nATOMIC_POSITIONS (angstrom)\n")
        for atom in atoms:
            # Species label based on Tag
            if atom.symbol == metal:
                label = f"{atom.symbol}{int(atom.tag)}" # Metal1 or Metal2
            else:
                label = "N"
            
            # Constraint: Tag 2 is fixed (0 0 0), Tag 1 is relaxed (1 1 1)
            # Tag 2 was defined in your maker script as the interior layers
            if int(atom.tag) == 2:
                fix = "0 0 0"
            else:
                fix = "1 1 1"
            
            f.write(f"  {label:4} {atom.position[0]:14.9f} {atom.position[1]:14.9f} {atom.position[2]:14.9f} {fix}\n")

        # 8. K_POINTS
        f.write(f"\nK_POINTS (automatic)\n  4 4 1 0 0 0\n")

        # 9. HUBBARD
        f.write(f"\nHUBBARD {{ortho-atomic}}\n")
        f.write(f"U {metal}1-{manifold} {u_surf}\n")
        f.write(f"U {metal}2-{manifold} {u_bulk}\n")

    print(f"Launching srun pw.x for {name}...")
    
    # Running the process
    # Output file will be e.g., VN_slab.out
    # Final execution line
    cmd = f"srun pw.x -nk 6 < {pwi_file.name} > {pwo_file.name}"
    
    print(f"Launching: {cmd}")
    try:
        # Run without 'stdout=out_f' because '>' handles it in the cmd string
        subprocess.run(cmd, shell=True, cwd=str(run_dir), check=True)
        print(f"Successfully finished {name}.")
    except Exception as e:
        print(f"Error during {name} run: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--structure", choices=SLAB_NAMES, help="Which metal nitride to run?")
    args = parser.parse_args()

    if args.structure:
        generate_and_run(args.structure)
    else:
        # If no argument, loop through all of them
        for target in SLAB_NAMES:
            generate_and_run(target)