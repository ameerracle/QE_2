#!/usr/bin/env python3
## RELAX THE SLABS WITH A FIXED U VALUE
import os
import subprocess
import argparse
from pathlib import Path
from ase.io import read

# --- Settings & Constants ---
HUBBARD_MAP = {
    "Ti": 5.36,
    "V":  5.87,
    "Sc": 3.56,
    "Nb": 3.29,
    "Zr": 2.97,
}
PSEUDOS = {
    "Ti": "ti_pbe_v1.4.uspp.F.UPF",
    "V":  "v_pbe_v1.4.uspp.F.UPF",
    "Sc": "sc_pbe_v1.uspp.F.UPF",
    "Nb": "nb_pbe_v1.uspp.F.UPF",
    "Zr": "zr_pbe_v1.uspp.F.UPF",
    "N":  "n_pbe_v1.2.uspp.F.UPF",
}

# --- Smart Pathing ---
SCRIPT_DIR  = Path(__file__).resolve().parent
SLAB_DIR    = SCRIPT_DIR / "uma144_slab"
RUN_ROOT    = SCRIPT_DIR
PSEUDO_DIR  = "/scratch/anizami/QE_2/USPP/"

SLAB_NAMES = ["TiN", "VN", "ScN", "NbN", "ZrN"]


def generate_and_run(name: str):
    metal = name.rstrip("N")
    if metal not in HUBBARD_MAP:
        print(f"Error: Metal {metal} not in Hubbard Map.")
        return

    u_value  = HUBBARD_MAP[metal]
    manifold = "4d" if metal in ["Nb", "Zr"] else "3d"

    xyz_file = (SLAB_DIR / f"{name}_144_slab_uma.xyz").resolve() # Updated filename based on user breadcrumbs
    if not xyz_file.exists():
        # Fallback to general name if specific one doesn't exist
        xyz_file = (SLAB_DIR / f"{name}_144_slab.xyz").resolve()
        
    if not xyz_file.exists():
        print(f"--- FAILED TO FIND FILE: {xyz_file} ---")
        return

    atoms = read(str(xyz_file))
    cell  = atoms.get_cell()

    run_type = "magnetic" if name == "VN" else "non_mag"
    run_dir  = (RUN_ROOT / f"{name}_{run_type}_relax_144").resolve()
    run_dir.mkdir(exist_ok=True)

    pwi_file = run_dir / f"{name}_slab.in"
    pwo_file = run_dir / f"{name}_slab.out"

    print(f"--- Preparing {name} Slab ({len(atoms)} atoms) | U={u_value} | Smearing=0.02 ---")

    with open(pwi_file, "w") as f:
        # 1. CONTROL
        f.write("&CONTROL\n")
        f.write("  calculation = 'relax'\n")
        f.write(f"  prefix = '{name}_slab'\n")
        f.write(f"  pseudo_dir = '{PSEUDO_DIR}'\n")
        f.write("  outdir = './tmp/'\n")
        f.write("  disk_io = 'low'\n")
        f.write("  verbosity = 'low'\n")
        f.write(f"  forc_conv_thr = 0.000778\n")
        f.write("  nstep = 150\n")
        f.write("/\n")

        # 2. SYSTEM
        f.write("&SYSTEM\n")
        f.write(f"  ibrav = 0, nat = {len(atoms)}, ntyp = 2\n")
        f.write("  ecutwfc = 45.0, ecutrho = 450.0\n")
        f.write("  occupations = 'smearing', smearing = 'cold', degauss = 0.015\n")
        f.write("  vdw_corr = 'dft-d3'\n")
        f.write("  dftd3_version = 4\n")
        if name == "VN":
            f.write("  nspin = 2\n")
            f.write("  starting_magnetization(1) = 1.0\n")
        else:
            f.write("  nspin = 1\n")
        f.write("/\n")

        # 3. ELECTRONS
        f.write("&ELECTRONS\n")
        f.write("  conv_thr = 5.0d-6\n")
        if name in ["VN", "TiN"]:
            f.write("  mixing_beta = 0.20\n")
            f.write("  mixing_mode = 'local-TF'\n")
            f.write("  electron_maxstep = 270\n")
        else:
            f.write("  mixing_beta = 0.3\n")
            f.write("  electron_maxstep = 200\n")
        f.write("/\n")

        # 4. IONS
        f.write("&IONS\n  ion_dynamics = 'bfgs'\n/\n")

        # 5. CELL_PARAMETERS
        f.write("\nCELL_PARAMETERS (angstrom)\n")
        for vec in cell:
            f.write(f"  {vec[0]:15.10f} {vec[1]:15.10f} {vec[2]:15.10f}\n")

        # 6. ATOMIC_SPECIES
        f.write("\nATOMIC_SPECIES\n")
        f.write(f"  {metal}  1.0  {PSEUDOS[metal]}\n")
        f.write(f"  N   1.0  {PSEUDOS['N']}\n")

        # 7. ATOMIC_POSITIONS
        f.write("\nATOMIC_POSITIONS (angstrom)\n")
        for atom in atoms:
            label = atom.symbol
            # Fix atoms with z < 13A
            fix   = "0 0 0" if atom.position[2] < 13.0 else "1 1 1"
            f.write(f"  {label:4} {atom.position[0]:14.9f} {atom.position[1]:14.9f} {atom.position[2]:14.9f} {fix}\n")

        # 8. K_POINTS
        f.write("\nK_POINTS (automatic)\n  4 4 1 0 0 0\n")

        # 9. HUBBARD
        f.write(f"\nHUBBARD {{ortho-atomic}}\n")
        f.write(f"U {metal}-{manifold} {u_value}\n")


    # --- Execution ---
    if "SLURM_JOB_ID" in os.environ:
        exe_cmd = f"srun pw.x -nk 2 -in {pwi_file.name} > {pwo_file.name}"
    else:
        exe_cmd = f"mpirun -np 64 pw.x -nk 2 -in {pwi_file.name} > {pwo_file.name}"

    try:
        subprocess.run(exe_cmd, shell=True, cwd=str(run_dir), check=True)
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
        for target in SLAB_NAMES:
            generate_and_run(target)