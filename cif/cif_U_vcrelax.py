import os
import subprocess
import textwrap
from pathlib import Path
from ase.io import read

# --- Paths & Commands ---
BASE_DIR      = Path("/scratch/anizami/QE_2/cif/")
INPUT_DIR     = BASE_DIR / "input_file"
U_PRIM_DIR    = BASE_DIR / "final_U_cif"
PSEUDO_DIR    = "/scratch/anizami/QE_2/USPP/"
PW_CMD        = "srun pw.x"

U_PRIM_DIR.mkdir(parents=True, exist_ok=True)

PSEUDO_MAP = {
    "Ti": "ti_pbe_v1.4.uspp.F.UPF", "V": "v_pbe_v1.4.uspp.F.UPF",
    "Sc": "sc_pbe_v1.uspp.F.UPF",   "Nb": "nb_pbe_v1.uspp.F.UPF",
    "Zr": "zr_pbe_v1.uspp.F.UPF",   "N":  "n_pbe_v1.2.uspp.F.UPF"
}

metals = ["Ti", "V", "Sc", "Nb", "Zr"]

for metal in metals:
    name = f"{metal}N"
    # Unit cells location for each respective metal compound
    extxyz_file = INPUT_DIR / f"{name}_relaxed.extxyz"
    
    try:
        atoms = read(extxyz_file)
    except Exception as e:
        print(f"Skipping {extxyz_file}: {e}")
        continue

    symbols = atoms.get_chemical_symbols()
    
    metal_indices    = [i for i, s in enumerate(symbols) if s == metal]
    nitrogen_indices = [i for i, s in enumerate(symbols) if s == "N"]
    atoms = atoms[metal_indices + nitrogen_indices]

    U_VALUES = {
        "Ti": 5.36, "V": 5.87, "Sc": 3.56,
        "Nb": 3.29, "Zr": 2.97,
    }
    u_val = U_VALUES[metal]
    manifold    = "4d" if metal in ["Nb", "Zr"] else "3d"
    is_magnetic = (metal == "V")

    work_dir = U_PRIM_DIR / name
    work_dir.mkdir(exist_ok=True)
    os.chdir(work_dir)

    pwi_content = textwrap.dedent(f"""\
&CONTROL
  calculation = 'vc-relax'
  prefix = '{name.lower()}'
  pseudo_dir = '{PSEUDO_DIR}'
  outdir = './tmp/'
/
&SYSTEM
  ibrav = 0, nat = {len(atoms)}, ntyp = 2
  ecutwfc = 45, ecutrho = 450.0
  occupations = 'smearing', smearing = 'mv', degauss = 0.08
  nspin = {2 if is_magnetic else 1}
  vdw_corr = 'dft-d3'
  dftd3_version = 4
""")

    if is_magnetic:
        pwi_content += "  starting_magnetization(1) = 0.02\n"

    pwi_content += textwrap.dedent(f"""\
/
&ELECTRONS
  conv_thr = 1.0d-6
  mixing_beta = 0.3
  mixing_mode = 'local-TF'
  electron_maxstep = 150
/
&IONS
  ion_dynamics = 'bfgs'
/
&CELL
  cell_dynamics = 'bfgs'
/
ATOMIC_SPECIES
 {metal} 1.0 {PSEUDO_MAP[metal]}
 N 1.0 {PSEUDO_MAP['N']}

CELL_PARAMETERS (angstrom)
""")

    for vec in atoms.get_cell():
        pwi_content += f" {vec[0]:15.10f} {vec[1]:15.10f} {vec[2]:15.10f}\n"

    pwi_content += "\nATOMIC_POSITIONS (angstrom)\n"
    for atom in atoms:
        pwi_content += f" {atom.symbol} {atom.position[0]:14.9f} {atom.position[1]:14.9f} {atom.position[2]:14.9f}\n"

    # Changed k-points to 5 5 5 and formatted U-value using `.2f` (with values already given to 2 decimals)
    pwi_content += f"\nK_POINTS (automatic)\n 5 5 5 0 0 0\n"
    pwi_content += f"\nHUBBARD {{ortho-atomic}}\n U {metal}-{manifold} {u_val:.2f}\n"

    with open("vcrelax.in", "w") as f:
        f.write(pwi_content)

    print(f"--- Running VC-RELAX for {name} (Reordered) ---")
    subprocess.run(f"{PW_CMD} < vcrelax.in > vcrelax.out", shell=True)
    print(f"--- {name} Finished Relaxation ---")
