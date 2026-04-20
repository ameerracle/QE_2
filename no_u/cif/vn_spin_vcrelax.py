import subprocess
import textwrap
from pathlib import Path

from ase.io import read


BASE_DIR = Path("/scratch/anizami/QE_2/cif/")
INPUT_DIR = BASE_DIR / "input_file"
RUN_DIR = BASE_DIR / "VN_spin_vcrelax"
PSEUDO_DIR = "/scratch/anizami/QE_2/USPP/"
PW_CMD = "srun pw.x"

PSEUDO_MAP = {
    "V": "v_pbe_v1.4.uspp.F.UPF",
    "N": "n_pbe_v1.2.uspp.F.UPF",
}


def build_input(atoms) -> str:
    pwi_content = textwrap.dedent(f"""\
    &CONTROL
      calculation = 'vc-relax'
      prefix = 'vn'
      pseudo_dir = '{PSEUDO_DIR}'
      outdir = './tmp/'
    /
    &SYSTEM
      ibrav = 0, nat = {len(atoms)}, ntyp = 2
      ecutwfc = 45, ecutrho = 450.0
      occupations = 'smearing', smearing = 'mv', degauss = 0.02
      nspin = 2
      starting_magnetization(1) = 1.0
      vdw_corr = 'dft-d3'
      dftd3_version = 4
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
     V 1.0 {PSEUDO_MAP['V']}
     N 1.0 {PSEUDO_MAP['N']}

    CELL_PARAMETERS (angstrom)
    """)

    for vec in atoms.get_cell():
        pwi_content += f" {vec[0]:15.10f} {vec[1]:15.10f} {vec[2]:15.10f}\n"

    pwi_content += "\nATOMIC_POSITIONS (angstrom)\n"
    for atom in atoms:
        pwi_content += (
            f" {atom.symbol} {atom.position[0]:14.9f} "
            f"{atom.position[1]:14.9f} {atom.position[2]:14.9f}\n"
        )

    pwi_content += "\nK_POINTS (automatic)\n 5 5 5 0 0 0\n"

    return pwi_content


def main() -> None:
    structure_path = INPUT_DIR / "VN_relaxed.extxyz"
    try:
        atoms = read(structure_path)
    except Exception as exc:
        print(f"Skipping {structure_path}: {exc}")
        return

    symbols = atoms.get_chemical_symbols()
    v_indices = [index for index, symbol in enumerate(symbols) if symbol == "V"]
    n_indices = [index for index, symbol in enumerate(symbols) if symbol == "N"]

    if not v_indices or not n_indices:
        print(f"Skipping {structure_path}: missing V or N atoms")
        return

    atoms = atoms[v_indices + n_indices]

    RUN_DIR.mkdir(parents=True, exist_ok=True)
    input_path = RUN_DIR / "vcrelax.in"
    output_path = RUN_DIR / "vcrelax.out"

    input_path.write_text(build_input(atoms))

    print("--- Running VC-RELAX for VN (spin polarized, no +U) ---")
    subprocess.run(
        f"{PW_CMD} < {input_path.name} > {output_path.name}",
        shell=True,
        cwd=RUN_DIR,
        check=False,
    )
    print("--- VN Finished Relaxation ---")


if __name__ == "__main__":
    main()