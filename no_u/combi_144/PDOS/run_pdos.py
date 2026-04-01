#!/usr/bin/env python3
"""
Generate and optionally run QE PDOS workflow for relaxed slabs.

Workflow per slab:
1) Read relaxed structure from QE output with ASE
2) Write pw.x SCF input
3) Write pw.x NSCF input
4) Write pp.x input (for Bader charge density extraction)
5) Write projwfc.x input
6) Optionally run: pw.x (scf) -> pp.x -> pw.x (nscf) -> projwfc.x

Example:
  python run_slab_pdos.py --structure TiN --run
  python run_slab_pdos.py --run
"""

import argparse
import os
import shutil
import subprocess
from pathlib import Path

from ase.io import read


PSEUDOS = {
    "Ti": "ti_pbe_v1.4.uspp.F.UPF",
    "V":  "v_pbe_v1.4.uspp.F.UPF",
    "Sc": "sc_pbe_v1.uspp.F.UPF",
    "Nb": "nb_pbe_v1.uspp.F.UPF",
    "Zr": "zr_pbe_v1.uspp.F.UPF",
    "N":  "n_pbe_v1.2.uspp.F.UPF",
    "Li": "li_pbe_v1.4.uspp.F.UPF",
    "S":  "s_pbe_v1.4.uspp.F.UPF"
}

SLABS = ["TiN", "VN", "ScN", "NbN", "ZrN"]


def parse_kgrid(text: str) -> list[int]:
    vals = [int(x) for x in text.split()]
    if len(vals) != 6:
        raise ValueError("K grid must contain 6 integers: kx ky kz sx sy sz")
    return vals


def read_relaxed_atoms(qe_out: Path):
    """Read the last structure from a QE pw.x output file."""
    try:
        return read(str(qe_out), format="espresso-out", index=-1)
    except Exception:
        # Fallback to ASE auto-detection if explicit parser fails.
        return read(str(qe_out), index=-1)


def write_scf_input(
    file_path: Path,
    atoms,
    slab_name: str,
    pseudo_dir: str,
    ecutwfc: float,
    ecutrho: float,
    degauss: float,
    k_scf: list[int],
):
    metal = slab_name.rstrip("N")
    nat = len(atoms)
    ntyp = len(set(atoms.get_chemical_symbols()))

    with open(file_path, "w", encoding="utf-8") as f:
        f.write("&CONTROL\n")
        f.write("  calculation = 'scf'\n")
        f.write(f"  prefix = '{slab_name}_pdos'\n")
        f.write(f"  pseudo_dir = '{pseudo_dir}'\n")
        f.write("  outdir = './tmp/'\n")
        f.write("  disk_io = 'medium'\n")  # Keeps charge density for pp.x
        f.write("  verbosity = 'low'\n")
        f.write("/\n")

        f.write("&SYSTEM\n")
        f.write(f"  ibrav = 0, nat = {nat}, ntyp = {ntyp}\n")
        f.write(f"  ecutwfc = {ecutwfc:.1f}, ecutrho = {ecutrho:.1f}\n")
        f.write(f"  occupations = 'smearing', smearing = 'cold', degauss = {degauss:.5f}\n")
        f.write("  nspin = 1\n")
        f.write("/\n")

        f.write("&ELECTRONS\n")
        f.write("  conv_thr = 1.0d-8\n")
        f.write("  mixing_beta = 0.25\n")
        f.write("  electron_maxstep = 250\n")
        f.write("/\n")

        f.write("\nCELL_PARAMETERS (angstrom)\n")
        for vec in atoms.get_cell():
            f.write(f"  {vec[0]:15.10f} {vec[1]:15.10f} {vec[2]:15.10f}\n")

        f.write("\nATOMIC_SPECIES\n")
        # Keep metal first, then nitrogen for consistency with existing scripts.
        if metal in PSEUDOS:
            f.write(f"  {metal:<3} 1.0  {PSEUDOS[metal]}\n")
        f.write(f"  N   1.0  {PSEUDOS['N']}\n")

        f.write("\nATOMIC_POSITIONS (angstrom)\n")
        for atom in atoms:
            f.write(
                f"  {atom.symbol:<3} "
                f"{atom.position[0]:14.9f} {atom.position[1]:14.9f} {atom.position[2]:14.9f}\n"
            )

        f.write("\nK_POINTS (automatic)\n")
        f.write(
            f"  {k_scf[0]} {k_scf[1]} {k_scf[2]} "
            f"{k_scf[3]} {k_scf[4]} {k_scf[5]}\n"
        )


def write_nscf_input(
    file_path: Path,
    atoms,
    slab_name: str,
    pseudo_dir: str,
    ecutwfc: float,
    ecutrho: float,
    k_nscf: list[int],
):
    metal = slab_name.rstrip("N")
    nat = len(atoms)
    ntyp = len(set(atoms.get_chemical_symbols()))

    with open(file_path, "w", encoding="utf-8") as f:
        f.write("&CONTROL\n")
        f.write("  calculation = 'nscf'\n")
        f.write(f"  prefix = '{slab_name}_pdos'\n")
        f.write(f"  pseudo_dir = '{pseudo_dir}'\n")
        f.write("  outdir = './tmp/'\n")
        f.write("  disk_io = 'low'\n")
        f.write("  verbosity = 'low'\n")
        f.write("/\n")

        f.write("&SYSTEM\n")
        f.write(f"  ibrav = 0, nat = {nat}, ntyp = {ntyp}\n")
        f.write(f"  ecutwfc = {ecutwfc:.1f}, ecutrho = {ecutrho:.1f}\n")
        f.write("  occupations = 'tetrahedra_opt'\n") # Crucial for clean PDOS!
        f.write("  nspin = 1\n")
        f.write("/\n")

        f.write("&ELECTRONS\n")
        f.write("  conv_thr = 1.0d-8\n")
        f.write("  mixing_beta = 0.25\n")
        f.write("  electron_maxstep = 250\n")
        f.write("/\n")

        f.write("\nCELL_PARAMETERS (angstrom)\n")
        for vec in atoms.get_cell():
            f.write(f"  {vec[0]:15.10f} {vec[1]:15.10f} {vec[2]:15.10f}\n")

        f.write("\nATOMIC_SPECIES\n")
        if metal in PSEUDOS:
            f.write(f"  {metal:<3} 1.0  {PSEUDOS[metal]}\n")
        f.write(f"  N   1.0  {PSEUDOS['N']}\n")

        f.write("\nATOMIC_POSITIONS (angstrom)\n")
        for atom in atoms:
            f.write(
                f"  {atom.symbol:<3} "
                f"{atom.position[0]:14.9f} {atom.position[1]:14.9f} {atom.position[2]:14.9f}\n"
            )

        f.write("\nK_POINTS (automatic)\n")
        f.write(
            f"  {k_nscf[0]} {k_nscf[1]} {k_nscf[2]} "
            f"{k_nscf[3]} {k_nscf[4]} {k_nscf[5]}\n"
        )


def write_projwfc_input(
    file_path: Path,
    slab_name: str,
    filpdos: str,
    delta_e: float,
    e_min: float,
    e_max: float,
):
    with open(file_path, "w", encoding="utf-8") as f:
        f.write("&PROJWFC\n")
        f.write(f"  prefix = '{slab_name}_pdos'\n")
        f.write("  outdir = './tmp/'\n")
        f.write("  ngauss = -1\n") # Note: -1 is ignored if tetrahedra are used, but good to have
        f.write(f"  DeltaE = {delta_e:.4f}\n")
        f.write(f"  Emin = {e_min:.2f}\n")
        f.write(f"  Emax = {e_max:.2f}\n")
        f.write(f"  filpdos = '{filpdos}'\n")
        f.write("/\n")


def write_pp_input(
    file_path: Path,
    slab_name: str,
    filplot: str,
    fileout: str,
):
    with open(file_path, "w", encoding="utf-8") as f:
        f.write("&INPUTPP\n")
        f.write(f"  prefix = '{slab_name}_pdos'\n")
        f.write("  outdir = './tmp/'\n")
        f.write("  plot_num = 0\n")
        f.write(f"  filplot = '{filplot}'\n")
        f.write("/\n")
        f.write("&PLOT\n")
        f.write("  iflag = 3\n")
        f.write("  output_format = 6\n")
        f.write(f"  fileout = '{fileout}'\n")
        f.write("/\n")


def run_chain(run_dir: Path, slab_name: str, np: int, nk: int):
    scf_in = f"{slab_name}_scf.in"
    scf_out = f"{slab_name}_scf.out"
    pp_in = f"{slab_name}_pp.in"
    pp_out = f"{slab_name}_pp.out"
    nscf_in = f"{slab_name}_nscf.in"
    nscf_out = f"{slab_name}_nscf.out"
    proj_in = f"{slab_name}_projwfc.in"
    proj_out = f"{slab_name}_projwfc.out"

    # Consolidate job execution commands
    if "SLURM_JOB_ID" in os.environ:
        scf_cmd = f"srun pw.x -nk {nk} < {scf_in} > {scf_out}"
        pp_cmd = f"srun pp.x < {pp_in} > {pp_out}"
        nscf_cmd = f"srun pw.x -nk {nk} < {nscf_in} > {nscf_out}"
        proj_cmd = f"srun projwfc.x < {proj_in} > {proj_out}"
    else:
        scf_cmd = f"mpirun -np {np} pw.x -nk {nk} < {scf_in} > {scf_out}"
        pp_cmd = f"mpirun -np {np} pp.x < {pp_in} > {pp_out}"
        nscf_cmd = f"mpirun -np {np} pw.x -nk {nk} < {nscf_in} > {nscf_out}"
        proj_cmd = f"mpirun -np {np} projwfc.x < {proj_in} > {proj_out}"

    # Execution order matters! PP directly after SCF guarantees pristine rho extraction.
    for label, cmd in [
        ("SCF", scf_cmd),
        ("PP", pp_cmd),
        ("NSCF", nscf_cmd),
        ("PROJWFC", proj_cmd),
    ]:
        print(f"  [{label}] {cmd}")
        subprocess.run(cmd, shell=True, cwd=str(run_dir), check=True)

        if label == "PP":
            cube_file = run_dir / f"{slab_name}_charge.cube"
            if cube_file.exists():
                print(f"  [PP] Wrote {cube_file.name}")


def prepare_single_slab(
    slab_name: str,
    output_root: Path,
    pseudo_dir: str,
    ecutwfc: float,
    ecutrho: float,
    degauss: float,
    k_scf: list[int],
    k_nscf: list[int],
    filpdos: str,
    delta_e: float,
    e_min: float,
    e_max: float,
    run_jobs: bool,
    np: int,
    nk: int,
):
    slab_dir = output_root / f"{slab_name}_slab_relax"
    qe_out = slab_dir / f"{slab_name}_slab.out"
    if not qe_out.exists():
        print(f"[SKIP] Missing relax output: {qe_out}")
        return

    try:
        atoms = read_relaxed_atoms(qe_out)
    except Exception as exc:
        print(f"[FAIL] Could not parse {qe_out}: {exc}")
        return

    run_dir = slab_dir / "pdos"
    run_dir.mkdir(exist_ok=True)
    (run_dir / "tmp").mkdir(exist_ok=True)

    scf_in = run_dir / f"{slab_name}_scf.in"
    nscf_in = run_dir / f"{slab_name}_nscf.in"
    pp_in = run_dir / f"{slab_name}_pp.in"
    proj_in = run_dir / f"{slab_name}_projwfc.in"

    print(f"\n--- {slab_name} | nat={len(atoms)} ---")
    print(f"Input structure: {qe_out}")
    print(f"PDOS directory:  {run_dir}")

    write_scf_input(
        file_path=scf_in,
        atoms=atoms,
        slab_name=slab_name,
        pseudo_dir=pseudo_dir,
        ecutwfc=ecutwfc,
        ecutrho=ecutrho,
        degauss=degauss,
        k_scf=k_scf,
    )
    write_nscf_input(
        file_path=nscf_in,
        atoms=atoms,
        slab_name=slab_name,
        pseudo_dir=pseudo_dir,
        ecutwfc=ecutwfc,
        ecutrho=ecutrho,
        k_nscf=k_nscf,
    )
    write_projwfc_input(
        file_path=proj_in,
        slab_name=slab_name,
        filpdos=filpdos,
        delta_e=delta_e,
        e_min=e_min,
        e_max=e_max,
    )
    write_pp_input(
        file_path=pp_in,
        slab_name=slab_name,
        filplot=f"{slab_name}_rho",
        fileout=f"{slab_name}_charge.cube",
    )

    print("Wrote input files:")
    print(f"  - {scf_in.name}")
    print(f"  - {pp_in.name}")
    print(f"  - {nscf_in.name}")
    print(f"  - {proj_in.name}")

    if run_jobs:
        try:
            run_chain(run_dir, slab_name, np=np, nk=nk)
            print(f"[OK] Finished PDOS chain for {slab_name}")
        except subprocess.CalledProcessError as exc:
            print(f"[FAIL] Run failed for {slab_name}: {exc}")


def build_parser():
    parser = argparse.ArgumentParser(
        description="Generate/run QE SCF+PP+NSCF+PROJWFC workflow for relaxed slabs"
    )
    parser.add_argument("--structure", choices=SLABS, help="Run only one slab")
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path(__file__).resolve().parent / "output_slab_144",
        help="Folder containing *_slab_relax subdirectories",
    )
    parser.add_argument("--pseudo-dir", default="/scratch/anizami/QE_2/USPP/")
    parser.add_argument("--ecutwfc", type=float, default=45.0)
    parser.add_argument("--ecutrho", type=float, default=450.0)
    parser.add_argument("--degauss", type=float, default=0.015)
    parser.add_argument("--k-scf", default="4 4 1 0 0 0")
    parser.add_argument("--k-nscf", default="8 8 1 0 0 0")
    parser.add_argument("--filpdos", default="slab_pdos")
    parser.add_argument("--deltae", type=float, default=0.05)
    parser.add_argument("--emin", type=float, default=-12.0)
    parser.add_argument("--emax", type=float, default=8.0)
    parser.add_argument("--np", type=int, default=64, help="MPI ranks (non-SLURM mode)")
    parser.add_argument("--nk", type=int, default=2, help="pw.x k-point pools")
    parser.add_argument(
        "--run",
        action="store_true",
        help="If set, run pw.x SCF/NSCF then projwfc.x after writing inputs",
    )
    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    try:
        k_scf = parse_kgrid(args.k_scf)
        k_nscf = parse_kgrid(args.k_nscf)
    except ValueError as exc:
        raise SystemExit(f"Invalid k-grid: {exc}") from exc

    targets = [args.structure] if args.structure else SLABS

    for slab_name in targets:
        prepare_single_slab(
            slab_name=slab_name,
            output_root=args.output_root,
            pseudo_dir=args.pseudo_dir,
            ecutwfc=args.ecutwfc,
            ecutrho=args.ecutrho,
            degauss=args.degauss,
            k_scf=k_scf,
            k_nscf=k_nscf,
            filpdos=args.filpdos,
            delta_e=args.deltae,
            e_min=args.emin,
            e_max=args.emax,
            run_jobs=args.run,
            np=args.np,
            nk=args.nk,
        )


if __name__ == "__main__":
    main()