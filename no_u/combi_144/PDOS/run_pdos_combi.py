#!/usr/bin/env python3
"""
Generate and optionally run QE PDOS workflow for relaxed slab+adsorbate combinations.

Workflow per combination:
1) Read final structure from XYZ (produced from QE relax outputs)
2) Write pw.x SCF input
3) Write pw.x NSCF input
4) Write pp.x input (for Bader charge density extraction)
5) Write projwfc.x input
6) Optionally run: pw.x (scf) -> pp.x -> pw.x (nscf) -> projwfc.x

Example:
    python run_pdos_combi.py --slab TiN --ads Li2S4
    python run_pdos_combi.py --ads S8 Li2S8 --mode bader-only --run
  python run_pdos_combi.py --run
"""

from __future__ import annotations

import argparse
import os
import subprocess
from pathlib import Path

from ase.io import read


PSEUDOS = {
    "Ti": "ti_pbe_v1.4.uspp.F.UPF",
    "V": "v_pbe_v1.4.uspp.F.UPF",
    "Sc": "sc_pbe_v1.uspp.F.UPF",
    "Nb": "nb_pbe_v1.uspp.F.UPF",
    "Zr": "zr_pbe_v1.uspp.F.UPF",
    "N": "n_pbe_v1.2.uspp.F.UPF",
    "Li": "li_pbe_v1.4.uspp.F.UPF",
    "S": "s_pbe_v1.4.uspp.F.UPF",
}

SLABS = ["TiN", "VN", "ScN", "NbN", "ZrN"]
FULL_PDOS_ADS = {"Li2S4"}
BADER_ONLY_ADS = {"S8", "Li2S8"}
FIXED_NK = 2

# Keep convergence at 1e-8 (do not tighten beyond this value).
SCF_CONV_THR = 1.0e-8
NSCF_CONV_THR = 1.0e-8


def parse_kgrid(text: str) -> list[int]:
    vals = [int(x) for x in text.split()]
    if len(vals) != 6:
        raise ValueError("K grid must contain 6 integers: kx ky kz sx sy sz")
    return vals


def ordered_species(tag: str, atoms) -> list[str]:
    slab = tag.split("_", 1)[0]
    metal = slab.rstrip("N")
    present = set(atoms.get_chemical_symbols())

    species = []
    for sym in [metal, "N", "Li", "S"]:
        if sym in present:
            species.append(sym)

    for sym in atoms.get_chemical_symbols():
        if sym not in species:
            species.append(sym)

    return species


def format_conv_thr(value: float) -> str:
    return f"{value:.1e}".replace("e", "d")


def write_scf_input(
    file_path: Path,
    atoms,
    tag: str,
    pseudo_dir: str,
    ecutwfc: float,
    ecutrho: float,
    degauss: float,
    k_scf: list[int],
):
    nat = len(atoms)
    species = ordered_species(tag, atoms)
    ntyp = len(species)

    with open(file_path, "w", encoding="utf-8") as f:
        f.write("&CONTROL\n")
        f.write("  calculation = 'scf'\n")
        f.write(f"  prefix = '{tag}_pdos'\n")
        f.write(f"  pseudo_dir = '{pseudo_dir}'\n")
        f.write("  outdir = './tmp/'\n")
        f.write("  disk_io = 'medium'\n")
        f.write("  verbosity = 'low'\n")
        f.write("/\n")

        f.write("&SYSTEM\n")
        f.write(f"  ibrav = 0, nat = {nat}, ntyp = {ntyp}\n")
        f.write(f"  ecutwfc = {ecutwfc:.1f}, ecutrho = {ecutrho:.1f}\n")
        f.write(f"  occupations = 'smearing', smearing = 'cold', degauss = {degauss:.5f}\n")
        f.write("  nspin = 1\n")
        f.write("/\n")

        f.write("&ELECTRONS\n")
        f.write(f"  conv_thr = {format_conv_thr(SCF_CONV_THR)}\n")
        f.write("  mixing_beta = 0.25\n")
        f.write("  electron_maxstep = 250\n")
        f.write("/\n")

        f.write("\nCELL_PARAMETERS (angstrom)\n")
        for vec in atoms.get_cell():
            f.write(f"  {vec[0]:15.10f} {vec[1]:15.10f} {vec[2]:15.10f}\n")

        f.write("\nATOMIC_SPECIES\n")
        for sym in species:
            pseudo = PSEUDOS.get(sym)
            if pseudo is None:
                raise ValueError(f"Missing pseudopotential mapping for element: {sym}")
            f.write(f"  {sym:<3} 1.0  {pseudo}\n")

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
    tag: str,
    pseudo_dir: str,
    ecutwfc: float,
    ecutrho: float,
    k_nscf: list[int],
):
    nat = len(atoms)
    species = ordered_species(tag, atoms)
    ntyp = len(species)

    with open(file_path, "w", encoding="utf-8") as f:
        f.write("&CONTROL\n")
        f.write("  calculation = 'nscf'\n")
        f.write(f"  prefix = '{tag}_pdos'\n")
        f.write(f"  pseudo_dir = '{pseudo_dir}'\n")
        f.write("  outdir = './tmp/'\n")
        f.write("  disk_io = 'low'\n")
        f.write("  verbosity = 'low'\n")
        f.write("/\n")

        f.write("&SYSTEM\n")
        f.write(f"  ibrav = 0, nat = {nat}, ntyp = {ntyp}\n")
        f.write(f"  ecutwfc = {ecutwfc:.1f}, ecutrho = {ecutrho:.1f}\n")
        f.write("  occupations = 'tetrahedra_opt'\n")
        f.write("  nspin = 1\n")
        f.write("/\n")

        f.write("&ELECTRONS\n")
        f.write(f"  conv_thr = {format_conv_thr(NSCF_CONV_THR)}\n")
        f.write("  mixing_beta = 0.25\n")
        f.write("  electron_maxstep = 250\n")
        f.write("/\n")

        f.write("\nCELL_PARAMETERS (angstrom)\n")
        for vec in atoms.get_cell():
            f.write(f"  {vec[0]:15.10f} {vec[1]:15.10f} {vec[2]:15.10f}\n")

        f.write("\nATOMIC_SPECIES\n")
        for sym in species:
            pseudo = PSEUDOS.get(sym)
            if pseudo is None:
                raise ValueError(f"Missing pseudopotential mapping for element: {sym}")
            f.write(f"  {sym:<3} 1.0  {pseudo}\n")

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
    tag: str,
    filpdos: str,
    delta_e: float,
    e_min: float,
    e_max: float,
):
    with open(file_path, "w", encoding="utf-8") as f:
        f.write("&PROJWFC\n")
        f.write(f"  prefix = '{tag}_pdos'\n")
        f.write("  outdir = './tmp/'\n")
        f.write("  ngauss = -1\n")
        f.write(f"  DeltaE = {delta_e:.4f}\n")
        f.write(f"  Emin = {e_min:.2f}\n")
        f.write(f"  Emax = {e_max:.2f}\n")
        f.write(f"  filpdos = '{filpdos}'\n")
        f.write("/\n")


def write_pp_input(file_path: Path, tag: str, filplot: str, fileout: str):
    with open(file_path, "w", encoding="utf-8") as f:
        f.write("&INPUTPP\n")
        f.write(f"  prefix = '{tag}_pdos'\n")
        f.write("  outdir = './tmp/'\n")
        f.write("  plot_num = 0\n")
        f.write(f"  filplot = '{filplot}'\n")
        f.write("/\n")
        f.write("&PLOT\n")
        f.write("  iflag = 3\n")
        f.write("  output_format = 6\n")
        f.write(f"  fileout = '{fileout}'\n")
        f.write("/\n")


def run_chain(run_dir: Path, tag: str, np: int, do_pdos: bool):
    scf_in = f"{tag}_scf.in"
    scf_out = f"{tag}_scf.out"
    pp_in = f"{tag}_pp.in"
    pp_out = f"{tag}_pp.out"
    nscf_in = f"{tag}_nscf.in"
    nscf_out = f"{tag}_nscf.out"
    proj_in = f"{tag}_projwfc.in"
    proj_out = f"{tag}_projwfc.out"

    if "SLURM_JOB_ID" in os.environ:
        scf_cmd = f"srun pw.x -nk {FIXED_NK} < {scf_in} > {scf_out}"
        pp_cmd = f"srun pp.x < {pp_in} > {pp_out}"
        nscf_cmd = f"srun pw.x -nk {FIXED_NK} < {nscf_in} > {nscf_out}"
        proj_cmd = f"srun projwfc.x < {proj_in} > {proj_out}"
    else:
        scf_cmd = f"mpirun -np {np} pw.x -nk {FIXED_NK} < {scf_in} > {scf_out}"
        pp_cmd = f"mpirun -np {np} pp.x < {pp_in} > {pp_out}"
        nscf_cmd = f"mpirun -np {np} pw.x -nk {FIXED_NK} < {nscf_in} > {nscf_out}"
        proj_cmd = f"mpirun -np {np} projwfc.x < {proj_in} > {proj_out}"

    steps = [("SCF", scf_cmd), ("PP", pp_cmd)]
    if do_pdos:
        steps.extend([("NSCF", nscf_cmd), ("PROJWFC", proj_cmd)])

    for label, cmd in steps:
        print(f"  [{label}] {cmd}")
        subprocess.run(cmd, shell=True, cwd=str(run_dir), check=True)

        if label == "PP":
            cube_file = run_dir / f"{tag}_charge.cube"
            if cube_file.exists():
                print(f"  [PP] Wrote {cube_file.name}")


def discover_tags(xyz_dir: Path, slab: str | None, ads_list: list[str]) -> list[str]:
    tags = []
    for ads in ads_list:
        pattern = f"*_{ads}_combi_final.xyz"
        for xyz in sorted(xyz_dir.glob(pattern)):
            name = xyz.name
            tag = name.removesuffix("_combi_final.xyz")
            slab_name = tag.split("_", 1)[0]
            if slab is not None and slab_name != slab:
                continue
            tags.append(tag)
    return tags


def resolve_workflow_mode(mode: str, ads: str) -> str:
    if mode == "full":
        return "full"
    if mode == "bader-only":
        return "bader-only"

    if ads in FULL_PDOS_ADS:
        return "full"
    if ads in BADER_ONLY_ADS:
        return "bader-only"
    return "full"


def prepare_single_combi(
    tag: str,
    xyz_dir: Path,
    run_root: Path,
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
    mode: str,
    run_jobs: bool,
    np: int,
):
    xyz_file = xyz_dir / f"{tag}_combi_final.xyz"
    if not xyz_file.exists():
        print(f"[SKIP] Missing final xyz: {xyz_file}")
        return

    try:
        atoms = read(str(xyz_file))
    except Exception as exc:
        print(f"[FAIL] Could not read {xyz_file}: {exc}")
        return

    run_dir = run_root / f"{tag}_pdos"
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "tmp").mkdir(exist_ok=True)

    scf_in = run_dir / f"{tag}_scf.in"
    pp_in = run_dir / f"{tag}_pp.in"

    ads = tag.split("_", 1)[1]
    workflow = resolve_workflow_mode(mode, ads)
    do_pdos = workflow == "full"

    nscf_in = run_dir / f"{tag}_nscf.in"
    proj_in = run_dir / f"{tag}_projwfc.in"

    print(f"\n--- {tag} | nat={len(atoms)} ---")
    print(f"Input structure: {xyz_file}")
    print(f"PDOS directory:  {run_dir}")

    try:
        write_scf_input(
            file_path=scf_in,
            atoms=atoms,
            tag=tag,
            pseudo_dir=pseudo_dir,
            ecutwfc=ecutwfc,
            ecutrho=ecutrho,
            degauss=degauss,
            k_scf=k_scf,
        )
        if do_pdos:
            write_nscf_input(
                file_path=nscf_in,
                atoms=atoms,
                tag=tag,
                pseudo_dir=pseudo_dir,
                ecutwfc=ecutwfc,
                ecutrho=ecutrho,
                k_nscf=k_nscf,
            )
    except ValueError as exc:
        print(f"[FAIL] {tag}: {exc}")
        return

    write_pp_input(
        file_path=pp_in,
        tag=tag,
        filplot=f"{tag}_rho",
        fileout=f"{tag}_charge.cube",
    )

    if do_pdos:
        write_projwfc_input(
            file_path=proj_in,
            tag=tag,
            filpdos=filpdos,
            delta_e=delta_e,
            e_min=e_min,
            e_max=e_max,
        )

    print("Wrote input files:")
    print(f"  - {scf_in.name}")
    print(f"  - {pp_in.name}")
    if do_pdos:
        print(f"  - {nscf_in.name}")
        print(f"  - {proj_in.name}")
    print(f"Workflow mode:   {workflow}")

    if run_jobs:
        try:
            run_chain(run_dir, tag=tag, np=np, do_pdos=do_pdos)
            print(f"[OK] Finished PDOS chain for {tag}")
        except subprocess.CalledProcessError as exc:
            print(f"[FAIL] Run failed for {tag}: {exc}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate/run QE SCF+PP+NSCF+PROJWFC workflow for final combi XYZ"
    )
    parser.add_argument("--slab", choices=SLABS, help="Run only one slab family (e.g., TiN)")
    parser.add_argument(
        "--ads",
        nargs="+",
        default=["Li2S4", "S8", "Li2S8"],
        help="Adsorbate labels in tag (default: Li2S4 S8 Li2S8)",
    )
    parser.add_argument(
        "--xyz-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "final_xyz",
        help="Directory containing *_combi_final.xyz",
    )
    parser.add_argument(
        "--run-root",
        type=Path,
        default=Path(__file__).resolve().parent / "runs",
        help="Directory where PDOS run folders are created",
    )
    parser.add_argument("--pseudo-dir", default="/scratch/anizami/QE_2/USPP/")
    parser.add_argument("--ecutwfc", type=float, default=45.0)
    parser.add_argument("--ecutrho", type=float, default=450.0)
    parser.add_argument("--degauss", type=float, default=0.015)
    parser.add_argument("--k-scf", default="4 4 1 0 0 0")
    parser.add_argument("--k-nscf", default="8 8 1 0 0 0")
    parser.add_argument("--filpdos", default="combi_pdos")
    parser.add_argument("--deltae", type=float, default=0.05)
    parser.add_argument("--emin", type=float, default=-12.0)
    parser.add_argument("--emax", type=float, default=8.0)
    parser.add_argument("--np", type=int, default=64, help="MPI ranks (non-SLURM mode)")
    parser.add_argument(
        "--mode",
        choices=["auto", "full", "bader-only"],
        default="auto",
        help="auto: Li2S4 full PDOS, S8/Li2S8 Bader-only; or force full/bader-only",
    )
    parser.add_argument(
        "--run",
        action="store_true",
        help="If set, run selected QE steps after writing inputs",
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    try:
        k_scf = parse_kgrid(args.k_scf)
        k_nscf = parse_kgrid(args.k_nscf)
    except ValueError as exc:
        raise SystemExit(f"Invalid k-grid: {exc}") from exc

    tags = discover_tags(args.xyz_dir, args.slab, args.ads)
    if not tags:
        print(
            f"[INFO] No matching xyz files in {args.xyz_dir} "
            f"for slab={args.slab or 'ALL'}, ads={','.join(args.ads)}"
        )
        return

    print(f"SCF conv_thr fixed to {SCF_CONV_THR:.1e}")
    print(f"NSCF conv_thr fixed to {NSCF_CONV_THR:.1e}")
    print(f"pw.x -nk fixed to {FIXED_NK}")

    for tag in tags:
        prepare_single_combi(
            tag=tag,
            xyz_dir=args.xyz_dir,
            run_root=args.run_root,
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
            mode=args.mode,
            run_jobs=args.run,
            np=args.np,
        )


if __name__ == "__main__":
    main()
