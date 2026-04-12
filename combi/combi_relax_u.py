#!/usr/bin/env python3
"""
Relax combined slab+adsorbate structures using Quantum ESPRESSO (pw.x) with DFT-D3 + Hubbard U.
Fixes bottom slab layers (z < 13 A), relaxes top slab layers and adsorbate.

Input files are searched in this folder with either naming style:
- {Metal}_{Adsorbate}_combi_144.xyz
- {Metal}_{Adsorbate}_combiU_144.xyz

Run:
python combi_relax_u.py [--structure TiN|VN|ScN|NbN|ZrN] [--adsorbate Li2S|Li2S2|Li2S4|Li2S6|Li2S8|S8] [--dry-run]
"""

import argparse
import os
import subprocess
from pathlib import Path

from ase.io import read


# --- Settings & Constants ---
HUBBARD_MAP = {
    "Ti": 4.44,
    "V": 2.50,
}

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

SCRIPT_DIR = Path(__file__).resolve().parent
COMBI_DIR = SCRIPT_DIR
RUN_ROOT = SCRIPT_DIR
PSEUDO_DIR = "/scratch/anizami/QE_2/USPP/"

SLAB_NAMES = ["TiN", "VN"]
ADS_NAMES = ["Li2S", "Li2S2", "Li2S4", "Li2S6", "Li2S8", "S8"]

ECUTWFC = 45.0
ECUTRHO = 450.0
Z_FIX = 13.0


def find_xyz_file(slab: str, ads: str) -> Path | None:
    """Return the first matching xyz path for supported naming variants."""
    tag = f"{slab}_{ads}"
    candidates = [
        COMBI_DIR / f"{tag}_combi_144.xyz",
        COMBI_DIR / f"{tag}_combiU_144.xyz",
    ]
    for xyz_file in candidates:
        if xyz_file.exists():
            return xyz_file
    return None


def generate_and_run(slab: str, ads: str, dry_run: bool = False) -> None:
    """Generate QE input and optionally run one slab+adsorbate relaxation."""
    metal = slab.rstrip("N")
    if metal not in HUBBARD_MAP:
        print(f"Skipping {slab}_{ads}: no Hubbard U value configured for {metal}.")
        return

    u_value = HUBBARD_MAP[metal]
    manifold = "4d" if metal in ["Nb", "Zr"] else "3d"

    tag = f"{slab}_{ads}"
    xyz_file = find_xyz_file(slab, ads)
    if xyz_file is None:
        print(f"Skipping {tag}: no matching xyz file found in {COMBI_DIR}.")
        return

    atoms = read(str(xyz_file))
    cell = atoms.get_cell()

    # Preserve first-occurrence order for species indexing in QE cards.
    seen = set()
    unique_elements = []
    for sym in atoms.get_chemical_symbols():
        if sym not in seen:
            seen.add(sym)
            unique_elements.append(sym)

    if metal not in unique_elements:
        print(f"Skipping {tag}: metal {metal} not present in structure.")
        return

    ntyp = len(unique_elements)

    run_dir = RUN_ROOT / f"{tag}_combi_relax_u"
    run_dir.mkdir(exist_ok=True)
    pwi_file = run_dir / f"{tag}_combi_u.in"
    pwo_file = run_dir / f"{tag}_combi_u.out"

    print(f"\n--- Preparing {tag} (+U) ---")
    print(f"Input: {xyz_file} | nat = {len(atoms)}, ntyp = {ntyp}")
    print(f"Elements: {unique_elements}")
    print(f"U({metal}-{manifold}) = {u_value}")
    print(f"Output dir: {run_dir}")

    with open(pwi_file, "w", encoding="utf-8") as f:
        # CONTROL
        f.write("&CONTROL\n")
        f.write("  calculation = 'relax'\n")
        f.write(f"  prefix = '{tag}_combi_u'\n")
        f.write(f"  pseudo_dir = '{PSEUDO_DIR}'\n")
        f.write("  outdir = './tmp/'\n")
        f.write("  disk_io = 'low'\n")
        f.write("  verbosity = 'low'\n")
        f.write("  forc_conv_thr = 0.00156\n")
        f.write("  nstep = 160\n")
        f.write("/\n")

        # SYSTEM
        f.write("&SYSTEM\n")
        f.write(f"  ibrav = 0, nat = {len(atoms)}, ntyp = {ntyp}\n")
        f.write(f"  ecutwfc = {ECUTWFC}, ecutrho = {ECUTRHO}\n")
        f.write("  occupations = 'smearing', smearing = 'cold', degauss = 0.025\n")
        f.write("  vdw_corr = 'dft-d3', dftd3_version = 4\n")
        if slab == "VN":
            metal_index = unique_elements.index(metal) + 1
            f.write("  nspin = 2\n")
            f.write(f"  starting_magnetization({metal_index}) = 1.0\n")
        else:
            f.write("  nspin = 1\n")
        f.write("/\n")

        # ELECTRONS
        f.write("&ELECTRONS\n")
        f.write("  conv_thr = 1.0d-6\n")
        f.write("  mixing_beta = 0.25\n")
        f.write("  electron_maxstep = 250\n")
        f.write("/\n")

        # IONS
        f.write("&IONS\n")
        f.write("  ion_dynamics = 'bfgs'\n")
        f.write("  upscale = 50\n")
        f.write("/\n")

        # CELL_PARAMETERS
        f.write("\nCELL_PARAMETERS (angstrom)\n")
        for vec in cell:
            f.write(f"  {vec[0]:15.10f} {vec[1]:15.10f} {vec[2]:15.10f}\n")

        # ATOMIC_SPECIES
        f.write("\nATOMIC_SPECIES\n")
        for sym in unique_elements:
            if sym not in PSEUDOS:
                raise KeyError(f"No pseudopotential defined for element '{sym}'")
            f.write(f"  {sym:<4} 1.0  {PSEUDOS[sym]}\n")

        # ATOMIC_POSITIONS
        f.write("\nATOMIC_POSITIONS (angstrom)\n")
        n_fixed = 0
        for atom in atoms:
            if atom.z < Z_FIX:
                fix = "0 0 0"
                n_fixed += 1
            else:
                fix = "1 1 1"
            f.write(
                f"  {atom.symbol:<4} "
                f"{atom.position[0]:14.9f} {atom.position[1]:14.9f} {atom.position[2]:14.9f}"
                f"  {fix}\n"
            )

        print(f"Fixed {n_fixed} atoms (z < {Z_FIX:.1f} A), {len(atoms) - n_fixed} atoms free.")

        # K_POINTS
        f.write("\nK_POINTS (automatic)\n")
        f.write("  4 4 1 0 0 0\n")

        # HUBBARD: requested atomic (not ortho-atomic)
        f.write("\nHUBBARD {atomic}\n")
        f.write(f"U {metal}-{manifold} {u_value}\n")

    if "SLURM_JOB_ID" in os.environ:
        cmd = f"srun pw.x -nk 4 -in {pwi_file.name} > {pwo_file.name}"
    else:
        cmd = f"mpirun -np 64 pw.x -nk 4 -in {pwi_file.name} > {pwo_file.name}"

    if dry_run:
        print(f"[dry-run] Input written to {pwi_file}")
        print(f"[dry-run] Would run: {cmd}")
        return

    print(f"Launching for {tag} ...")
    print(f"Command: {cmd}\n")

    try:
        subprocess.run(cmd, shell=True, cwd=str(run_dir), check=True)
        print(f"Success: {tag}")
    except subprocess.CalledProcessError as exc:
        print(f"Failed: {tag} (exit code {exc.returncode})")
    except Exception as exc:
        print(f"Error: {tag}: {exc}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Relax slab+adsorbate combined structures with QE DFT-D3 + Hubbard U (atomic)."
    )
    parser.add_argument(
        "--structure",
        choices=SLAB_NAMES,
        default=None,
        help="Metal nitride slab (default: run all).",
    )
    parser.add_argument(
        "--adsorbate",
        choices=ADS_NAMES,
        default=None,
        help="Adsorbate species (default: run all).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Write input files only; do not launch pw.x.",
    )
    args = parser.parse_args()

    slabs = [args.structure] if args.structure else SLAB_NAMES
    adsorbates = [args.adsorbate] if args.adsorbate else ADS_NAMES

    for slab_name in slabs:
        for ads_name in adsorbates:
            generate_and_run(slab_name, ads_name, dry_run=args.dry_run)

    print("\nDone!")
