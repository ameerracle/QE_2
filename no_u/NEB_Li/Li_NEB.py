#!/usr/bin/env python3
"""
Generate a Quantum ESPRESSO `neb.x` input file with unique subdirectories
to prevent HDF5 file-locking errors during parallel NEB runs.
"""

from __future__ import annotations
import argparse
from pathlib import Path
from ase.io import read

WORK_DIR = Path(__file__).resolve().parent
PSEUDO_DIR = "/scratch/anizami/QE_2/USPP/"

# Global QE Settings
ECUTWFC = 45.0
ECUTRHO = 450.0
Z_FIX = 13.0
NUM_IMAGES = 5
OPT_SCHEME = "broyden"
CI_SCHEME = "auto"
NSTEP_PATH = 120

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

def load_qe_output(path: Path):
    return read(str(path), format="espresso-out", index=-1)

def get_metal_label(atoms):
    """Extracts the metal symbol (Ti, V, Sc, etc.) from the structure."""
    metals = {"Ti", "V", "Sc", "Nb", "Zr"}
    symbols = set(atoms.get_chemical_symbols())
    found = metals.intersection(symbols)
    return list(found)[0] if found else "Unknown"

def atom_fix_flags(atoms):
    flags = []
    for atom in atoms:
        if atom.position[2] < Z_FIX:
            flags.append("0 0 0")
        else:
            flags.append("1 1 1")
    return flags

def build_neb_input(initial_atoms, final_atoms, prefix, outdir, pseudo_dir):
    nat = len(initial_atoms)
    ntyp = len(set(initial_atoms.get_chemical_symbols()))
    flags = atom_fix_flags(initial_atoms)

    lines = [
        "BEGIN",
        "BEGIN_PATH_INPUT",
        "&PATH",
        "  restart_mode = 'restart',",
        "  string_method = 'neb',",
        f"  nstep_path = {NSTEP_PATH},",
        f"  opt_scheme = '{OPT_SCHEME}',",
        f"  num_of_images = {NUM_IMAGES},",
        f"  CI_scheme = '{CI_SCHEME}',",
        "/",
        "END_PATH_INPUT",
        "BEGIN_ENGINE_INPUT",
        "&CONTROL",
        "  calculation = 'scf',",
        f"  prefix = '{prefix}',",
        f"  outdir = '{outdir}',",
        f"  pseudo_dir = '{pseudo_dir}',",
        "  disk_io = 'low',",
        "  verbosity = 'low',",
        "/",
        "&SYSTEM",
        "  ibrav = 0,",
        f"  nat = {nat},",
        f"  ntyp = {ntyp},",
        f"  ecutwfc = {ECUTWFC},",
        f"  ecutrho = {ECUTRHO},",
        "  occupations = 'smearing',",
        "  smearing = 'cold',",
        "  degauss = 0.025,",
        "  vdw_corr = 'dft-d3',",
        "  dftd3_version = 4,",
        "/",
        "&ELECTRONS",
        "  conv_thr = 1.0d-6,",
        "  mixing_beta = 0.25,",
        "  electron_maxstep = 250,",
        "/",
    ]
    
    # Species lines
    lines.append("ATOMIC_SPECIES")
    seen = []
    for s in initial_atoms.get_chemical_symbols():
        if s not in seen:
            lines.append(f"  {s:<4} 1.0  {PSEUDOS[s]}")
            seen.append(s)

    # Cell
    lines.append("CELL_PARAMETERS (angstrom)")
    for vec in initial_atoms.cell:
        lines.append(f"  {vec[0]:15.10f} {vec[1]:15.10f} {vec[2]:15.10f}")

    lines.append("BEGIN_POSITIONS")
    lines.append("FIRST_IMAGE")
    lines.append("ATOMIC_POSITIONS (angstrom)")
    for atom, fix in zip(initial_atoms, flags):
        lines.append(f"  {atom.symbol:<4} {atom.position[0]:14.9f} {atom.position[1]:14.9f} {atom.position[2]:14.9f}  {fix}")
    
    lines.append("LAST_IMAGE")
    lines.append("ATOMIC_POSITIONS (angstrom)")
    for atom, fix in zip(final_atoms, flags):
        lines.append(f"  {atom.symbol:<4} {atom.position[0]:14.9f} {atom.position[1]:14.9f} {atom.position[2]:14.9f}  {fix}")
    
    lines.append("END_POSITIONS")
    lines.append("K_POINTS (automatic)")
    lines.append("  4 4 1 0 0 0")
    lines.append("END_ENGINE_INPUT")
    lines.append("END")

    return "\n".join(lines)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--initial", required=True, type=Path)
    parser.add_argument("--final", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path, help="Explicit path to write the neb.in file")
    parser.add_argument("--prefix", required=True, type=str, help="Prefix for QE")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    init_atoms = load_qe_output(args.initial)
    final_atoms = load_qe_output(args.final)
    
    # outdir is forced to be the 'tmp' folder next to wherever the output file goes
    outdir_path = args.output.parent / "tmp"

    neb_text = build_neb_input(init_atoms, final_atoms, args.prefix, str(outdir_path), PSEUDO_DIR)

    if args.dry_run:
        print(f"\n--- DRY RUN: {args.output.name} ---\n")
        print(neb_text)
        return

    # Create directories and write file strictly based on --output
    args.output.parent.mkdir(parents=True, exist_ok=True)
    outdir_path.mkdir(exist_ok=True)
    args.output.write_text(neb_text)
    
    print(f"Input written to: {args.output}")

if __name__ == "__main__":
    main()