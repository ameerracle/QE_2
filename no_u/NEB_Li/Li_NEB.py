#!/usr/bin/env python3
"""Generate a Quantum ESPRESSO `neb.x` input file from two relaxed QE outputs.

The script reads the initial and final structures with ASE from QE `.out` files,
then writes a complete NEB input file with:
- `num_of_images = 5`
- `opt_scheme = 'broyden'`
- `CI_scheme = 'auto'`

Use `--dry-run` to print the exact input text instead of writing it.

Example:
    python Li_NEB.py --initial ScN_Li_combi2_144_relax/ScN_Li_combi2_144.out \
                     --final VN_Li_combi1_144_relax/VN_Li_combi1_144.out \
                     --output neb.in
"""

from __future__ import annotations

import argparse
from pathlib import Path

from ase.io import read


WORK_DIR = Path(__file__).resolve().parent
PSEUDO_DIR = "/scratch/anizami/QE_2/USPP/"

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
	"""Read the final structure from a Quantum ESPRESSO output file."""
	return read(str(path), format="espresso-out", index=-1)


def get_unique_elements(atoms):
	"""Return element symbols preserving their first occurrence order."""
	seen = set()
	unique = []
	for symbol in atoms.get_chemical_symbols():
		if symbol not in seen:
			seen.add(symbol)
			unique.append(symbol)
	return unique


def atom_fix_flags(atoms):
	"""Build QE positional constraints using the same z-threshold logic as relax."""
	flags = []
	for atom in atoms:
		if atom.position[2] < Z_FIX:
			flags.append("0 0 0")
		else:
			flags.append("1 1 1")
	return flags


def compare_structures(initial_atoms, final_atoms):
	"""Make sure the two endpoint images are compatible for NEB."""
	if len(initial_atoms) != len(final_atoms):
		raise ValueError(
			f"Initial and final structures have different nat: {len(initial_atoms)} vs {len(final_atoms)}"
		)

	if initial_atoms.get_chemical_symbols() != final_atoms.get_chemical_symbols():
		raise ValueError(
			"Initial and final structures do not have the same atom order or species. "
			"NEB requires matching endpoints."
		)

	initial_cell = initial_atoms.cell.array
	final_cell = final_atoms.cell.array
	max_diff = 0.0
	for row_a, row_b in zip(initial_cell, final_cell):
		for value_a, value_b in zip(row_a, row_b):
			max_diff = max(max_diff, abs(value_a - value_b))

	if max_diff > 1.0e-6:
		raise ValueError(
			f"Initial and final cells differ by up to {max_diff:.3e}; NEB expects the same cell."
		)


def format_species_lines(atoms):
	lines = ["ATOMIC_SPECIES"]
	for symbol in get_unique_elements(atoms):
		if symbol not in PSEUDOS:
			raise KeyError(f"No pseudopotential defined for element '{symbol}'")
		lines.append(f"  {symbol:<4} 1.0  {PSEUDOS[symbol]}")
	return lines


def format_cell_lines(atoms):
	lines = ["CELL_PARAMETERS (angstrom)"]
	for vector in atoms.cell:
		lines.append(f"  {vector[0]:15.10f} {vector[1]:15.10f} {vector[2]:15.10f}")
	return lines


def format_atomic_positions(atoms, flags):
	lines = ["ATOMIC_POSITIONS (angstrom)"]
	for atom, fix in zip(atoms, flags):
		lines.append(
			f"  {atom.symbol:<4} {atom.position[0]:14.9f} {atom.position[1]:14.9f} "
			f"{atom.position[2]:14.9f}  {fix}"
		)
	return lines


def build_neb_input(initial_atoms, final_atoms, prefix, outdir, pseudo_dir):
	compare_structures(initial_atoms, final_atoms)
	unique_elements = get_unique_elements(initial_atoms)
	nat = len(initial_atoms)
	nyp = len(unique_elements)
	flags = atom_fix_flags(initial_atoms)

	lines = [
		"BEGIN",
		"BEGIN_PATH_INPUT",
		"&PATH",
		"  restart_mode = 'from_scratch',",
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
		f"  ntyp = {nyp},",
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

	lines.extend(format_species_lines(initial_atoms))
	lines.extend(format_cell_lines(initial_atoms))
	lines.append("BEGIN_POSITIONS")
	lines.append("FIRST_IMAGE")
	lines.extend(format_atomic_positions(initial_atoms, flags))
	lines.append("LAST_IMAGE")
	lines.extend(format_atomic_positions(final_atoms, flags))
	lines.append("END_POSITIONS")
	lines.append("K_POINTS (automatic)")
	lines.append("  4 4 1 0 0 0")
	lines.append("END_ENGINE_INPUT")
	lines.append("END")

	return "\n".join(lines) + "\n"


def build_default_prefix(initial_path: Path, final_path: Path) -> str:
	"""Create a readable default prefix if the user does not supply one."""
	initial_name = initial_path.stem.replace("_relax", "")
	final_name = final_path.stem.replace("_relax", "")
	if initial_name == final_name:
		return initial_name
	return "neb_li"


def main():
	parser = argparse.ArgumentParser(
		description="Generate a Quantum ESPRESSO neb.x input file from two QE .out files."
	)
	parser.add_argument(
		"--initial",
		required=True,
		type=Path,
		help="QE output file for the initial image (typically in a *_relax folder).",
	)
	parser.add_argument(
		"--final",
		required=True,
		type=Path,
		help="QE output file for the final image (typically in a *_relax folder).",
	)
	parser.add_argument(
		"--output",
		type=Path,
		default=WORK_DIR / "neb.in",
		help="Output neb.x input file to write.",
	)
	parser.add_argument(
		"--prefix",
		default=None,
		help="PW prefix written in the CONTROL namelist. Defaults to a simple name derived from the inputs.",
	)
	parser.add_argument(
		"--outdir",
		default="./tmp/",
		help="PW outdir written in the CONTROL namelist.",
	)
	parser.add_argument(
		"--pseudo-dir",
		default=PSEUDO_DIR,
		help="QE pseudopotential directory written in the CONTROL namelist.",
	)
	parser.add_argument(
		"--dry-run",
		action="store_true",
		help="Print the exact neb.x input file instead of writing it.",
	)
	args = parser.parse_args()

	initial_path = args.initial if args.initial.is_absolute() else (Path.cwd() / args.initial)
	final_path = args.final if args.final.is_absolute() else (Path.cwd() / args.final)

	if not initial_path.exists():
		raise FileNotFoundError(f"Initial QE output not found: {initial_path}")
	if not final_path.exists():
		raise FileNotFoundError(f"Final QE output not found: {final_path}")

	initial_atoms = load_qe_output(initial_path)
	final_atoms = load_qe_output(final_path)
	prefix = args.prefix or build_default_prefix(initial_path, final_path)
	neb_text = build_neb_input(initial_atoms, final_atoms, prefix, args.outdir, args.pseudo_dir)

	print(f"Initial: {initial_path}")
	print(f"Final:   {final_path}")
	print(f"nat = {len(initial_atoms)}")
	print(f"num_of_images = {NUM_IMAGES}, opt_scheme = '{OPT_SCHEME}', CI_scheme = '{CI_SCHEME}'")

	if args.dry_run:
		print("\n--- neb.x input ---\n")
		print(neb_text, end="")
		return

	args.output.parent.mkdir(parents=True, exist_ok=True)
	args.output.write_text(neb_text, encoding="utf-8")
	print(f"Wrote NEB input to {args.output}")


if __name__ == "__main__":
	main()
