#!/usr/bin/env python3
import argparse
import os
import re
import subprocess
from pathlib import Path


# --- Hardware & Account Settings ---
ACCOUNT = "def-peslherb"
TASKS_PER_NODE = 64
MEM = "145G"
TIME = "23:50:00"

# --- Logging Settings ---
LOG_FILE = "completed_neb_jobs.txt"

# --- Paths ---
WORK_DIR = Path(__file__).resolve().parent
NEB_SCRIPT = WORK_DIR / "Li_NEB.py"
NEB_EXE = "neb.x"

INITIAL_SUFFIX = "combi1_144_relax"
FINAL_SUFFIX = "combi2_144_relax"


def normalize_metal_names(raw_names):
	"""Normalize argparse metal names and allow comma-separated input."""
	names = []
	for value in raw_names or []:
		for item in value.split(","):
			item = item.strip()
			if item:
				names.append(item)
	return names


def discover_pairs(selected_metals=None):
	"""Find matching combi1/combi2 relaxation folders for the requested metals."""
	pattern = re.compile(r"^(?P<metal>.+)_Li_" + re.escape(INITIAL_SUFFIX) + r"$")
	pairs = []

	for initial_dir in sorted(WORK_DIR.glob(f"*_Li_{INITIAL_SUFFIX}")):
		match = pattern.match(initial_dir.name)
		if not match:
			continue

		metal = match.group("metal")
		if selected_metals and metal not in selected_metals:
			continue

		final_dir = WORK_DIR / f"{metal}_Li_{FINAL_SUFFIX}"
		if not final_dir.exists():
			print(f"Skipping {metal}: missing final folder {final_dir.name}")
			continue

		initial_out = initial_dir / f"{metal}_Li_{INITIAL_SUFFIX.replace('_relax', '')}.out"
		final_out = final_dir / f"{metal}_Li_{FINAL_SUFFIX.replace('_relax', '')}.out"

		if not initial_out.exists():
			print(f"Skipping {metal}: missing initial output {initial_out.name}")
			continue
		if not final_out.exists():
			print(f"Skipping {metal}: missing final output {final_out.name}")
			continue

		pairs.append((metal, initial_out, final_out))

	return pairs


def build_job_script(metal, initial_out, final_out):
	"""Create a SLURM script that generates the NEB input and launches neb.x."""
	job_name = f"nebLi_{metal}"
	neb_input = WORK_DIR / f"{metal}_neb.in"
	neb_log = WORK_DIR / f"{metal}_neb.out"
	initial_arg = initial_out.relative_to(WORK_DIR).as_posix()
	final_arg = final_out.relative_to(WORK_DIR).as_posix()

	header = [
		"#!/bin/bash\n",
		f"#SBATCH --job-name={job_name}\n",
		f"#SBATCH --output={job_name}_%j.out\n",
		f"#SBATCH --error={job_name}_%j.err\n",
		f"#SBATCH --account={ACCOUNT}\n",
		"#SBATCH --nodes=1\n",
		f"#SBATCH --ntasks-per-node={TASKS_PER_NODE}\n",
		f"#SBATCH --mem={MEM}\n",
		f"#SBATCH --time={TIME}\n",
		"\n",
		"set -euo pipefail\n",
		"module load quantumespresso\n",
		"source ~/ase/bin/activate\n",
		"export OMP_NUM_THREADS=1\n",
		f"cd {WORK_DIR}\n\n",
	]

	command = (
		f"python {NEB_SCRIPT.name} --initial {initial_arg} --final {final_arg} --output {neb_input.name} && "
		f"srun {NEB_EXE} -input {neb_input.name} > {neb_log.name} 2>&1 && "
		f"echo \"$(date): {metal} SUCCESS\" >> {LOG_FILE} || "
		f"echo \"$(date): {metal} FAILED/TIMEOUT\" >> {LOG_FILE}\n"
	)

	return job_name, header, command


def submit_jobs(pairs, dry_run=False):
	"""Submit one SLURM job per matched metal pair."""
	if not pairs:
		print(f"No matching combi1/combi2 pairs found in {WORK_DIR}")
		return

	for metal, initial_out, final_out in pairs:
		job_name, header, command = build_job_script(metal, initial_out, final_out)
		submit_script = WORK_DIR / f"temp_submit_{metal}.sh"

		with open(submit_script, "w", encoding="utf-8") as handle:
			handle.writelines(header)
			handle.write(command)

		print(f"Preparing NEB job: {job_name}")
		print(f"  initial: {initial_out.name}")
		print(f"  final:   {final_out.name}")

		if dry_run:
			print(f"  [dry-run] would submit {submit_script.name}")
			continue

		try:
			result = subprocess.run(
				["sbatch", str(submit_script)],
				capture_output=True,
				text=True,
				check=True,
			)
			print(f"  -> {result.stdout.strip()}")
		except subprocess.CalledProcessError as exc:
			print(f"  !! Error submitting {job_name}: {exc.stderr}")
		finally:
			if submit_script.exists():
				os.remove(submit_script)


def main():
	parser = argparse.ArgumentParser(
		description="Submit one NEB SLURM job per matched Li metal pair (combi1 -> combi2)."
	)
	parser.add_argument(
		"--metals",
		nargs="*",
		default=None,
		help="Optional list of metal names to submit, e.g. --metals ScN VN. If omitted, all matching pairs are submitted.",
	)
	parser.add_argument(
		"--all",
		action="store_true",
		help="Submit every matching metal pair found in this folder.",
	)
	parser.add_argument(
		"--dry-run",
		action="store_true",
		help="Print the jobs that would be submitted, but do not call sbatch.",
	)
	args = parser.parse_args()

	selected_metals = normalize_metal_names(args.metals)
	if args.all:
		selected_metals = []

	pairs = discover_pairs(selected_metals or None)
	submit_jobs(pairs, dry_run=args.dry_run)


if __name__ == "__main__":
	main()