#!/usr/bin/env python3
import argparse
import os
import re
import subprocess
from pathlib import Path


ACCOUNT = "def-peslherb"
TASKS_PER_NODE = 64
MEM = "145G"
TIME = "23:50:00"

LOG_FILE = "completed_neb_jobs.txt"

WORK_DIR = Path(__file__).resolve().parent
NEB_SCRIPT = WORK_DIR / "Li_NEB.py"
NEB_EXE = "neb.x"

INITIAL_PATTERN = "*_combi1u_144_final.xyz"
FINAL_PATTERN = "*_combi2u_144_final.xyz"


def normalize_metal_names(raw_names):
	names = []
	for value in raw_names or []:
		for item in value.split(","):
			item = item.strip()
			if item:
				names.append(item)
	return names


def discover_pairs(selected_metals=None):
	pattern = re.compile(r"^(?P<metal>.+)_Li_combi1u_144_final\.xyz$")
	pairs = []

	for initial_xyz in sorted(WORK_DIR.glob(INITIAL_PATTERN)):
		match = pattern.match(initial_xyz.name)
		if not match:
			continue

		metal = match.group("metal")
		if selected_metals and metal not in selected_metals:
			continue

		final_xyz = WORK_DIR / f"{metal}_Li_combi2u_144_final.xyz"
		if not final_xyz.exists():
			print(f"Skipping {metal}: missing final file {final_xyz.name}")
			continue

		pairs.append((metal, initial_xyz, final_xyz))

	return pairs


def build_job_script(metal, initial_xyz, final_xyz):
	job_name = f"nebLi_{metal}"
	neb_input = WORK_DIR / f"{metal}_neb.in"
	neb_log = WORK_DIR / f"{metal}_neb.out"
	initial_arg = initial_xyz.name
	final_arg = final_xyz.name

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
		"module load quantumespresso/7.5 || true\n",
		"source ~/ase/bin/activate\n",
		"export OMP_NUM_THREADS=1\n",
		f"cd {WORK_DIR}\n\n",
	]

	command = (
		f"python {NEB_SCRIPT.name} --initial {initial_arg} --final {final_arg} --output {neb_input.name} --prefix {metal}_Li && "
		f"srun {NEB_EXE} -input {neb_input.name} > {neb_log.name} 2>&1 && "
		f"echo \"$(date): {metal} SUCCESS\" >> {LOG_FILE} || "
		f"echo \"$(date): {metal} FAILED/TIMEOUT\" >> {LOG_FILE}\n"
	)

	return job_name, header, command


def submit_jobs(pairs, dry_run=False):
	if not pairs:
		print(f"No matching combi1u/combi2u pairs found in {WORK_DIR}")
		return

	for metal, initial_xyz, final_xyz in pairs:
		job_name, header, command = build_job_script(metal, initial_xyz, final_xyz)
		submit_script = WORK_DIR / f"temp_submit_{metal}.sh"

		with open(submit_script, "w", encoding="utf-8") as handle:
			handle.writelines(header)
			handle.write(command)

		print(f"Preparing NEB job: {job_name}")
		print(f"  initial: {initial_xyz.name}")
		print(f"  final:   {final_xyz.name}")

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
		description="Submit one NEB SLURM job per matched Li metal pair (combi1u -> combi2u) with Hubbard +U."
	)
	parser.add_argument(
		"--metals",
		nargs="*",
		default=None,
		help="Optional list of metal names to submit, e.g. --metals VN. If omitted, all matching pairs are submitted.",
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
