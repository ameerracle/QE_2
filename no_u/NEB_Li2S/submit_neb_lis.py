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
TIME = "43:50:00"

# --- Logging Settings ---
LOG_FILE = "completed_neb_lis_jobs.txt"

# --- Paths ---
WORK_DIR = Path(__file__).resolve().parent
INPUT_DIR = WORK_DIR / "input_files"
NEB_SCRIPT = WORK_DIR / "LiS_NEB.py"
NEB_EXE = "neb.x"


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
    """Find matching *_Li2S_1.xyz and *_LiS_2.xyz endpoint files in input_files."""
    pattern = re.compile(r"^(?P<metal>.+)_Li2S_1\.xyz$")
    pairs = []

    for initial_file in sorted(INPUT_DIR.glob("*_Li2S_1.xyz")):
        match = pattern.match(initial_file.name)
        if not match:
            continue

        metal = match.group("metal")
        if selected_metals and metal not in selected_metals:
            continue

        final_file = INPUT_DIR / f"{metal}_LiS_2.xyz"
        if not final_file.exists():
            print(f"Skipping {metal}: missing final file {final_file.name}")
            continue

        pairs.append((metal, initial_file, final_file))

    return pairs


def build_job_script(metal, initial_xyz, final_xyz):
    """Create a SLURM script that builds NEB input and launches neb.x."""
    job_name = f"nebLiS_{metal}"
    run_dir = WORK_DIR / f"{metal}_neb_run"
    neb_input = run_dir / "neb.in"
    neb_log = run_dir / "neb.out"
    prefix = f"{metal}_neb"

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
        f"mkdir -p {run_dir}\n",
        f"cd {WORK_DIR}\n\n",
    ]

    command = (
        f"python {NEB_SCRIPT.name} --metal {metal} --input-dir {INPUT_DIR.name} "
        f"--output {neb_input} --prefix {prefix} && "
        f"srun {NEB_EXE} -input {neb_input} > {neb_log} 2>&1 && "
        f"echo \"$(date): {metal} SUCCESS\" >> {LOG_FILE} || "
        f"echo \"$(date): {metal} FAILED/TIMEOUT\" >> {LOG_FILE}\n"
    )

    return job_name, header, command


def submit_jobs(pairs, dry_run=False):
    """Submit one SLURM job per matched metal pair."""
    if not pairs:
        print(f"No matching Li2S/LiS pairs found in {INPUT_DIR}")
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
        description="Submit one NEB SLURM job per matched Li2S(1) -> LiS(2) metal pair."
    )
    parser.add_argument(
        "--metals",
        nargs="*",
        default=None,
        help="Optional metal list, e.g. --metals ScN VN. If omitted, all pairs are submitted.",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Submit every matching metal pair found in input_files.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print jobs that would be submitted without calling sbatch.",
    )
    args = parser.parse_args()

    selected_metals = normalize_metal_names(args.metals)
    if args.all:
        selected_metals = []

    pairs = discover_pairs(selected_metals or None)
    submit_jobs(pairs, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
