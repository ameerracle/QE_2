#!/usr/bin/env python3
import os
import subprocess
from pathlib import Path

# --- Hardware & Account Settings ---
ACCOUNT = "def-peslherb"
TASKS_PER_NODE = 64
MEM = "135G"
TIME = "23:50:00"

# --- Logging Settings ---
LOG_FILE = "completed_jobs.txt"

# --- Paths ---
WORK_DIR = Path(__file__).resolve().parent
RELAX_SCRIPT = WORK_DIR / "Li_relax.py"


def submit_individual_jobs():
    """Submit one SLURM job per .xyz structure in this folder."""
    xyz_files = sorted(WORK_DIR.glob("*.xyz"))

    if not xyz_files:
        print(f"No .xyz files found in {WORK_DIR}")
        return

    if not RELAX_SCRIPT.exists():
        print(f"Missing relax script: {RELAX_SCRIPT}")
        return

    for xyz_path in xyz_files:
        xyz_name = xyz_path.name
        base = xyz_path.stem

        job_name = f"relaxLi_{base}"
        submit_script = WORK_DIR / f"temp_submit_{base}.sh"

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
            f"python {RELAX_SCRIPT.name} --xyz {xyz_name} && "
            f"echo \"$(date): {xyz_name} SUCCESS\" >> {LOG_FILE} || "
            f"echo \"$(date): {xyz_name} FAILED/TIMEOUT\" >> {LOG_FILE}\n"
        )

        with open(submit_script, "w", encoding="utf-8") as f:
            f.writelines(header)
            f.write(command)

        print(f"Submitting individual job: {job_name}")
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
            if os.path.exists(submit_script):
                os.remove(submit_script)


if __name__ == "__main__":
    submit_individual_jobs()
