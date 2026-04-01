#!/usr/bin/env python3
import subprocess
import os

ACCOUNT = "def-peslherb"
TASKS_PER_NODE = 64
MEM = "140G"
TIME = "40:50:00"

# Only the missing jobs -- S8 running or done, excluded here
TARGETS = [
    ("ScN", "Li2S6"),
]

def submit_restarts():
    for metal, ads in TARGETS:
        job_tag = f"{metal}_{ads}"
        job_name = f"rst144_{job_tag}"
        submit_script = f"temp_restart_{job_tag}.sh"

        header = [
            f"#!/bin/bash\n",
            f"#SBATCH --job-name={job_name}\n",
            f"#SBATCH --output={job_name}_%j.out\n",
            f"#SBATCH --error={job_name}_%j.err\n",
            f"#SBATCH --account={ACCOUNT}\n",
            f"#SBATCH --nodes=1\n",
            f"#SBATCH --ntasks-per-node={TASKS_PER_NODE}\n",
            f"#SBATCH --mem={MEM}\n",
            f"#SBATCH --time={TIME}\n",
            "\n",
            "set -euo pipefail\n",
            "module load quantumespresso/7.5 || true\n",
            "source ~/ase/bin/activate\n",
            "export OMP_NUM_THREADS=1\n\n",
        ]
        command = f"python combi_restart.py --structure {metal} --adsorbate {ads}\n"

        with open(submit_script, "w") as f:
            f.writelines(header)
            f.write(command)

        print(f"Submitting: {job_name}")
        try:
            result = subprocess.run(["sbatch", submit_script], capture_output=True, text=True, check=True)
            print(f"  -> {result.stdout.strip()}")
        except subprocess.CalledProcessError as e:
            print(f"  !! Error submitting {job_name}: {e.stderr}")
        finally:
            if os.path.exists(submit_script):
                os.remove(submit_script)

if __name__ == "__main__":
    submit_restarts()