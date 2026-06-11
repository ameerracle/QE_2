#!/usr/bin/env python3
"""
Generate and submit sbatch scripts for adsorbate-only PDOS (Li/S).
One SLURM job per combi system.
"""

import os
import subprocess
from pathlib import Path

ACCOUNT = "def-peslherb"
TASKS_PER_NODE = 8
MEM = "16G"
TIME = "1:00:00"

SYSTEMS = [
    "TiN_Li2S4", "TiN_Li2S8", "TiN_S8",
    "ScN_Li2S4", "ScN_Li2S8", "ScN_S8",
    "VN_Li2S4", "VN_Li2S8", "VN_S8",
    "NbN_Li2S4", "NbN_Li2S8", "NbN_S8",
    "ZrN_Li2S4", "ZrN_Li2S8", "ZrN_S8"
]

LOG_FILE = "completed_adsorbate_pdos_jobs.txt"


def submit_individual_jobs() -> None:
    script_dir = Path(__file__).resolve().parent

    for system in SYSTEMS:
        job_name = f"ads_pdos_{system}"
        submit_script = script_dir / f"temp_submit_{system}.sh"

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
            f"cd {script_dir}\n\n",
        ]

        command = (
            f"python adsorbate_pdos.py --system {system} --run --np {TASKS_PER_NODE} && "
            f"echo \"$(date): {system} SUCCESS\" >> {LOG_FILE} || "
            f"echo \"$(date): {system} FAILED/TIMEOUT\" >> {LOG_FILE}\n"
        )

        with open(submit_script, "w", encoding="utf-8") as f:
            f.writelines(header)
            f.write(command)

        print(f"Submitting: {job_name}")
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
