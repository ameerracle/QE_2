#!/usr/bin/env python3
import subprocess
from pathlib import Path

# --- Hardware & Account Settings ---
ACCOUNT = "def-peslherb"
TASKS_PER_NODE = 40
MEM = "200G"
TIME = "12:00:00"

LOG_FILE = "completed_slab_pdos_jobs.txt"
SCRIPT_DIR = Path(__file__).resolve().parent
SLAB_SCRIPT = SCRIPT_DIR / "slab_pdos_bader.py"


def submit_slab_job():
    job_tag = "VN_slab"
    job_name = f"slab_pdos_{job_tag}"
    submit_script = SCRIPT_DIR / "temp_submit_VN_slab.sh"

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
        f"cd {SCRIPT_DIR}\n\n",
    ]

    command = (
        f"python {SLAB_SCRIPT} --structure VN --run && "
        f"echo \"$(date): {job_tag} SUCCESS\" >> {LOG_FILE} || "
        f"echo \"$(date): {job_tag} FAILED/TIMEOUT\" >> {LOG_FILE}\n"
    )

    with open(submit_script, "w", encoding="utf-8") as f:
        f.writelines(header)
        f.write(command)

    print(f"Submitting job: {job_name}")
    try:
        result = subprocess.run(["sbatch", str(submit_script)], capture_output=True, text=True, check=True)
        print(f"  -> {result.stdout.strip()}")
    except subprocess.CalledProcessError as e:
        print(f"  !! Error submitting {job_name}: {e.stderr}")
    finally:
        if submit_script.exists():
            submit_script.unlink()


if __name__ == "__main__":
    log_path = SCRIPT_DIR / LOG_FILE
    if not log_path.exists():
        log_path.touch()

    submit_slab_job()
