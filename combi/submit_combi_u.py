#!/usr/bin/env python3
import os
import subprocess

# --- Hardware & Account Settings ---
ACCOUNT = "def-peslherb"
TASKS_PER_NODE = 64
MEM = "145G"
TIME = "43:50:00"

# --- Target Systems ---
METALS = ["VN", "TiN"]
ADSORBATES = ["S8", "Li2S8", "Li2S6", "Li2S4", "Li2S2", "Li2S"]

# --- Logging Settings ---
LOG_FILE = "completed_jobs_u.txt"


def submit_individual_jobs():
    """Submits a unique SLURM job for every single Metal-Adsorbate combination."""

    for ads in ADSORBATES:
        for metal in METALS:
            job_tag = f"{metal}_{ads}"
            job_name = f"relax144u_{job_tag}"
            submit_script = f"temp_submit_u_{job_tag}.sh"

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
                "export OMP_NUM_THREADS=1\n\n",
            ]

            command = (
                f"python combi_relax_u.py --structure {metal} --adsorbate {ads} && "
                f"echo \"$(date): {job_tag} SUCCESS\" >> {LOG_FILE} || "
                f"echo \"$(date): {job_tag} FAILED/TIMEOUT\" >> {LOG_FILE}\n"
            )

            with open(submit_script, "w", encoding="utf-8") as f:
                f.writelines(header)
                f.write(command)

            print(f"Submitting individual job: {job_name}")
            try:
                result = subprocess.run(["sbatch", submit_script], capture_output=True, text=True, check=True)
                print(f"  -> {result.stdout.strip()}")
            except subprocess.CalledProcessError as e:
                print(f"  !! Error submitting {job_name}: {e.stderr}")
            finally:
                if os.path.exists(submit_script):
                    os.remove(submit_script)


if __name__ == "__main__":
    submit_individual_jobs()
