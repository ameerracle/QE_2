#!/usr/bin/env python3
import os
import subprocess
from pathlib import Path


# --- Hardware & Account Settings ---
ACCOUNT = "rrg-peslherb-ac"
TASKS_PER_NODE = 48
MEM = "96G"
TIME = "12:00:00"

# --- Target Systems ---
METALS = ["VN", "TiN", "ScN", "NbN", "ZrN"]
ADSORBATES = ["Li2S4", "S8", "Li2S8"]

# --- Logging Settings ---
LOG_FILE = "completed_pdos_jobs.txt"


def submit_individual_jobs() -> None:
    """Submit one SLURM job per metal-adsorbate combination."""
    script_dir = Path(__file__).resolve().parent

    for ads in ADSORBATES:
        for metal in METALS:
            job_tag = f"{metal}_{ads}"
            job_name = f"pdos144_{job_tag}"
            submit_script = script_dir / f"temp_submit_{job_tag}.sh"

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

            # run_pdos_combi.py has -nk hard-fixed to 2.
            command = (
                f"python run_pdos_combi.py --slab {metal} --ads {ads} --mode auto --run --np {TASKS_PER_NODE} && "
                f"echo \"$(date): {job_tag} SUCCESS\" >> {LOG_FILE} || "
                f"echo \"$(date): {job_tag} FAILED/TIMEOUT\" >> {LOG_FILE}\n"
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
