#!/usr/bin/env python3
import subprocess
import os

# --- Hardware & Account Settings ---
ACCOUNT = "rrg-peslherb-ac"
TASKS_PER_NODE = 64
MEM = "135G" 
TIME = "23:50:00"

# --- Target Systems ---
METALS = ["VN", "TiN", "ScN", "NbN", "ZrN"]
ADSORBATES = ["S8", "Li2S8", "Li2S6", "Li2S4", "Li2S2", "Li2S"]

# --- Logging Settings ---
LOG_FILE = "completed_jobs.txt"

def submit_individual_jobs():
    """Submits a unique SLURM job for every single Metal-Adsorbate combination."""
    
    for ads in ADSORBATES:
        for metal in METALS:
            # Unique identifiers for this specific pair
            job_tag = f"{metal}_{ads}"
            job_name = f"relax144_{job_tag}"
            submit_script = f"temp_submit_{job_tag}.sh"
            
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
                "export OMP_NUM_THREADS=1\n\n"
            ]

            # Standardized 12-space indentation for the command block
            command = (
                f"python combi_relax.py --structure {metal} --adsorbate {ads} && "
                f"echo \"$(date): {job_tag} SUCCESS\" >> {LOG_FILE} || "
                f"echo \"$(date): {job_tag} FAILED/TIMEOUT\" >> {LOG_FILE}\n"
            )

            # Write the individual temporary script
            with open(submit_script, "w") as f:
                f.writelines(header)
                f.write(command)

            # Submit to the queue
            print(f"Submitting individual job: {job_name}")
            try:
                result = subprocess.run(["sbatch", submit_script], capture_output=True, text=True, check=True)
                print(f"  -> {result.stdout.strip()}")
            except subprocess.CalledProcessError as e:
                print(f"  !! Error submitting {job_name}: {e.stderr}")
            finally:
                # Cleanup temporary submission file
                if os.path.exists(submit_script):
                    os.remove(submit_script)

if __name__ == "__main__":
    submit_individual_jobs()