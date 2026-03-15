import subprocess
from pathlib import Path

# --- Settings ---
METALS = ["TiN", "VN",
"ScN", "NbN", "ZrN"
]
ACCOUNT = "def-peslherb"
NTASKS = 64
MEM = "135G"
TIME = "23:00:00"

def create_and_submit():
    for metal in METALS:
        job_name = f"relax144_{metal}"
        sh_filename = f"submit_{metal}.sh"
        
        # Build the SLURM script string
        script_content = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output={job_name}_%j.out
#SBATCH --error={job_name}_%j.err
#SBATCH --account={ACCOUNT}
#SBATCH --ntasks-per-node={NTASKS}
#SBATCH --mem={MEM}
#SBATCH --time={TIME}

set -euo pipefail

# Load modules and environment
module load quantumespresso/7.5 || true
source ~/ase/bin/activate
export OMP_NUM_THREADS=1

# Run the python script for this specific metal
python slab_relax.py --structure {metal}
"""
        
        # Write the .sh file
        with open(sh_filename, "w") as f:
            f.write(script_content)
        
        # Submit the .sh file via sbatch
        print(f"Submitting job for {metal}...")
        subprocess.run(["sbatch", sh_filename])

if __name__ == "__main__":
    create_and_submit()