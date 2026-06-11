#!/usr/bin/env python3
"""
Generate sbatch scripts for Bader charge analysis on combi systems.
Creates 5 separate job scripts, one per metal, processing all adsorbates.
"""

import os
import subprocess
from pathlib import Path

ACCOUNT = "def-peslherb"
TASKS_PER_NODE = 4
MEM = "16G"
TIME = "1:30:00"

METALS = ["ScN", "TiN", "VN", "NbN", "ZrN"]
ADSORBATES = ["Li2S4", "Li2S8", "S8"]
BASE_DIR = Path("/lustre10/scratch/anizami/QE_2")
BADER_SCRIPT = BASE_DIR / "bader_combi.py"
OUTPUT_DIR = BASE_DIR / "submit_scripts"

def write_sbatch(metal):
    output_dir = OUTPUT_DIR / metal
    output_dir.mkdir(parents=True, exist_ok=True)
    
    sbatch_path = output_dir / f"bader_{metal}.sh"
    
    adsorbate_list = " ".join(ADSORBATES)
    
    content = f"""#!/bin/bash
#SBATCH --job-name=bader_{metal}
#SBATCH --output={output_dir}/bader_{metal}.out
#SBATCH --error={output_dir}/bader_{metal}.err
#SBATCH --account={ACCOUNT}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node={TASKS_PER_NODE}
#SBATCH --mem={MEM}
#SBATCH --time={TIME}

set -euo pipefail
module load quantumespresso/7.5 || true
source ~/ase/bin/activate
export OMP_NUM_THREADS=1
cd {BASE_DIR}

python {BADER_SCRIPT} --metals {metal} --adsorbates {adsorbate_list}

echo "Done for {metal}"
"""
    
    with open(sbatch_path, "w") as f:
        f.write(content)
    
    print(f"Created: {sbatch_path}")
    
    print(f"Submitting: {sbatch_path}")
    try:
        result = subprocess.run(
            ["sbatch", str(sbatch_path)],
            capture_output=True,
            text=True,
            check=True,
        )
        print(f"  -> {result.stdout.strip()}")
    except subprocess.CalledProcessError as exc:
        print(f"  !! Error submitting: {exc.stderr}")


if __name__ == "__main__":
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    print("Generating and submitting sbatch scripts for Bader analysis (combi only)...\n")
    
    for metal in METALS:
        write_sbatch(metal)
    
    print(f"\nAll scripts generated and submitted from {OUTPUT_DIR}")
