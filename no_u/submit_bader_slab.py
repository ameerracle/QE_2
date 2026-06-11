#!/usr/bin/env python3
"""
Generate and submit sbatch script for Bader slab analysis.
Single job to process all 5 metal slabs.
"""

import subprocess
from pathlib import Path

ACCOUNT = "def-peslherb"
TASKS_PER_NODE = 4
MEM = "16G"
TIME = "1:30:00"

METALS = ["ScN", "TiN", "VN", "NbN", "ZrN"]
BASE_DIR = Path("/lustre10/scratch/anizami/QE_2")
BADER_SCRIPT = BASE_DIR / "bader_slab.py"
OUTPUT_DIR = BASE_DIR / "submit_scripts" / "slab"


def write_sbatch():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    sbatch_path = OUTPUT_DIR / "bader_slab.sh"
    
    metals_list = " ".join(METALS)
    
    content = f"""#!/bin/bash
#SBATCH --job-name=bader_slab
#SBATCH --output={OUTPUT_DIR}/bader_slab.out
#SBATCH --error={OUTPUT_DIR}/bader_slab.err
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

python {BADER_SCRIPT} --metals {metals_list}

echo "Done for all slabs"
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
    print("Generating and submitting sbatch script for Bader slab analysis...\n")
    write_sbatch()
    print("\nDone.")
