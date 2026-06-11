#!/usr/bin/env python3
"""
Generate and submit sbatch script for adsorbate PDOS calculations.
Single job to process all 6 adsorbates.
"""

import argparse
import subprocess
from pathlib import Path

ACCOUNT = "def-peslherb"
TASKS_PER_NODE = 1
BASE_DIR = Path("/lustre10/scratch/anizami/QE_2")
SCRIPT = BASE_DIR / "adsorbates" / "adsorbate_pdos.py"
OUTPUT_DIR = BASE_DIR / "submit_scripts" / "adsorbates"


def write_sbatch(nspin: int = 1):
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    tag = f"nspin{nspin}" if nspin == 2 else "nspin1"
    mem = "4G" if nspin == 1 else "8G"
    time = "0:20:00" if nspin == 1 else "5:00:00"
    sbatch_path = OUTPUT_DIR / f"adsorbate_pdos_{tag}.sh"

    nspin_flag = f"--nspin {nspin}" if nspin == 2 else ""

    content = f"""#!/bin/bash
#SBATCH --job-name=ads_pdos_{tag}
#SBATCH --output={OUTPUT_DIR}/adsorbate_pdos_{tag}.out
#SBATCH --error={OUTPUT_DIR}/adsorbate_pdos_{tag}.err
#SBATCH --account={ACCOUNT}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node={TASKS_PER_NODE}
#SBATCH --mem={mem}
#SBATCH --time={time}

set -eo pipefail
module load quantumespresso/7.5 || true
source ~/ase/bin/activate
export OMP_NUM_THREADS=1
cd {BASE_DIR}

python {SCRIPT} --mode full --step full {nspin_flag} --run

echo "Done for all adsorbates ({tag})"
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
    parser = argparse.ArgumentParser()
    parser.add_argument("--nspin", type=int, choices=[1, 2], default=1)
    args = parser.parse_args()

    print(f"Generating and submitting sbatch script for adsorbate PDOS (nspin={args.nspin})...\n")
    write_sbatch(nspin=args.nspin)
    print("\nDone.")
