#!/usr/bin/env python3
"""
Submit SLURM jobs for adsorbate relaxations.

Creates and submits a separate job for each adsorbate in the hard-coded list.
Simple / minimal: fixed 3h walltime, 30 tasks, 30GB RAM, sets OMP_NUM_THREADS=1,
and uses SBATCH --account=def-peslherb. Each job runs `ads_relax.py <adsorbate>`
from the submission directory.

Usage:
    python submit_ads.py

"""
from pathlib import Path
import subprocess
import textwrap

# keep the adsorbate list in sync with ads_relax.py's defaults
ADSORBATES = ["Li2S", "Li2S2", "Li2S4", "Li2S6", "Li2S8", "S8"]

# cluster resources (minimal, as requested)
WALLTIME = "03:00:00"
NTASKS = 24
MEM = "24GB"
ACCOUNT = "def-peslherb"

SCRIPT_NAME = "ads_relax.py"


def create_slurm_script(adsorbate: str) -> str:
    job_name = f"ads_relax_{adsorbate}"
    return textwrap.dedent(f"""\
    #!/bin/bash
    #SBATCH --job-name={job_name}
    #SBATCH --output={job_name}_%j.out
    #SBATCH --error={job_name}_%j.err
    #SBATCH --account={ACCOUNT}
    #SBATCH --ntasks-per-node={NTASKS}
    #SBATCH --mem={MEM}
    #SBATCH --time={WALLTIME}

    set -euo pipefail

    # load modules / activate env if needed (adjust for your cluster)
    module load quantumespresso/7.5 || true
    # source ~/qe/bin/activate || true

    export OMP_NUM_THREADS=1

    # run from the submission directory so relative paths in the repo work
    cd "$SLURM_SUBMIT_DIR"

    python {SCRIPT_NAME} {adsorbate}
    """)


def submit_all(dst_dir: Path) -> None:
    for ads in ADSORBATES:
        script_path = dst_dir / f"submit_{ads}.sh"
        script_text = create_slurm_script(ads)
        script_path.write_text(script_text)
        script_path.chmod(0o755)

        try:
            res = subprocess.run(["sbatch", str(script_path)], capture_output=True, text=True, cwd=str(dst_dir))
        except FileNotFoundError as e:
            print(f"✗ sbatch not found: {e}")
            return
        except Exception as e:
            print(f"✗ error while submitting {ads}: {e}")
            continue

        if res.returncode == 0:
            print(f"✓ Submitted {ads}: {res.stdout.strip()}")
        else:
            print(f"✗ Failed to submit {ads}: {res.stderr.strip() or res.stdout.strip()}")


if __name__ == "__main__":
    repo_dir = Path(__file__).parent
    print("Submitting one SLURM job per adsorbate (3:00:00, 30 tasks, 30GB)...")
    submit_all(repo_dir)
    print("Done.")
