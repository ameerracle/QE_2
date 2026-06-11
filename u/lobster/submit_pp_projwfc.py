#!/usr/bin/env python3
import argparse
import subprocess
from pathlib import Path

ACCOUNT = "def-peslherb"
NTASKS = 32
MEM = "200G"
TIME = "02:00:00"
BASE = Path("/scratch/anizami/QE_2/u/lobster")

SYSTEMS = ["VN_Li2S4", "VN_Li2S8", "VN_S8"]
TASKS = ["pp_plot0", "projwfc"]


def main():
    parser = argparse.ArgumentParser(description="Submit pp_plot0 and/or projwfc for VN+U lobster systems")
    parser.add_argument("--systems", nargs="+", choices=SYSTEMS, default=SYSTEMS)
    parser.add_argument("--tasks", nargs="+", choices=TASKS, default=TASKS)
    parser.add_argument("--mem", default=MEM)
    parser.add_argument("--time", default=TIME)
    parser.add_argument("--ntasks", type=int, default=NTASKS)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    steps = []
    for sys in args.systems:
        for task in args.tasks:
            steps.append((sys, task))

    if not steps:
        print("No tasks selected.")
        return

    script_lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name=u_pp_proj",
        f"#SBATCH --output=u_pp_proj_%j.out",
        f"#SBATCH --error=u_pp_proj_%j.err",
        f"#SBATCH --account={ACCOUNT}",
        f"#SBATCH --nodes=1",
        f"#SBATCH --ntasks-per-node={args.ntasks}",
        f"#SBATCH --mem={args.mem}",
        f"#SBATCH --time={args.time}",
        "",
        "source ~/ase/bin/activate",
        "module load quantumespresso",
        "set -eo pipefail",
        "",
    ]

    for sys, task in steps:
        run_dir = BASE / f"{sys}_combi"
        if task == "pp_plot0":
            inp = f"{sys}_pp_plot0.in"
            out = f"{sys}_pp_plot0.out"
            cmd = f"srun pp.x < {inp} > {out}"
        elif task == "projwfc":
            inp = f"{sys}_projwfc.in"
            out = f"{sys}_projwfc.out"
            cmd = f"srun projwfc.x < {inp} > {out}"
        else:
            continue

        script_lines.append(f"cd {run_dir}")
        script_lines.append(f"echo \"=== {sys} {task} ===\"")
        script_lines.append(f"{cmd}")
        script_lines.append("")

    script_lines.append('echo "ALL DONE"')

    script_text = "\n".join(script_lines) + "\n"

    if args.dry_run:
        print(script_text)
        return

    script_path = BASE / "submit_pp_projwfc_auto.sh"
    with open(script_path, "w") as f:
        f.write(script_text)

    result = subprocess.run(["sbatch", str(script_path)], capture_output=True, text=True)
    print(result.stdout.strip())
    if result.stderr:
        print(result.stderr.strip())


if __name__ == "__main__":
    main()
