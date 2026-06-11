#!/usr/bin/env python3
import subprocess
import re
from pathlib import Path

# --- Cluster & Account Settings ---
ACCOUNT = "def-peslherb"
TASKS_PER_NODE = 32
MEM = "50G"
TIME = "43:00:00"
ADSORBATE_SPECIES = ['Li', 'S']
SLAB_METALS = ['Ti', 'Nb', 'Zr', 'Sc', 'V']

# --- Path Configuration ---
LOBSTER_EXE = "/home/anizami/lobster/lobster-5.1.1"
LOG_FILE = "automated_lobster_submissions.txt"
RUN_ROOT = Path("/lustre10/scratch/anizami/QE_2/no_u/lobster")


def get_prefix(file_path):
    prefix = "pwscf"
    try:
        with open(file_path, 'r') as f:
            content = f.read()
            match = re.search(r"prefix\s*=\s*['\"]([^'\"]+)['\"]", content, re.IGNORECASE)
            if match:
                prefix = match.group(1)
    except Exception:
        pass
    return prefix


def find_save_folder(run_dir, prefix):
    paths_to_check = [run_dir, run_dir / "tmp"]
    for p in paths_to_check:
        candidate = p / f"{prefix}.save"
        if candidate.exists():
            return candidate
    return None


def prepare_lobster_files(run_dir, nscf_filename):
    src_input = run_dir / nscf_filename
    prefix = get_prefix(src_input)

    lobster_cwd = run_dir / "lobster_resolved_600steps"
    lobster_cwd.mkdir(parents=True, exist_ok=True)

    tmp_dir = lobster_cwd / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    with open(src_input, 'r') as f:
        input_content = f.read()

    input_content = re.sub(
        r"outdir\s*=\s*['\"]?[^,\n'\"]+['\"]?",
        "outdir = './tmp'",
        input_content,
        flags=re.IGNORECASE
    )

    dst_input = lobster_cwd / f"{prefix}.scf.in"
    with open(dst_input, 'w') as f:
        f.write(input_content)

    real_save_folder = find_save_folder(run_dir, prefix)
    if not real_save_folder:
        raise FileNotFoundError(f"Could not find {prefix}.save in {run_dir}")

    save_link = tmp_dir / f"{prefix}.save"
    if save_link.is_symlink() or save_link.exists():
        save_link.unlink()
    save_link.symlink_to(real_save_folder)

    lobsterin_lines = [
        "basisSet pbeVaspFit2015",
        "useRecommendedBasisFunctions",
        "COHPStartEnergy -12",
        "COHPEndEnergy 8",
        "COHPsteps 600",
        "gaussianSmearingWidth 0.27",
        "skipPopulationAnalysis",
        "skipGrossPopulation",
        "skipCOOP",
        "skipDOS",
        "skipCOBI",
        "saveProjectionToFile",
        "",
        "# Distance-limited adsorbate-metal pair generation",
    ]

    for ads_species in ADSORBATE_SPECIES:
        for metal_species in SLAB_METALS:
            lobsterin_lines.append(
                f"cohpGenerator from 1.0 to 3.5 type {ads_species} type {metal_species}"
            )

    with open(lobster_cwd / "lobsterin", "w") as f:
        f.write("\n".join(lobsterin_lines))

    return lobster_cwd, prefix


def submit_all_jobs():
    nscf_files = list(RUN_ROOT.glob("**/*nscf.in"))

    if not nscf_files:
        print(f"No *nscf.in files found under {RUN_ROOT}")
        return

    print(f"Found {len(nscf_files)} NSCF input(s) to process.")

    with open(LOG_FILE, "w") as log:
        for nscf_path in nscf_files:
            run_dir = nscf_path.parent
            job_tag = run_dir.name
            print(f"\nProcessing: {job_tag}")

            try:
                lobster_cwd, prefix = prepare_lobster_files(
                    run_dir, nscf_path.name
                )
            except Exception as e:
                print(f"  Skipping {job_tag}: {e}")
                continue

            job_name = f"LOB_{job_tag}"
            submit_script = lobster_cwd / "submit_lobster.sh"

            sbatch_lines = [
                "#!/bin/bash",
                f"#SBATCH --job-name={job_name}",
                f"#SBATCH --account={ACCOUNT}",
                f"#SBATCH --cpus-per-task={TASKS_PER_NODE}",
                f"#SBATCH --mem={MEM}",
                f"#SBATCH --time={TIME}",
                "",
                "source ~/ase/bin/activate",
                "module load quantumespresso",
                f"export OMP_NUM_THREADS={TASKS_PER_NODE}",
                f"cd {lobster_cwd}",
                f"{LOBSTER_EXE} > lobster.out 2>&1",
            ]

            with open(submit_script, "w") as f:
                f.write("\n".join(sbatch_lines))

            result = subprocess.run(
                ["sbatch", str(submit_script)],
                capture_output=True, text=True
            )

            if result.returncode == 0:
                job_id = result.stdout.strip()
                print(f"  Submitted: {job_name} -- {job_id}")
                log.write(f"{job_name}\t{job_id}\t{lobster_cwd}\n")
            else:
                print(f"  sbatch FAILED for {job_name}:\n{result.stderr.strip()}")
                log.write(f"{job_name}\tFAILED\t{lobster_cwd}\n")


if __name__ == "__main__":
    submit_all_jobs()