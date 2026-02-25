import os
from pathlib import Path
from ase import io
from ase.build import add_vacuum
from ase.calculators.espresso import Espresso, EspressoProfile
import argparse

# default list of adsorbates
DEFAULT_ADSORBATES = ["Li2S", "Li2S2", "Li2S4", "Li2S6", "Li2S8", "S8"]

# simple CLI: positional single adsorbate (keeps original behaviour when omitted)
parser = argparse.ArgumentParser(
    description="Relax a single adsorbate (positional). If omitted, all default adsorbates are run."
)
parser.add_argument(
    "adsorbate",
    nargs="?",
    choices=DEFAULT_ADSORBATES,
    help="Name of adsorbate to run (e.g. Li2S)."
)
args = parser.parse_args()
if args.adsorbate:
    adsorbates = [args.adsorbate]
else:
    adsorbates = DEFAULT_ADSORBATES

run_root = Path(".")
run_root.mkdir(exist_ok=True)
pseudo_dir = Path("/scratch/anizami/QE_2/USPP/")
vacuum = 10
ecutwfc = 45.0
ecutrho = 10 * ecutwfc
kpts = (1, 1, 1)
# Convergence thresholds (in Rydberg and Ry/Bohr)


# use pseudopotential filenames from `dir_notes.txt` (case-sensitive)
# Nibi / server pseudos: li_pbe_v1.4.uspp.F.UPF, s_pbe_v1.4.uspp.F.UPF
pseudos = {"Li": "li_pbe_v1.4.uspp.F.UPF", "S": "s_pbe_v1.4.uspp.F.UPF"}

for ads in adsorbates:
    structure_path = Path(f"input_ads/{ads}_final.extxyz")
    run_dir = run_root / ads 
    run_dir.mkdir(parents=True, exist_ok=True)

    atoms = io.read(structure_path)
    add_vacuum(atoms, vacuum)

    os.environ.setdefault("OMP_NUM_THREADS", "1")
    # Prefer the scheduler launcher (srun) so MPI ranks match the allocation.
    # Allow override via QE_PW_COMMAND for testing or non-SLURM systems.
    pw_cmd = os.environ.get("QE_PW_COMMAND")
    if not pw_cmd:
        if "SLURM_NTASKS" in os.environ:
            pw_cmd = "srun pw.x"
        else:
            pw_cmd = "pw.x"

    profile = EspressoProfile(command=pw_cmd, pseudo_dir=str(pseudo_dir))
    
    calc = Espresso(
        profile=profile,
        pseudopotentials=pseudos,
        input_data={
            "control": {
                "calculation": "relax",      # Uses QE native optimizer
                "prefix": ads,
            },
            "system": {
                "ecutwfc": ecutwfc, 
                "ecutrho": ecutrho,
                "input_dft": "pbe",    
                "occupations": "smearing",
                "smearing": "gaussian",
                "degauss": 0.001,
                "vdw_corr": 'dft-d3',
                "dftd3_version": 4,
            },
            "electrons": {
                "conv_thr": 1.0e-6,
                "mixing_beta": 0.5,
                "electron_maxstep": 150,
                "mixing_mode": "local-TF",
            },
            "ions": {
                "ion_dynamics": "bfgs",      # QE default for relax
            }
        },
        kpts=kpts,
        pseudo_dir=str(pseudo_dir),
        directory=str(run_dir),
    )
    
    # Run calculation
    atoms.calc = calc
    atoms.get_potential_energy()

    # Save final optimized structure to a clean .extxyz
    output_path = run_dir / f"{ads}.relaxed.extxyz"
    io.write(output_path, atoms, format="extxyz")
    
    print(f"Relaxation complete for {ads}. Structure -> {output_path}")