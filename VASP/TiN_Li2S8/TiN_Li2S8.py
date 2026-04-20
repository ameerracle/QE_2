#!/usr/bin/env python3

from pathlib import Path

from ase.io import read, Trajectory
from ase.constraints import FixAtoms
from ase.calculators.vasp import Vasp
from ase.optimize import BFGS

# ---------------------------------------------------------------------------
# Input / output
# ---------------------------------------------------------------------------

XYZ_FILE = Path("TiN_Li2S8_combi_144.xyz")
RUN_DIR = "TiN_Li2S8_vasp_relax"
TRAJ_FILE = "TiN_Li2S8_relax.traj"

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------

Z_FIX = 13.0
FMAX = 0.04
ENCUT = 600
SIGMA = 0.2
ISMEAR = 1
AMIX = 0.20
KPTS = (4, 4, 1)
GAMMA = False


def main():
    if not XYZ_FILE.exists():
        raise FileNotFoundError(f"Input file not found: {XYZ_FILE}")

    atoms = read(str(XYZ_FILE))

    # --- Constraints: Fix bottom slab layers -------------------------------
    indices_fixed = [i for i, a in enumerate(atoms) if a.position[2] < Z_FIX]
    atoms.set_constraint(FixAtoms(indices=indices_fixed))

    # --- VASP Calculator (Internal VASP Tags) ------------------------------
    calc = Vasp(
        xc='PBE',
        prec='Normal',
        encut=ENCUT,
        ismear=ISMEAR,
        sigma=SIGMA,
        amix=AMIX,
        ivdw=12,         # VASP internal DFT-D3(BJ)
        kpts=KPTS,
        gamma=GAMMA,
        nsw=0,           # Required: ASE handles ionic steps
        ibrion=-1,       # Required: Disable VASP internal relaxation
        lwave=False,
        lcharg=False,
        directory=RUN_DIR,
        txt="TiN_Li2S8_vasp.out"
    )

    atoms.calc = calc

    # --- BFGS Relaxation ---------------------------------------------------
    # Logfile now named TiN_Li2S8_vasp.log as requested
    qn = BFGS(atoms, logfile="TiN_Li2S8_vasp.log")
    traj = Trajectory(TRAJ_FILE, 'w', atoms)
    qn.attach(traj)

    print(f"Starting BFGS relaxation (fmax={FMAX})...")
    try:
        qn.run(fmax=FMAX)
        print("Done.")
    except Exception as e:
        print(f"Run failed: {e}")


if __name__ == "__main__":
    main()

 