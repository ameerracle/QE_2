#!/usr/bin/env python3

from pathlib import Path

from ase.calculators.vasp import Vasp
from ase.constraints import FixAtoms
from ase.io import Trajectory, read
from ase.optimize import BFGS

# ---------------------------------------------------------------------------
# Input / output
# ---------------------------------------------------------------------------

XYZ_FILE = Path("TiN_Li2S8_combi_144.xyz")
RUN_DIR = "TiN_Li2S8_vasp_relax"
TRAJ_FILE = Path("TiN_Li2S8_relax.traj")

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
KPAR = 4
NCORE = 4


def resolve_traj_path() -> Path:
    """Support both layouts: trajectory in current folder or one level up."""
    candidates = [TRAJ_FILE, Path("..") / TRAJ_FILE]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError(
        f"Restart trajectory not found. Tried: {', '.join(str(p) for p in candidates)}"
    )


def main():
    traj_path = resolve_traj_path()

    # Continue from the last geometry reached in the previous run.
    atoms = read(str(traj_path), index=-1)

    # Keep fixed atoms consistent with the original slab setup.
    if XYZ_FILE.exists():
        ref_atoms = read(str(XYZ_FILE))
    else:
        ref_atoms = read(str(traj_path), index=0)

    indices_fixed = [i for i, a in enumerate(ref_atoms) if a.position[2] < Z_FIX]
    atoms.set_constraint(FixAtoms(indices=indices_fixed))

    calc = Vasp(
        xc="PBE",
        prec="Normal",
        encut=ENCUT,
        ismear=ISMEAR,
        sigma=SIGMA,
        amix=AMIX,
        ivdw=12,
        kpts=KPTS,
        gamma=GAMMA,
        kpar=KPAR,
        ncore=NCORE,
        nsw=0,
        ibrion=-1,
        lwave=False,
        lcharg=False,
        directory=RUN_DIR,
        txt="TiN_Li2S8_vasp_restart.out",
    )

    atoms.calc = calc

    with Trajectory(str(traj_path), "r") as history:
        nframes = len(history)

    qn = BFGS(atoms, logfile="TiN_Li2S8_vasp_restart.log")

    # Rebuild BFGS curvature information from previous steps.
    if nframes > 1:
        qn.replay_trajectory(str(traj_path))

    # Append new structures to the existing trajectory.
    traj = Trajectory(str(traj_path), "a", atoms)
    qn.attach(traj)

    print(
        f"Restarting BFGS from frame {nframes - 1} "
        f"(fmax={FMAX}, KPAR={KPAR}, NCORE={NCORE})..."
    )
    try:
        qn.run(fmax=FMAX)
        print("Restart finished.")
    except Exception as exc:
        print(f"Restart run failed: {exc}")


if __name__ == "__main__":
    main()
