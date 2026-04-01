#!/usr/bin/env python3
"""
Relax TiN_S8 combined slab+adsorbate structure using VASP via ASE.

Conversion notes (all Ry->eV via 1 Ry = 13.6057 eV):
  ecutwfc  45.0 Ry        ->  encut  = 600 eV
  conv_thr               ->  ediff  = VASP default (1e-4 eV); not set explicitly
  forc_conv_thr 0.00156   ->  ediffg = -0.04 eV/A  (0.00156 Ry/Bohr * 25.7110)
  degauss  0.025 Ry       ->  sigma  = 0.2 eV
  smearing = 'cold'       ->  ismear = 1  (Methfessel-Paxton order 1)
  mixing_beta 0.25        ->  amix   = 0.25  (VASP default is 0.4)
  electron_maxstep 250    ->  nelm   = 250
  nstep 160               ->  nsw    = 160
  ion_dynamics = 'bfgs'   ->  ibrion = 2  (CG; VASP has no native BFGS)
  upscale 750             ->  (no VASP equivalent; omitted)
  dftd3_version = 4       ->  ivdw   = 12 (DFT-D3 with BJ damping)
  kpts 4 4 1 0 0 0        ->  kpts=(4,4,1), gamma=False (Monkhorst-Pack, no shift)

Pseudopotentials: VASP default PAW-PBE (no setups keyword = ASE 'minimal',
one PAW folder per element: Ti, N, S).

ecutrho has no VASP equivalent; PAW handles the augmentation cutoff
automatically via prec='Normal' (~4*encut).

Environment variables required on the cluster:
  export ASE_VASP_COMMAND="mpirun vasp_std"  # adjust to your scheduler
  export VASP_PP_PATH=/path/to/vasp/potentials

Run:
  python combi_relax_vasp.py [--dry-run]
"""

import argparse
from pathlib import Path

from ase.io import read
from ase.constraints import FixAtoms
from ase.calculators.vasp import Vasp

# ---------------------------------------------------------------------------
# Input / output
# ---------------------------------------------------------------------------
XYZ_FILE = Path("TiN_S8_combi_144.xyz")
RUN_DIR  = Path("TiN_S8_vasp_relax")

# ---------------------------------------------------------------------------
# DFT parameters
# ---------------------------------------------------------------------------
Z_FIX  = 13.0   # Angstrom -- atoms below this z are frozen (bottom 2 slab layers)

ENCUT  = 600    # eV
EDIFFG = -0.04  # eV/A  (negative = force criterion)
SIGMA  = 0.2    # eV
ISMEAR = 1      # Methfessel-Paxton order 1
IBRION = 2      # conjugate gradient
NSW    = 160    # max ionic steps
NELM   = 250    # max SCF iterations per ionic step
AMIX   = 0.25   # charge density mixing (QE used 0.25; VASP default is 0.4)
IVDW   = 12     # DFT-D3 with BJ damping
KPTS   = (4, 4, 1)
GAMMA  = False  # Monkhorst-Pack, no shift

# ---------------------------------------------------------------------------
def main(dry_run: bool = False):
    if not XYZ_FILE.exists():
        raise FileNotFoundError(f"Input file not found: {XYZ_FILE}")

    atoms = read(str(XYZ_FILE))
    nat   = len(atoms)
    print(f"Read {XYZ_FILE}  |  nat = {nat}")

    # --- Selective dynamics: fix all atoms below Z_FIX ----------------------
    indices_fixed = [i for i, a in enumerate(atoms) if a.position[2] < Z_FIX]
    n_fixed = len(indices_fixed)
    atoms.set_constraint(FixAtoms(indices=indices_fixed))
    print(f"Fixed {n_fixed} atoms (z < {Z_FIX:.1f} A), "
          f"{nat - n_fixed} atoms free to relax.")

    # --- Output directory ---------------------------------------------------
    RUN_DIR.mkdir(exist_ok=True)
    print(f"Output dir: {RUN_DIR}")

    # --- VASP calculator ----------------------------------------------------
    # ASE writes INCAR, POSCAR, POTCAR, and KPOINTS into RUN_DIR.
    calc = Vasp(
        # XC / pseudopotentials
        # No 'setups' keyword: ASE resolves to default PAW-PBE folders per element.
        xc     = "PBE",

        # Plane-wave cutoff
        encut  = ENCUT,
        prec   = "Normal",   # augmentation cutoff ~4*encut; standard PAW choice

        # SCF convergence (ediff not set -> VASP default 1e-4 eV)
        nelm   = NELM,

        # Smearing
        ismear = ISMEAR,
        sigma  = SIGMA,

        # Charge density mixing
        amix   = AMIX,

        # Ionic relaxation
        ibrion = IBRION,     # CG (no native BFGS in VASP)
        nsw    = NSW,
        ediffg = EDIFFG,
        isif   = 2,          # relax ions only, keep cell fixed

        # DFT-D3(BJ)
        ivdw   = IVDW,

        # k-points
        kpts   = KPTS,
        gamma  = GAMMA,

        # I/O
        nwrite = 1,          # reduce OUTCAR verbosity
        lwave  = False,      # skip WAVECAR (set True if restarts are needed)
        lcharg = False,      # skip CHGCAR  (set True for charge analysis)

        # Parallelisation -- uncomment and tune for your cluster:
        # ncore = 16,        # cores per band; ~sqrt(total cores) is a good start
        # kpar  = 1,         # k-point parallelism; 1 is fine for a 4x4x1 mesh

        directory = str(RUN_DIR),
        txt       = "TiN_S8_vasp.out",
    )

    atoms.calc = calc

    if dry_run:
        # Writes INCAR / POSCAR / POTCAR / KPOINTS without running VASP.
        calc.write_input(atoms)
        print(f"[dry-run] Input files written to {RUN_DIR}/")
        return

    print("Launching VASP ...")
    try:
        atoms.get_potential_energy()
        print("Done.")
    except Exception as e:
        print(f"VASP run failed: {e}")


# ---------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Relax TiN_S8 slab+adsorbate with VASP DFT-D3 via ASE."
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Write input files only; do not launch VASP.",
    )
    args = parser.parse_args()
    main(dry_run=args.dry_run)