#!/usr/bin/env python3
"""Relax each ordered metal nitride slab with a fixed single setup."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

from ase import io
from ase.calculators.espresso import Espresso


PSEUDOS = {
    "Ti": "ti_pbe_v1.4.uspp.F.UPF",
    "V": "v_pbe_v1.4.uspp.F.UPF",
    "Sc": "sc_pbe_v1.uspp.F.UPF",
    "Nb": "nb_pbe_v1.uspp.F.UPF",
    "Zr": "zr_pbe_v1.uspp.F.UPF",
    "N": "n_pbe_v1.2.uspp.F.UPF",
}

SLAB_NAMES = ["TiN", "VN", "ScN", "NbN", "ZrN"]
SLAB_DIR = Path("slab/input_slab")
RUN_ROOT = Path(".")
PSEUDO_DIR = Path("/scratch/anizami/QE_2/USPP/")
SLAB_SUFFIX = "_144_slab.xyz"
KPTS = (4, 4, 1)
ECUT_WFC = 45.0
ECUT_RHO = 8 * ECUT_WFC


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Relax the selected metal nitride slab with a fixed setup."
    )
    parser.add_argument(
        "--structure",
        choices=SLAB_NAMES,
        help="Name of the slab to relax (TiN, VN, ScN, NbN, ZrN)."
    )
    return parser.parse_args()


def get_structure_path(name: str) -> Path:
    path = SLAB_DIR / f"{name}{SLAB_SUFFIX}"
    if not path.exists():
        raise SystemExit(f"Missing slab file for {name}: {path}")
    return path


def relax_slab(name: str) -> None:
    metal = name.rstrip("N")
    pseudo_map = {metal: PSEUDOS[metal], "N": PSEUDOS["N"]}
    structure = get_structure_path(name)
    run_dir = RUN_ROOT / structure.stem
    run_dir.mkdir(parents=True, exist_ok=True)

    os.environ.setdefault("OMP_NUM_THREADS", "1")

    atoms = io.read(structure)
    calc = Espresso(
        command="pw.x",
        directory=str(run_dir),
        pseudopotentials=pseudo_map,
        pseudo_dir=str(PSEUDO_DIR),
        kpts=KPTS,
        input_data={
            "control": {
                "calculation": "relax",
                "prefix": "slab",
                "disk_io": "low",
            },
            "system": {
                "ecutwfc": ECUT_WFC,
                "ecutrho": ECUT_RHO,
            },
            "input_dft": "pbe",
            "vdw_corr": "dft-d3",
            "dftd3_version": 4,
            "occupations": "smearing",
            "smearing": "cold",
            "degauss": 0.005,
        },
    )

    atoms.calc = calc
    energy = atoms.get_potential_energy()
    output_file = run_dir / f"{structure.stem}_relaxed.xyz"
    io.write(output_file, atoms, format="xyz")

    print(f"Relaxed {structure.name} -> {output_file} | energy = {energy:.6f} eV")


def main() -> None:
    args = parse_arguments()
    RUN_ROOT.mkdir(parents=True, exist_ok=True)
    SLAB_DIR.mkdir(parents=True, exist_ok=True)

    targets = [args.structure] if args.structure else SLAB_NAMES
    for name in targets:
        relax_slab(name)


if __name__ == "__main__":
    main()
