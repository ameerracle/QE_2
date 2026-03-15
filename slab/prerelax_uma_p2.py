#!/usr/bin/env python3
"""Pre-relax slab structures with FairChem UMA-S-P2.

Default behavior:
- Input folder:  slab/input_slab144
- Input files:   *_slab.xyz
- Output files:  *_slab_uma.xyz
- Log files:     *_slab_uma.log

Example:
    python slab/prerelax_uma_p2.py
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Iterable, List

import torch
from ase.constraints import FixAtoms
from ase.io import read, write
from ase.optimize import BFGS
from fairchem.core import FAIRChemCalculator, pretrained_mlip


def resolve_device(device: str) -> str:
    """Resolve runtime device with ORB_DEVICE override support."""
    env_device = os.getenv("ORB_DEVICE")
    if env_device:
        device = env_device

    if device == "gpu":
        device = "cuda"
    if device == "auto":
        device = "cuda" if torch.cuda.is_available() else "cpu"
    if device.startswith("cuda") and not torch.cuda.is_available():
        raise RuntimeError("CUDA requested but not available on this machine")
    return device


def get_uma_calculator(task: str = "omat", device: str = "auto") -> FAIRChemCalculator:
    """Load FairChem UMA-S-P2 model for solid-state pre-relaxation."""
    resolved = resolve_device(device)
    model = pretrained_mlip.get_predict_unit("uma-s-1p2", device=resolved)
    return FAIRChemCalculator(model, task_name=task)


def fixed_indices_below_z(atoms, z_threshold: float = 13.0) -> List[int]:
    """Return indices for atoms with Cartesian z strictly below z_threshold (Angstrom)."""
    z = atoms.get_positions()[:, 2]
    return [i for i, zi in enumerate(z) if zi < z_threshold]


def iter_input_files(input_dir: Path) -> Iterable[Path]:
    """Yield slab input files in deterministic order."""
    yield from sorted(input_dir.glob("*_slab.xyz"))


def relax_one(atoms, calc: FAIRChemCalculator, fmax: float, max_steps: int, logfile: Path):
    """Run BFGS geometry optimization for a single slab."""
    atoms.pbc = True
    atoms.calc = calc
    dyn = BFGS(atoms, logfile=str(logfile), trajectory=None)
    dyn.run(fmax=fmax, steps=max_steps)
    return atoms


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    input_dir = (script_dir / "input_slab144").resolve()
    output_dir = (script_dir / "uma144_slab").resolve()
    fmax = 0.006
    max_steps = 400
    freeze_z_below = 13.0
    device = "auto"
    task = "omat"

    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")

    files = list(iter_input_files(input_dir))
    if not files:
        raise FileNotFoundError(f"No *_slab.xyz files found in {input_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {output_dir}")

    print(f"Loading UMA-S-P2 on device={resolve_device(device)} ...")
    calc = get_uma_calculator(task=task, device=device)
    print("Model loaded.")

    ok = 0
    fail = 0

    for src in files:
        stem = src.stem
        out_path = output_dir / f"{stem}_uma.xyz"
        log_path = output_dir / f"{stem}_uma.log"

        print(f"\nProcessing {src.name}")
        try:
            atoms = read(str(src))
            frozen_idx = fixed_indices_below_z(
                atoms,
                z_threshold=freeze_z_below,
            )
            if frozen_idx:
                atoms.set_constraint(FixAtoms(indices=frozen_idx))
                print(
                    f"  freezing {len(frozen_idx)} atoms with z < {freeze_z_below:.2f} A"
                )
                if len(frozen_idx) != 72:
                    print("  WARNING: expected 72 fixed atoms for this slab")
            else:
                print(f"  freezing 0 atoms with z < {freeze_z_below:.2f} A")

            relaxed = relax_one(
                atoms,
                calc=calc,
                fmax=fmax,
                max_steps=max_steps,
                logfile=log_path,
            )
            write(str(out_path), relaxed)
            print(f"  wrote {out_path.name}")
            ok += 1
        except Exception as exc:
            print(f"  FAILED: {src.name} -> {exc}")
            fail += 1

    print("\nDone.")
    print(f"Success: {ok}")
    print(f"Failed:  {fail}")


if __name__ == "__main__":
    main()
