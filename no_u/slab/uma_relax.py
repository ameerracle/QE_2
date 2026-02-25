#!/usr/bin/env python3
"""
Relax slabs from input_slab_144 using the UMA model (fairchem),
fixing the bottom 4 layers.

Run from: ~/Git_projects/QE_2/no_u/slab/

Output written as: {metal}_{num}_slab_uma.xyz
"""
import os
from pathlib import Path

import numpy as np
import torch
from ase.constraints import FixAtoms
from ase.optimize import BFGS
import ase.io

from fairchem.core import FAIRChemCalculator, pretrained_mlip


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def get_uma_calculator(task: str = "omat", device: str = "auto") -> FAIRChemCalculator:
    """Load the UMA-SM model for solid-state materials (omat task)."""
    env_device = os.getenv("ORB_DEVICE")
    if env_device:
        device = env_device
    if device == "gpu":
        device = "cuda"
    if device == "auto":
        device = "cuda" if torch.cuda.is_available() else "cpu"
    if device.startswith("cuda") and not torch.cuda.is_available():
        raise RuntimeError("CUDA requested but not available")

    model = pretrained_mlip.get_predict_unit("uma-s-1p1", device=device)
    return FAIRChemCalculator(model, task_name=task)


def bottom_layer_indices(atoms, layer_count: int = 4, axis: int = 2, tolerance: float = 1e-3):
    """Return indices of atoms in the bottom *layer_count* layers (by scaled coords)."""
    scaled = atoms.get_scaled_positions()
    if scaled.size == 0:
        return [], [], []
    coords = np.mod(scaled[:, axis], 1.0)
    order = np.argsort(coords)
    layers, current = [], [order[0]]
    reference = coords[order[0]]
    for idx in order[1:]:
        if abs(coords[idx] - reference) <= tolerance:
            current.append(idx)
        else:
            layers.append(current)
            if len(layers) >= layer_count:
                break
            current, reference = [idx], coords[idx]
    if len(layers) < layer_count:
        layers.append(current)
    selected = layers[:layer_count]
    return (
        [i for layer in selected for i in layer],
        [np.mean(coords[layer]) for layer in selected],
        [len(layer) for layer in selected],
    )


def relax(atoms, calc, fmax: float = 0.05, logfile=None):
    """Attach *calc*, apply constraints and run BFGS until convergence."""
    atoms.pbc = True
    atoms.calc = calc
    dyn = BFGS(atoms, logfile=logfile, trajectory=None)
    dyn.run(fmax=fmax)
    return atoms


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    metals = ["TiN", "VN", "ScN", "ZrN", "NbN"]
    fmax   = 0.008  # eV/Å

    # Load UMA once and reuse across all slabs
    print("Loading UMA model ...")
    calc = get_uma_calculator(task="omat")
    print("Model ready.")

    # Find input_slab_144 relative to this script's location
    script_dir = Path(__file__).resolve().parent
    slab_dir = script_dir / "input_slab_144"
    if not slab_dir.exists():
        print(f"Directory {slab_dir} not found.")
        return

    for metal in metals:
        print(f"\n{metal}:")
        for src in sorted(slab_dir.glob(f"{metal}_*_slab.xyz")):
            match = src.stem.split("_")
            if len(match) < 2:
                print(f"  skip  {src.name} (not found)")
                continue

            num = match[1]
            dst = slab_dir / f"{metal}_{num}_slab_uma.xyz"
            logfile = slab_dir / f"{metal}_{num}_slab_uma.log"

            try:
                atoms = ase.io.read(str(src))

                # Fix bottom 1 layer so only the top/middle surface relaxes
                frozen, heights, counts = bottom_layer_indices(atoms, layer_count=2)
                if frozen:
                    atoms.set_constraint(FixAtoms(indices=frozen))
                    print(f"  {metal} ({num}): fixing {len(frozen)} atoms in bottom 2 layers")
                else:
                    mask = atoms.positions[:, 2] <= 13
                    n_fixed = np.sum(mask)
                    atoms.set_constraint(FixAtoms(mask=mask))
                    print(f"  {metal} ({num}): layer detection failed — fixing {n_fixed} atoms with z <= 13 A")

                atoms = relax(atoms, calc, fmax=fmax, logfile=str(logfile))
                ase.io.write(str(dst), atoms)
                print(f"  OK  {dst.name}")

            except Exception as exc:
                print(f"  FAIL  {metal}_{num}_slab.xyz : {exc}")

    print("\nAll done.")


if __name__ == "__main__":
    main()