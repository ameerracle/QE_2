#!/usr/bin/env python3
"""
Fix adsorbate XYZ files by creating smaller cubic cells with 10Å vacuum.
"""

from pathlib import Path
from ase.io import read, write
from ase.build import molecule

INPUT_DIR = Path("/lustre10/scratch/anizami/QE_2/adsorbates/final_xyz")
OUTPUT_DIR = Path("/lustre10/scratch/anizami/QE_2/adsorbates/final_xyz_fixed")

ADSORBATES = ["Li2S", "Li2S2", "Li2S4", "Li2S6", "Li2S8", "S8"]

def fix_cell(adsorbate):
    input_path = INPUT_DIR / f"{adsorbate}.xyz"
    
    if not input_path.exists():
        print(f"  ERROR: {input_path} not found")
        return
    
    print(f"Reading {input_path}...")
    atoms = read(str(input_path))
    
    # Get current positions
    positions = atoms.get_positions()
    
    # Calculate center
    center = positions.mean(axis=0)
    
    # Center the molecule
    atoms.set_positions(positions - center)
    
    # Create a cubic cell with 10Å vacuum in each direction
    max_extent = (positions.max(axis=0) - positions.min(axis=0)).max()
    cell_size = max_extent + 10.0  # 10Å vacuum
    
    atoms.set_cell([cell_size, cell_size, cell_size])
    atoms.set_pbc([True, True, True])
    
    output_path = OUTPUT_DIR / f"{adsorbate}.xyz"
    write(str(output_path), atoms, format="extxyz")
    
    print(f"  Saved: {output_path}")
    print(f"  Cell: {cell_size:.2f} x {cell_size:.2f} x {cell_size:.2f} Å")
    print(f"  Atoms: {len(atoms)}, Formula: {atoms.get_chemical_formula()}")


if __name__ == "__main__":
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    print("Fixing adsorbate cell sizes...\n")
    
    for ads in ADSORBATES:
        fix_cell(ads)
    
    print(f"\nDone. All structures saved to {OUTPUT_DIR}")
