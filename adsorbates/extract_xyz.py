#!/usr/bin/env python3
"""
Extract final structures from espresso.pwo files and save as XYZ in final_xyz folder.
"""

from pathlib import Path
from ase.io import read, write

BASE_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = BASE_DIR / "final_xyz"

ADSORBATES = ["Li2S4", "Li2S8", "S8", "Li2S6", "Li2S2", "Li2S"]


def extract_structure(adsorbate):
    pwo_path = BASE_DIR / adsorbate / "espresso.pwo"
    
    if not pwo_path.exists():
        print(f"  ERROR: {pwo_path} not found")
        return
    
    print(f"Reading {pwo_path}...")
    
    try:
        atoms = read(str(pwo_path), index="-1")
        
        output_path = OUTPUT_DIR / f"{adsorbate}.xyz"
        
        write(str(output_path), atoms)
        
        print(f"  Saved: {output_path}")
        print(f"  Atoms: {len(atoms)}, Formula: {atoms.get_chemical_formula()}")
        
    except Exception as e:
        print(f"  ERROR reading {pwo_path}: {e}")


if __name__ == "__main__":
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    print("Extracting final structures from espresso.pwo files...\n")
    
    for ads in ADSORBATES:
        extract_structure(ads)
    
    print(f"\nDone. All structures saved to {OUTPUT_DIR}")
