#!/usr/bin/env python3
"""
Extract final structures from slab .out files and save as XYZ in final_xyz folder.
"""

from pathlib import Path
from ase.io import read, write

BASE_DIR = Path(__file__).resolve().parent.parent
SLAB_DIR = BASE_DIR / "slab" / "output_slab_144"
OUTPUT_DIR = SLAB_DIR / "final_xyz"

SLAB_FOLDERS = ["ScN_slab_relax", "TiN_slab_relax", "VN_slab_relax", "NbN_slab_relax", "ZrN_slab_relax"]


def extract_structure(folder_name):
    folder_path = SLAB_DIR / folder_name
    out_files = list(folder_path.glob("*.out"))
    
    if not out_files:
        print(f"  ERROR: No .out files found in {folder_path}")
        return
    
    out_file = out_files[0]
    print(f"Reading {out_file}...")
    
    try:
        atoms = read(str(out_file), index="-1")
        
        metal = folder_name.replace("_slab_relax", "")
        output_path = OUTPUT_DIR / f"{metal}.xyz"
        
        write(str(output_path), atoms)
        
        print(f"  Saved: {output_path}")
        print(f"  Atoms: {len(atoms)}, Formula: {atoms.get_chemical_formula()}")
        
    except Exception as e:
        print(f"  ERROR reading {out_file}: {e}")


if __name__ == "__main__":
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    print("Extracting final structures from slab .out files...\n")
    
    for folder in SLAB_FOLDERS:
        extract_structure(folder)
    
    print(f"\nDone. All structures saved to {OUTPUT_DIR}")
