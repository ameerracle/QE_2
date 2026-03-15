import os
import numpy as np
from pathlib import Path

def get_a_from_file(filepath):
    """Extraction of the cubic lattice constant 'a' from a QE file (in or out).
    Assuming the cell is cubic/near-cubic, using (Volume)^(1/3) since nat=8.
    """
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    cell_matrix = None
    in_cell_params = False
    
    # Iterate backwards for output files (to get final params), forwards for input files
    # However, for consistency we'll just find the *last* instance of CELL_PARAMETERS
    found_idx = -1
    for i, line in enumerate(lines):
        if 'CELL_PARAMETERS' in line:
            found_idx = i
            
    if found_idx != -1:
        # Extract the next 3 lines as the matrix
        matrix_lines = lines[found_idx+1:found_idx+4]
        try:
            cell_matrix = []
            for mline in matrix_lines:
                cell_matrix.append([float(x) for x in mline.split()])
            cell_matrix = np.array(cell_matrix)
            vol = np.abs(np.linalg.det(cell_matrix))
            # Since nat=8, it's the 2x2x2 supercell (cubic unit cell).
            # So a = vol^(1/3)
            return vol**(1/3.0)
        except:
            return None
    return None

def main():
    base_dir = Path("/home/ameer_ubuntu/Git_projects/QE_2/cif/final_U_cif")
    compounds = ["TiN", "VN", "ScN", "NbN", "ZrN"]
    
    print(f"{'Compound':<10} | {'Initial a (Å)':<15} | {'Final a (Å)':<15} | {'Expansion (%)':<15}")
    print("-" * 65)
    
    for c in compounds:
        c_dir = base_dir / c
        in_file = c_dir / "vcrelax.in"
        out_file = c_dir / "vcrelax.out"
        
        if not (in_file.exists() and out_file.exists()):
            print(f"{c:<10} | {'Files Missing':<15} | {'--':<15} | {'--':<15}")
            continue
            
        a_init = get_a_from_file(in_file)
        a_final = get_a_from_file(out_file)
        
        if a_init and a_final:
            expansion = (a_final - a_init) / a_init * 100
            print(f"{c:<10} | {a_init:<15.4f} | {a_final:<15.4f} | {expansion:<15.2f}")
        else:
            print(f"{c:<10} | {'Parse Error':<15} | {'--':<15} | {'--':<15}")

if __name__ == "__main__":
    main()
