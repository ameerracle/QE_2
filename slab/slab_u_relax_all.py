#!/usr/bin/env python3

import os
import subprocess
import re
import shutil
import time
import argparse
from pathlib import Path
from ase import io
from ase import Atoms
from ase.constraints import FixAtoms
import numpy as np

# --- Settings ---
PSEUDO_DIR = "/scratch/anizami/QE_2/USPP"
QE_BIN = "/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/MPI/gcc12/openmpi4/quantumespresso/7.5/bin/"
PW_CMD = "srun " + QE_BIN + "pw.x"
HP_CMD = "srun " + QE_BIN + "hp.x"
TOLERANCE = 0.01

START_U_VALUES = {
    "VN": 5.8675,
    "TiN": 5.5479,
    "NbN": 3.5140,
    "ZrN": 2.9683,
    "ScN": 3.5695
}

PSEUDO_MAP = {
    "V": "v_pbe_v1.4.uspp.F.UPF",
    "Ti": "ti_pbe_v1.4.uspp.F.UPF",
    "Sc": "sc_pbe_v1.uspp.F.UPF",
    "Nb": "nb_pbe_v1.uspp.F.UPF",
    "Zr": "zr_pbe_v1.uspp.F.UPF",
    "N": "n_pbe_v1.2.uspp.F.UPF"
}


def parse_cell_from_in(relax_in_file):
    """Extract CELL_PARAMETERS from relax input file."""
    if not os.path.exists(relax_in_file):
        return None
    with open(relax_in_file, 'r') as f:
        lines = f.readlines()
    cell = []
    cell_idx = -1
    
    for i, line in enumerate(lines):
        if "CELL_PARAMETERS" in line:
            cell_idx = i
            break
    
    if cell_idx == -1:
        return None
    
    try:
        for j in range(1, 4):
            cell.append([float(x) for x in lines[cell_idx + j].split()])
    except:
        return None
    return np.array(cell)


def parse_coords_from_out(relax_out_file):
    """Extract atomic coordinates from relax output file (Begin final coordinates)."""
    if not os.path.exists(relax_out_file):
        return None
    with open(relax_out_file, 'r') as f:
        lines = f.readlines()
    positions = []
    symbols = []
    
    begin_idx = -1
    end_idx = -1
    
    for i, line in enumerate(lines):
        if "Begin final coordinates" in line:
            begin_idx = i
        if begin_idx != -1 and "End final coordinates" in line:
            end_idx = i
            break
    
    if begin_idx == -1 or end_idx == -1:
        return None
    
    try:
        for i in range(begin_idx + 1, end_idx):
            line = lines[i].strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) >= 4:
                symbols.append(parts[0])
                positions.append([float(parts[1]), float(parts[2]), float(parts[3])])
    except:
        return None
    if len(positions) == 0:
        return None
    return symbols, np.array(positions)


def parse_relaxed_geometry(relax_in_file, relax_out_file):
    """Extract geometry: cell from .in file, coords from .out file."""
    cell = parse_cell_from_in(relax_in_file)
    coords_data = parse_coords_from_out(relax_out_file)
    
    if cell is None or coords_data is None:
        return None
    
    symbols, positions = coords_data
    return Atoms(symbols=symbols, positions=positions, cell=cell, pbc=True)


def extract_u_from_dat(folder, prefix):
    """Reads Hubbard U from .dat file."""
    dat_path = folder / (prefix + ".Hubbard_parameters.dat")
    if not dat_path.exists():
        return None
    
    try:
        with open(str(dat_path), 'r') as f:
            content = f.read()
            matches = re.findall(r"[34]d\s+([\d\.]+)", content)
            if matches:
                return float(matches[-1])
            return None
    except:
        return None


def get_element_symbol(metal):
    """Extract element symbol from metal compound (e.g., 'TiN' -> 'Ti')."""
    if metal[:2] in ["Ti", "Sc", "Nb", "Zr", "Va"]:
        return metal[:2]
    else:
        return metal[0]


def get_manifold(symbol):
    """Get the orbital manifold for the element."""
    return "4d" if symbol in ["Nb", "Zr"] else "3d"


def get_occupancy_options(symbol):
    """Get occupancy options based on element type."""
    if symbol == "Sc":
        # ScN is a semiconductor/insulator with bandgap
        return {'occupations': 'fixed'}
    else:
        # Others are metallic
        return {
            'occupations': 'smearing',
            'smearing': 'cold',
            'degauss': 0.005
        }


def setup_relax_input(iter_dir, prefix, atoms, symbol, manifold, current_u):
    """Create relaxation input file and return atoms with reordered indices."""
    relax_in = iter_dir / (prefix + ".relax.in")
    tmp_dir = iter_dir / "tmp"
    
    # Reorder atoms so Hubbard atoms come first
    relax_atoms = atoms.copy()
    symbols_list = relax_atoms.get_chemical_symbols()
    hubbard_indices = [i for i, sym in enumerate(symbols_list) if sym == symbol]
    other_indices = [i for i, sym in enumerate(symbols_list) if sym != symbol]
    new_order = hubbard_indices + other_indices
    relax_atoms = relax_atoms[new_order]
    
    relax_data = {
        'control': {
            'calculation': 'relax',
            'prefix': prefix,
            'pseudo_dir': PSEUDO_DIR,
            'outdir': str(tmp_dir),
            'tprnfor': True,
            'forc_conv_thr': 0.001,
            'nstep': 150
        },
        'system': {
            'ecutwfc': 45.0,
            'ecutrho': 450.0
        },
        'electrons': {
            'conv_thr': 1.0e-6,
            'mixing_beta': 0.5,
            'electron_maxstep': 150
        },
        'ions': {
            'ion_dynamics': 'bfgs'
        }
    }
    relax_data['system'].update(get_occupancy_options(symbol))
    
    io.write(
        str(relax_in),
        relax_atoms,
        format='espresso-in',
        input_data=relax_data,
        pseudopotentials={symbol: PSEUDO_MAP[symbol], "N": PSEUDO_MAP["N"]},
        kpts=(5, 5, 1)
    )
    
    with open(str(relax_in), 'a') as f:
        f.write("\nHUBBARD {ortho-atomic}\nU " + symbol + "-" + manifold + " " + str(current_u) + "\n")
    
    return relax_in


def setup_scf_input(iter_dir, prefix, atoms, symbol, manifold, current_u):
    """Create SCF input file and return atoms with reordered indices."""
    scf_in = iter_dir / (prefix + ".scf.in")
    tmp_dir = iter_dir / "tmp"
    
    # Reorder atoms so Hubbard atoms come first
    scf_atoms = atoms.copy()
    symbols_list = scf_atoms.get_chemical_symbols()
    hubbard_indices = [i for i, sym in enumerate(symbols_list) if sym == symbol]
    other_indices = [i for i, sym in enumerate(symbols_list) if sym != symbol]
    new_order = hubbard_indices + other_indices
    scf_atoms = scf_atoms[new_order]
    
    scf_data = {
        'control': {
            'calculation': 'scf',
            'prefix': prefix,
            'pseudo_dir': PSEUDO_DIR,
            'outdir': str(tmp_dir),
        },
        'system': {
            'ecutwfc': 45.0,
            'ecutrho': 450.0
        },
        'electrons': {
            'conv_thr': 1.0e-7,
            'mixing_beta': 0.5,
            'electron_maxstep': 150
        }
    }
    scf_data['system'].update(get_occupancy_options(symbol))
    
    io.write(
        str(scf_in),
        scf_atoms,
        format='espresso-in',
        input_data=scf_data,
        pseudopotentials={symbol: PSEUDO_MAP[symbol], "N": PSEUDO_MAP["N"]},
        kpts=(5, 5, 1)
    )
    
    with open(str(scf_in), 'a') as f:
        f.write("\nHUBBARD {ortho-atomic}\nU " + symbol + "-" + manifold + " " + str(current_u) + "\n")
    
    return scf_in


def run_relax(iter_dir, prefix, symbol, manifold, current_atoms, current_u):
    """Run relaxation step."""
    relax_in = setup_relax_input(iter_dir, prefix, current_atoms, symbol, manifold, current_u)
    relax_out = iter_dir / (prefix + ".relax.out")
    
    cmd = PW_CMD + " -in " + relax_in.name + " > " + relax_out.name
    result = subprocess.run(cmd, shell=True, cwd=str(iter_dir))
    
    if result.returncode != 0:
        print("ERROR: relax failed")
        return None
    
    relaxed_atoms = parse_relaxed_geometry(str(relax_in), str(relax_out))
    if relaxed_atoms is None:
        print("ERROR: Could not extract final geometry from relax")
        return None
    
    # Reapply constraint to new geometry
    mask = relaxed_atoms.positions[:, 2] < 13.0
    constraint = FixAtoms(mask=mask)
    relaxed_atoms.set_constraint(constraint)
    
    return relaxed_atoms


def run_scf(iter_dir, prefix, symbol, manifold, current_atoms, current_u):
    """Run SCF step."""
    scf_in = setup_scf_input(iter_dir, prefix, current_atoms, symbol, manifold, current_u)
    scf_out = iter_dir / (prefix + ".scf.out")
    
    cmd = PW_CMD + " -in " + scf_in.name + " > " + scf_out.name
    result = subprocess.run(cmd, shell=True, cwd=str(iter_dir))
    
    if result.returncode != 0:
        print("ERROR: SCF failed")
        return False
    
    return True


def run_hp(iter_dir, prefix):
    """Run HP.X step."""
    hp_in = iter_dir / (prefix + ".hp.in")
    hp_out = iter_dir / (prefix + ".hp.out")
    tmp_dir = iter_dir / "tmp"
    
    hp_content = "&inputhp\n"
    hp_content += "   prefix = '" + prefix + "',\n"
    hp_content += "   outdir = '" + str(tmp_dir) + "/',\n"
    hp_content += "   nq1 = 2, nq2 = 2, nq3 = 1,\n"
    hp_content += "   conv_thr_chi = 1.0d-8,\n"
    hp_content += "   iverbosity = 2,\n"
    hp_content += "   compute_hp = .false.,\n"
    hp_content += "/\n"
    
    with open(str(hp_in), 'w') as f:
        f.write(hp_content)
    
    cmd = HP_CMD + " -in " + hp_in.name + " > " + hp_out.name
    result = subprocess.run(cmd, shell=True, cwd=str(iter_dir))
    
    if result.returncode != 0:
        print("ERROR: HP.X failed")
        return None
    
    new_u = extract_u_from_dat(iter_dir, prefix)
    return new_u


def process_metal(metal):
    """Process a single metal compound through the consistency loop."""
    print("\n" + "="*70)
    print("SLAB: " + metal)
    print("="*70)
    
    symbol = get_element_symbol(metal)
    manifold = get_manifold(symbol)
    
    slab_file = Path("./slim") / (metal + "_slim_slab.xyz")
    if not slab_file.exists():
        print("ERROR: Slab file not found: " + str(slab_file))
        return False
    
    # Load slab structure
    try:
        slab_atoms = io.read(str(slab_file))
    except Exception as e:
        print("ERROR: Could not read slab file: " + str(e))
        return False
    
    # Freeze atoms with z < 13 (bulk layer)
    mask = slab_atoms.positions[:, 2] < 13.0
    constraint = FixAtoms(mask=mask)
    slab_atoms.set_constraint(constraint)
    
    n_frozen = sum(mask)
    n_free = len(slab_atoms) - n_frozen
    print("Frozen atoms (z < 13): " + str(n_frozen))
    print("Free atoms:            " + str(n_free))
    
    run_root = Path("./" + metal + "_slab_consistency").resolve()
    run_root.mkdir(exist_ok=True)
    
    current_atoms = slab_atoms
    current_u = START_U_VALUES[metal]
    prefix = metal.lower() + "_slab"
    
    # --- Consistency loop ---
    for i in range(1, 6):
        print("\n" + "="*70)
        print("ITERATION " + str(i) + " | U = {:.4f}".format(current_u))
        print("="*70)
        
        iter_dir = run_root / ("iter_" + str(i))
        iter_dir.mkdir(exist_ok=True)
        tmp_dir = iter_dir / "tmp"
        tmp_dir.mkdir(exist_ok=True)
        
        # ====== STEP 1: RELAX (skip on first iteration) ======
        if i == 1:
            print("[1/3] RELAX: Skipping on first iteration to save time")
        else:
            print("[1/3] RELAX: Relaxing atomic positions with U={:.4f}".format(current_u))
            relaxed_atoms = run_relax(iter_dir, prefix, symbol, manifold, current_atoms, current_u)
            if relaxed_atoms is None:
                break
            current_atoms = relaxed_atoms
            print("  OK: Relaxation complete")
        
        # ====== STEP 2: SCF ======
        print("[2/3] SCF: Self-consistent calculation")
        scf_success = run_scf(iter_dir, prefix, symbol, manifold, current_atoms, current_u)
        if not scf_success:
            break
        print("  OK: SCF converged")
        
        # ====== STEP 3: HP.X ======
        print("[3/3] HP.X: Computing Hubbard U")
        new_u = run_hp(iter_dir, prefix)
        
        if new_u is not None:
            diff = abs(new_u - current_u)
            print("  OK: HP converged")
            print("  Previous U: {:.6f}".format(current_u))
            print("  New U:      {:.6f}".format(new_u))
            print("  Shift:      {:.6f}".format(diff))
            
            if diff < TOLERANCE:
                print("\n" + "="*70)
                print("CONVERGED: " + metal)
                print("Final U = {:.6f}".format(new_u))
                print("="*70)
                break
            
            current_u = new_u
        else:
            print("ERROR: Could not extract U")
            break
        
        if tmp_dir.exists():
            shutil.rmtree(str(tmp_dir))
        
        time.sleep(2)
    
    print()
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Compute Hubbard U parameters for metal nitride slabs using QE.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s TiN VN
  %(prog)s -m TiN -m NbN -m ZrN
  %(prog)s TiN -t 0.005
        """
    )
    
    parser.add_argument(
        'metals',
        nargs='*',
        help='Metal compound(s) to process (e.g., TiN, VN, NbN)'
    )
    
    parser.add_argument(
        '-m', '--metal',
        action='append',
        dest='metal_list',
        help='Metal compound to process (can be used multiple times)'
    )
    
    parser.add_argument(
        '-t', '--tolerance',
        type=float,
        default=TOLERANCE,
        help='Convergence tolerance for U (default: %(default)s)'
    )
    
    parser.add_argument(
        '-p', '--pseudo-dir',
        default=PSEUDO_DIR,
        help='Path to pseudopotential directory (default: %(default)s)'
    )
    
    parser.add_argument(
        '-u', '--start-u',
        type=float,
        help='Override starting U value for all metals (use carefully)'
    )
    
    args = parser.parse_args()
    
    # Combine positional and optional metal arguments
    metals_to_process = args.metals.copy() if args.metals else []
    if args.metal_list:
        metals_to_process.extend(args.metal_list)
    
    if not metals_to_process:
        parser.print_help()
        print("\nERROR: At least one metal compound must be specified")
        return 1
    
    # Validate metals
    valid_metals = set(START_U_VALUES.keys())
    for metal in metals_to_process:
        if metal not in valid_metals:
            print(f"ERROR: Unknown metal '{metal}'. Valid options: {', '.join(sorted(valid_metals))}")
            return 1
    
    # Update global settings if provided
    global TOLERANCE
    TOLERANCE = args.tolerance
    
    if args.start_u is not None:
        for metal in metals_to_process:
            START_U_VALUES[metal] = args.start_u
    
    # Process each metal
    print("\n" + "="*70)
    print("SLAB CONSISTENCY WORKFLOW")
    print("="*70)
    print(f"Metals to process: {', '.join(metals_to_process)}")
    print(f"Convergence tolerance: {TOLERANCE}")
    print()
    
    for metal in metals_to_process:
        success = process_metal(metal)
        if not success:
            print(f"WARNING: Processing of {metal} did not complete successfully")
    
    print("\n" + "="*70)
    print("SLAB WORKFLOW COMPLETE")
    print("="*70)
    
    return 0


if __name__ == "__main__":
    exit(main())