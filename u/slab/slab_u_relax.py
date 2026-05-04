import os
import subprocess
import re
import shutil
import time
from pathlib import Path
from ase import io
from ase import Atoms
from ase.constraints import FixAtoms
import numpy as np

# --- Settings ---
metals = ["ZrN"]
pseudo_dir = "/scratch/anizami/QE_2/USPP"
qe_bin = "/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/MPI/gcc12/openmpi4/quantumespresso/7.5/bin/"
pw_cmd = "srun " + qe_bin + "pw.x"
hp_cmd = "srun " + qe_bin + "hp.x"

tolerance = 0.01

start_u_values = {
    "VN": 5.87,
    "TiN": 5.55,
    "NbN": 3.52,
    "ZrN": 2.97,
    "ScN": 3.57
}

pseudo_map = {
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
    if metal[:2] in ["Ti", "Sc", "Nb", "Zr"]:
        return metal[:2]
    else:
        return metal[0]


def get_manifold(symbol):
    return "4d" if symbol in ["Nb", "Zr"] else "3d"


# --- Main ---
for metal in metals:
    print("\n" + "="*70)
    print("SLAB: " + metal)
    print("="*70)
    
    symbol = get_element_symbol(metal)
    manifold = get_manifold(symbol)
    
    slab_file = Path("./slim") / (metal + "_slim_slab.xyz")
    if not slab_file.exists():
        print("ERROR: Slab file not found: " + str(slab_file))
        continue
    
    # Load slab structure
    try:
        slab_atoms = io.read(str(slab_file))
    except Exception as e:
        print("ERROR: Could not read slab file: " + str(e))
        continue
    
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
    current_u = start_u_values[metal]
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
        
        relax_in = iter_dir / (prefix + ".relax.in")
        relax_out = iter_dir / (prefix + ".relax.out")
        scf_in = iter_dir / (prefix + ".scf.in")
        scf_out = iter_dir / (prefix + ".scf.out")
        hp_in = iter_dir / (prefix + ".hp.in")
        hp_out = iter_dir / (prefix + ".hp.out")
        
        # ====== STEP 1: RELAX ======
        if relax_out.exists():
            print("[1/3] RELAX: Already completed, reading geometry...")
            relaxed_atoms = parse_relaxed_geometry(str(relax_in), str(relax_out))
            if relaxed_atoms is None:
                print("ERROR: Could not parse existing relax output")
                break
            mask = relaxed_atoms.positions[:, 2] < 13.0
            constraint = FixAtoms(mask=mask)
            relaxed_atoms.set_constraint(constraint)
            current_atoms = relaxed_atoms
            print("  OK: Geometry loaded from previous run")
        else:
            print("[1/3] RELAX: Relaxing atomic positions with U={:.4f}".format(current_u))
        
            # ScN is a semiconductor/insulator with bandgap, others are metallic
            if symbol == "Sc":
                occ_opts = {'occupations': 'fixed'}
            else:
                occ_opts = {
                    'occupations': 'smearing',
                    'smearing': 'cold',
                    'degauss': 0.005
                }
            
            relax_data = {
                'control': {
                    'calculation': 'relax',
                    'prefix': prefix,
                    'pseudo_dir': pseudo_dir,
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
            relax_data['system'].update(occ_opts)
            
            io.write(
                str(relax_in),
                current_atoms,
                format='espresso-in',
                input_data=relax_data,
                pseudopotentials={symbol: pseudo_map[symbol], "N": pseudo_map["N"]},
                kpts=(5, 5, 1)
            )
            
            with open(str(relax_in), 'a') as f:
                f.write("\nHUBBARD {ortho-atomic}\nU " + symbol + "-" + manifold + " " + str(current_u) + "\n")
            
            cmd = pw_cmd + " -in " + relax_in.name + " > " + relax_out.name
            result = subprocess.run(cmd, shell=True, cwd=str(iter_dir))
            
            if result.returncode != 0:
                print("ERROR: relax failed")
                break
            
            relaxed_atoms = parse_relaxed_geometry(str(relax_in), str(relax_out))
            if relaxed_atoms is None:
                print("ERROR: Could not extract final geometry from relax")
                break
            
            # Reapply constraint to new geometry
            mask = relaxed_atoms.positions[:, 2] < 13.0
            constraint = FixAtoms(mask=mask)
            relaxed_atoms.set_constraint(constraint)
            
            current_atoms = relaxed_atoms
            print("  OK: Relaxation complete")
        
        # ====== STEP 2: SCF (ALWAYS RUN) ======
        print("[2/3] SCF: Self-consistent calculation")
        
        if scf_out.exists():
            print("  OK: SCF already completed, skipping...")
        else:
            # Reorder atoms so Hubbard atoms (Ti) come first
            scf_atoms = current_atoms.copy()
            symbols_list = scf_atoms.get_chemical_symbols()
            ti_indices = [i for i, sym in enumerate(symbols_list) if sym == symbol]
            n_indices = [i for i, sym in enumerate(symbols_list) if sym != symbol]
            new_order = ti_indices + n_indices
            scf_atoms = scf_atoms[new_order]
            
            scf_data = {
                'control': {
                    'calculation': 'scf',
                    'prefix': prefix,
                    'pseudo_dir': pseudo_dir,
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
            
            # ScN is a semiconductor/insulator with bandgap, others are metallic
            if symbol == "Sc":
                occ_opts = {'occupations': 'fixed'}
            else:
                occ_opts = {
                    'occupations': 'smearing',
                    'smearing': 'cold',
                    'degauss': 0.005
                }
            scf_data['system'].update(occ_opts)
            
            io.write(
                str(scf_in),
                scf_atoms,
                format='espresso-in',
                input_data=scf_data,
                pseudopotentials={symbol: pseudo_map[symbol], "N": pseudo_map["N"]},
                kpts=(5, 5, 1)
            )
            
            with open(str(scf_in), 'a') as f:
                f.write("\nHUBBARD {ortho-atomic}\nU " + symbol + "-" + manifold + " " + str(current_u) + "\n")
            
            cmd = pw_cmd + " -in " + scf_in.name + " > " + scf_out.name
            result = subprocess.run(cmd, shell=True, cwd=str(iter_dir))
            
            if result.returncode != 0:
                print("ERROR: SCF failed")
                break
            
            print("  OK: SCF converged")
        
        # ====== STEP 3: HP.X (ALWAYS RUN) ======
        print("[3/3] HP.X: Computing Hubbard U")
        
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
        
        cmd = hp_cmd + " -in " + hp_in.name + " > " + hp_out.name
        result = subprocess.run(cmd, shell=True, cwd=str(iter_dir))
        
        if result.returncode != 0:
            print("ERROR: HP.X failed")
            break
        
        new_u = extract_u_from_dat(iter_dir, prefix)
        
        if new_u is not None:
            diff = abs(new_u - current_u)
            print("  OK: HP converged")
            print("  Previous U: {:.6f}".format(current_u))
            print("  New U:      {:.6f}".format(new_u))
            print("  Shift:      {:.6f}".format(diff))
            
            if diff < tolerance:
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

print("\n" + "="*70)
print("SLAB WORKFLOW COMPLETE")
print("="*70)