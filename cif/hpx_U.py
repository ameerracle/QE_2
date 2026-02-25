import os
import re
from pathlib import Path
from ase import io

# --- Setup ---
structure_path = Path("/lustre07/scratch/anizami/QE_2/cif/VN_vcrelax/VN_relaxed.extxyz")
pseudo_dir = "/scratch/anizami/QE_2/USPP"
run_root = Path("./VN_self_consistent_U")
run_root.mkdir(exist_ok=True)

qe_path = "/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/MPI/gcc12/openmpi4/quantumespresso/7.5/bin/"
pw_cmd = f"srun {qe_path}pw.x"
hp_cmd = f"srun {qe_path}hp.x"

pseudos = {"V": "v_pbe_v1.4.uspp.F.UPF", "N": "n_pbe_v1.2.uspp.F.UPF"}
current_u = 0.0

# Parse structure once with ASE
atoms = io.read(structure_path)

def write_pw_input(filename, atoms, input_data, pseudos, kpts, hubbard_u):
    """Manually write QE 7.5 input file with new HUBBARD card syntax"""
    
    # Get cell and atomic info
    cell = atoms.cell
    pos = atoms.get_positions()
    symbols = atoms.get_chemical_symbols()
    
    # Build input string
    input_str = "&control\n"
    for key, val in input_data['control'].items():
        if isinstance(val, str):
            input_str += f"  {key} = '{val}'\n"
        else:
            input_str += f"  {key} = {val}\n"
    input_str += "/\n\n"
    
    # System card - NO lda_plus_u or lda_plus_u_kind for QE 7.5!
    input_str += "&system\n"
    input_str += f"  ibrav = 0\n"
    input_str += f"  nat = {len(atoms)}\n"
    input_str += f"  ntyp = {len(set(symbols))}\n"
    for key, val in input_data['system'].items():
        if isinstance(val, str):
            input_str += f"  {key} = '{val}'\n"
        else:
            input_str += f"  {key} = {val}\n"
    input_str += "/\n\n"
    
    # Electrons card
    input_str += "&electrons\n"
    for key, val in input_data['electrons'].items():
        if isinstance(val, str):
            input_str += f"  {key} = '{val}'\n"
        else:
            input_str += f"  {key} = {val}\n"
    input_str += "/\n\n"
    
    # Ions card (only if vc-relax)
    if input_data['control']['calculation'] == 'vc-relax':
        input_str += "&ions\n"
        for key, val in input_data['ions'].items():
            if isinstance(val, str):
                input_str += f"  {key} = '{val}'\n"
            else:
                input_str += f"  {key} = {val}\n"
        input_str += "/\n\n"
        
        input_str += "&cell\n"
        for key, val in input_data['cell'].items():
            if isinstance(val, str):
                input_str += f"  {key} = '{val}'\n"
            else:
                input_str += f"  {key} = {val}\n"
        input_str += "/\n\n"
    
    # Atomic species
    input_str += "ATOMIC_SPECIES\n"
    unique_symbols = sorted(set(symbols))
    for symbol in unique_symbols:
        atomic_mass = {"V": 50.9415, "N": 14.0067}
        input_str += f"  {symbol}  {atomic_mass[symbol]:.4f}  {pseudos[symbol]}\n"
    input_str += "\n"
    
    # Cell parameters
    input_str += "CELL_PARAMETERS (angstrom)\n"
    for i in range(3):
        input_str += f"  {cell[i, 0]:15.10f}  {cell[i, 1]:15.10f}  {cell[i, 2]:15.10f}\n"
    input_str += "\n"
    
    # Atomic positions
    input_str += "ATOMIC_POSITIONS (angstrom)\n"
    for symbol, p in zip(symbols, pos):
        input_str += f"  {symbol}  {p[0]:15.10f}  {p[1]:15.10f}  {p[2]:15.10f}\n"
    input_str += "\n"
    
    # K-points
    input_str += "K_POINTS (automatic)\n"
    input_str += f"  {kpts[0]}  {kpts[1]}  {kpts[2]}  0  0  0\n"
    input_str += "\n"
    
    # NEW QE 7.5 HUBBARD SYNTAX (ortho-atomic)
    # Format: U element-manifold value
    input_str += "HUBBARD (ortho-atomic)\n"
    input_str += f"  U V-3d {hubbard_u:.4f}\n"
    
    with open(filename, 'w') as f:
        f.write(input_str)


def extract_hubbard_u(hp_output_file):
    """Extract converged Hubbard U value from hp.out"""
    try:
        with open(hp_output_file, 'r') as f:
            lines = f.readlines()
        
        # Look for the Hubbard U value in hp output
        # QE hp.x writes something like: "Hubbard U value = X.XXXX"
        for line in reversed(lines):  # Start from end, most recent values last
            if 'Hubbard' in line and 'U' in line and '=' in line:
                try:
                    # Extract number after =
                    val_str = line.split('=')[-1].strip().split()[0]
                    return float(val_str)
                except (ValueError, IndexError):
                    continue
        
        print(f"Warning: Could not extract U value from {hp_output_file}")
        return None
        
    except FileNotFoundError:
        print(f"Error: {hp_output_file} not found")
        return None


# Main iteration loop
for i in range(5):
    iter_dir = run_root / f"iteration_{i}"
    iter_dir.mkdir(exist_ok=True)
    os.chdir(iter_dir)
    
    print(f"\n{'='*60}")
    print(f"Iteration {i}: U = {current_u:.4f}")
    print(f"{'='*60}\n")
    
    # Input data template (NO lda_plus_u!)
    input_data = {
        'control': {'calculation': 'vc-relax', 'prefix': 'vn_calc', 'pseudo_dir': pseudo_dir},
        'system': {
            'ecutwfc': 45.0, 'ecutrho': 360.0,
            'occupations': 'smearing', 'smearing': 'cold', 'degauss': 0.005,
        },
        'electrons': {'conv_thr': 1e-7, 'mixing_beta': 0.35},
        'ions': {'ion_dynamics': 'bfgs'},
        'cell': {'cell_dynamics': 'bfgs'}
    }
    
    # 1. VC-Relax with PW.X
    print("Running pw.x (vc-relax)...")
    write_pw_input('pw.in', atoms, input_data, pseudos, (4, 4, 4), current_u)
    ret = os.system(f"{pw_cmd} < pw.in > pw.out 2>&1")
    if ret != 0:
        print(f"ERROR: pw.x failed with code {ret}")
        break
    
    # Update structure from relaxation
    xml_file = 'vn_calc.save/data-file-schema.xml'
    try:
        atoms = io.read(xml_file)
        print(f"? Loaded relaxed structure from {xml_file}")
    except FileNotFoundError:
        print(f"ERROR: Could not find {xml_file}")
        print("Check vn_calc.save/ directory contents")
        break
    
    # 2. SCF calculation
    print("Running pw.x (scf)...")
    input_data['control']['calculation'] = 'scf'
    write_pw_input('scf.in', atoms, input_data, pseudos, (4, 4, 4), current_u)
    ret = os.system(f"{pw_cmd} < scf.in > scf.out 2>&1")
    if ret != 0:
        print(f"ERROR: pw.x (scf) failed with code {ret}")
        break
    
    # 3. HP.X calculation
    print("Running hp.x...")
    hp_input = f"""&inputhp
    prefix='vn_calc'
    outdir='./'
    nq1=2, nq2=2, nq3=2
    conv_thr_chi=1.0d-8
    iverbosity=2
/
"""
    with open('hp.in', 'w') as f:
        f.write(hp_input)
    
    ret = os.system(f"{hp_cmd} < hp.in > hp.out 2>&1")
    if ret != 0:
        print(f"ERROR: hp.x failed with code {ret}")
        break
    
    # 4. Extract new U and check convergence
    print("Parsing hp.out...")
    new_u = extract_hubbard_u('hp.out')
    
    if new_u is not None:
        print(f"Previous U: {current_u:.4f}")
        print(f"Calculated U: {new_u:.4f}")
        print(f"Change: {abs(new_u - current_u):.4f}\n")
        
        if abs(new_u - current_u) < 0.01:  # Convergence threshold
            print("? Hubbard U converged!")
            current_u = new_u
            break
        
        current_u = new_u
    else:
        print("Warning: Could not extract U value from hp.out")
        print("Check hp.out manually for errors")
        break
    
    os.chdir("../../")

print(f"\nFinal U value: {current_u:.4f}")