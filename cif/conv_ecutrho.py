#!/usr/bin/env python3
"""
QE ecutrho Convergence Test (ASE + Espresso)

Sweep ecutrho = multiplier * ecutwfc for a fixed ecutwfc (45 Ry).
Saves per-run directories and a summary plot (PNG) in `conv_ecutrho/`.

Usage: adjust MATERIALS, pseudos or run with MPI (recommended):
    mpirun -n 16 python conv_ecutrho.py

Notes:
- All cutoffs are in Rydberg (Ry).
- This follows the style of `conv_ecut.py` in this repo.
"""

import os
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from ase.io import read
from ase.calculators.espresso import Espresso, EspressoProfile

try:
    from mpi4py import MPI
    rank = MPI.COMM_WORLD.Get_rank()
except Exception:
    rank = 0
    os.environ.setdefault('OMP_NUM_THREADS', '1')

# --- User settings ---
ROOT = Path(__file__).resolve().parent
output_dir = ROOT / 'conv_ecutrho'
output_dir.mkdir(exist_ok=True)

# Materials to test (adjust as needed)
MATERIALS = ['TiN', 'NbN', 'ScN', 'VN', 'ZrN']

# Fixed ecutwfc (Ry) and ecutrho multipliers to test
ECUTWFC_RY = 45.0
ECUTRHO_MULTIPLIERS = [5, 6, 7, 8]
CUTOFFRHO_LIST_RY = [int(m * ECUTWFC_RY) for m in ECUTRHO_MULTIPLIERS]

# Use local USPP directory in this workspace (adjust if running on cluster)
pseudo_dir = ROOT / 'USPP'

# Map element -> pseudo filename (keep in sync with your repo's USPP/)
pseudos = {
    'Ti': 'ti_pbe_v1.4.uspp.F.UPF',
    'V':  'v_pbe_v1.4.uspp.F.UPF',
    'Sc': 'sc_pbe_v1.uspp.F.UPF',
    'Nb': 'nb_pbe_v1.uspp.F.UPF',
    'Zr': 'zr_pbe_v1.uspp.F.UPF',
    'N':  'n_pbe_v1.2.uspp.F.UPF',
}

# Keep launcher consistent with conv_ecut.py (adjust for your cluster)
profile = EspressoProfile(command='mpirun -n 16 --oversubscribe env OMP_NUM_THREADS=1 pw.x', pseudo_dir=str(pseudo_dir))

# k-point grid used for SCF runs
kpts = (5, 5, 5)


def find_structure(material):
    search_dir = ROOT / 'CIF' / 'input_cif'
    if not search_dir.exists():
        return None, None

    for p in sorted(search_dir.glob(f'{material}*')):
        if not p.is_file():
            continue
        try:
            atoms = read(str(p))
            return atoms, p
        except Exception:
            continue

    return None, None


def test_ecutrho_convergence(material):
    if rank == 0:
        print('\n' + '='*60)
        print(f'Testing ecutrho convergence for {material} (ecutwfc={ECUTWFC_RY} Ry)')
        print('='*60)

    atoms, path = find_structure(material)
    if atoms is None:
        if rank == 0:
            print(f'✗ Could not find a structure for {material}')
        return None

    if rank == 0:
        print(f'✓ Loaded structure from {path}')

    energies = []
    cutoffs = []
    prev_energy = None

    for ecutrho in CUTOFFRHO_LIST_RY:
        if rank == 0:
            print(f'\nCalculating energy for ecutrho = {ecutrho} Ry (ecutwfc = {ECUTWFC_RY} Ry)')

        run_dir = output_dir / f'{material}_rho_{int(ecutrho)}Ry'
        run_dir.mkdir(exist_ok=True)

        input_data = {
            'control': {
                'calculation': 'scf',
                'prefix': f'{material}_rho_{int(ecutrho)}',
                'pseudo_dir': str(pseudo_dir),
                'verbosity': 'high',
                'disk_io': 'low',
            },
            'system': {
                'ecutwfc': float(ECUTWFC_RY),
                'ecutrho': float(ecutrho),
            },
            'electrons': {
                'conv_thr': 1e-6,
            },
            'input_dft': 'PBE',
            'vdw_corr': 'dft-d3',
            'dftd3_version': 4,
        }

        calc = Espresso(
            pseudopotentials=pseudos,
            pseudo_dir=str(pseudo_dir),
            input_data=input_data,
            kpts=kpts,
            profile=profile,
            directory=str(run_dir),
            label=f'{material}_rho_{int(ecutrho)}',
        )

        atoms.calc = calc

        energy = atoms.get_potential_energy()
        energies.append(energy)
        cutoffs.append(ecutrho)

        if rank == 0:
            if prev_energy is not None:
                delta_total = energy - prev_energy
                delta_per_atom = delta_total / len(atoms)
                print(f'  Energy: {energy:.6f} eV | ΔE: {delta_total:.6f} eV | ΔE/atom: {delta_per_atom:.6f} eV/atom')
            else:
                print(f'  Energy: {energy:.6f} eV')

        prev_energy = energy

    if rank == 0:
        print(f'\n✓ ecutrho convergence test completed for {material}')

    return {'material': material, 'ecutrho': cutoffs, 'energies': energies}


def main():
    if rank == 0:
        print('\n' + '='*60)
        print('QE ecutrho Convergence Tests')
        print(f'Fixed ecutwfc = {ECUTWFC_RY} Ry; testing ecutrho multipliers: {ECUTRHO_MULTIPLIERS}')
        print('='*60)

    all_results = {}
    for material in MATERIALS:
        try:
            res = test_ecutrho_convergence(material)
            if res:
                all_results[material] = res
        except Exception as e:
            if rank == 0:
                print(f'✗ Error during {material} ecutrho test: {e}')

    if rank == 0 and all_results:
        plt.figure(figsize=(10, 6))
        for material, results in all_results.items():
            plt.plot(results['ecutrho'], results['energies'], marker='o', label=material, linewidth=2)
        plt.xlabel('ecutrho (Ry)')
        plt.ylabel('Total Energy (eV)')
        plt.title(f'ecutrho Convergence (ecutwfc={ECUTWFC_RY} Ry, PBE+D3)')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        out_png = output_dir / 'qe_ecutrho_convergence.png'
        plt.savefig(out_png, dpi=300)
        print(f'✓ Plot saved to {out_png}')


if __name__ == '__main__':
    main()
