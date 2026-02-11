#!/usr/bin/env python3
"""
QE Energy Cutoff Convergence Test (ASE + Espresso)

This script performs energy cutoff convergence tests for a list of materials
using Quantum ESPRESSO (`pw.x`) via ASE's `Espresso` calculator. It reads
structures from `input_cif` (or `input_cif_prim`) and sweeps a set of
energy cutoffs (provided in eV; converted to Ry for QE).

Usage: adjust `MATERIALS`, `CUTOFF_RANGE_EEV`, `pseudo_dir` and `pseudos`.
Run with MPI (recommended): `mpirun -n 16 python conv_check_qe.py`
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

# Force single OpenMP thread per MPI rank to avoid oversubscription issues


# --- User settings ---
ROOT = Path(__file__).resolve().parent
output_dir = ROOT / 'conv'
output_dir.mkdir(exist_ok=True)

# Materials to test
MATERIALS = ['TiN',
    'NbN', 'ScN',  'VN', 'ZrN'
    ]
# Energy cutoffs in Rydberg (Ry)
# Sweep from 35 Ry to 55 Ry inclusive in steps of 5 Ry
CUTOFF_RANGE_RY = np.arange(35, 56, 5)

# Use local USPP directory in this workspace
pseudo_dir = ROOT / 'USPP'

# Map element -> pseudo filename (adjust to your pseudo names)
pseudos = {
    'Ti': 'ti_pbe_v1.4.uspp.F.UPF',
    'V':  'v_pbe_v1.4.uspp.F.UPF',
    'Sc': 'sc_pbe_v1.uspp.F.UPF',
    'Nb': 'nb_pbe_v1.uspp.F.UPF',
    'Zr': 'zr_pbe_v1.uspp.F.UPF',
    'N':  'n_pbe_v1.2.uspp.F.UPF',
}

# Use oversubscribe and enforce single-threaded ranks to avoid MPI/OpenMP conflicts
profile = EspressoProfile(command='mpirun -n 16 --oversubscribe env OMP_NUM_THREADS=1 pw.x', pseudo_dir=str(pseudo_dir))

# k-point grid used for SCF runs
kpts = (5, 5, 5)


def find_structure(material):
    # Restrict search to CIF/input_cif only
    search_dir = ROOT / 'CIF' / 'input_cif'
    if not search_dir.exists():
        return None, None

    # look for any file in CIF/input_cif starting with the material name
    for p in sorted(search_dir.glob(f'{material}*')):
        if not p.is_file():
            continue
        try:
            atoms = read(str(p))
            return atoms, p
        except Exception:
            continue

    return None, None


def test_cutoff_convergence(material):
    if rank == 0:
        print('\n' + '='*60)
        print(f'Testing cutoff convergence for {material}')
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

    for ecutwfc in CUTOFF_RANGE_RY:
        if rank == 0:
            print(f'\nCalculating energy for cutoff = {ecutwfc} Ry')

        # create a per-run output directory to keep pw.x files separate
        run_dir = output_dir / f'{material}_cutoff_{int(ecutwfc)}Ry'
        run_dir.mkdir(exist_ok=True)

        input_data = {
            'control': {
                'calculation': 'scf',
                'prefix': f'{material}_cut_{int(ecutwfc)}',
                'pseudo_dir': str(pseudo_dir),
                'verbosity': 'high',
                'disk_io': 'low',
            },
            'system': {
                'ecutwfc': float(ecutwfc),
                # ecutrho will be set to 4 * ecutwfc (in Ry) unless your pseudos
                # require a different multiplier.
                'ecutrho': 5 * float(ecutwfc),
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
            label=f'{material}_cut_{int(ecutwfc)}',
        )

        atoms.calc = calc

        # run single-point SCF and collect energy
        energy = atoms.get_potential_energy()
        energies.append(energy)
        cutoffs.append(ecutwfc)

        if rank == 0:
            if prev_energy is not None:
                delta_total = energy - prev_energy
                delta_per_atom = delta_total / len(atoms)
                print(f'  Energy: {energy:.6f} eV | ΔE: {delta_total:.6f} eV | ΔE/atom: {delta_per_atom:.6f} eV/atom')
            else:
                print(f'  Energy: {energy:.6f} eV')

        prev_energy = energy

    if rank == 0:
        print(f'\n✓ Convergence test completed for {material}')

    return {'material': material, 'cutoffs': cutoffs, 'energies': energies}


def main():
    if rank == 0:
        print('\n' + '='*60)
        print('QE Energy Cutoff Convergence Tests')
        print('Method: PBE + D3 (via QE)')
        print('='*60)

    all_results = {}
    for material in MATERIALS:
        try:
            res = test_cutoff_convergence(material)
            if res:
                all_results[material] = res
        except Exception as e:
            if rank == 0:
                print(f'✗ Error during {material} convergence test: {e}')

    # plotting (rank 0)
    if rank == 0 and all_results:
        plt.figure(figsize=(10, 6))
        for material, results in all_results.items():
            plt.plot(results['cutoffs'], results['energies'], marker='o', label=material, linewidth=2)
        plt.xlabel('Energy Cutoff (eV)')
        plt.ylabel('Total Energy (eV)')
        plt.title('Energy Cutoff Convergence (QE RPBE+D3)')
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        out_png = output_dir / 'qe_cutoff_convergence.png'
        plt.savefig(out_png, dpi=300)
        print(f'✓ Plot saved to {out_png}')


if __name__ == '__main__':
    main()
