"""
QE hp.x bulk Hubbard U calculation for transition-metal nitrides.
Paired counterpart to gpaw_U_clda.py — both target the same primitive cells
so their U values can be directly compared.

Workflow per compound:
  1. SCF (pw.x) on primitive cell — Hubbard block with U=0 first pass
  2. hp.x linear-response calculation on the SCF charge density
  3. Extract U from .Hubbard_parameters.dat

Using USPP (USPP/ dir) — required for hp.x (ONCV/NC PPs lack <PP_PSWFC>
atomic wavefunctions that hp.x uses to define the Hubbard projector manifold).

Usage:
    python qe_U_bulk.py                    # run all compounds
    python qe_U_bulk.py TiN VN             # run specific compounds
    python qe_U_bulk.py --compare-only     # print comparison table only
"""

import os
import re
import subprocess
import argparse
import numpy as np
from pathlib import Path
from ase.io import read

# ============================================================
# CONFIGURATION
# ============================================================

# Resolved relative to this script so it works on any machine
_HERE       = Path(__file__).parent
POSCAR_DIR  = _HERE / 'primitive'
PSEUDO_DIR  = _HERE.parent / 'USPP'       # USPP — required for hp.x
WORK_DIR    = _HERE / 'qe_bulk_U'         # output tree

# QE binaries — override with env var QE_BIN or edit here
QE_BIN  = os.environ.get('QE_BIN', '')
PW_CMD  = (QE_BIN + '/pw.x') if QE_BIN else 'pw.x'
HP_CMD  = (QE_BIN + '/hp.x') if QE_BIN else 'hp.x'

# Compound database  (same as gpaw_U_clda.py)
COMPOUND_DB = {
    'NbN': ('NbN_prim.poscar', 'Nb'),
    'ScN': ('ScN_prim.poscar', 'Sc'),
    'TiN': ('TiN_prim.poscar', 'Ti'),
    'VN':  ('VN_prim.poscar',  'V'),
    'ZrN': ('ZrN_prim.poscar', 'Zr'),
}

# USPP pseudopotential filenames  (these include <PP_PSWFC>, required by hp.x)
PSEUDO_MAP = {
    'Ti': 'ti_pbe_v1.4.uspp.F.UPF',
    'V':  'v_pbe_v1.4.uspp.F.UPF',
    'Sc': 'sc_pbe_v1.uspp.F.UPF',
    'Nb': 'nb_pbe_v1.uspp.F.UPF',
    'Zr': 'zr_pbe_v1.uspp.F.UPF',
    'N':  'n_pbe_v1.2.uspp.F.UPF',
}

# DFT settings
ECUTWFC  = 45.0    # Ry
ECUTRHO  = 450.0   # Ry
KPTS     = (8, 8, 8)   # primitive-cell k-mesh (= 4x4x4 for 2x2x2 supercell)
QPTS     = (4, 4, 4)   # hp.x q-mesh
CONV_SCF = 1.0e-10     # tighter than slab: bulk has no vacuum
CONV_HP  = 1.0e-8      # conv_thr_chi

# Orbital manifold map
MANIFOLD = {
    'Nb': '4d', 'Zr': '4d',
    'Ti': '3d', 'V':  '3d', 'Sc': '3d',
}


# ============================================================
# HELPERS
# ============================================================

def occupancy_block(metal_sym):
    """Return dict of smearing/occupations settings for this metal."""
    if metal_sym == 'Sc':
        # ScN is a semiconductor — no smearing
        return {'occupations': "'fixed'"}
    return {
        'occupations': "'smearing'",
        'smearing':    "'cold'",
        'degauss':     0.005,
    }


def write_scf_input(scf_in, atoms, metal_sym, prefix, tmp_dir):
    """Write a pw.x SCF input file from an ASE Atoms object."""
    cell = atoms.get_cell()          # Angstrom
    positions = atoms.get_positions()
    symbols   = atoms.get_chemical_symbols()

    # Reorder: metal atoms first (required so Hubbard atom index = 1 in QE)
    order = ([i for i, s in enumerate(symbols) if s == metal_sym] +
             [i for i, s in enumerate(symbols) if s != metal_sym])
    symbols   = [symbols[i]   for i in order]
    positions = positions[order]

    unique_syms = list(dict.fromkeys(symbols))   # preserve order, deduplicate
    nat   = len(symbols)
    ntyp  = len(unique_syms)

    occ = occupancy_block(metal_sym)
    occ_lines = '\n'.join(f'   {k} = {v},' for k, v in occ.items())

    cell_bohr = cell * 1.8897259886    # Å → Bohr
    pos_bohr  = positions * 1.8897259886

    lines = [
        "&CONTROL",
        f"   calculation = 'scf',",
        f"   prefix      = '{prefix}',",
        f"   pseudo_dir  = '{PSEUDO_DIR}',",
        f"   outdir      = '{tmp_dir}/',",
        f"   tprnfor     = .true.,",
        "/",
        "&SYSTEM",
        f"   ibrav       = 0,",
        f"   nat         = {nat},",
        f"   ntyp        = {ntyp},",
        f"   ecutwfc     = {ECUTWFC},",
        f"   ecutrho     = {ECUTRHO},",
        f"   lda_plus_u  = .false.,",
        occ_lines,
        "/",
        "&ELECTRONS",
        f"   conv_thr        = {CONV_SCF},",
        f"   mixing_beta     = 0.4,",
        f"   electron_maxstep = 200,",
        "/",
        "",
        "ATOMIC_SPECIES",
    ]

    for sym in unique_syms:
        # Approximate mass (doesn't affect SCF)
        mass_map = {'Ti': 47.87, 'V': 50.94, 'Sc': 44.96,
                    'Nb': 92.91, 'Zr': 91.22, 'N': 14.01}
        lines.append(f"   {sym}  {mass_map.get(sym, 1.0)}  {PSEUDO_MAP[sym]}")

    lines += [
        "",
        "CELL_PARAMETERS {bohr}",
    ]
    for row in cell_bohr:
        lines.append(f"   {row[0]:.12f}  {row[1]:.12f}  {row[2]:.12f}")

    lines += [
        "",
        "ATOMIC_POSITIONS {bohr}",
    ]
    for sym, pos in zip(symbols, pos_bohr):
        lines.append(f"   {sym}  {pos[0]:.12f}  {pos[1]:.12f}  {pos[2]:.12f}")

    lines += [
        "",
        "K_POINTS {automatic}",
        f"   {KPTS[0]} {KPTS[1]} {KPTS[2]}  0 0 0",
    ]

    scf_in.write_text('\n'.join(lines) + '\n')


def write_hp_input(hp_in, prefix, tmp_dir):
    """Write an hp.x input file."""
    content = (
        "&inputhp\n"
        f"   prefix       = '{prefix}',\n"
        f"   outdir       = '{tmp_dir}/',\n"
        f"   nq1 = {QPTS[0]}, nq2 = {QPTS[1]}, nq3 = {QPTS[2]},\n"
        f"   conv_thr_chi = {CONV_HP},\n"
        f"   iverbosity   = 2,\n"
        "/\n"
    )
    hp_in.write_text(content)


def extract_u_from_dat(folder, prefix):
    """Parse Hubbard U from the .Hubbard_parameters.dat file."""
    dat = folder / f"{prefix}.Hubbard_parameters.dat"
    if not dat.exists():
        return None
    text = dat.read_text()
    # Matches lines like:  Ti  3d  0.000  5.5479
    matches = re.findall(r'[34]d\s+([\d.]+)', text)
    return float(matches[-1]) if matches else None


def run_cmd(cmd, cwd, log):
    """Run a shell command, stream stdout to log file."""
    print(f"    $ {cmd}")
    with open(log, 'w') as flog:
        result = subprocess.run(
            cmd, shell=True, cwd=str(cwd),
            stdout=flog, stderr=subprocess.STDOUT
        )
    return result.returncode == 0


# ============================================================
# PER-COMPOUND CALCULATION
# ============================================================

def run_compound_qe(compound):
    """Run SCF + hp.x for one compound. Returns U in eV or None."""
    poscar_file, metal = COMPOUND_DB[compound]
    manifold = MANIFOLD[metal]
    prefix   = compound.lower() + '_bulk'

    comp_dir = WORK_DIR / compound
    tmp_dir  = comp_dir / 'tmp'
    comp_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir.mkdir(exist_ok=True)

    print()
    print('#' * 60)
    print(f'#  {compound}  (metal={metal}, manifold={manifold})')
    print('#' * 60)

    # Read primitive cell POSCAR
    poscar_path = POSCAR_DIR / poscar_file
    atoms = read(str(poscar_path))

    # --------------------------------------------------------
    # STEP 1 — SCF
    # --------------------------------------------------------
    print('  [1/2] SCF  (pw.x)')
    scf_in  = comp_dir / f'{prefix}.scf.in'
    scf_log = comp_dir / f'{prefix}.scf.out'

    write_scf_input(scf_in, atoms, metal, prefix, tmp_dir)

    ok = run_cmd(f'pw.x -in {scf_in.name} > {scf_log.name}', comp_dir, scf_log)
    if not ok:
        print(f'  ERROR: pw.x failed — see {scf_log}')
        return None
    print(f'  OK: SCF done  -> {scf_log.name}')

    # --------------------------------------------------------
    # STEP 2 — hp.x
    # --------------------------------------------------------
    print('  [2/2] Linear response  (hp.x)')
    hp_in  = comp_dir / f'{prefix}.hp.in'
    hp_log = comp_dir / f'{prefix}.hp.out'

    write_hp_input(hp_in, prefix, tmp_dir)

    ok = run_cmd(f'hp.x -in {hp_in.name} > {hp_log.name}', comp_dir, hp_log)
    if not ok:
        print(f'  ERROR: hp.x failed — see {hp_log}')
        return None

    U = extract_u_from_dat(comp_dir, prefix)
    if U is None:
        print(f'  ERROR: Could not parse U from .dat file')
        return None

    print(f'  OK: U({metal} {manifold}) = {U:.4f} eV')
    return U


# ============================================================
# COMPARISON TABLE
# ============================================================

def print_comparison(results_qe):
    """
    Load GPAW cLDA results (.npy files) and print a side-by-side table.
    results_qe: dict compound -> U (eV) or None
    """
    print()
    print('=' * 70)
    print('COMPARISON TABLE  —  Bulk Hubbard U (eV)')
    print('=' * 70)
    header = f"{'Compound':<8}  {'QE hp.x':>10}  {'GPAW cLDA':>10}  {'Diff':>8}  {'Diff%':>7}"
    print(header)
    print('-' * 70)

    for compound in COMPOUND_DB:
        u_qe = results_qe.get(compound)

        # Try to read GPAW result from saved .npy
        npy = _HERE / f'{compound}_cLDA_results.npy'
        u_gpaw = None
        if npy.exists():
            try:
                data = np.load(str(npy), allow_pickle=True).item()
                u_gpaw = float(data['U_eV'])
            except Exception:
                pass

        u_qe_s   = f'{u_qe:.4f}'   if u_qe   is not None else '    —   '
        u_gpaw_s = f'{u_gpaw:.4f}' if u_gpaw is not None else '    —   '

        if u_qe is not None and u_gpaw is not None:
            diff    = u_gpaw - u_qe
            diff_s  = f'{diff:+.4f}'
            diffp_s = f'{100*diff/u_qe:+.1f}%'
        else:
            diff_s  = '    —'
            diffp_s = '    —'

        print(f"  {compound:<6}  {u_qe_s:>10}  {u_gpaw_s:>10}  {diff_s:>8}  {diffp_s:>7}")

    print('-' * 70)
    print('  Diff = GPAW cLDA − QE hp.x')
    print()


# ============================================================
# MAIN
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description='QE hp.x bulk U calculation — pair to gpaw_U_clda.py',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python qe_U_bulk.py                  # run all compounds
  python qe_U_bulk.py TiN VN           # run only TiN and VN
  python qe_U_bulk.py --compare-only   # print comparison table (no QE run)
        """
    )
    parser.add_argument('compounds', nargs='*',
                        help='Compounds to run (default: all)')
    parser.add_argument('--compare-only', action='store_true',
                        help='Skip QE runs, just print comparison table')
    parser.add_argument('--qe-bin', default='',
                        help='Path to QE bin directory (overrides QE_BIN env var)')
    args = parser.parse_args()

    if args.qe_bin:
        global PW_CMD, HP_CMD
        PW_CMD = args.qe_bin.rstrip('/') + '/pw.x'
        HP_CMD = args.qe_bin.rstrip('/') + '/hp.x'

    to_run = args.compounds if args.compounds else list(COMPOUND_DB.keys())

    # Validate
    for c in to_run:
        if c not in COMPOUND_DB:
            print(f"ERROR: Unknown compound '{c}'. Valid: {list(COMPOUND_DB)}")
            return 1

    results_qe = {}

    if not args.compare_only:
        WORK_DIR.mkdir(parents=True, exist_ok=True)
        print('=' * 60)
        print('QE hp.x BULK U CALCULATIONS')
        print('=' * 60)
        print(f'  Compounds   : {to_run}')
        print(f'  POSCAR dir  : {POSCAR_DIR}')
        print(f'  Pseudo dir  : {PSEUDO_DIR}')
        print(f'  Work dir    : {WORK_DIR}')
        print(f'  k-mesh      : {KPTS}')
        print(f'  q-mesh      : {QPTS}')
        print(f'  Cutoff (Ry) : {ECUTWFC} / {ECUTRHO}')

        for compound in to_run:
            U = run_compound_qe(compound)
            results_qe[compound] = U

        # Save QE results
        out_npy = _HERE / 'qe_bulk_U_results.npy'
        np.save(str(out_npy), results_qe)
        print(f'\n  QE results saved -> {out_npy.name}')

    else:
        # Load previously saved QE results if available
        out_npy = _HERE / 'qe_bulk_U_results.npy'
        if out_npy.exists():
            results_qe = np.load(str(out_npy), allow_pickle=True).item()
            print(f'Loaded QE results from {out_npy.name}')
        else:
            print('No saved QE results found. Run without --compare-only first.')

    print_comparison(results_qe)
    return 0


if __name__ == '__main__':
    exit(main())
