"""
cLDA Hubbard U calculation for transition-metal nitrides using GPAW.
Equivalent to QE hp.x linear response approach (finite difference method).
All compounds non-magnetic (spin-unpolarized).

Uses DEFAULT GPAW PAW datasets (built-in to GPAW installation).

Why a supercell?
  setups={metal: ':d,alpha'} on primitive perturbs ALL metal atoms (periodic
  images). A 2x2x2 supercell with integer atom index perturbs ONE site only -
  the correct Cococcioni-de Gironcoli procedure.

Strategy
  chi  (screened): full SCF at each alpha, D_asp from the one metal atom
  chi0 (bare)    : analytical from GS wavefunctions
                   chi0 = sum_{n,k} (df/de)|_{e_nk} * [d-char_{nk}]^2

Directory layout (all relative to this script):
  ./primitive/        POSCARs
  ./GPAW_U/           per-alpha SCF logs
  ./U_GPAW/           final .npy results

Usage:
    python gpaw_U_clda_default_paw.py             # serial
    gpaw -P 16 python gpaw_U_clda_default_paw.py  # parallel (recommended)
"""

import numpy as np
from pathlib import Path
from ase.io import read
from gpaw import GPAW, PW, Mixer, FermiDirac
from gpaw.utilities import unpack_density
from gpaw.mpi import world

# ============================================================
# BASE DIRECTORIES
# ============================================================
BASE       = Path(__file__).resolve().parent
DIR_PRIM   = BASE / 'primitive'
DIR_LOGS   = BASE / 'GPAW_U'
DIR_RESULT = BASE / 'U_GPAW'

# ============================================================
# CONFIGURATION
# ============================================================

COMPOUND_DB = {
    'NbN': ('NbN_prim.poscar', 'Nb'),
    'ScN': ('ScN_prim.poscar', 'Sc'),
    'TiN': ('TiN_prim.poscar', 'Ti'),
    'VN':  ('VN_prim.poscar',  'V'),
    'ZrN': ('ZrN_prim.poscar', 'Zr'),
}

# Set to a single string to run one compound, or None to run all.
COMPOUND = None   # e.g. 'TiN'

# Supercell
SUPER = (2, 2, 2)

# Perturbation strengths (eV)
ALPHAS = [-0.20, -0.15, -0.10, -0.05, 0.00, 0.05, 0.10, 0.15, 0.20]

# DFT settings  (8x8x8 on primitive = 4x4x4 on 2x2x2 supercell)
KPTS  = (4, 4, 4)
ECUT  = 600         # eV
SMEAR = 0.1         # eV  Fermi-Dirac smearing

# Reference U values from QE hp.x — USPP pseudopotentials
QE_U_REF = {
    'VN':  5.8675,
    'TiN': 5.5479,
    'NbN': 3.5140,
    'ZrN': 2.9683,
    'ScN': 3.5695,
}

# ============================================================
# HELPERS
# ============================================================

def ensure_dirs():
    for d in [DIR_PRIM, DIR_LOGS, DIR_RESULT]:
        d.mkdir(parents=True, exist_ok=True)
        if not d.exists():
            raise RuntimeError(f"Failed to create directory: {d}")


def d_indices_for_atom(calc, atom_index):
    """
    Projector indices for the VALENCE l=2 (d) channel only.
    Some PAW datasets (e.g. Nb, Zr) have two d blocks (semi-core + valence).
    We always take the LAST d block, which is the valence one.
    """
    setup = calc.wfs.setups[atom_index]
    idx, d_blocks = 0, []
    for l in setup.l_j:
        block = list(range(idx, idx + 2 * l + 1))
        if l == 2:
            d_blocks.append(block)
        idx += 2 * l + 1
    if world.rank == 0:
        print(f"  [d_idx debug] atom={atom_index}  "
              f"d_blocks={len(d_blocks)}  "
              f"using last block (valence): {d_blocks[-1] if d_blocks else []}")
    return d_blocks[-1] if d_blocks else []


def get_d_occupation_Dasp(calc, atom_index):
    """Total d occupation from PAW density matrix D_asp (screened response)."""
    try:
        D_asp = calc.density.D_asp
        if atom_index not in D_asp:
            return None
        D_sp  = D_asp[atom_index]
        d_idx = d_indices_for_atom(calc, atom_index)
        total = 0.0
        for s in range(D_sp.shape[0]):
            D = unpack_density(D_sp[s])
            total += float(sum(D[i, i] for i in d_idx))
        return total
    except Exception as ex:
        if world.rank == 0:
            print(f"  Warning (D_asp): {ex}")
        return None


def compute_chi0_analytical(calc, atom_index, kT=SMEAR):
    """
    Bare response chi0 from GS wavefunctions. MPI-safe.

    chi0 = sum_{n,k} w_k * (df/de)|_{e_nk} * [d-char_{nk}]^2

    Spin-unpolarised FermiDirac: df/de = -f * (2-f) / (2*kT)
    """
    d_idx = d_indices_for_atom(calc, atom_index)
    if world.rank == 0:
        print(f"  [chi0 debug] atom_index={atom_index}  d_idx={d_idx}")

    chi0          = 0.0
    found_kpts    = 0
    missing_pani  = 0
    for kpt in calc.wfs.kpt_u:
        if kpt.s != 0:
            continue
        found_kpts += 1
        wk  = kpt.weight
        f_n = kpt.f_n
        if atom_index not in kpt.P_ani:
            missing_pani += 1
            continue
        P = kpt.P_ani[atom_index]
        for n in range(len(f_n)):
            d_char = sum(float(abs(P[n, i])**2) for i in d_idx)
            fn     = float(f_n[n])
            dfdeps = -fn * (wk - fn) / (wk * kT)
            chi0  += dfdeps * d_char**2

    chi0 = world.sum(chi0)
    if world.rank == 0:
        if missing_pani == found_kpts and found_kpts > 0:
            print(f"  [chi0 debug] WARNING: atom {atom_index} missing from "
                  f"P_ani at ALL {found_kpts} local k-points - "
                  f"chi0 will be 0. Check metal_index is correct.")
        else:
            print(f"  [chi0 debug] k-points processed={found_kpts}  "
                  f"missing P_ani={missing_pani}  chi0={chi0:.6f}")
    return chi0


def extract_U(chi0, occ_screened):
    """Linear fit for chi (screened slope), return U = 1/chi0 - 1/chi."""
    alphas_arr = np.array(ALPHAS)
    occ_s_arr  = np.array(occ_screened)
    chi = np.polyfit(alphas_arr, occ_s_arr, 1)[0]
    if world.rank == 0:
        print(f"  [debug] alpha vs occ_screened:")
        for a, o in zip(ALPHAS, occ_screened):
            print(f"    alpha={a:+.2f}  occ={o:.6f}")
    if abs(chi0) < 1e-10:
        if world.rank == 0:
            print("  WARNING: chi0 ~ 0 - bare response failed. "
                  "Check d_idx and P_ani. U undefined.")
        return float('nan'), float(chi)
    if abs(chi) < 1e-10:
        if world.rank == 0:
            print("  WARNING: chi ~ 0 - screened occ flat. U undefined.")
        return float('nan'), float(chi)
    U = (1.0 / chi0) - (1.0 / chi)
    return float(U), float(chi)

# ============================================================
# CORE cLDA RUN
# ============================================================

def clda_run(compound, atoms, metal_index, setups_gs, setups_label):
    """
    Full cLDA cycle: GS + alpha loop + U extraction.
    setups_label : 'PAW' (used for file naming)
    Returns result dict or None on failure.
    """
    metal   = COMPOUND_DB[compound][1]
    prim    = read(DIR_PRIM / COMPOUND_DB[compound][0])
    gs_gpw  = DIR_LOGS / f'{compound}_{setups_label}_gs.gpw'
    gs_txt  = DIR_LOGS / f'{compound}_{setups_label}_gs.txt'

    # ---- Ground state ----
    if world.rank == 0:
        print(f"\n  [{setups_label}] Ground state ...")

    calc_gs = GPAW(
        mode=PW(ECUT),
        xc='PBE',
        kpts=KPTS,
        setups=setups_gs,
        occupations=FermiDirac(SMEAR),
        mixer=Mixer(beta=0.05, nmaxold=5, weight=50),
        convergence={'energy': 1e-6, 'density': 1e-6},
        txt=str(gs_txt),
    )
    atoms.calc = calc_gs
    try:
        e_gs = atoms.get_potential_energy()
    except Exception as ex:
        if world.rank == 0:
            print(f"  [{setups_label}] GS FAILED: {ex}")
        return None

    # Extract chi0/occ BEFORE write — P_ani only lives in memory here
    occ_gs = get_d_occupation_Dasp(calc_gs, metal_index)
    chi0   = compute_chi0_analytical(calc_gs, metal_index)
    calc_gs.write(str(gs_gpw))

    if world.rank == 0:
        print(f"  [{setups_label}] E_gs={e_gs:.4f} eV  "
              f"{metal} d-occ={occ_gs:.4f}  chi0={chi0:.4f} eV^-1")

    # ---- Alpha loop ----
    occ_screened = []
    for alpha in ALPHAS:
        if world.rank == 0:
            print(f"  [{setups_label}] alpha={alpha:+.2f} eV ...", end='', flush=True)

        atoms_pert = prim.repeat(SUPER)
        log_file   = DIR_LOGS / f'{compound}_{setups_label}_alpha_{alpha:+.2f}.txt'

        setups_pert = {metal_index: f':d,{alpha}', 'default': 'paw'}

        calc_pert = GPAW(
            mode=PW(ECUT),
            xc='PBE',
            kpts=KPTS,
            setups=setups_pert,
            occupations=FermiDirac(SMEAR),
            mixer=Mixer(beta=0.05, nmaxold=5, weight=50),
            convergence={'energy': 1e-7, 'density': 1e-7},
            txt=str(log_file),
        )
        atoms_pert.calc = calc_pert
        try:
            atoms_pert.get_potential_energy()
            occ_s = get_d_occupation_Dasp(calc_pert, metal_index)
        except Exception as ex:
            if world.rank == 0:
                print(f" FAILED ({ex})")
            occ_s = None

        occ_screened.append(occ_s if occ_s is not None else 0.0)
        if world.rank == 0:
            print(f" d-occ={occ_s:.4f}" if occ_s else " d-occ=None")

    # ---- Extract U ----
    U, chi = extract_U(chi0, occ_screened)
    if world.rank == 0:
        print(f"  [{setups_label}] chi={chi:.4f} eV^-1  =>  U={U:.3f} eV")

    return {
        'compound':     compound,
        'metal':        metal,
        'label':        setups_label,
        'supercell':    SUPER,
        'kpts':         KPTS,
        'ecut_eV':      ECUT,
        'smear_eV':     SMEAR,
        'alphas':       ALPHAS,
        'occ_screened': occ_screened,
        'chi0':         chi0,
        'chi':          chi,
        'U_eV':         U,
    }

# ============================================================
# PER-COMPOUND WRAPPER
# ============================================================

def run_compound(compound):
    poscar_file, metal = COMPOUND_DB[compound]
    poscar_path = DIR_PRIM / poscar_file

    if not poscar_path.exists():
        if world.rank == 0:
            print(f"  SKIP {compound}: POSCAR not found at {poscar_path}")
        return {}

    if world.rank == 0:
        print()
        print("#" * 60)
        print(f"#  {compound}  ({metal}N rocksalt)")
        print("#" * 60)

    prim        = read(poscar_path)
    atoms       = prim.repeat(SUPER)
    metal_index = next(
        i for i, s in enumerate(atoms.get_chemical_symbols()) if s == metal
    )
    if world.rank == 0:
        print(f"  Supercell: {len(atoms)} atoms  |  {metal} supercell index: {metal_index}")

    results = {}

    # ---- PAW Run ----
    data = clda_run(compound, atoms.copy(), metal_index,
                    setups_gs='paw', setups_label='PAW')
    if data:
        data['U_QE_ref'] = QE_U_REF.get(compound)
        results['PAW']   = data
        if world.rank == 0:
            np.save(str(DIR_RESULT / f'{compound}_PAW_cLDA.npy'), data)

    return results

# ============================================================
# FINAL COMPARISON TABLE
# ============================================================

def print_table(all_results):
    """
    Simple table: Compound | GPAW PAW | QE USPP | Diff
    """
    sep = '  ' + '-' * 50
    hdr = (f"  {'Compound':<8}  {'GPAW PAW':>10}  {'QE USPP':>10}  {'Diff':>8}")

    print()
    print('=' * 55)
    print('FINAL COMPARISON TABLE - Hubbard U (eV)')
    print('=' * 55)
    print(hdr)
    print(sep)

    for compound in COMPOUND_DB:
        res = all_results.get(compound, {})

        # PAW columns
        paw = res.get('PAW')
        if paw:
            u_paw    = f"{paw['U_eV']:.3f}"
            u_qe     = paw.get('U_QE_ref')
            u_qe_s   = f"{u_qe:.3f}" if u_qe else '   -   '
            diff_paw = f"{paw['U_eV']-u_qe:+.3f}" if u_qe else '   -   '
        else:
            u_paw = u_qe_s = diff_paw = '   -   '

        print(f"  {compound:<8}  {u_paw:>10}  {u_qe_s:>10}  {diff_paw:>8}")

    print(sep)
    print('  Diff = GPAW cLDA - QE hp.x  (positive = GPAW higher)')
    print()

# ============================================================
# MAIN
# ============================================================

if __name__ == '__main__':
    ensure_dirs()

    to_run = [COMPOUND] if COMPOUND else list(COMPOUND_DB.keys())

    if world.rank == 0:
        print(f"Base dir   : {BASE}")
        print(f"Compounds  : {to_run}")
        print(f"Supercell  : {SUPER}")
        print(f"k-mesh     : {KPTS}")
        print(f"Cutoff     : {ECUT} eV")
        print(f"Setups     : GPAW default PAW (built-in)")

    all_results = {}
    for _compound in to_run:
        all_results[_compound] = run_compound(_compound)

    if world.rank == 0:
        print_table(all_results)
        print("All done.")
