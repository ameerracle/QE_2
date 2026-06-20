"""
PDOS comparison plots: slab vs combi for each metal + adsorbate.
Handles spin-polarized VN+U separately.
Produces per-metal PNGs in PDOS_Plots/comparison/, then 3x2 subplots per adsorbate.
"""

from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
sys.path.insert(0, str(Path(__file__).parent))
from plot_db import write_table

mpl.rcParams.update({'font.size': 11, 'axes.linewidth': 1.2})

BASE = Path('/home/ameer_ubuntu/Git_projects/QE_2')
OUT_DIR = BASE / 'Final_plots/PDOS_Plots/comparison'
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── d-band center: highest value per compound = surface atom ────────────────
_dband_df = pd.read_csv(BASE / 'Final_plots/d_band_values.csv')
DBAND = _dband_df.groupby('compound')['d_band_center_eV'].max().to_dict()
# map keys to metal code
DBAND_MAP = {
    'Sc': DBAND.get('ScN'), 'Ti': DBAND.get('TiN'), 'V': DBAND.get('VN'),
    'Nb': DBAND.get('NbN'), 'Zr': DBAND.get('ZrN'), 'V_U': DBAND.get('VN (+U)'),
}

# ── directories ─────────────────────────────────────────────────────────────
SLAB_DIRS = {
    'Sc':  BASE / 'no_u/PDOS/ScN_slab',
    'Ti':  BASE / 'no_u/PDOS/TiN_slab',
    'V':   BASE / 'no_u/PDOS/VN_slab',
    'Nb':  BASE / 'no_u/PDOS/NbN_slab',
    'Zr':  BASE / 'no_u/PDOS/ZrN_slab',
    'V_U': BASE / 'u/PDOS/VN_slab',
}
COMBI_DIRS = {
    'Sc':  {'Li2S4': BASE/'no_u/PDOS/ScN_Li2S4',  'Li2S8': BASE/'no_u/PDOS/ScN_Li2S8',  'S8': BASE/'no_u/PDOS/ScN_S8'},
    'Ti':  {'Li2S4': BASE/'no_u/PDOS/TiN_Li2S4',  'Li2S8': BASE/'no_u/PDOS/TiN_Li2S8',  'S8': BASE/'no_u/PDOS/TiN_S8'},
    'V':   {'Li2S4': BASE/'no_u/PDOS/VN_Li2S4',   'Li2S8': BASE/'no_u/PDOS/VN_Li2S8',   'S8': BASE/'no_u/PDOS/VN_S8'},
    'Nb':  {'Li2S4': BASE/'no_u/PDOS/NbN_Li2S4',  'Li2S8': BASE/'no_u/PDOS/NbN_Li2S8',  'S8': BASE/'no_u/PDOS/NbN_S8'},
    'Zr':  {'Li2S4': BASE/'no_u/PDOS/ZrN_Li2S4',  'Li2S8': BASE/'no_u/PDOS/ZrN_Li2S8',  'S8': BASE/'no_u/PDOS/ZrN_S8'},
    'V_U': {'Li2S4': BASE/'u/lobster/VN_Li2S4_combi','Li2S8': BASE/'u/lobster/VN_Li2S8_combi','S8': BASE/'u/lobster/VN_S8_combi'},
}

METAL_ELEM  = {'Sc':'Sc','Ti':'Ti','V':'V','Nb':'Nb','Zr':'Zr','V_U':'V'}
METAL_LABEL = {'Sc':'ScN','Ti':'TiN','V':'VN','Nb':'NbN','Zr':'ZrN','V_U':'VN (+U)'}

# Fermi energies from nscf.out on rorqual/nibi (eV, absolute)
FERMI = {
    'Sc': {'Li2S4': 2.0664, 'Li2S8': 1.9508, 'S8': 2.3108},
    'Ti': {'Li2S4': 3.9018, 'Li2S8': 3.6401, 'S8': 3.6653},
    'V':  {'Li2S4': 3.6604, 'Li2S8': 3.9278, 'S8': 3.7265},
    'Nb': {'Li2S4': 4.7543, 'Li2S8': 4.5951, 'S8': 4.5377},
    'Zr': {'Li2S4': 4.4116, 'Li2S8': 4.4226, 'S8': 4.4330},
    'V_U':{'Li2S4': 3.4832, 'Li2S8': 3.2664, 'S8': 3.2146},
}

EMIN, EMAX = -6.5, 2.0


# ── PDOS reader ──────────────────────────────────────────────────────────────
def _parse_fname(fname):
    m = re.search(r'atm#\d+\(([A-Za-z]+)\)_wfc#\d+\(([spdf])\)', fname)
    return (m.group(1), m.group(2)) if m else (None, None)


def load_pdos(directory, spin=False, fermi=0.0):
    """
    Accumulate pdos_atm files by (element, orbital), shifting E by -fermi.
    Non-spin: returns {key: (E, dos)}
    Spin:     returns {key: (E, dos_up, dos_dw)}   (down stored as positive)
    """
    directory = Path(directory)
    acc = {}
    e_ref = None

    for f in sorted(directory.glob('*.pdos_atm*')):
        elem, orb = _parse_fname(f.name)
        if elem is None:
            continue
        try:
            data = np.loadtxt(f, comments='#')
        except Exception:
            continue
        if data.ndim < 2 or len(data) == 0:
            continue

        E = data[:, 0] - fermi
        if e_ref is None:
            e_ref = E
        elif len(E) != len(e_ref):
            continue

        key = (elem, orb)
        if spin:
            # header: E ldosup ldosdw pdosup pdosdw [pdosup pdosdw ...]
            ncols = data.shape[1]
            up_cols = list(range(3, ncols, 2))
            dw_cols = list(range(4, ncols, 2))
            dos_up = data[:, up_cols].sum(axis=1) if up_cols else data[:, 3]
            dos_dw = data[:, dw_cols].sum(axis=1) if dw_cols else data[:, 4]
            if key in acc:
                acc[key] = (E, acc[key][1] + dos_up, acc[key][2] + dos_dw)
            else:
                acc[key] = (E, dos_up, dos_dw)
        else:
            # header: E ldos pdos [pdos ...]
            dos = data[:, 2:].sum(axis=1)
            if key in acc:
                acc[key] = (E, acc[key][1] + dos)
            else:
                acc[key] = (E, dos)

    return acc


def _get(pdos, elem, orb, spin=False, component='up'):
    key = (elem, orb)
    if key not in pdos:
        return None, None
    v = pdos[key]
    E = v[0]
    if spin:
        return E, (v[1] if component == 'up' else v[2])
    return E, v[1]


# ── single panel ─────────────────────────────────────────────────────────────
def plot_panel(metal, adsorbate, ax):
    spin = (metal == 'V_U')
    elem = METAL_ELEM[metal]
    eps_d = DBAND_MAP[metal]

    ef = FERMI[metal][adsorbate]
    combi = load_pdos(COMBI_DIRS[metal][adsorbate], spin=spin, fermi=ef)

    ax.axhline(0, color='k', lw=0.5)
    ax.axvline(0, color='gray', lw=0.8, ls='--')

    if spin:
        e = next((v[0] for v in combi.values()), None)
        _, d_up  = _get(combi, elem, 'd',  spin=True, component='up')
        _, d_dw  = _get(combi, elem, 'd',  spin=True, component='dw')
        _, p_up  = _get(combi, 'N',  'p',  spin=True, component='up')
        _, p_dw  = _get(combi, 'N',  'p',  spin=True, component='dw')
        _, li_up = _get(combi, 'Li', 's',  spin=True, component='up')
        _, li_dw = _get(combi, 'Li', 's',  spin=True, component='dw')
        _, s_up  = _get(combi, 'S',  'p',  spin=True, component='up')
        _, s_dw  = _get(combi, 'S',  'p',  spin=True, component='dw')

        if e is not None:
            if d_up is not None:
                ax.plot(e,  d_up, lw=1.2, color='steelblue', label=f'{elem}-d')
                ax.plot(e, -d_dw, lw=1.2, color='steelblue')
            if p_up is not None:
                ax.plot(e,  p_up, lw=1.0, color='salmon',    label='N-p')
                ax.plot(e, -p_dw, lw=1.0, color='salmon')
            if s_up is not None:
                ax.plot(e,  s_up, lw=1.2, color='darkorange', label='S-p')
                ax.plot(e, -s_dw, lw=1.2, color='darkorange')
            if li_up is not None:
                ax.plot(e,  li_up, lw=1.2, color='mediumseagreen', label='Li-s')
                ax.plot(e, -li_dw, lw=1.2, color='mediumseagreen')

        if eps_d is not None:
            ax.axvline(eps_d, color='k', lw=1.0, ls=':')
            ax.text(eps_d - 0.1, 0.97, rf'$\varepsilon_{{d,+U}}={eps_d:.2f}$',
                    transform=ax.get_xaxis_transform(), ha='right', va='top', fontsize=9)

    else:
        e, d  = _get(combi, elem, 'd')
        _, p  = _get(combi, 'N',  'p')
        _, s  = _get(combi, 'S',  'p')
        _, li = _get(combi, 'Li', 's')

        if e is not None:
            if d  is not None: ax.plot(e, d,  lw=1.2, color='steelblue',      label=f'{elem}-d')
            if p  is not None: ax.plot(e, p,  lw=1.0, color='salmon',          label='N-p')
            if s  is not None: ax.plot(e, s,  lw=1.2, color='darkorange',      label='S-p')
            if li is not None: ax.plot(e, li, lw=1.2, color='mediumseagreen',  label='Li-s')

        sym_map = {'Sc': r'\varepsilon_{Sc}', 'Ti': r'\varepsilon_D',
                   'V': r'\varepsilon_D', 'Nb': r'\varepsilon_{Nb}', 'Zr': r'\varepsilon_D'}
        if eps_d is not None:
            ax.axvline(eps_d, color='k', lw=1.0, ls=':')
            sym = sym_map.get(metal, r'\varepsilon_D')
            ax.text(eps_d - 0.1, 0.97, rf'${sym}={eps_d:.2f}$',
                    transform=ax.get_xaxis_transform(), ha='right', va='top', fontsize=9)

    ax.set_xlim(EMIN, EMAX)
    ax.set_ylim(-200, 200) if spin else ax.set_ylim(0, 200)
    ax.set_xlabel(r'$E - E_F$ (eV)', fontsize=10, fontweight='bold')
    ax.set_ylabel('PDOS (states/eV)', fontsize=10, fontweight='bold')
    ax.legend(fontsize=7, loc='upper left', frameon=False)


# ── generate subplots ────────────────────────────────────────────────────────
METAL_ORDER = ['Sc', 'Ti', 'V', 'V_U', 'Nb', 'Zr']
ADSORBATES  = ['Li2S4', 'Li2S8', 'S8']
ADS_LABELS  = {'Li2S4': r'Li$_2$S$_4$', 'Li2S8': r'Li$_2$S$_8$', 'S8': r'S$_8$'}
SUBPLOT_OUT = {
    'Li2S4': BASE / 'Final_plots/Figure_3_PDOS_subplot.png',
    'Li2S8': BASE / 'Final_plots/Figure_S2_PDOS_subplot.png',
    'S8':    BASE / 'Final_plots/Figure_S3_PDOS_subplot.png',
}

for ads in ADSORBATES:
    fig, axes = plt.subplots(3, 2, figsize=(12, 15))
    axes = axes.ravel()
    panel_labels = [f'({chr(97+i)})' for i in range(6)]

    for i, metal in enumerate(METAL_ORDER):
        ax = axes[i]
        plot_panel(metal, ads, ax)
        ax.text(-0.18, 1.03, panel_labels[i], transform=ax.transAxes,
                fontsize=14, fontweight='bold', ha='left', va='top')

        # also save individual PNG
        fig_s, ax_s = plt.subplots(figsize=(5, 4))
        plot_panel(metal, ads, ax_s)
        fig_s.tight_layout()
        png = OUT_DIR / f'{metal}_{ads}_comparison.png'
        fig_s.savefig(png, dpi=150, bbox_inches='tight')
        plt.close(fig_s)
        print(f'  saved {png.name}')

    fig.tight_layout(pad=1.5)
    fig.savefig(SUBPLOT_OUT[ads], dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved subplot → {SUBPLOT_OUT[ads].name}')

    # Extra 2-row × 3-col layout for Li2S4 (Figure 3 main)
    if ads == 'Li2S4':
        fig2, axes2 = plt.subplots(2, 3, figsize=(18, 10))
        axes2 = axes2.ravel()
        for i, metal in enumerate(METAL_ORDER):
            plot_panel(metal, ads, axes2[i])
            axes2[i].text(-0.18, 1.03, panel_labels[i], transform=axes2[i].transAxes,
                          fontsize=14, fontweight='bold', ha='left', va='top')
        fig2.tight_layout(pad=1.5)
        out2 = BASE / 'Final_plots/Figure_3_2row_PDOS_subplot.png'
        fig2.savefig(out2, dpi=300, bbox_inches='tight')
        plt.close(fig2)
        print(f'Saved 2-row layout → {out2.name}')

# ── archive plotted PDOS curves to SQLite ────────────────────────────────────
# Long-form: one row per (metal, adsorbate, element, orbital, spin_channel, E point)
_db_rows = []
_curves = [('metal_d', 'd'), ('N', 'p'), ('S', 'p'), ('Li', 's')]
for metal in METAL_ORDER:
    spin = (metal == 'V_U')
    melem = METAL_ELEM[metal]
    for ads in ADSORBATES:
        ef = FERMI[metal][ads]
        combi = load_pdos(COMBI_DIRS[metal][ads], spin=spin, fermi=ef)
        for label, orb in _curves:
            elem = melem if label == 'metal_d' else label
            if spin:
                for comp in ('up', 'dw'):
                    E, dos = _get(combi, elem, orb, spin=True, component=comp)
                    if E is None:
                        continue
                    for e_val, d_val in zip(E, dos):
                        _db_rows.append({'metal': METAL_LABEL[metal], 'adsorbate': ads,
                                         'element': elem, 'orbital': orb, 'spin': comp,
                                         'E_minus_Ef': float(e_val), 'pdos': float(d_val)})
            else:
                E, dos = _get(combi, elem, orb)
                if E is None:
                    continue
                for e_val, d_val in zip(E, dos):
                    _db_rows.append({'metal': METAL_LABEL[metal], 'adsorbate': ads,
                                     'element': elem, 'orbital': orb, 'spin': 'none',
                                     'E_minus_Ef': float(e_val), 'pdos': float(d_val)})
write_table(pd.DataFrame(_db_rows), 'figure_3_pdos')

print('Done.')
