# Recovered from IPython history - session cell ending line 153045
import matplotlib
matplotlib.use('Agg')
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import subplot_mosaic
from matplotlib.patches import ConnectionPatch, Rectangle
import seaborn as sns

# > Figure: ICOHP Grouped Bar Plot ? Mean |ICOHP| by Metal / Bond Type / Adsorbate
# --------------------------------------------------------------------------
plt.rcParams.update({'font.family': 'Liberation Sans', 'mathtext.fontset': 'stix',
                     'font.weight': 'bold', 'axes.labelweight': 'bold', 'axes.titleweight': 'bold'})

base_dir = Path('/home/ameer_ubuntu/Git_projects/QE_2')
lobster_dir = base_dir / 'Final_plots/lobster_out'

metal_order = ['ScN', 'VN', 'VN_U', 'TiN', 'NbN', 'ZrN']
metal_labels = {m: m.replace('VN_U', 'VN (+U)') for m in metal_order}
ads_order = ['Li2S4', 'Li2S8', 'S8']
bond_cats = ['M\u2013Li', 'M\u2013S']

ADSORBATE_SPECIES = {'S', 'Li'}
INTRA_ADSORBATE = {'S-S', 'Li-S', 'S-Li', 'Li-Li'}
METALS = {'Sc', 'V', 'Ti', 'Nb', 'Zr'}

def map_bond_cat(row):
    a1, a2 = row['atom1'], row['atom2']
    pair = {a1, a2}
    if 'S' in pair and pair.intersection(METALS):
        return 'M\u2013S'
    if 'Li' in pair and pair.intersection(METALS):
        return 'M\u2013Li'
    return 'Other'

icohp_data = []
for d in sorted(lobster_dir.glob('*_combi')):
    stem = d.name.replace('_combi', '')
    if stem.startswith('VN_U'):
        metal, ads = 'VN_U', stem.replace('VN_U_', '')
    else:
        metal, ads = stem.split('_', 1)
    icohp_file = d / 'ICOHPLIST.lobster'
    if not icohp_file.exists():
        continue
    with open(icohp_file) as f:
        lines = f.readlines()
    for line in lines:
        parts = line.split()
        if len(parts) < 7 or not parts[0].isdigit():
            continue
        el1 = ''.join(c for c in parts[1] if c.isalpha())
        el2 = ''.join(c for c in parts[2] if c.isalpha())
        icohp_data.append({
            'metal': metal, 'adsorbate': ads, 'atom1': el1, 'atom2': el2,
            'bond_type': f'{el1}-{el2}', 'distance': float(parts[3]),
            'ICOHP': float(parts[7]),
        })

icohp_df = pd.DataFrame(icohp_data)
icohp_df = icohp_df[icohp_df['atom1'].isin(ADSORBATE_SPECIES) | icohp_df['atom2'].isin(ADSORBATE_SPECIES)]
icohp_df = icohp_df[~icohp_df['bond_type'].isin(INTRA_ADSORBATE)]
icohp_df['bond_cat'] = icohp_df.apply(map_bond_cat, axis=1)
icohp_df = icohp_df[icohp_df['bond_cat'] != 'Other']

summary = icohp_df.groupby(['metal', 'adsorbate', 'bond_cat']).agg(
    ICOHP_mean=('ICOHP', 'mean'),
).reset_index()
summary['ICOHP_abs'] = summary['ICOHP_mean'].abs()

ads_colors = dict(zip(ads_order, sns.color_palette('pastel', len(ads_order))))
ads_colors['S8'] = '#FFD65C'

mosaic = [metal_labels[m] for m in metal_order]
fig, axd = subplot_mosaic([mosaic], figsize=(16, 6), sharey=True)

for mi, metal in enumerate(metal_order):
    ax = axd[metal_labels[metal]]
    m_data = summary[summary['metal'] == metal]
    x_pos = 0
    for bi, bcat in enumerate(bond_cats):
        bm = m_data[m_data['bond_cat'] == bcat]
        if bm.empty:
            continue
        for ai, ads in enumerate(ads_order):
            row = bm[bm['adsorbate'] == ads]
            if row.empty:
                continue
            val = row['ICOHP_abs'].values[0]
            ax.bar(x_pos, val, color=ads_colors[ads], edgecolor='white', lw=0.5, width=0.7)
            ax.text(x_pos, val + 0.05, f'{val:.2f}', ha='center', va='bottom',
                    fontsize=7, fontweight='bold')
            x_pos += 1
        x_pos += 0.4

    # Bond-type labels higher up (just below bars)
    sub_x = 0
    for bi, bcat in enumerate(bond_cats):
        bm = m_data[m_data['bond_cat'] == bcat]
        n_bars = len(bm['adsorbate'].unique())
        if n_bars > 0:
            centre = sub_x + (n_bars - 1) / 2
            ax.text(centre, -0.06, bcat, transform=ax.get_xaxis_transform(),
                    ha='center', va='top', fontsize=9, fontweight='bold')
        sub_x += n_bars + 0.4

    # Metal label at the very bottom
    ax.text(0.5, -0.16, metal_labels[metal], transform=ax.transAxes,
            ha='center', va='top', fontsize=12, fontweight='bold')

    # Clean up axes ? match reference style
    if mi == 0:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    else:
        ax.spines[:].set_visible(False)
    if mi != 0:
        ax.yaxis.set_tick_params(left=False)
        ax.yaxis.set_tick_params(left=False)
    ax.xaxis.set_tick_params(bottom=False)
    ax.set_xticks([])
    ax.margins(y=0.25)

# Y-axis label only on leftmost
axd[metal_labels[metal_order[0]]].set_ylabel('Mean |ICOHP| (eV)', fontsize=13)

# Legend on last axis, upper right
handles = [Rectangle((0, 0), 1, 1, facecolor=ads_colors[a], edgecolor='white') for a in ads_order]
axd[metal_labels[metal_order[-1]]].legend(
    handles, ads_order,
    title='Adsorbate', ncol=3, loc='upper right',
    bbox_to_anchor=(1, 1), frameon=False, fontsize=10)

# Connecting line across bottom
conn = ConnectionPatch(
    xyA=(0, 0), coordsA=axd[metal_labels[metal_order[0]]].transAxes,
    xyB=(1, 0), coordsB=axd[metal_labels[metal_order[-1]]].transAxes,
    lw=1.0
)
fig.add_artist(conn)

fig.tight_layout(pad=0.3)
fig.subplots_adjust(wspace=0.20)
plt.savefig(base_dir / 'Final_plots' / 'Figure_6_ICOHP_Stacked_Bar.png',
            dpi=600, bbox_inches='tight')
