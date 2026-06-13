# Recovered from IPython history - session cell ending line 130328
# --- Figure 6: ICOHP at metal-adsorbate interface bonds ---
# Lobster ICOHP data (no_u only); more negative = stronger covalent interaction
import matplotlib
matplotlib.use('Agg')
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams.update({'font.family': 'Liberation Sans', 'mathtext.fontset': 'stix',
                     'font.weight': 'bold', 'axes.labelweight': 'bold', 'axes.titleweight': 'bold'})

base_dir = Path('/home/ameer_ubuntu/Git_projects/QE_2')
lobster_dir = base_dir / 'Final_plots/lobster_out'

metal_order_icohp = ['ScN', 'VN', 'TiN', 'NbN', 'ZrN']
ads_order_icohp = ['Li2S4', 'Li2S8', 'S8']
ads_colors_icohp = dict(zip(ads_order_icohp, sns.color_palette('plasma', len(ads_order_icohp))))
ADSORBATE_SPECIES = {'S', 'Li'}
INTRA_ADSORBATE = {'S-S', 'Li-S', 'S-Li', 'Li-Li'}

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
            'metal': metal,
            'adsorbate': ads,
            'atom1': el1,
            'atom2': el2,
            'bond_type': f'{el1}-{el2}',
            'distance': float(parts[3]),
            'ICOHP': float(parts[7]),
        })

icohp_df = pd.DataFrame(icohp_data)
icohp_df = icohp_df[icohp_df['atom1'].isin(ADSORBATE_SPECIES) | icohp_df['atom2'].isin(ADSORBATE_SPECIES)]
icohp_df = icohp_df[~icohp_df['bond_type'].isin(INTRA_ADSORBATE)]

# --- Plot 1: ICOHP by metal (3 rows x 2 columns... but 5 metals, so 3x2 with empty slot) ---
fig, axes = plt.subplots(3, 2, figsize=(14, 16), sharey=True, sharex=True)
axes_flat = axes.ravel()

for idx, metal in enumerate(metal_order_icohp):
    ax = axes_flat[idx]
    for ads in ads_order_icohp:
        subset = icohp_df[(icohp_df['metal'] == metal) & (icohp_df['adsorbate'] == ads)]
        if subset.empty:
            continue
        ax.scatter(
            subset['distance'],
            subset['ICOHP'],
            s=30,
            alpha=0.7,
            color=ads_colors_icohp[ads],
            label=ads,
            edgecolors='none',
        )
    ax.axhline(0, color='grey', linewidth=0.8, linestyle='--')
    ax.axhline(-1, color='grey', linewidth=0.6, linestyle=':')
    ax.axhline(-2, color='grey', linewidth=0.6, linestyle=':')
    ax.annotate(metal, xy=(0.95, 0.95), xycoords='axes fraction',
                fontsize=16, fontweight='bold', ha='right', va='top')
    ax.tick_params(axis='both', labelsize=12, labelbottom=True)
    ax.set_xlabel('Bond length (Angstrom)', fontsize=13)
    if idx % 2 == 0:
        ax.set_ylabel('ICOHP (eV)', fontsize=13)

# Hide empty 6th subplot
axes_flat[5].set_visible(False)

axes_flat[0].legend(title='Adsorbate', loc='lower right',
                  frameon=False, fontsize=10)

fig.suptitle('ICOHP at Metal-Adsorbate Interface', fontsize=16, fontweight='bold', y=0.995)
plt.tight_layout()
plt.subplots_adjust(hspace=0.3)
plt.savefig(base_dir / 'Final_plots' / 'Figure_ICOHP.png', dpi=600, bbox_inches='tight')

# --- Plot 2: ICOHP by adsorbate (comparing all metals) ---
metal_colors_icohp = dict(zip(metal_order_icohp, sns.color_palette('tab10', len(metal_order_icohp))))
metal_markers = dict(zip(metal_order_icohp, ['o', 's', '^', 'D', 'v']))

fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True, sharex=True)

for idx, ads in enumerate(ads_order_icohp):
    ax = axes[idx]
    for metal in metal_order_icohp:
        subset = icohp_df[(icohp_df['metal'] == metal) & (icohp_df['adsorbate'] == ads)]
        if subset.empty:
            continue
        ax.scatter(
            subset['distance'],
            subset['ICOHP'],
            s=30,
            alpha=0.7,
            color=metal_colors_icohp[metal],
            marker=metal_markers[metal],
            label=metal,
            edgecolors='none',
        )
    ax.axhline(0, color='grey', linewidth=0.8, linestyle='--')
    ax.axhline(-1, color='grey', linewidth=0.6, linestyle=':')
    ax.axhline(-2, color='grey', linewidth=0.6, linestyle=':')
    ax.annotate(ads, xy=(0.95, 0.95), xycoords='axes fraction',
                fontsize=16, fontweight='bold', ha='right', va='top')
    ax.tick_params(axis='both', labelsize=12, labelbottom=True)
    ax.set_xlabel('Bond length (Angstrom)', fontsize=13)
    if idx == 0:
        ax.set_ylabel('ICOHP (eV)', fontsize=13)

handles, labels = axes[0].get_legend_handles_labels()
axes[0].legend(handles, labels, title='Metal', loc='lower right',
               frameon=False, fontsize=9, ncol=2)

fig.suptitle('ICOHP at BCPs: All Metals', fontsize=16, fontweight='bold', y=1.01)
plt.tight_layout()
plt.savefig(base_dir / 'Final_plots' / 'Figure_ICOHP_by_ads.png', dpi=600, bbox_inches='tight')
