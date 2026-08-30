# Recovered from IPython history - session cell ending line 145843
# > Figure 4c: QTAIM 2D Bond Classification Map (|V|/G vs ???)
import matplotlib
matplotlib.use('Agg')
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.lines import Line2D

plt.rcParams.update({'font.family': 'Liberation Sans', 'mathtext.fontset': 'stix', 'font.weight': 'bold', 'axes.labelweight': 'bold', 'axes.titleweight': 'bold'})

base_dir = Path('/home/ameer_ubuntu/Git_projects/QE_2')
bcp_dir = base_dir / 'no_u/bcp_filtered'
bcp_dir_u = base_dir / 'Final_plots/bcp_filtered_u'

metal_order = ['ScN', 'VN', 'TiN', 'NbN', 'ZrN', 'VN_U']
ads_order = ['Li2S4', 'Li2S8', 'S8']
ADSORBATE_SPECIES = {'S', 'Li', 'O'}
INTRA_ADSORBATE = {'S-S', 'Li-S', 'S-Li'}

bcp_data = []

for csv_path in sorted(bcp_dir.glob('*.csv')):
    stem = csv_path.stem.replace('_bcp_filtered', '')
    metal, ads = stem.split('_', 1)
    df_tmp = pd.read_csv(csv_path)
    df_tmp['metal'] = metal
    df_tmp['adsorbate'] = ads
    bcp_data.append(df_tmp)

for csv_path in sorted(bcp_dir_u.glob('*.csv')):
    stem = csv_path.stem.replace('_bcp_filtered', '')
    if stem.startswith('VN_'):
        metal = 'VN_U'
        ads = stem.replace('VN_', '')
        df_tmp = pd.read_csv(csv_path)
        df_tmp['metal'] = metal
        df_tmp['adsorbate'] = ads
        bcp_data.append(df_tmp)

bcp_df = pd.concat(bcp_data, ignore_index=True)
bcp_df = bcp_df[bcp_df['atom1'].isin(ADSORBATE_SPECIES) | bcp_df['atom2'].isin(ADSORBATE_SPECIES)]
bcp_df = bcp_df[~bcp_df['bond_type'].isin(INTRA_ADSORBATE)]

metal_colors = dict(zip(metal_order, sns.color_palette('tab10', len(metal_order))))
ads_markers = {'Li2S4': 'o', 'Li2S8': 's', 'S8': '^'}

fig, ax = plt.subplots(figsize=(10, 8))

for metal in metal_order:
    for ads in ads_order:
        subset = bcp_df[(bcp_df['metal'] == metal) & (bcp_df['adsorbate'] == ads)]
        if subset.empty:
            continue
        sizes = (subset['bond_length_ang'] ** 2) * 8
        ax.scatter(
            subset['laplacian'],
            subset['abs_V_over_G'],
            s=sizes,
            alpha=0.7,
            color=metal_colors[metal],
            marker=ads_markers[ads],
            edgecolors='black',
            linewidths=0.3,
        )

# QTAIM reference lines
ax.axvline(0, color='black', linewidth=1.0, linestyle='-')
ax.axhline(1.0, color='grey', linewidth=0.8, linestyle='--')
ax.axhline(2.0, color='grey', linewidth=0.8, linestyle='--')

# Zone labels
ax.text(0.02, 0.22, 'Weak / vdW\n(|V|/G < 1)', transform=ax.transAxes, fontsize=10, color='grey', ha='left', va='bottom')
ax.text(0.02, 0.55, 'Intermediate\n(1 < |V|/G < 2)', transform=ax.transAxes, fontsize=10, color='grey', ha='left', va='bottom')
ax.text(0.02, 0.88, 'Shared-shell\n(|V|/G > 2)', transform=ax.transAxes, fontsize=10, color='grey', ha='left', va='top')
ax.text(0.70, 0.95, '??? > 0\n(closed-shell)', transform=ax.transAxes, fontsize=10, color='grey', ha='center', va='top')
ax.text(0.25, 0.95, '??? < 0\n(covalent)', transform=ax.transAxes, fontsize=10, color='grey', ha='center', va='top')

ax.set_xlabel('??? (a.u.)', fontsize=14)
ax.set_ylabel('|V|/G (a.u.)', fontsize=14)
ax.set_title('QTAIM Bond Classification Map', fontsize=16, fontweight='bold')

# Legends
metal_handles = [Line2D([0], [0], marker='o', color='w', markerfacecolor=metal_colors[m], markeredgecolor='black', markeredgewidth=0.8, markersize=11, label=m.replace('VN_U', 'VN (+U)')) for m in metal_order]
legend1 = ax.legend(handles=metal_handles, title='Metal', loc='upper left', frameon=False, fontsize=9, bbox_to_anchor=(1.01, 1))
ax.add_artist(legend1)

ads_handles = [Line2D([0], [0], marker=ads_markers[a], color='w', markeredgecolor='black', markerfacecolor='grey', markersize=10, label=a) for a in ads_order]
ax.legend(handles=ads_handles, title='Adsorbate', loc='lower right', frameon=False, fontsize=9)

for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontweight('bold')

plt.tight_layout()
plt.savefig(base_dir / 'Final_plots' / 'Figure_4c.png', dpi=600, bbox_inches='tight')
plt.close(fig)
