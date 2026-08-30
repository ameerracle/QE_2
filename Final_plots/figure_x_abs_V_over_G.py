# Recovered from IPython history - session cell ending line 141236
#
# > Figure 5: |V|/G at BCPs (bonding character) ---
# Metal-adsorbate interface bonds only; |V|/G < 1 = closed-shell (ionic/vdW), > 2 = strongly covalent
import matplotlib
matplotlib.use('Agg')
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, NullFormatter
import seaborn as sns

plt.rcParams.update({'font.family': 'Liberation Sans', 'mathtext.fontset': 'stix',
                     'font.weight': 'bold', 'axes.labelweight': 'bold', 'axes.titleweight': 'bold'})

base_dir = Path('/home/ameer_ubuntu/Git_projects/QE_2')
bcp_dir = base_dir / 'no_u/bcp_filtered'
bcp_dir_u = base_dir / 'Final_plots/bcp_filtered_u'

metal_order_vg = ['ScN', 'VN', 'TiN', 'NbN', 'ZrN', 'VN_U']
ads_order_vg = ['Li2S4', 'Li2S8', 'S8']
ads_colors_vg = dict(zip(ads_order_vg, sns.color_palette('plasma', len(ads_order_vg))))
ADSORBATE_SPECIES = {'S', 'Li', 'O'}
INTRA_ADSORBATE = {'S-S', 'Li-S', 'S-Li'}

vg_data = []

for csv_path in sorted(bcp_dir.glob('*.csv')):
    stem = csv_path.stem.replace('_bcp_filtered', '')
    metal, ads = stem.split('_', 1)
    df_tmp = pd.read_csv(csv_path)
    df_tmp['metal'] = metal
    df_tmp['adsorbate'] = ads
    vg_data.append(df_tmp)

for csv_path in sorted(bcp_dir_u.glob('*.csv')):
    stem = csv_path.stem.replace('_bcp_filtered', '')
    if stem.startswith('VN_'):
        metal = 'VN_U'
        ads = stem.replace('VN_', '')
        df_tmp = pd.read_csv(csv_path)
        df_tmp['metal'] = metal
        df_tmp['adsorbate'] = ads
        vg_data.append(df_tmp)

vg_df = pd.concat(vg_data, ignore_index=True)
# Keep only adsorbate-involved BCPs, exclude intra-adsorbate bonds (S-S, Li-S) ? only metal-adsorbate interface
vg_df = vg_df[vg_df['atom1'].isin(ADSORBATE_SPECIES) | vg_df['atom2'].isin(ADSORBATE_SPECIES)]
vg_df = vg_df[~vg_df['bond_type'].isin(INTRA_ADSORBATE)]

# --- Plot 1: |V|/G by metal (3x2 grid) ---
fig, axes = plt.subplots(2, 3, figsize=(18, 10), sharey=True, sharex=True)
axes_flat = axes.ravel()

def set_bond_length_ticks(ax):
    ax.xaxis.set_major_locator(MultipleLocator(0.4))
    ax.xaxis.set_minor_locator(MultipleLocator(0.2))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.tick_params(axis='x', which='major', length=6)
    ax.tick_params(axis='x', which='minor', length=3)
    for label in ax.get_xticklabels(minor=False) + ax.get_yticklabels(minor=False):
        label.set_fontweight('bold')
    ax.xaxis.label.set_fontweight('bold')
    ax.yaxis.label.set_fontweight('bold')

for idx, metal in enumerate(metal_order_vg):
    ax = axes_flat[idx]
    display_name = metal.replace('VN_U', 'VN (+U)')
    for ads in ads_order_vg:
        subset = vg_df[(vg_df['metal'] == metal) & (vg_df['adsorbate'] == ads)]
        if subset.empty:
            continue
        ax.scatter(
            subset['bond_length_ang'],
            subset['abs_V_over_G'],
            s=30,
            alpha=0.7,
            color=ads_colors_vg[ads],
            label=ads,
            edgecolors='none',
        )
    ax.axhline(1.0, color='grey', linewidth=0.8, linestyle='--')
    ax.axhline(2.0, color='black', linewidth=0.8, linestyle='--')
    ax.annotate(display_name, xy=(0.95, 0.95), xycoords='axes fraction',
                fontsize=16, fontweight='bold', ha='right', va='top')
    ax.tick_params(axis='both', labelsize=12, labelbottom=True)
    ax.set_xlabel('Bond length (Angstrom)', fontsize=13)
    set_bond_length_ticks(ax)
    if idx % 3 == 0:
        ax.set_ylabel('|V|/G (a.u.)', fontsize=13)

axes_flat[0].legend(title='Adsorbate', loc='lower left', markerscale=1.25,
                  frameon=False, fontsize=10)

plt.tight_layout()
plt.subplots_adjust(hspace=0.15)
plt.savefig(base_dir / 'Final_plots' / 'Figure_abs_V_over_G.png', dpi=600, bbox_inches='tight')

# --- Plot 2: |V|/G by adsorbate (comparing all metals) ---
metal_order_plot2 = ['ScN', 'VN', 'TiN', 'NbN', 'ZrN', 'VN_U']
metal_colors_vg = dict(zip(metal_order_plot2, sns.color_palette('tab10', len(metal_order_plot2))))
metal_markers = dict(zip(metal_order_plot2, ['o', 's', '^', 'D', 'v', 'P']))

fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True, sharex=True)

for idx, ads in enumerate(ads_order_vg):
    ax = axes[idx]
    for metal in metal_order_plot2:
        subset = vg_df[(vg_df['metal'] == metal) & (vg_df['adsorbate'] == ads)]
        if subset.empty:
            continue
        ax.scatter(
            subset['bond_length_ang'],
            subset['abs_V_over_G'],
            s=30,
            alpha=0.7,
            color=metal_colors_vg[metal],
            marker=metal_markers[metal],
            label=metal.replace('VN_U', 'VN (+U)'),
            edgecolors='none',
        )
    ax.axhline(1.0, color='grey', linewidth=0.8, linestyle='--')
    ax.axhline(2.0, color='black', linewidth=0.8, linestyle='--')
    ax.annotate(ads, xy=(0.95, 0.95), xycoords='axes fraction',
                fontsize=16, fontweight='bold', ha='right', va='top')
    ax.tick_params(axis='both', labelsize=12, labelbottom=True)
    ax.set_xlabel('Bond length (Angstrom)', fontsize=13)
    set_bond_length_ticks(ax)
    if idx == 0:
        ax.set_ylabel('|V|/G (a.u.)', fontsize=13)

handles, labels = axes[0].get_legend_handles_labels()
axes[0].legend(handles, labels, title='Metal', loc='lower left',
               frameon=False, fontsize=10, ncol=2)

plt.tight_layout()
plt.savefig(base_dir / 'Final_plots' / 'Figure_abs_V_over_G_by_ads.png', dpi=600, bbox_inches='tight')
