# Recovered from IPython history - session cell ending line 129866
# --- Figure 4: BCP Laplacian scatter, metal-adsorbate interface bonds only (including VN +U) ---
from pathlib import Path
base_dir = Path("/home/ameer_ubuntu/Git_projects/QE_2")
bcp_dir = base_dir / 'no_u/bcp_filtered'
bcp_dir_u = base_dir / 'Final_plots/bcp_filtered_u'

metal_order_bcp = ['ScN', 'VN', 'TiN', 'NbN', 'ZrN', 'VN_U']
ads_order_bcp = ['Li2S4', 'Li2S8', 'S8']
ads_colors_bcp = dict(zip(ads_order_bcp, sns.color_palette('plasma', len(ads_order_bcp))))
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
# Keep only adsorbate-involved BCPs, exclude intra-adsorbate bonds (S-S, Li-S) ? only metal-adsorbate interface
bcp_df = bcp_df[bcp_df['atom1'].isin(ADSORBATE_SPECIES) | bcp_df['atom2'].isin(ADSORBATE_SPECIES)]
bcp_df = bcp_df[~bcp_df['bond_type'].isin(INTRA_ADSORBATE)]

# --- Plot 1: One subplot per metal (3 rows x 2 columns for 6 metals) ---
fig, axes = plt.subplots(2, 3, figsize=(18, 10), sharey=True, sharex=True)
axes_flat = axes.ravel()

for idx, metal in enumerate(metal_order_bcp):
    ax = axes_flat[idx]
    display_name = metal.replace('VN_U', 'VN (+U)')
    for ads in ads_order_bcp:
        subset = bcp_df[(bcp_df['metal'] == metal) & (bcp_df['adsorbate'] == ads)]
        if subset.empty:
            continue
        ax.scatter(
            subset['bond_length_ang'],
            subset['laplacian'],
            s=30,
            alpha=0.7,
            color=ads_colors_bcp[ads],
            label=ads,
            edgecolors='none',
        )
    ax.axhline(0, color='grey', linewidth=0.8, linestyle='--')
    ax.annotate(display_name, xy=(0.95, 0.95), xycoords='axes fraction',
                fontsize=16, fontweight='bold', ha='right', va='top')
    ax.tick_params(axis='both', labelsize=12, labelbottom=True)
    ax.set_xlabel('Bond length (Angstrom)', fontsize=13)
    if idx % 3 == 0:
        ax.set_ylabel('Laplacian (a.u.)', fontsize=13)

axes_flat[0].legend(title='Adsorbate', loc='upper left',
                  frameon=False, fontsize=10)

fig.suptitle('BCP Laplacian at Metal-Adsorbate Interface', fontsize=16, fontweight='bold', y=0.995)
plt.tight_layout()
plt.subplots_adjust(hspace=0.3)
plt.savefig(base_dir / 'Final_plots' / 'Figure_bcp_laplacian.png', dpi=600, bbox_inches='tight')
plt.show()

# --- Plot 2: BCP Laplacian by adsorbate (comparing all metals including VN +U) ---
metal_order_plot2 = ['ScN', 'VN', 'TiN', 'NbN', 'ZrN', 'VN_U']
metal_colors_bcp = dict(zip(metal_order_plot2, sns.color_palette('tab10', len(metal_order_plot2))))
metal_markers = dict(zip(metal_order_plot2, ['o', 's', '^', 'D', 'v', '*']))

fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True, sharex=True)

for idx, ads in enumerate(ads_order_bcp):
    ax = axes[idx]
    for metal in metal_order_plot2:
        subset = bcp_df[(bcp_df['metal'] == metal) & (bcp_df['adsorbate'] == ads)]
        if subset.empty:
            continue
        ax.scatter(
            subset['bond_length_ang'],
            subset['laplacian'],
            s=30,
            alpha=0.7,
            color=metal_colors_bcp[metal],
            marker=metal_markers[metal],
            label=metal.replace('VN_U', 'VN (+U)'),
            edgecolors='none',
        )
    ax.axhline(0, color='grey', linewidth=0.8, linestyle='--')
    ax.annotate(ads, xy=(0.95, 0.95), xycoords='axes fraction',
                fontsize=16, fontweight='bold', ha='right', va='top')
    ax.tick_params(axis='both', labelsize=12, labelbottom=True)
    ax.set_xlabel('Bond length (Angstrom)', fontsize=13)
    if idx == 0:
        ax.set_ylabel('Laplacian (a.u.)', fontsize=13)

handles, labels = axes[0].get_legend_handles_labels()
axes[0].legend(handles, labels, title='Metal', loc='lower left',
               frameon=False, fontsize=9, ncol=2)

fig.suptitle('BCP Laplacian: All Metals (including VN +U)', fontsize=16, fontweight='bold', y=1.01)
plt.tight_layout()
plt.savefig(base_dir / 'Final_plots' / 'Figure_bcp_laplacian_by_ads.png', dpi=600, bbox_inches='tight')
plt.show()
