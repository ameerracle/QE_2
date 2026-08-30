# Recovered from IPython history - session cell ending line 152327
import matplotlib
matplotlib.use('Agg')
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D
from plot_db import write_table

# > Figure: QTAIM Bond Classification Diagram (|V|/G vs ???) ? Interface Bonds Only
# --------------------------------------------------------------------------
plt.rcParams.update({'font.family': 'Liberation Sans', 'mathtext.fontset': 'stix', 
                     'font.weight': 'bold', 'axes.labelweight': 'bold', 'axes.titleweight': 'bold'})

base_dir = Path("/home/ameer_ubuntu/Git_projects/QE_2")
bcp_dir = base_dir / 'Final_plots/bcp_filtered'
bcp_dir_u = base_dir / 'Final_plots/bcp_filtered_u'

metal_order_bcp = ['ScN', 'VN', 'TiN', 'NbN', 'ZrN', 'VN_U']
metal_plot_labels = {m: m.replace('VN_U', 'VN (+U)') for m in metal_order_bcp}
ads_order_bcp = ['Li2S4', 'Li2S8', 'S8']

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

# Keep only interface bonds: adsorbate-involved, exclude intra-adsorbate bonds
bcp_df = bcp_df[bcp_df['atom1'].isin(ADSORBATE_SPECIES) | bcp_df['atom2'].isin(ADSORBATE_SPECIES)]
bcp_df = bcp_df[~bcp_df['bond_type'].isin(INTRA_ADSORBATE)]

# Filter out physically unreasonable outliers that ruin plot scaling
bcp_df = bcp_df[bcp_df['abs_V_over_G'] < 120]

# Drop weak/noise-floor critical points: require meaningful density at the BCP
bcp_df = bcp_df[bcp_df['rho'] > 0.02]

metal_colors_bcp = dict(zip(metal_order_bcp, sns.color_palette('tab10', len(metal_order_bcp))))
ads_markers_bcp = {'Li2S4': 'o', 'Li2S8': 's', 'S8': '^'}

def bold_axis_text(ax):
    ax.xaxis.label.set_fontweight('bold')
    ax.yaxis.label.set_fontweight('bold')
    ax.title.set_fontweight('bold')
    for tick_label in ax.get_xticklabels() + ax.get_yticklabels():
        tick_label.set_fontweight('bold')

fig, ax = plt.subplots(figsize=(10, 8))
ax.text(-0.04, 1.03, '(a)', transform=ax.transAxes,
        fontsize=16, fontweight='bold', ha='left', va='top')

# Plot each point: color = metal, marker = adsorbate, uniform size
for metal in metal_order_bcp:
    for ads in ads_order_bcp:
        subset = bcp_df[(bcp_df['metal'] == metal) & (bcp_df['adsorbate'] == ads)]
        if subset.empty:
            continue
        ax.scatter(
            subset['laplacian'],
            subset['abs_V_over_G'],
            s=50,
            alpha=0.78,
            color=metal_colors_bcp[metal],
            marker=ads_markers_bcp[ads],
            edgecolors='black',
            linewidths=0.35,
        )

# --- Custom legend items ---
metal_handles = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor=metal_colors_bcp[m],
           markersize=10, label=metal_plot_labels[m])
    for m in metal_order_bcp
]
ads_handles = [
    Line2D([0], [0], marker=ads_markers_bcp[ads], color='k',
           markersize=10, label=ads, linestyle='None')
    for ads in ads_order_bcp
]

# Zone boundary lines
ax.axvline(0, color='grey', linewidth=1.0, linestyle='--', zorder=0)
ax.axhline(1, color='grey', linewidth=1.0, linestyle='--', zorder=0)
ax.axhline(2, color='grey', linewidth=1.0, linestyle='--', zorder=0)

# Zone annotations (filled white background, no outline, on top of dots)
# Laplacian-sign labels (x-axis concept) sit just above the x-axis itself,
# straddling the x=0 line (data coords for x, axes-fraction for y) - the
# nearest data to x=0 is at |V|/G >= 1.21, well above this band, so they
# never collide with points and stay inside the plot frame.
_xaxis_t = ax.get_xaxis_transform()
ax.text(-0.005, 0.04, 'Covalent\n(shared-shell)', transform=_xaxis_t,
        fontsize=11, fontweight='bold', ha='right', va='bottom', color='darkgreen',
        clip_on=False, zorder=5)
ax.text(0.005, 0.04, 'Ionic /\nClosed-shell', transform=_xaxis_t,
        fontsize=11, fontweight='bold', ha='left', va='bottom', color='darkred',
        clip_on=False, zorder=5)
# y-axis (|V|/G) zone labels - three bands matching the Espinosa-Molins-
# Lecomte thresholds at 1 and 2. Top band placed top-left since the only
# points above |V|/G=2 sit on the right side of the plot (high Laplacian).
ax.text(0.02, 0.98, 'Strongly covalent\n(shared-shell, |V|/G > 2)', transform=ax.transAxes,
        fontsize=10, fontweight='bold', ha='left', va='top', color='purple',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='none', alpha=0.85), zorder=5)
ax.text(0.02, 0.72, 'Intermediate\n(1 < |V|/G < 2)', transform=ax.transAxes,
        fontsize=10, fontweight='bold', ha='left', va='bottom', color='darkgoldenrod',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='none', alpha=0.85), zorder=5)
ax.text(0.02, 0.18, 'Weak / vdW\n(|V|/G < 1)', transform=ax.transAxes,
        fontsize=10, fontweight='bold', ha='left', va='bottom', color='steelblue',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='none', alpha=0.85), zorder=5)

ax.set_xlabel(r'Laplacian $\nabla^2\rho$ (a.u.)', fontsize=13)
ax.set_ylabel(r'$|V|/G$', fontsize=13)

# Single combined legend at bottom right
all_handles = metal_handles + ads_handles
ax.legend(handles=all_handles, loc='lower right', frameon=False,
          fontsize=10, ncol=3, handletextpad=0.4, columnspacing=0.8)

bold_axis_text(ax)

plt.tight_layout()
plt.savefig(base_dir / 'Final_plots' / 'Figure_5_QTAIM_Bond_Classification.png',
            dpi=600, bbox_inches='tight')
plt.close()

# --- Plot 2: same dots, colored by bond group ---
METALS = {'Sc', 'V', 'Ti', 'Nb', 'Zr'}

def classify_bond_group(row):
    a1, a2 = row['atom1'], row['atom2']
    pair = {a1, a2}
    if pair.intersection(METALS) and 'S' in pair:
        return 'M–S'
    if pair.intersection(METALS) and 'Li' in pair:
        return 'M–Li'
    if 'N' in pair and 'S' in pair:
        return 'N–S'
    if 'N' in pair and 'Li' in pair:
        return 'N–Li'
    return 'Other'

bcp_df['bond_group'] = bcp_df.apply(classify_bond_group, axis=1)
bcp_df_bg = bcp_df[bcp_df['bond_group'] != 'Other']

bond_group_order = ['M–S', 'N–Li', 'N–S']
bond_group_markers = {'M–S': 'o', 'N–Li': 'D', 'N–S': '^'}

# Archive plotted data: BCP descriptors for every plotted interface bond
_db_cols = ['metal', 'adsorbate', 'atom1', 'atom2', 'bond_type',
            'rho', 'laplacian', 'abs_V_over_G', 'bond_group']
write_table(bcp_df[[c for c in _db_cols if c in bcp_df.columns]], 'figure_5_qtaim_bcp')

fig, ax = plt.subplots(figsize=(10, 8))
ax.text(-0.04, 1.03, '(b)', transform=ax.transAxes,
        fontsize=16, fontweight='bold', ha='left', va='top')

for metal in metal_order_bcp:
    for bg in bond_group_order:
        subset = bcp_df_bg[(bcp_df_bg['metal'] == metal) & (bcp_df_bg['bond_group'] == bg)]
        if subset.empty:
            continue
        ax.scatter(
            subset['laplacian'],
            subset['abs_V_over_G'],
            s=50,
            alpha=0.78,
            color=metal_colors_bcp[metal],
            marker=bond_group_markers[bg],
            edgecolors='black',
            linewidths=0.35,
        )

# Zone boundary lines
ax.axvline(0, color='grey', linewidth=1.0, linestyle='--', zorder=0)
ax.axhline(1, color='grey', linewidth=1.0, linestyle='--', zorder=0)
ax.axhline(2, color='grey', linewidth=1.0, linestyle='--', zorder=0)

# Zone annotations
# Laplacian-sign labels (x-axis concept) sit just above the x-axis itself,
# straddling the x=0 line (data coords for x, axes-fraction for y) - the
# nearest data to x=0 is at |V|/G >= 1.21, well above this band, so they
# never collide with points and stay inside the plot frame.
_xaxis_t = ax.get_xaxis_transform()
ax.text(-0.005, 0.04, 'Covalent\n(shared-shell)', transform=_xaxis_t,
        fontsize=11, fontweight='bold', ha='right', va='bottom', color='darkgreen',
        clip_on=False, zorder=5)
ax.text(0.005, 0.04, 'Ionic /\nClosed-shell', transform=_xaxis_t,
        fontsize=11, fontweight='bold', ha='left', va='bottom', color='darkred',
        clip_on=False, zorder=5)
# y-axis (|V|/G) zone labels - three bands matching the Espinosa-Molins-
# Lecomte thresholds at 1 and 2. Top band placed top-left since the only
# points above |V|/G=2 sit on the right side of the plot (high Laplacian).
ax.text(0.02, 0.98, 'Strongly covalent\n(shared-shell, |V|/G > 2)', transform=ax.transAxes,
        fontsize=10, fontweight='bold', ha='left', va='top', color='purple',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='none', alpha=0.85), zorder=5)
ax.text(0.02, 0.72, 'Intermediate\n(1 < |V|/G < 2)', transform=ax.transAxes,
        fontsize=10, fontweight='bold', ha='left', va='bottom', color='darkgoldenrod',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='none', alpha=0.85), zorder=5)
ax.text(0.02, 0.18, 'Weak / vdW\n(|V|/G < 1)', transform=ax.transAxes,
        fontsize=10, fontweight='bold', ha='left', va='bottom', color='steelblue',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='none', alpha=0.85), zorder=5)

ax.set_xlabel(r'Laplacian $\nabla^2\rho$ (a.u.)', fontsize=13)
ax.set_ylabel(r'$|V|/G$', fontsize=13)

metal_handles2 = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor=metal_colors_bcp[m],
           markersize=10, label=metal_plot_labels[m])
    for m in metal_order_bcp
]
bg_handles = [
    Line2D([0], [0], marker=bond_group_markers[bg], color='k',
           markersize=10, label=bg, linestyle='None')
    for bg in bond_group_order
]
ax.legend(handles=metal_handles2 + bg_handles, loc='lower right', frameon=False,
          fontsize=10, ncol=3, handletextpad=0.4, columnspacing=0.8)

bold_axis_text(ax)
plt.tight_layout()
plt.savefig(base_dir / 'Final_plots' / 'Figure_5_QTAIM_Bond_Group.png',
            dpi=600, bbox_inches='tight')
plt.close()

# --- Combined Figure 5: 2 columns, 1 row ---
import matplotlib.image as mpimg

fig, axes = plt.subplots(1, 2, figsize=(20, 8))
panel_files = [base_dir / 'Final_plots' / 'Figure_5_QTAIM_Bond_Classification.png',
               base_dir / 'Final_plots' / 'Figure_5_QTAIM_Bond_Group.png']
panel_labels = ['(a)', '(b)']

for ax, img_path in zip(axes, panel_files):
    img = mpimg.imread(img_path)
    ax.imshow(img)
    ax.set_aspect('auto')
    ax.axis('off')

fig.tight_layout(pad=0.2)
plt.savefig(base_dir / 'Final_plots' / 'Figure_5_combined.png', dpi=600, bbox_inches='tight')
plt.close()
