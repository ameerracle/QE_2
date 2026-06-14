# Recovered from IPython history - session cell ending line 143763
#
# > Figure 1 a and b, seperate
import matplotlib
matplotlib.use('Agg')
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams.update({'font.family': 'Liberation Sans', 'mathtext.fontset': 'stix', 'font.weight': 'bold', 'axes.labelweight': 'bold', 'axes.titleweight': 'bold'})

base_dir = Path('/home/ameer_ubuntu/Git_projects/QE_2')
primary_csv = base_dir / 'no_u/combi_144/adsorption_energies_from_outs.csv'
plus_u_csv = base_dir / 'u/combi/final_xyz/adsorption_energies.csv'
dband_csv = base_dir / 'Final_plots/d_band_values.csv'

ads_order = ['Li2S', 'Li2S2', 'Li2S4', 'Li2S6', 'Li2S8', 'S8']

# Create a softer pastel color scheme for the bars, with S8 in honey yellow.
colors = list(sns.color_palette('pastel', len(ads_order)))
colors[-1] = "#FFD65C"
adsorbate_colors = dict(zip(ads_order, colors))

# Load and normalize both datasets.
df_primary = pd.read_csv(primary_csv)
df_plus_u = pd.read_csv(plus_u_csv)

# Keep only VN from the second CSV, then relabel it as VN (+U).
df_plus_u = df_plus_u[df_plus_u['surface'] == 'VN'].copy()

# Keep only the columns we need, and standardize names.
df_primary = df_primary.rename(columns={'slab': 'metal'})[['metal', 'adsorbate', 'E_adsorption_eV']].copy()
df_plus_u = df_plus_u.rename(columns={'surface': 'metal'})[['metal', 'adsorbate', 'E_adsorption_eV']].copy()
df_plus_u['metal'] = df_plus_u['metal'].replace({'VN': 'VN (+U)'})

# Make adsorption energies positive for plotting.
df_primary['E_adsorption_eV'] = df_primary['E_adsorption_eV'].abs()
df_plus_u['E_adsorption_eV'] = df_plus_u['E_adsorption_eV'].abs()

df = pd.concat([df_primary, df_plus_u], ignore_index=True)

# Order metals by lowest average adsorption energy (left to right).
metal_order = df.groupby('metal')['E_adsorption_eV'].mean().sort_values().index.tolist()

df['metal'] = pd.Categorical(df['metal'], categories=metal_order, ordered=True)
df['adsorbate'] = pd.Categorical(df['adsorbate'], categories=ads_order, ordered=True)
df = df.sort_values(['metal', 'adsorbate']).reset_index(drop=True)

#display(df)

# Load d-band center values and compute max per metal.
dband_df = pd.read_csv(dband_csv)
dband_max = dband_df.groupby('compound')['d_band_center_eV'].max().to_dict()

def bold_axis_text(ax):
    ax.xaxis.label.set_fontweight('bold')
    ax.yaxis.label.set_fontweight('bold')
    ax.title.set_fontweight('bold')
    for tick_label in ax.get_xticklabels() + ax.get_yticklabels():
        tick_label.set_fontweight('bold')

# --- Plot 1: grouped bars, metals on x-axis and adsorbates as separate bars ---
fig, ax = plt.subplots(figsize=(14, 6))

metal_positions = range(len(metal_order))
bar_width = 0.12
offsets = [
    (i - (len(ads_order) - 1) / 2) * bar_width
    for i in range(len(ads_order))
]

for ads_idx, ads in enumerate(ads_order):
    heights = []
    for metal in metal_order:
        subset = df[(df['metal'] == metal) & (df['adsorbate'] == ads)]
        heights.append(subset['E_adsorption_eV'].iloc[0] if not subset.empty else float('nan'))

    ax.bar(
        [pos + offsets[ads_idx] for pos in metal_positions],
        heights,
        width=bar_width,
        label=ads,
        color=adsorbate_colors[ads],
    )

ax.set_xticks(list(metal_positions))
ax.set_xticklabels(metal_order)
ax.set_xlabel('Metal')
ax.set_ylabel('Adsorption energy (eV)')

# Annotate each metal with its highest d-band center value.
for pos, metal in enumerate(metal_order):
    val = dband_max.get(metal, None)
    if val is not None:
        ax.annotate(
            f"d-band: {val:.2f} eV",
            xy=(pos, 0),
            xytext=(0, 10),
            textcoords='offset points',
            ha='center',
            va='top',
            fontsize=11,
            fontweight='bold',
        )

bold_axis_text(ax)
ax.legend(title='Adsorbate', ncol=3, frameon=False, loc='upper left', prop={'weight': 'bold'})
plt.tight_layout()
plt.savefig(base_dir / 'Final_plots' / 'Figure_2a.png', dpi=600, bbox_inches='tight')
plt.close(fig)

# --- Plot 2: one subplot per metal, with the six adsorbates on the x-axis ---
fig, axes = plt.subplots(2, 3, figsize=(16, 8), sharey=True)
axes = axes.ravel()

for idx, (ax, metal) in enumerate(zip(axes, metal_order)):
    subset = df[df['metal'] == metal].set_index('adsorbate').reindex(ads_order)
    values = subset['E_adsorption_eV'].tolist()

    ax.bar(ads_order, values, color=[adsorbate_colors[a] for a in ads_order])
    ax.text(0.02, 0.95, metal, transform=ax.transAxes, fontsize=13, fontweight='bold', va='top', ha='left')
    ax.set_xlabel('Adsorbate', labelpad=8)
    ax.tick_params(axis='x', rotation=0)

axes[0].set_ylabel('Adsorption energy (eV)')
axes[3].set_ylabel('Adsorption energy (eV)')
for a in axes:
    bold_axis_text(a)
plt.tight_layout()
plt.savefig(base_dir / 'Final_plots' / 'Figure_2b.png', dpi=600, bbox_inches='tight')
plt.close(fig)

# --- Combined Figure 1: 2 columns, 1 row ---
import matplotlib.image as mpimg

fig, axes = plt.subplots(2, 1, figsize=(14, 12))
panel_labels = ['(a)', '(b)']
panel_files = [base_dir / 'Final_plots' / 'Figure_2a.png',
               base_dir / 'Final_plots' / 'Figure_2b.png']

for ax, label, img_path in zip(axes, panel_labels, panel_files):
    img = mpimg.imread(img_path)
    ax.imshow(img)
    ax.set_aspect('auto')
    ax.axis('off')
    ax.text(-0.04, 1.03, label, transform=ax.transAxes,
            fontsize=16, fontweight='bold', ha='left', va='top')

fig.tight_layout(pad=0.2)
plt.savefig(base_dir / 'Final_plots' / 'Figure_2_combined.png', dpi=600, bbox_inches='tight')
plt.close(fig)
