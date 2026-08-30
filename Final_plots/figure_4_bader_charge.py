# > Figure 4: Bader Charge Redistribution (stacked dQ by species, per metal)
import matplotlib
matplotlib.use('Agg')
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from plot_db import write_table

plt.rcParams.update({'font.family': 'Liberation Sans', 'mathtext.fontset': 'stix',
                     'font.weight': 'bold', 'axes.labelweight': 'bold', 'axes.titleweight': 'bold'})

base_dir = Path('/home/ameer_ubuntu/Git_projects/QE_2')
bader_dir = base_dir / 'Final_plots/bader_output'

metal_order = ['ScN', 'TiN', 'VN', 'VN_U', 'NbN', 'ZrN']
metal_labels = {m: m.replace('VN_U', 'VN (+U)').replace('N', '') if m != 'VN_U' else 'VN (+U)'
                for m in metal_order}
# Clean metal axis labels: Sc, Ti, V, VN (+U), Nb, Zr
metal_labels = {'ScN': 'Sc', 'TiN': 'Ti', 'VN': 'VN', 'VN_U': 'VN (+U)', 'NbN': 'Nb', 'ZrN': 'Zr'}
ads_order = ['Li2S4', 'Li2S8', 'S8']

# Species (stacking order) and colors
species_order = ['S', 'Li', 'N', 'Metal']
species_colors = {'S': '#FFD65C', 'Li': '#A1D99B', 'N': '#FBB4AE', 'Metal': '#A6CEE3'}

# Map region -> species
region_species = {
    'adsorbate_S': 'S',
    'adsorbate_Li': 'Li',
    'surface_N': 'N',
    'surface_metal': 'Metal',
}

# Load and aggregate delta_q per (metal, adsorbate, species)
records = {}
for metal in metal_order:
    for ads in ads_order:
        csv_path = bader_dir / f'{metal}_{ads}_bader.csv'
        if not csv_path.exists():
            continue
        df = pd.read_csv(csv_path)
        df['species'] = df['region'].map(region_species)
        sums = df.groupby('species')['delta_q'].sum()
        records[(metal, ads)] = sums

# Archive plotted data: summed delta_q per metal/adsorbate/species
_db_rows = []
for (metal, ads), sums in records.items():
    for sp in species_order:
        _db_rows.append({'metal': metal, 'adsorbate': ads, 'species': sp,
                         'delta_q': float(sums.get(sp, 0.0))})
write_table(pd.DataFrame(_db_rows), 'figure_4_bader_charge')

# --- Plot ---
fig, ax = plt.subplots(figsize=(13, 7))

bar_width = 0.8
group_gap = 1.0
x_pos = 0.0
xtick_centres = []
ads_labels = {'Li2S4': r'Li$_2$S$_4$', 'Li2S8': r'Li$_2$S$_8$', 'S8': r'S$_8$'}
bar_labels = []  # (x_pos, ads)

for metal in metal_order:
    group_start = x_pos
    n_bars = 0
    for ads in ads_order:
        if (metal, ads) not in records:
            continue
        sums = records[(metal, ads)]
        pos_base = 0.0
        neg_base = 0.0
        for sp in species_order:
            val = sums.get(sp, 0.0)
            if val == 0.0:
                continue
            base = pos_base if val >= 0 else neg_base
            ax.bar(x_pos, val, bottom=base, width=bar_width,
                   color=species_colors[sp], edgecolor='white', linewidth=0.4)
            if val >= 0:
                pos_base += val
            else:
                neg_base += val
        bar_labels.append((x_pos, ads))
        x_pos += 1.0
        n_bars += 1
    if n_bars:
        xtick_centres.append((group_start + x_pos - 1.0) / 2)
    x_pos += group_gap

ax.axhline(0, color='black', linewidth=1.0)
ax.set_xticks([])
ax.set_xlabel('')
ax.set_ylabel(r'Charge Redistribution $\Delta Q$ ($e^-$)', fontsize=13)

# Per-bar adsorbate labels, just below the axis (figure-6 style).
ax.set_xlim(-group_gap / 2, x_pos - group_gap / 2)
for bx, ads in bar_labels:
    ax.text(bx, -0.02, ads_labels[ads], transform=ax.get_xaxis_transform(),
             ha='center', va='top', fontsize=10, fontweight='bold')

# Metal label further down, below the adsorbate labels.
for centre, metal in zip(xtick_centres, metal_order):
    ax.text(centre, -0.10, metal_labels[metal], transform=ax.get_xaxis_transform(),
             ha='center', va='top', fontsize=12, fontweight='bold')

ax.xaxis.set_tick_params(bottom=False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_position(('axes', 0))
for tl in ax.get_yticklabels():
    tl.set_fontweight('bold')

# Legend
handles = [Rectangle((0, 0), 1, 1, facecolor=species_colors[sp], edgecolor='white')
           for sp in species_order]
ax.legend(handles, species_order, ncol=4, loc='upper left', frameon=False,
          fontsize=11, prop={'weight': 'bold'})

plt.tight_layout()
plt.savefig(base_dir / 'Final_plots' / 'Figure_4_Bader_Charge.png', dpi=600, bbox_inches='tight')
plt.close(fig)
