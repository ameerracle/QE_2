# Recovered from IPython history - session cell ending line 153045
import matplotlib
matplotlib.use('Agg')
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import subplot_mosaic
from matplotlib.patches import ConnectionPatch, Rectangle
from matplotlib.ticker import AutoMinorLocator
import seaborn as sns
from plot_db import write_table

# > Figure: ICOHP Grouped Bar Plot ? Mean |ICOHP| by Metal / Bond Type / Adsorbate
# --------------------------------------------------------------------------
plt.rcParams.update({'font.family': 'Liberation Sans', 'mathtext.fontset': 'stix',
                     'font.weight': 'bold', 'axes.labelweight': 'bold', 'axes.titleweight': 'bold'})

base_dir = Path('/home/ameer_ubuntu/Git_projects/QE_2')
lobster_dir = base_dir / 'Final_plots/lobster_out'

metal_order = ['ScN', 'VN', 'VN_U', 'TiN', 'NbN', 'ZrN']
metal_labels = {m: m.replace('VN_U', 'VN (+U)') for m in metal_order}
ads_order = ['Li2S4', 'Li2S8', 'S8']
bond_cats = ['Li\u2013M', 'S\u2013M', 'Li\u2013N', 'S\u2013N']

ADSORBATE_SPECIES = {'S', 'Li'}
INTRA_ADSORBATE = {'S-S', 'Li-S', 'S-Li', 'Li-Li'}
METALS = {'Sc', 'V', 'Ti', 'Nb', 'Zr'}

def map_bond_cat(row):
    a1, a2 = row['atom1'], row['atom2']
    pair = {a1, a2}
    if 'S' in pair and pair.intersection(METALS):
        return 'S\u2013M'
    if 'Li' in pair and pair.intersection(METALS):
        return 'Li\u2013M'
    if 'N' in pair and 'Li' in pair:
        return 'Li\u2013N'
    if 'N' in pair and 'S' in pair:
        return 'S\u2013N'
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

# Drop weak/far interactions: noise-floor ICOHP and non-bonding distances
icohp_df = icohp_df[icohp_df['ICOHP'].abs() >= 0.05]
icohp_df = icohp_df[icohp_df['distance'] <= 3.0]

summary = icohp_df.groupby(['metal', 'adsorbate', 'bond_cat']).agg(
    ICOHP_mean=('ICOHP', 'mean'),
    n=('ICOHP', 'size'),
).reset_index()
summary['ICOHP_abs'] = summary['ICOHP_mean'].abs()

# Short-contact (<= 2.5 A) mean: drops weak long-range tails that drag the
# full (<= 3.0 A) mean down. Restricting to shorter bonds can only raise or
# hold the mean, so the hollow bar always sits at/above the colored bar.
short_df = icohp_df[icohp_df['distance'] <= 2.5]
short = short_df.groupby(['metal', 'adsorbate', 'bond_cat']).agg(
    ICOHP_mean_short=('ICOHP', 'mean'),
    n_short=('ICOHP', 'size'),
).reset_index()
summary = summary.merge(short, on=['metal', 'adsorbate', 'bond_cat'], how='left')
summary['ICOHP_abs_short'] = summary['ICOHP_mean_short'].abs()
summary['n_short'] = summary['n_short'].fillna(0).astype(int)

# Archive plotted data: per-bond ICOHP and the mean summary that is plotted
write_table(icohp_df[['metal', 'adsorbate', 'atom1', 'atom2', 'bond_type',
                      'bond_cat', 'distance', 'ICOHP']], 'figure_6_icohp_bonds')
write_table(summary, 'figure_6_icohp_summary')

ads_colors = dict(zip(ads_order, sns.color_palette('pastel', len(ads_order))))
ads_colors['S8'] = '#FFD65C'

# Selectively show BOTH the colored (<=3.0 A) and hollow (<=2.5 A) values only
# for the groups with large noise-driven deviations between the two.
DUAL_LABEL = {
    ('ScN', 'Li2S8', 'S–N'),
    ('NbN', 'S8', 'S–M'),
    ('ZrN', 'Li2S4', 'S–M'),
    ('ZrN', 'Li2S8', 'S–M'),
}

# 2-row x 3-column layout for legibility
row0 = [metal_labels[m] for m in metal_order[:3]]
row1 = [metal_labels[m] for m in metal_order[3:]]
fig, axd = subplot_mosaic([row0, row1], figsize=(16, 11), sharey=True)
left_col = {metal_order[0], metal_order[3]}  # leftmost axis in each row

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
            val = row['ICOHP_abs'].values[0]            # mean over d <= 3.0 A
            short = row['ICOHP_abs_short'].values[0]     # mean over d <= 2.5 A (may be NaN)
            n_short = row['n_short'].values[0]
            ax.bar(x_pos, val, color=ads_colors[ads], edgecolor='white', lw=0.5, width=0.7)

            if pd.notna(short) and short > val + 1e-6:
                # hollow outline = denoised short-contact mean, rising above the bar
                ax.bar(x_pos, short, facecolor='none', edgecolor='black',
                       lw=1.2, width=0.7, zorder=3)
                star = r'$^*$' if n_short == 1 else ''
                ax.text(x_pos, short + 0.08, f'{short:.2f}{star}', ha='center',
                        va='bottom', fontsize=10, fontweight='bold')
                # for the big-deviation groups, also label the diluted (<=3.0 A)
                # mean inside the colored fill so both values are visible
                if (metal, ads, bcat) in DUAL_LABEL:
                    ax.text(x_pos, val - 0.08, f'{val:.2f}', ha='center', va='top',
                            fontsize=8, fontweight='bold', color='black',
                            bbox=dict(boxstyle='round,pad=0.12', facecolor='white',
                                      edgecolor='none', alpha=0.7), zorder=4)
            else:
                ax.text(x_pos, val + 0.08, f'{val:.2f}', ha='center', va='bottom',
                        fontsize=10, fontweight='bold')
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

    # Clean up axes ? match reference style: keep left spine + ticks only on the
    # leftmost column of each row
    if metal in left_col:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.tick_params(axis='y', which='major', length=5, left=True)
        ax.tick_params(axis='y', which='minor', length=3, left=True)
    else:
        ax.spines[:].set_visible(False)
        ax.yaxis.set_tick_params(left=False, which='both')
    ax.xaxis.set_tick_params(bottom=False)
    ax.set_xticks([])
    ax.margins(y=0.25)

# Y-axis label on leftmost of each row
for m in left_col:
    axd[metal_labels[m]].set_ylabel('Mean |ICOHP| (eV)', fontsize=13)

# Legend on row 1, column 3 (top-right panel), upper right
handles = [Rectangle((0, 0), 1, 1, facecolor=ads_colors[a], edgecolor='white') for a in ads_order]
axd[metal_labels[metal_order[2]]].legend(
    handles, ads_order,
    title='Adsorbate', ncol=3, loc='upper right',
    bbox_to_anchor=(1, 1), frameon=False, fontsize=10)

# Connecting line across bottom of each row
for left, right in [(metal_order[0], metal_order[2]), (metal_order[3], metal_order[5])]:
    conn = ConnectionPatch(
        xyA=(0, 0), coordsA=axd[metal_labels[left]].transAxes,
        xyB=(1, 0), coordsB=axd[metal_labels[right]].transAxes,
        lw=1.0
    )
    fig.add_artist(conn)

fig.tight_layout(pad=0.3)
fig.subplots_adjust(wspace=0.18, hspace=0.40)
plt.savefig(base_dir / 'Final_plots' / 'Figure_6_ICOHP_Stacked_Bar.png',
            dpi=600, bbox_inches='tight')
