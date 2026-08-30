"""
R2_check.py — cross-descriptor correlation scan over plot_data.db.

Builds one row per (metal, adsorbate) from the figure tables, then computes
pairwise R^2 (squared Pearson r) between every physical descriptor to surface
expected and unexpected correlations across adsorption / charge / QTAIM / ICOHP.
"""
import sqlite3
from pathlib import Path
from itertools import combinations

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

base = Path(__file__).parent
con = sqlite3.connect(base / 'plot_data.db')

COMMON_ADS = ['Li2S4', 'Li2S8', 'S8']


def norm_metal(s):
    return s.replace('VN (+U)', 'VN_U')


# ── adsorption energy + d-band (per metal/adsorbate) ─────────────────────────
ads = pd.read_sql('select * from figure_2_adsorption_energy', con)
ads['metal'] = ads['metal'].map(norm_metal)
ads = ads[ads['adsorbate'].isin(COMMON_ADS)]
ads = ads.rename(columns={'E_adsorption_eV': 'E_ads', 'd_band_center_eV': 'd_band'})
base_tbl = ads[['metal', 'adsorbate', 'E_ads', 'd_band']]

# ── Bader charge: signed dQ per region, plus total transferred to adsorbate ───
bader = pd.read_sql('select * from figure_4_bader_charge', con)
bader['metal'] = bader['metal'].map(norm_metal)
bq = bader.pivot_table(index=['metal', 'adsorbate'], columns='species',
                       values='delta_q', aggfunc='sum').reset_index()
bq = bq.rename(columns={'S': 'dQ_S', 'Li': 'dQ_Li', 'N': 'dQ_N', 'Metal': 'dQ_metal'})
bq['dQ_adsorbate'] = bq.get('dQ_S', 0) + bq.get('dQ_Li', 0)

# ── QTAIM BCP descriptors: mean over M–S interface bonds (the real anchor) ────
bcp = pd.read_sql('select * from figure_5_qtaim_bcp', con)
bcp['metal'] = bcp['metal'].map(norm_metal)
ms = bcp[bcp['bond_group'] == 'M–S']
bcp_ms = ms.groupby(['metal', 'adsorbate']).agg(
    rho_MS=('rho', 'mean'),
    lap_MS=('laplacian', 'mean'),
    VG_MS=('abs_V_over_G', 'mean'),
    n_MS=('rho', 'size'),
    rho_MS_sum=('rho', 'sum'),
).reset_index()

# ── ICOHP: mean |ICOHP| of M–S bonds (bond-strength descriptor) ──────────────
icohp = pd.read_sql('select * from figure_6_icohp_summary', con)
icohp['metal'] = icohp['metal'].map(norm_metal)
ic_ms = icohp[icohp['bond_cat'] == 'M–S'][['metal', 'adsorbate', 'ICOHP_abs']]
ic_ms = ic_ms.rename(columns={'ICOHP_abs': 'ICOHP_MS'})

# ── merge into one descriptor table ──────────────────────────────────────────
df = (base_tbl
      .merge(bq, on=['metal', 'adsorbate'], how='left')
      .merge(bcp_ms, on=['metal', 'adsorbate'], how='left')
      .merge(ic_ms, on=['metal', 'adsorbate'], how='left'))

descriptors = ['E_ads', 'd_band', 'dQ_adsorbate', 'dQ_S', 'dQ_metal', 'dQ_N',
               'rho_MS', 'lap_MS', 'VG_MS', 'rho_MS_sum', 'ICOHP_MS']
descriptors = [d for d in descriptors if d in df.columns]

print(f'Descriptor table: {len(df)} (metal, adsorbate) rows\n')
print(df[['metal', 'adsorbate'] + descriptors].round(3).to_string(index=False))

# ── pairwise R^2 ─────────────────────────────────────────────────────────────
results = []
for a, b in combinations(descriptors, 2):
    sub = df[[a, b]].dropna()
    if len(sub) < 4:
        continue
    r, p = stats.pearsonr(sub[a], sub[b])
    results.append({'x': a, 'y': b, 'n': len(sub), 'r': r, 'R2': r**2, 'p': p})

res = pd.DataFrame(results).sort_values('R2', ascending=False).reset_index(drop=True)

print('\n' + '=' * 70)
print('PAIRWISE CORRELATIONS (sorted by R^2)')
print('=' * 70)
print(res.assign(r=res['r'].round(3), R2=res['R2'].round(3),
                 p=res['p'].round(4)).to_string(index=False))

print('\nNotable (R^2 >= 0.4):')
for _, row in res[res['R2'] >= 0.4].iterrows():
    sign = 'positive' if row['r'] > 0 else 'negative'
    print(f"  {row['x']:>12} vs {row['y']:<12}  R2={row['R2']:.3f}  "
          f"({sign}, r={row['r']:+.2f}, p={row['p']:.3f}, n={int(row['n'])})")

# ── correlation heatmap (signed r) ───────────────────────────────────────────
corr = df[descriptors].corr()
fig, ax = plt.subplots(figsize=(9, 7.5))
im = ax.imshow(corr, cmap='RdBu_r', vmin=-1, vmax=1)
ax.set_xticks(range(len(descriptors)))
ax.set_yticks(range(len(descriptors)))
ax.set_xticklabels(descriptors, rotation=45, ha='right')
ax.set_yticklabels(descriptors)
for i in range(len(descriptors)):
    for j in range(len(descriptors)):
        ax.text(j, i, f'{corr.iloc[i, j]:.2f}', ha='center', va='center',
                fontsize=8, color='black')
fig.colorbar(im, ax=ax, label='Pearson r')
ax.set_title('Descriptor correlation (signed r)')
fig.tight_layout()
fig.savefig(base / 'R2_correlation_heatmap.png', dpi=200, bbox_inches='tight')
print('\nSaved heatmap -> R2_correlation_heatmap.png')

res.to_csv(base / 'R2_correlation_results.csv', index=False)
print('Saved results  -> R2_correlation_results.csv')
