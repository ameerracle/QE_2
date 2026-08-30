# Recovered from IPython history - session cell ending line 144493
# --- Figure 7: COHP (Crystal Orbital Hamilton Population) Analysis ---
# Lobster COHP data: -pCOHP shows bonding/antibonding character vs energy
# Convention: bonding states appear on the right (positive -pCOHP), antibonding on the left (negative)
# X-axis: -pCOHP (states/eV), Y-axis: Energy (eV)
from pathlib import Path
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, NullFormatter

base_dir = Path('/home/ameer_ubuntu/Git_projects/QE_2')
lobster_dir = base_dir / 'Final_plots/lobster_out'

metal_order_cohp = ['ScN', 'VN', 'TiN', 'NbN', 'ZrN', 'VN_U']
ads_order_cohp = ['Li2S4', 'Li2S8', 'S8']
ads_colors_cohp = dict(zip(ads_order_cohp, sns.color_palette('plasma', len(ads_order_cohp))))

def parse_cohpcar(cohpcar_path):
    """Parse COHPCAR.lobster file.
    Returns: energy array, average pCOHP array, and list of individual pCOHP arrays
    """
    with open(cohpcar_path, 'r') as f:
        lines = f.readlines()

    # Line 0: header
    # Line 1: metadata (N_interactions, N_points, E_min, E_max, E_step)
    metadata = lines[1].split()
    n_interactions = int(metadata[0])
    n_points = int(metadata[1])

    # Lines 2 to 2+n_interactions: labels
    # Data starts after labels
    data_start = 2 + n_interactions

    energy = []
    pcohp_avg = []
    ipCohp_avg = []

    for i in range(data_start, len(lines)):
        parts = lines[i].split()
        if len(parts) < 3:
            continue
        energy.append(float(parts[0]))
        pcohp_avg.append(float(parts[1]))
        ipCohp_avg.append(float(parts[2]))

    return np.array(energy), np.array(pcohp_avg), np.array(ipCohp_avg)

# Collect COHP data
cohp_data_dict = {}

for d in sorted(lobster_dir.glob('*_combi')):
    stem = d.name.replace('_combi', '')
    # Handle VN_U folders
    if stem.startswith('VN_U_'):
        metal = 'VN_U'
        ads = stem.replace('VN_U_', '')
    else:
        parts = stem.split('_', 1)
        if len(parts) == 2:
            metal, ads = parts
        else:
            continue

    cohpcar_file = d / 'COHPCAR.lobster'
    if not cohpcar_file.exists():
        continue

    try:
        energy, pcohp_avg, ipCohp_avg = parse_cohpcar(cohpcar_file)
        cohp_data_dict[f'{metal}_{ads}'] = {
            'metal': metal,
            'adsorbate': ads,
            'energy': energy,
            'pcohp_avg': pcohp_avg,
            'ipCohp_avg': ipCohp_avg,
        }
    except Exception as e:
        print(f"Error parsing {cohpcar_file}: {e}")

# Calculate global x-axis limits (using negated values for -pCOHP convention)
all_pcohp_values = []
for data in cohp_data_dict.values():
    all_pcohp_values.extend(-data['pcohp_avg'])  # Negate for -COHP plotting convention
pcohp_min = np.min(all_pcohp_values)
pcohp_max = np.max(all_pcohp_values)
# Add 10% padding
x_padding = (pcohp_max - pcohp_min) * 0.1
x_lim_min = pcohp_min - x_padding
x_lim_max = pcohp_max + x_padding

# --- Plot 1: -pCOHP by metal (2 rows x 3 columns) ---
fig, axes = plt.subplots(2, 3, figsize=(16, 10))
axes_flat = axes.ravel()

for idx, metal in enumerate(metal_order_cohp):
    ax = axes_flat[idx]
    display_name = metal.replace('VN_U', 'VN (+U)')

    for ads in ads_order_cohp:
        key = f'{metal}_{ads}'
        if key not in cohp_data_dict:
            continue

        data = cohp_data_dict[key]
        ax.plot(
            -data['pcohp_avg'],  # X-axis: negated pCOHP (bonding to the right)
            data['energy'],      # Y-axis: Energy
            linewidth=2.0,
            label=ads,
            color=ads_colors_cohp[ads],
            alpha=0.8,
        )

    ax.axvline(0, color='grey', linewidth=0.8, linestyle='--')
    ax.axhline(0, color='black', linewidth=1.0, linestyle='-', alpha=0.3)
    ax.annotate(display_name, xy=(0.03, 0.95), xycoords='axes fraction',
                fontsize=16, fontweight='bold', ha='left', va='top')
    ax.tick_params(axis='both', labelsize=11)
    ax.set_ylabel('Energy (eV)', fontsize=12)
    ax.set_ylim(-12, 8)
    ax.set_xlim(x_lim_min, x_lim_max)

    # Set up y-axis ticks: 4 eV major (labeled), 2 eV minor (unlabeled)
    ax.yaxis.set_major_locator(MultipleLocator(4))
    ax.yaxis.set_minor_locator(MultipleLocator(2))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.tick_params(axis='y', which='major', length=6)
    ax.tick_params(axis='y', which='minor', length=3)

    if idx % 3 == 0:
        ax.set_xlabel('-pCOHP (states/eV)', fontsize=12)

    # Make labels bold
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('bold')

axes_flat[0].legend(title='Adsorbate', loc='lower right', frameon=False, fontsize=10)

plt.tight_layout()
plt.subplots_adjust(hspace=0.25, wspace=0.3)
plt.savefig(base_dir / 'Final_plots' / 'Figure_COHP.png', dpi=600, bbox_inches='tight')
plt.show()

# --- Plot 2: -pCOHP by adsorbate (comparing all metals) ---
metal_colors_cohp = dict(zip(metal_order_cohp, sns.color_palette('tab10', len(metal_order_cohp))))
metal_styles = {
    'ScN': '-',
    'VN': '--',
    'TiN': '-.',
    'NbN': ':',
    'ZrN': '-',
    'VN_U': '--',
}

fig, axes = plt.subplots(1, 3, figsize=(18, 6))

for idx, ads in enumerate(ads_order_cohp):
    ax = axes[idx]

    for metal in metal_order_cohp:
        key = f'{metal}_{ads}'
        if key not in cohp_data_dict:
            continue

        data = cohp_data_dict[key]
        ax.plot(
            -data['pcohp_avg'],  # X-axis: negated pCOHP (bonding to the right)
            data['energy'],      # Y-axis: Energy
            linewidth=2.0,
            label=metal.replace('VN_U', 'VN (+U)'),
            color=metal_colors_cohp[metal],
            linestyle=metal_styles.get(metal, '-'),
            alpha=0.8,
        )

    ax.axvline(0, color='grey', linewidth=0.8, linestyle='--')
    ax.axhline(0, color='black', linewidth=1.0, linestyle='-', alpha=0.3)
    ax.annotate(ads, xy=(0.03, 0.95), xycoords='axes fraction',
                fontsize=16, fontweight='bold', ha='left', va='top')
    ax.tick_params(axis='both', labelsize=11)
    ax.set_ylabel('Energy (eV)', fontsize=12)
    ax.set_ylim(-12, 8)
    ax.set_xlim(x_lim_min, x_lim_max)

    # Set up y-axis ticks: 4 eV major (labeled), 2 eV minor (unlabeled)
    ax.yaxis.set_major_locator(MultipleLocator(4))
    ax.yaxis.set_minor_locator(MultipleLocator(2))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.tick_params(axis='y', which='major', length=6)
    ax.tick_params(axis='y', which='minor', length=3)

    if idx == 0:
        ax.set_xlabel('-pCOHP (states/eV)', fontsize=12)

    # Make labels bold
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('bold')

handles, labels = axes[0].get_legend_handles_labels()
axes[0].legend(handles, labels, title='Metal', loc='lower right', frameon=False, fontsize=10, ncol=1)

plt.tight_layout()
plt.subplots_adjust(hspace=0.25, wspace=0.3)
plt.savefig(base_dir / 'Final_plots' / 'Figure_COHP_by_ads.png', dpi=600, bbox_inches='tight')
plt.show()

# --- Plot 3: ScN only, separated by adsorbate ------------------------------
# A dedicated ScN figure keeps the three adsorbate curves readable without
# competing with the other metal systems.
fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True)

for idx, ads in enumerate(ads_order_cohp):
    ax = axes[idx]
    key = f'ScN_{ads}'

    if key in cohp_data_dict:
        data = cohp_data_dict[key]
        ax.plot(
            -data['pcohp_avg'],
            data['energy'],
            linewidth=2.0,
            color=ads_colors_cohp[ads],
            label=ads,
            alpha=0.8,
        )

    ax.axvline(0, color='grey', linewidth=0.8, linestyle='--')
    ax.axhline(0, color='black', linewidth=1.0, linestyle='-', alpha=0.3)
    ax.annotate(ads, xy=(0.03, 0.95), xycoords='axes fraction',
                fontsize=16, fontweight='bold', ha='left', va='top')
    ax.tick_params(axis='both', labelsize=11)
    ax.set_xlabel('-pCOHP (states/eV)', fontsize=12)
    ax.set_ylim(-12, 8)
    ax.set_xlim(x_lim_min, x_lim_max)
    ax.yaxis.set_major_locator(MultipleLocator(4))
    ax.yaxis.set_minor_locator(MultipleLocator(2))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.tick_params(axis='y', which='major', length=6)
    ax.tick_params(axis='y', which='minor', length=3)

    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('bold')

axes[0].set_ylabel('Energy (eV)', fontsize=12)
fig.suptitle('COHP at ScN–Adsorbate Interfaces', fontsize=16,
             fontweight='bold')
fig.tight_layout()
fig.subplots_adjust(top=0.88, wspace=0.25)
plt.savefig(base_dir / 'Final_plots' / 'Figure_COHP_ScN.png',
            dpi=600, bbox_inches='tight')
plt.show()

print(f"COHP plots saved (using -pCOHP convention, rotated orientation). Data retrieved for {len(cohp_data_dict)} systems.")
print(f"X-axis range: [{x_lim_min:.4f}, {x_lim_max:.4f}] states/eV")
print(f"Bonding states appear to the right (positive -pCOHP), antibonding to the left (negative)")
print(f"Y-axis: Energy from -12 to 8 eV (Fermi level at 0)")
