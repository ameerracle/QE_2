# Recovered from IPython history - session cell ending line 152840
from pathlib import Path
import matplotlib.image as mpimg
import matplotlib.pyplot as plt

# > Figure S3: PDOS S8 comparison plots (all metals) ? Supplemental
base_dir = Path('/home/ameer_ubuntu/Git_projects/QE_2/Final_plots/PDOS_Plots/comparison')
out_file = Path('/home/ameer_ubuntu/Git_projects/QE_2/Final_plots/Figure_S3_PDOS_S8.png')

metal_order = ['Sc', 'Ti', 'V', 'V_U', 'Nb', 'Zr']

plot_files = [base_dir / f'{m}_S8_comparison.png' for m in metal_order]

n_plots = len(plot_files)
n_cols = 2
n_rows = (n_plots + n_cols - 1) // n_cols

fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 6 * n_rows))
axes = axes.ravel()

panel_labels = [f'({chr(97 + i)})' for i in range(n_plots)]

for ax, panel_label, img_path in zip(axes, panel_labels, plot_files):
    if img_path.exists():
        image = mpimg.imread(img_path)
        ax.imshow(image)
        ax.set_aspect('auto') 
    else:
        ax.text(0.5, 0.5, f'Missing\n{img_path.name}', ha='center', va='center', fontsize=12)
    ax.axis('off')
    ax.text(-0.04, 1.03, panel_label, transform=ax.transAxes,
            fontsize=16, fontweight='bold', ha='left', va='top')

for ax in axes[n_plots:]:
    ax.axis('off')

fig.tight_layout(pad=0.8)
fig.savefig(out_file, dpi=600, bbox_inches='tight')
plt.show()
plt.close(fig)
