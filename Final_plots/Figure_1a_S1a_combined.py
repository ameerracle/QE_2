"""
Combine the top-down (001) renders of VN+S8 (Figure 1a) and ScN+Li2S8
(Figure S1a) into a single 1-row, 2-column figure.
"""
import matplotlib
matplotlib.use('Agg')
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

base_dir = Path(__file__).parent

panel_files = [base_dir / 'Figure_1a_VN_S8_slab.png',
               base_dir / 'Figure_S1a_ScN_Li2S8_slab.png']
panel_labels = ['(a)', '(b)']

fig, axes = plt.subplots(1, 2, figsize=(16, 8))
for ax, label, img_path in zip(axes, panel_labels, panel_files):
    img = mpimg.imread(img_path)
    ax.imshow(img)
    ax.set_aspect('auto')
    ax.axis('off')
    ax.text(-0.02, 1.02, label, transform=ax.transAxes,
            fontsize=18, fontweight='bold', ha='left', va='top')

fig.tight_layout(pad=0.3)
fig.savefig(base_dir / 'Figure_1a_S1a_combined.png', dpi=600, bbox_inches='tight')
plt.close(fig)
print('Saved -> Figure_1a_S1a_combined.png')
