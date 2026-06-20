"""
Crop the top-down (001) renders of VN+S8 (Figure 1a) and ScN+Li2S8 (Figure S1a)
to just the unit-cell box.

Method: re-render a throwaway "locator" copy of each with a unique magenta cell
color (same geometry, only the cell color differs), find the magenta box's pixel
bounds, then crop the real default-colored renders to that box.
"""
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from PIL import Image

base_dir = Path(__file__).parent
xyz_dir = base_dir.parent / 'no_u/combi_144/final_combi'

SURFACE = "6,7,14,15,22,23,30,31,38,39,46,47,54,55,62,63,70,71,79,80,87,88,95,96,103,104,111,112,119,120,127,128,135,136,143,144"

LOCATOR = '#ff00ff'  # unique magenta, absent from the real palette

panels = [
    # (real_png, xyz, only_atoms)
    (base_dir / 'Figure_1a_VN_S8_slab.png',
     xyz_dir / 'VN_S8_combi_final.xyz',
     SURFACE + ",145,146,147,148,149,150,151,152"),
    (base_dir / 'Figure_S1a_ScN_Li2S8_slab.png',
     xyz_dir / 'ScN_Li2S8_combi_final.xyz',
     SURFACE + ",145,146,147,148,149,150,151,152,153,154"),
]


PAD = 36  # px of breathing room around the cell + axis arrows (~2% more)


def cell_bbox(xyz, only):
    """Render a locator; return (left, top, right, bottom) covering the magenta
    cell box AND the red/green a/b axis arrows."""
    with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tf:
        loc = Path(tf.name)
    cmd = [
        "xyzrender", str(xyz), "--only", only,
        '-b', '0', "-B", "white", "--no-fog", "--cell",
        "--cell-color", LOCATOR, "--cell-width", "3",
        "-S", "400", "--axis", "001", "-o", str(loc),
    ]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    arr = np.asarray(Image.open(loc).convert('RGB'))
    r, g, b = arr[..., 0], arr[..., 1], arr[..., 2]
    cell = (r > 180) & (g < 100) & (b > 180)           # magenta cell box
    red_arrow = (r > 150) & (g < 100) & (b < 100)      # 'a' arrow
    green_arrow = (g > 100) & (r < 100) & (b < 100)    # 'b' arrow
    mask = cell | red_arrow | green_arrow
    ys, xs = np.where(mask)
    H, W = arr.shape[:2]
    loc.unlink(missing_ok=True)
    return (max(xs.min() - PAD, 0), max(ys.min() - PAD, 0),
            min(xs.max() + 1 + PAD, W), min(ys.max() + 1 + PAD, H))


cropped = []
for real_png, xyz, only in panels:
    box = cell_bbox(xyz, only)
    img = Image.open(real_png).convert('RGBA').crop(box)
    out = real_png.with_name(real_png.stem + '_cropped.png')
    img.save(out)
    cropped.append(out)
    print(f'  cropped {real_png.name} -> {out.name}  box={box}')

# Combined 1-row, 2-column of the cropped cells
fig, axes = plt.subplots(1, 2, figsize=(14, 7.6))
for ax, img_path in zip(axes, cropped):
    ax.imshow(mpimg.imread(img_path))
    ax.set_aspect('auto')
    ax.axis('off')

# (a) top-left of left panel, (b) top-left of right panel — same y, each on their side
axes[0].text(-0.05, 0.97, '(a)', transform=axes[0].transAxes,
             fontsize=18, fontweight='bold', ha='left', va='top')
axes[1].text(-0.05, 0.97, '(b)', transform=axes[1].transAxes,
             fontsize=18, fontweight='bold', ha='left', va='top')

fig.tight_layout(pad=0.6)
fig.savefig(base_dir / 'Figure_1a_S1a_cropped_combined.png', dpi=600, bbox_inches='tight')
plt.close(fig)
print('Saved -> Figure_1a_S1a_cropped_combined.png')
