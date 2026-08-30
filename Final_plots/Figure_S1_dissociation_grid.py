"""
Supplementary Figure S1: top-down xyzrender views for all 18 metal-nitride /
adsorbate combinations (6 metals x 3 adsorbates), same rendering style as
Figure_1_dissociation.py.
"""

import subprocess
import tempfile
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from PIL import Image

BASE_DIR = Path(__file__).parent.parent
OUT_DIR = Path(__file__).parent
PANEL_DIR = OUT_DIR / 'Figure_S1_panels'
PANEL_DIR.mkdir(exist_ok=True)

NO_U_DIR = BASE_DIR / 'no_u/combi_144/final_combi'
U_DIR = BASE_DIR / 'u/combi/final_xyz'

metal_rows = ['ScN', 'VN', 'TiN', 'NbN', 'ZrN', 'VN_U']
ads_cols = ['Li2S4', 'Li2S8', 'S8']

metal_labels = {m: m.replace('VN_U', 'VN (+U)') for m in metal_rows}


def xyz_path(metal, ads):
    if metal == 'VN_U':
        return U_DIR / f'VN_{ads}_final.xyz'
    return NO_U_DIR / f'{metal}_{ads}_combi_final.xyz'


def read_xyz(path):
    with open(path) as f:
        lines = f.readlines()
    n = int(lines[0])
    atoms = []
    for line in lines[2:2 + n]:
        parts = line.split()
        x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
        atoms.append((x, y, z))
    return atoms


def z_gt_16_indices(path):
    atoms = read_xyz(path)
    return [i + 1 for i, (x, y, z) in enumerate(atoms) if z > 16]


XYZRENDER = '/home/ameer_ubuntu/miniforge3/envs/qe/bin/xyzrender'
cmd_base_tail = ['-b', '0', '-B', 'white', '--fog', '-F', '1.7', '--cell', '-S', '400']

LOCATOR = '#ff00ff'  # unique magenta, absent from the real palette
PAD = 36  # px of breathing room around the cell + axis arrows


def cell_bbox(xyz, only):
    """Render a throwaway locator (magenta cell box, same geometry) and
    return the pixel bbox covering the cell box AND the a/b axis arrows,
    the same trick used in Figure_1a_S1a_crop.py."""
    with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tf:
        loc = Path(tf.name)
    cmd = [
        XYZRENDER, str(xyz), '--only', only,
        '-b', '0', '-B', 'white', '--no-fog', '--cell',
        '--cell-color', LOCATOR, '--cell-width', '3',
        '-S', '400', '--axis', '001', '-o', str(loc),
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


panel_files = {}
for metal in metal_rows:
    for ads in ads_cols:
        xyz = xyz_path(metal, ads)
        if not xyz.exists():
            print(f'  MISSING: {xyz}')
            continue
        only = ','.join(str(i) for i in z_gt_16_indices(xyz))
        out_png = PANEL_DIR / f'{metal}_{ads}_top.png'
        cmd = [XYZRENDER, str(xyz), '--only', only, '--axis', '001',
               '-o', str(out_png)] + cmd_base_tail
        print(f'Rendering {metal} / {ads} ...')
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f'  ERROR: {result.stderr.strip()}')
            continue

        box = cell_bbox(xyz, only)
        cropped_png = PANEL_DIR / f'{metal}_{ads}_top_cropped.png'
        Image.open(out_png).convert('RGBA').crop(box).save(cropped_png)
        panel_files[(metal, ads)] = cropped_png
        print(f'  -> {cropped_png}  box={box}')

# --- Assemble 6-row x 3-column montage ---
fig, axes = plt.subplots(len(metal_rows), len(ads_cols),
                          figsize=(3 * len(ads_cols), 3 * len(metal_rows)))

for r, metal in enumerate(metal_rows):
    for c, ads in enumerate(ads_cols):
        ax = axes[r][c]
        png = panel_files.get((metal, ads))
        if png is not None:
            img = mpimg.imread(png)
            ax.imshow(img)
        ax.axis('off')
        if r == 0:
            ax.set_title(ads, fontsize=14, fontweight='bold')
        if c == 0:
            ax.text(-0.05, 0.5, metal_labels[metal], transform=ax.transAxes,
                    fontsize=14, fontweight='bold', ha='right', va='center',
                    rotation=90)

fig.tight_layout(pad=0.3)
fig.savefig(OUT_DIR / 'Figure_S1_dissociation_grid.png', dpi=300, bbox_inches='tight')
plt.close(fig)
print('saved Figure_S1_dissociation_grid.png')
