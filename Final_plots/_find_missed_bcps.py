import json
import math
import csv
from pathlib import Path

BOHR = 0.529177249
BASE = Path('/home/ameer_ubuntu/Git_projects/QE_2')

Z2EL = {21: 'Sc', 7: 'N', 16: 'S', 3: 'Li', 23: 'V', 22: 'Ti', 41: 'Nb', 40: 'Zr'}

with open(BASE / 'Final_plots/no_u_headers.json') as f:
    no_u_headers = json.load(f)
with open(BASE / 'Final_plots/u_headers.json') as f:
    u_headers = json.load(f)

metals = ['ScN', 'VN', 'TiN', 'NbN', 'ZrN']
ads_list = ['Li2S4', 'Li2S8', 'S8']

systems = []
for m in metals:
    for a in ads_list:
        path = f'/lustre10/scratch/anizami/QE_2/no_u/PDOS/{m}_{a}/{m}_{a}_charge.cube'
        systems.append(('no_u', m, a, path, no_u_headers[path]))
for a in ads_list:
    path = f'/scratch/anizami/QE_2/u/PDOS/VN_{a}/VN_{a}_charge.cube'
    systems.append(('u', 'VN_U', a, path, u_headers[path]))

CUTOFF = 3.5  # angstrom, generous

def load_existing_bcps(metal, ads, group):
    if group == 'u':
        p = BASE / 'Final_plots/bcp_filtered_u' / f'{metal.replace("VN_U", "VN")}_{ads}_bcp_filtered.csv'
    else:
        p = BASE / 'Final_plots/bcp_filtered' / f'{metal}_{ads}_bcp_filtered.csv'
    rows = []
    with open(p, newline='') as f:
        for row in csv.DictReader(f):
            rows.append(row)
    return rows, p

report = []

for group, metal, ads, path, hdr in systems:
    atoms = hdr['atoms_bohr']  # [Z, x, y, z] in bohr, index = order in file (1-based = idx+1)
    cell = hdr['cell_bohr']
    a_len, b_len = cell[0][0], cell[1][1]

    # classify atoms
    idx_el = []
    for i, (Z, x, y, z) in enumerate(atoms):
        el = Z2EL.get(Z, str(Z))
        idx_el.append((i + 1, el, x, y, z))

    substrate_els = {'Sc', 'V', 'Ti', 'Nb', 'Zr', 'N'}
    adsorbate_els = {'S', 'Li'}

    subs = [t for t in idx_el if t[1] in substrate_els]
    ads_atoms = [t for t in idx_el if t[1] in adsorbate_els]

    existing_rows, csv_path = load_existing_bcps(metal, ads, group)
    existing_pts = []
    for r in existing_rows:
        existing_pts.append((float(r['x_ang']), float(r['y_ang']), float(r['z_ang'])))

    def min_image_dist(x1, y1, z1, x2, y2, z2):
        dx = x1 - x2
        dy = y1 - y2
        dz = z1 - z2
        # minimum image in xy (periodic slab), not in z
        dx -= a_len * round(dx / a_len)
        dy -= b_len * round(dy / b_len)
        return math.sqrt(dx * dx + dy * dy + dz * dz), dx, dy, dz

    candidates = []
    for si, sel, sx, sy, sz in subs:
        for ai, ael, ax, ay, az in ads_atoms:
            d_bohr, dx, dy, dz = min_image_dist(sx, sy, sz, ax, ay, az)
            d_ang = d_bohr * BOHR
            if d_ang < CUTOFF / BOHR * BOHR and d_ang < CUTOFF:
                candidates.append((si, sel, sx, sy, sz, ai, ael, ax - round((ax-sx)/a_len)*0 - 0, ay, az, d_ang))

    # recompute properly with wrapped ads coords for seeding (apply same min-image shift to ads atom)
    final_candidates = []
    for si, sel, sx, sy, sz in subs:
        for ai, ael, ax, ay, az in ads_atoms:
            dx = ax - sx
            dy = ay - sy
            dz = az - sz
            shift_x = a_len * round(dx / a_len)
            shift_y = b_len * round(dy / b_len)
            ax_w = ax - shift_x
            ay_w = ay - shift_y
            d_bohr = math.sqrt((ax_w - sx) ** 2 + (ay_w - sy) ** 2 + (az - sz) ** 2)
            d_ang = d_bohr * BOHR
            if d_ang < CUTOFF:
                final_candidates.append({
                    'sub_idx': si, 'sub_el': sel, 'sub_xyz': (sx, sy, sz),
                    'ads_idx': ai, 'ads_el': ael, 'ads_xyz_wrapped': (ax_w, ay_w, az),
                    'dist_ang': d_ang,
                })

    # check against existing BCPs: is there an existing BCP near the midpoint / line of this pair?
    unmatched = []
    for c in final_candidates:
        sx, sy, sz = c['sub_xyz']
        ax, ay, az = c['ads_xyz_wrapped']
        mx, my, mz = (sx + ax) / 2 * BOHR, (sy + ay) / 2 * BOHR, (sz + az) / 2 * BOHR
        found = False
        for (ex, ey, ez) in existing_pts:
            dd = math.sqrt((ex - mx) ** 2 + (ey - my) ** 2 + (ez - mz) ** 2)
            if dd < 1.3:  # generous tolerance: BCP can be off-midpoint for polar bonds
                found = True
                break
        if not found:
            unmatched.append(c)

    if unmatched:
        report.append({
            'group': group, 'metal': metal, 'ads': ads, 'cube_path': path,
            'cell': cell, 'candidates': unmatched,
        })

with open(BASE / 'Final_plots/_missed_candidates.json', 'w') as f:
    json.dump(report, f, indent=1)

total = sum(len(r['candidates']) for r in report)
print(f'{len(report)} systems with unmatched candidates, {total} total candidate pairs')
for r in report:
    print(f"  {r['metal']}_{r['ads']}: {len(r['candidates'])} candidates")
    for c in r['candidates']:
        print(f"      {c['sub_el']}{c['sub_idx']} -- {c['ads_el']}{c['ads_idx']}  d={c['dist_ang']:.3f} A")
