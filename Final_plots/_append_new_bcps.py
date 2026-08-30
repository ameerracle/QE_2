import json
import math
import csv
from pathlib import Path

BASE = Path('/home/ameer_ubuntu/Git_projects/QE_2')
BOHR = 0.529177249

with open(BASE / 'Final_plots/_cri_jobs/_new_bonds.json') as f:
    new_bonds = json.load(f)
with open(BASE / 'Final_plots/no_u_headers.json') as f:
    no_u_headers = json.load(f)
with open(BASE / 'Final_plots/u_headers.json') as f:
    u_headers = json.load(f)

def get_cell(group, metal, ads):
    if group == 'u':
        path = f'/scratch/anizami/QE_2/u/PDOS/VN_{ads}/VN_{ads}_charge.cube'
        return u_headers[path]['cell_bohr']
    else:
        path = f'/lustre10/scratch/anizami/QE_2/no_u/PDOS/{metal}_{ads}/{metal}_{ads}_charge.cube'
        return no_u_headers[path]['cell_bohr']

def abramov(rho, lap):
    pref = (3 / 10) * (3 * math.pi ** 2) ** (2 / 3)
    G = pref * rho ** (5 / 3) + (1 / 6) * lap
    H = -(1 / 4) * lap
    V = H - G
    return G, V, H

# get next cp_no per file to avoid collisions (not critical, just informational)
rows_by_file = {}

for b in new_bonds:
    metal, ads, group = b['metal'], b['ads'], b['group']
    cell = get_cell(group, metal, ads)
    a_len, b_len, c_len = cell[0][0], cell[1][1], cell[2][2]
    cx, cy, cz = b['cryst_xyz']
    x_bohr, y_bohr, z_bohr = cx * a_len, cy * b_len, cz * c_len
    x_ang, y_ang, z_ang = x_bohr * BOHR, y_bohr * BOHR, z_bohr * BOHR
    rho, lap = b['rho'], b['laplacian']
    G, V, H = abramov(rho, lap)
    abs_v_over_g = abs(V) / G

    row = {
        'cp_no': 9000 + len(rows_by_file),
        'cp_type': '"(3,-1)"',
        'bond_type': f"{b['el1']}-{b['el2']}",
        'atom1': b['el1'], 'atom2': b['el2'],
        'x_bohr': x_bohr, 'y_bohr': y_bohr, 'z_bohr': z_bohr,
        'x_ang': x_ang, 'y_ang': y_ang, 'z_ang': z_ang,
        'rho': rho, 'laplacian': lap,
        'G_ked': G, 'V_ped': V, 'H_ted': H, 'abs_V_over_G': abs_v_over_g,
    }

    metal_dir_name = 'VN' if group == 'u' else metal
    targets = []
    if group == 'u':
        targets.append(BASE / 'Final_plots/bcp_filtered_u' / f'{metal_dir_name}_{ads}_bcp_filtered.csv')
    else:
        targets.append(BASE / 'Final_plots/bcp_filtered' / f'{metal}_{ads}_bcp_filtered.csv')
        targets.append(BASE / 'no_u/bcp_filtered' / f'{metal}_{ads}_bcp_filtered.csv')

    cols = ['cp_no', 'cp_type', 'bond_type', 'atom1', 'atom2', 'x_bohr', 'y_bohr', 'z_bohr',
            'x_ang', 'y_ang', 'z_ang', 'rho', 'laplacian', 'G_ked', 'V_ped', 'H_ted', 'abs_V_over_G']
    line = ','.join(str(row[c]) for c in cols)

    for t in targets:
        with open(t, 'rb') as f:
            data = f.read()
        ending = b'\r\n' if b'\r\n' in data else b'\n'
        if not data.endswith(b'\n'):
            data += ending
        data += line.encode() + ending
        with open(t, 'wb') as f:
            f.write(data)
        print(f'appended to {t.relative_to(BASE)}: {row["bond_type"]} rho={rho:.4f} lap={lap:.4f} |V|/G={abs_v_over_g:.3f}')

    rows_by_file[len(rows_by_file)] = row

print(f'\nTotal new rows appended: {len(new_bonds)} bonds x targets each')
