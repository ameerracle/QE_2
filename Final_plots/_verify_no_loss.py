import json
import csv
import math
from pathlib import Path

BASE = Path('/home/ameer_ubuntu/Git_projects/QE_2')
JOBS = BASE / 'Final_plots/_cri_jobs'
BOHR = 0.529177249

with open(JOBS / '_job_index.json') as f:
    job_index = json.load(f)

for job in job_index:
    metal, ads, group = job['metal'], job['ads'], job['group']
    cro_path = Path(job['cro'])
    lines = cro_path.read_text(errors='replace').splitlines()

    # parse all bond cps with positions from "Critical point list, final report"
    try:
        cpstart = next(i for i, l in enumerate(lines) if 'ncp   pg  type' in l)
    except StopIteration:
        print(f'{metal}_{ads}: no cp list found')
        continue
    found_bonds = []  # (x_cryst,y_cryst,z_cryst)
    j = cpstart + 1
    while j < len(lines) and lines[j].strip():
        parts = lines[j].split()
        if len(parts) >= 8 and 'bond' in lines[j]:
            cx, cy, cz = float(parts[4]), float(parts[5]), float(parts[6])
            found_bonds.append((cx, cy, cz))
        j += 1

    # cell for conversion
    with open(BASE / 'Final_plots/no_u_headers.json') as f:
        no_u_headers = json.load(f)
    with open(BASE / 'Final_plots/u_headers.json') as f:
        u_headers = json.load(f)
    if group == 'u':
        cube_path = f'/scratch/anizami/QE_2/u/PDOS/VN_{ads}/VN_{ads}_charge.cube'
        cell = u_headers[cube_path]['cell_bohr']
    else:
        cube_path = f'/lustre10/scratch/anizami/QE_2/no_u/PDOS/{metal}_{ads}/{metal}_{ads}_charge.cube'
        cell = no_u_headers[cube_path]['cell_bohr']
    a_len, b_len, c_len = cell[0][0], cell[1][1], cell[2][2]

    found_ang = []
    for cx, cy, cz in found_bonds:
        x_ang = cx * a_len * BOHR
        y_ang = cy * b_len * BOHR
        z_ang = cz * c_len * BOHR
        found_ang.append((x_ang, y_ang, z_ang))

    # load the ORIGINAL filtered csv as it existed before our run -- but we already
    # appended new rows to the live file, so just check all rows EXCEPT cp_no>=9000 / 155
    if group == 'u':
        csv_path = BASE / 'Final_plots/bcp_filtered_u' / f'VN_{ads}_bcp_filtered.csv'
    else:
        csv_path = BASE / 'Final_plots/bcp_filtered' / f'{metal}_{ads}_bcp_filtered.csv'

    missing = []
    with open(csv_path, newline='') as f:
        for row in csv.DictReader(f):
            cp_no = int(row['cp_no'])
            if cp_no >= 9000 or cp_no == 155:
                continue  # skip our own newly-added rows
            ex, ey, ez = float(row['x_ang']), float(row['y_ang']), float(row['z_ang'])
            best = min(math.sqrt((ex - fx) ** 2 + (ey - fy) ** 2 + (ez - fz) ** 2) for fx, fy, fz in found_ang) if found_ang else 999
            if best > 0.05:
                missing.append((cp_no, row['bond_type'], best))

    status = 'OK' if not missing else f'MISSING {len(missing)}'
    print(f'{metal}_{ads}: {status}  (existing={sum(1 for _ in csv.DictReader(open(csv_path)))}, found_in_rerun={len(found_ang)})')
    for m in missing:
        print(f'    cp{m[0]} {m[1]}  nearest_match_dist={m[2]:.3f}A')
