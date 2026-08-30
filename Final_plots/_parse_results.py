import json
import re
import math
from pathlib import Path

BASE = Path('/home/ameer_ubuntu/Git_projects/QE_2')
JOBS = BASE / 'Final_plots/_cri_jobs'
BOHR = 0.529177249

with open(JOBS / '_job_index.json') as f:
    job_index = json.load(f)

bond_re = re.compile(
    r'^\s*\d+\s+(\S+)\s*\((\d+)\)\s+(\S+)\s*\((\d+)\)\s+'
    r'([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*$'
)

cp_line_re = re.compile(
    r'^\s*(\d+)\s+\S+\s+\(3,-1\)\s+bond\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s+\d+\s+\S+\s+'
    r'([\d.eE+-]+)\s+([\d.eE+-]+)\s+([\d.eE+-]+)\s*$'
)

new_bonds = []

for job in job_index:
    cro_path = Path(job['cro'])
    text = cro_path.read_text(errors='replace')
    lines = text.splitlines()

    # find "Analysis of system bonds" table
    try:
        start = next(i for i, l in enumerate(lines) if 'ncp  End-1' in l)
    except StopIteration:
        print(f"WARN: no bond table for {job['metal']}_{job['ads']}")
        continue

    bond_rows = []  # (ncp, el1, idx1, el2, idx2)
    i = start + 1
    while i < len(lines) and lines[i].strip():
        l = lines[i]
        parts = l.split()
        if len(parts) >= 5 and parts[0].isdigit():
            ncp = int(parts[0])
            el1 = parts[1]
            idx1 = int(parts[2].strip('()'))
            el2 = parts[3]
            idx2 = int(parts[4].strip('()'))
            bond_rows.append((ncp, el1, idx1, el2, idx2))
        i += 1

    # candidate index pairs for this system
    cand_pairs = set()
    cand_lookup = {}
    for c in job['candidates']:
        si, ai = c['sub_idx'], c['ads_idx']
        cand_pairs.add((si, ai))
        cand_pairs.add((ai, si))
        cand_lookup[(si, ai)] = c
        cand_lookup[(ai, si)] = c

    matched_ncps = []
    for ncp, el1, idx1, el2, idx2 in bond_rows:
        if (idx1, idx2) in cand_pairs:
            matched_ncps.append((ncp, el1, idx1, el2, idx2))

    if not matched_ncps:
        continue

    # now find the cp data lines (rho, |grad|, laplacian) for these ncp numbers
    # in the "Critical point list, final report" section
    cpdata = {}
    try:
        cpstart = next(i for i, l in enumerate(lines) if 'ncp   pg  type' in l)
    except StopIteration:
        cpstart = None
    if cpstart is not None:
        j = cpstart + 1
        while j < len(lines) and lines[j].strip():
            l = lines[j]
            m = re.match(r'^\s*(\d+)\s+\S+\s+\(3,-1\)\s+bond\s+([\d.\-]+)\s+([\d.\-]+)\s+([\d.\-]+)\s+\d+\s+\S+\s+([\d.eE+\-]+)\s+([\d.eE+\-]+)\s+([\d.eE+\-]+)', l)
            if m:
                ncp = int(m.group(1))
                cx, cy, cz = float(m.group(2)), float(m.group(3)), float(m.group(4))
                rho = float(m.group(5))
                grad = float(m.group(6))
                lap = float(m.group(7))
                cpdata[ncp] = (cx, cy, cz, rho, grad, lap)
            j += 1

    cell = None
    for r_orig in []:
        pass

    for ncp, el1, idx1, el2, idx2 in matched_ncps:
        if ncp not in cpdata:
            print(f"WARN: ncp {ncp} bond row found but no cp data ({job['metal']}_{job['ads']})")
            continue
        cx, cy, cz, rho, grad, lap = cpdata[ncp]
        c = cand_lookup[(idx1, idx2)]
        new_bonds.append({
            'group': job['group'], 'metal': job['metal'], 'ads': job['ads'],
            'el1': el1, 'idx1': idx1, 'el2': el2, 'idx2': idx2,
            'cryst_xyz': (cx, cy, cz), 'rho': rho, 'grad': grad, 'laplacian': lap,
            'dist_ang': c['dist_ang'],
        })

with open(JOBS / '_new_bonds.json', 'w') as f:
    json.dump(new_bonds, f, indent=1)

print(f'Found {len(new_bonds)} new bond critical points among candidates:')
for b in new_bonds:
    print(f"  {b['metal']}_{b['ads']}: {b['el1']}{b['idx1']}-{b['el2']}{b['idx2']}  "
          f"d={b['dist_ang']:.3f}A  rho={b['rho']:.4f}  lap={b['laplacian']:.4f}  grad={b['grad']:.2e}")
