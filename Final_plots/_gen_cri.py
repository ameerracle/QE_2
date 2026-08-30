import json
from pathlib import Path

BASE = Path('/home/ameer_ubuntu/Git_projects/QE_2')
OUT = BASE / 'Final_plots/_cri_jobs'
OUT.mkdir(exist_ok=True)

with open(BASE / 'Final_plots/_missed_candidates.json') as f:
    report = json.load(f)

job_index = []

for r in report:
    metal, ads, group = r['metal'], r['ads'], r['group']
    cell = r['cell']
    a_len, b_len, c_len = cell[0][0], cell[1][1], cell[2][2]
    cube_path = r['cube_path']
    name = f"{metal}_{ads}"
    cri_lines = [
        f"crystal {cube_path}",
        f"load {cube_path} id rho_tot",
        "reference rho_tot",
    ]
    for c in r['candidates']:
        sx, sy, sz = c['sub_xyz']
        ax, ay, az = c['ads_xyz_wrapped']
        x0 = sx / a_len
        y0 = sy / b_len
        z0 = sz / c_len
        x1 = ax / a_len
        y1 = ay / b_len
        z1 = az / c_len
        cri_lines.append(
            f"auto seed line x0 {x0:.10f} {y0:.10f} {z0:.10f} "
            f"x1 {x1:.10f} {y1:.10f} {z1:.10f} npts 15"
        )
    cri_lines.append("cpreport")
    cri_path = OUT / f"{name}.cri"
    with open(cri_path, 'w') as f:
        f.write("\n".join(cri_lines) + "\n")
    job_index.append({
        'metal': metal, 'ads': ads, 'group': group,
        'cri': str(cri_path), 'cro': str(OUT / f"{name}.cro"),
        'candidates': r['candidates'],
    })

with open(OUT / '_job_index.json', 'w') as f:
    json.dump(job_index, f, indent=1)

print(f'wrote {len(job_index)} .cri files to {OUT}')
