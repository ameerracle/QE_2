import numpy as np
import os
import glob
import csv


def get_fermi_energy(nscf_file):
    if not os.path.exists(nscf_file):
        return None
    with open(nscf_file, 'r') as f:
        for line in f:
            if "the Fermi energy is" in line:
                parts = line.split()
                for i, p in enumerate(parts):
                    if p == "is":
                        return float(parts[i+1])
    return None


def calculate_dband_center(pdos_file, e_fermi):
    try:
        data = np.loadtxt(pdos_file)
    except Exception as e:
        print(f"Error reading {pdos_file}: {e}")
        return None

    energy = data[:, 0] - e_fermi

    with open(pdos_file, 'r') as f:
        header = f.readline()
        is_spin = "dw" in header or "spin" in header.lower()

    if is_spin:
        dos = data[:, 1] + data[:, 2]
    else:
        dos = data[:, 1]

    numerator = np.trapz(dos * energy, energy)
    denominator = np.trapz(dos, energy)

    if denominator == 0:
        return 0.0

    return numerator / denominator


# --- Configuration ---
metals = {
    "TiN": "Ti",
    "VN":  "V",
    "ScN": "Sc",
    "ZrN": "Zr",
    "NbN": "Nb",
}

for compound, element in metals.items():
    nscf_filename = f"{compound}_slab_nscf.out"
    slab_dir = f"{compound}_slab"
    pdos_pattern = f"slab_pdos.pdos_atm#*({element})_wfc#*(d)"

    if os.path.exists(nscf_filename):
        nscf_path = nscf_filename
        search_path = pdos_pattern
    elif os.path.exists(os.path.join(slab_dir, nscf_filename)):
        nscf_path = os.path.join(slab_dir, nscf_filename)
        search_path = os.path.join(slab_dir, pdos_pattern)
    else:
        print(f"Could not find {nscf_filename}, skipping {compound}")
        continue

    ef = get_fermi_energy(nscf_path)
    if ef is None:
        print(f"Could not parse Fermi Energy from {nscf_path}, skipping {compound}")
        continue

    print(f"\n--- {compound} d-Band Center Analysis ---")
    print(f"Source: {nscf_path}")
    print(f"Fermi Energy (Ef): {ef:.4f} eV")
    print("-" * 40)

    pdos_files = glob.glob(search_path)
    if not pdos_files:
        print(f"No PDOS files found for {compound} (searched: {search_path})")
        continue

    pdos_files.sort(key=lambda x: int(x.split('#')[1].split('(')[0]))

    rows = []
    for pf in pdos_files:
        fname = os.path.basename(pf)
        atom_id = int(fname.split('#')[1].split('(')[0])
        dbc = calculate_dband_center(pf, ef)
        if dbc is not None:
            print(f"  Atom {atom_id:3}: d-band center = {dbc:8.4f} eV")
            rows.append({"atom_index": atom_id, "d_band_center_eV": round(dbc, 6)})

    if rows:
        csv_path = f"d_band_{compound}.csv"
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["atom_index", "d_band_center_eV"])
            writer.writeheader()
            writer.writerows(rows)
        avg = np.mean([r["d_band_center_eV"] for r in rows])
        print(f"  Average: {avg:.4f} eV -> saved to {csv_path}")
