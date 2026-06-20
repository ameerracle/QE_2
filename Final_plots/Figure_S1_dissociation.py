"""
Artistic render of ScN slab + Li2S8 adsorbate using xyzrender.
Same settings as Figure_1, for ScN_Li2S8.
"""

import subprocess
from pathlib import Path

XYZ = Path(__file__).parent.parent / "no_u/combi_144/final_combi/ScN_Li2S8_combi_final.xyz"
OUT_DIR = Path(__file__).parent
OUT_DIR.mkdir(exist_ok=True)

OUT_PNG_A = OUT_DIR / "Figure_S1a_ScN_Li2S8_slab.png"
OUT_SVG_A = OUT_DIR / "Figure_S1a_ScN_Li2S8_slab.svg"
OUT_PNG_B = OUT_DIR / "Figure_S1b_ScN_Li2S8_slab.png"
OUT_SVG_B = OUT_DIR / "Figure_S1b_ScN_Li2S8_slab.svg"
OUT_PNG_C = OUT_DIR / "Figure_S1c_ScN_Li2S8_slab.png"
OUT_SVG_C = OUT_DIR / "Figure_S1c_ScN_Li2S8_slab.svg"

# Surface atoms with Z > 16 Å (1-indexed) + Li2S8 adsorbate (145-154)
ONLY = "6,7,14,15,22,23,30,31,38,39,46,47,54,55,62,63,70,71,79,80,87,88,95,96,103,104,111,112,119,120,127,128,135,136,143,144,145,146,147,148,149,150,151,152,153,154"

cmd_base = [
    "xyzrender", str(XYZ),
    "--only", ONLY,
    '-b','0',
    "-B", "white",
    "--fog", "-F", "0.85",
    "--cell",
    "-S", "400",
]

# Figure a: top-down, looking down Z (c-axis)
for label, out in [("PNG", OUT_PNG_A), ("SVG", OUT_SVG_A)]:
    print(f"Rendering Figure a {label}...")
    result = subprocess.run(cmd_base + ["--axis", "001", "-o", str(out)], capture_output=True, text=True)
    if result.returncode == 0:
        print(f"  -> {out}")
    else:
        print(f"  ERROR: {result.stderr.strip()}")

# Figure b: looking down Y (b-axis)
for label, out in [("PNG", OUT_PNG_B), ("SVG", OUT_SVG_B)]:
    print(f"Rendering Figure b {label}...")
    result = subprocess.run(cmd_base + ["--axis", "010", "-o", str(out)], capture_output=True, text=True)
    if result.returncode == 0:
        print(f"  -> {out}")
    else:
        print(f"  ERROR: {result.stderr.strip()}")

# Figure c: looking down X (a-axis)
for label, out in [("PNG", OUT_PNG_C), ("SVG", OUT_SVG_C)]:
    print(f"Rendering Figure c {label}...")
    result = subprocess.run(cmd_base + ["--axis", "100", "-o", str(out)], capture_output=True, text=True)
    if result.returncode == 0:
        print(f"  -> {out}")
    else:
        print(f"  ERROR: {result.stderr.strip()}")
