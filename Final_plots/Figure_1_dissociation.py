"""
Artistic render of VN slab + S8 adsorbates using xyzrender.
"""

import subprocess
from pathlib import Path

XYZ = Path(__file__).parent.parent / "no_u/combi_144/final_combi/VN_S8_combi_final.xyz"
OUT_DIR = Path(__file__).parent / "Figure_0"
OUT_DIR.mkdir(exist_ok=True)

OUT_PNG = OUT_DIR / "VN_S8_slab.png"
OUT_SVG = OUT_DIR / "VN_S8_slab.svg"

cmd_base = [
    "xyzrender", str(XYZ),
    "--vdw",
    "-B", "black",
    "--fog", "-F", "0.6",
    "--grad",
    "--atom-gradient-strength", "1.5",
    "--vdw-gradient-strength", "1.5",
    "--vdw-interlocking",
    "-S", "3000",
]

for label, out in [("PNG", OUT_PNG), ("SVG", OUT_SVG)]:
    print(f"Rendering {label}...")
    result = subprocess.run(cmd_base + ["-o", str(out)], capture_output=True, text=True)
    if result.returncode == 0:
        print(f"  -> {out}")
    else:
        print(f"  ERROR: {result.stderr.strip()}")
