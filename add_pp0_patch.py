#!/usr/bin/env python3
import re

filepath = "/lustre10/scratch/anizami/QE_2/adsorbates/adsorbate_pdos.py"
with open(filepath, "r") as f:
    content = f.read()

# Add write_pp_plot0_input function before run_chain
new_func = '''
def write_pp_plot0_input(file_path, tag, filplot, fileout):
    with open(file_path, "w", encoding="utf-8") as f:
        f.write("&INPUTPP\\n")
        f.write(f"  prefix = '{tag}'\\n")
        f.write("  outdir = './tmp/'\\n")
        f.write("  plot_num = 0\\n")
        f.write(f"  filplot = '{filplot}'\\n")
        f.write("/\\n")
        f.write("&PLOT\\n")
        f.write("  iflag = 3\\n")
        f.write("  output_format = 6\\n")
        f.write(f"  fileout = '{fileout}'\\n")
        f.write("/\\n")

'''

content = content.replace("def run_chain(run_dir, tag, np, step, do_pdos):", new_func + "def run_chain(run_dir, tag, np, step, do_pdos):")

# Add PP0 step to run_chain after PP step
old_pp = 'steps.append(("PP", f"srun pp.x < {tag}_pp.in > {tag}_pp.out" if "SLURM_JOB_ID" in os.environ else f"mpirun -np {np} pp.x < {tag}_pp.in > {tag}_pp.out"))'
new_pp = '''steps.append(("PP", f"srun pp.x < {tag}_pp.in > {tag}_pp.out" if "SLURM_JOB_ID" in os.environ else f"mpirun -np {np} pp.x < {tag}_pp.in > {tag}_pp.out"))
        steps.append(("PP0", f"srun pp.x < {tag}_pp_plot0.in > {tag}_pp_plot0.out" if "SLURM_JOB_ID" in os.environ else f"mpirun -np {np} pp.x < {tag}_pp_plot0.in > {tag}_pp_plot0.out"))'''

content = content.replace(old_pp, new_pp)

# Add write_pp_plot0_input call in main() after write_pp_input
old_call = 'write_pp_input(run_dir / f"{tag}_pp.in", tag, f"{tag}_rho", f"{tag}_charge.cube")'
new_call = '''write_pp_input(run_dir / f"{tag}_pp.in", tag, f"{tag}_rho", f"{tag}_charge.cube")
            write_pp_plot0_input(run_dir / f"{tag}_pp_plot0.in", tag, f"{tag}_rho_plot0", f"{tag}_charge_plot0.cube")'''

content = content.replace(old_call, new_call)

with open(filepath, "w") as f:
    f.write(content)

print("Patched successfully")
