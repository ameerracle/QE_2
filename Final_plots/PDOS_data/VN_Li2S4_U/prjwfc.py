#!/usr/bin/env python3
"""
Generate and run QE Projwfc post-processing for completed SCF/NSCF runs.
Usage: python run_projwfc.py --prefix VN_Li2S4 --run --np 4
"""

from __future__ import annotations
import argparse
import os
import subprocess
from pathlib import Path

def write_projwfc_input(
    file_path: Path, 
    tag: str, 
    outdir: str, 
    filpdos: str, 
    delta_e: float, 
    e_min: float, 
    e_max: float
) -> None:
    """Writes the input deck for projwfc.x."""
    with open(file_path, "w", encoding="utf-8") as f:
        f.write("&PROJWFC\n")
        f.write(f"  prefix = '{tag}'\n")
        f.write(f"  outdir = '{outdir}'\n")
        f.write("  ngauss = -1\n")
        f.write(f"  DeltaE = {delta_e:.4f}\n")
        f.write(f"  Emin = {e_min:.2f}\n")
        f.write(f"  Emax = {e_max:.2f}\n")
        f.write(f"  filpdos = '{filpdos}'\n")
        f.write("/\n")

def run_projwfc(run_dir: Path, tag: str, np: int) -> None:
    """Runs projwfc.x using Slurm srun or local mpirun."""
    if "SLURM_JOB_ID" in os.environ:
        cmd = f"srun projwfc.x < {tag}_projwfc.in > {tag}_projwfc.out"
    else:
        cmd = f"mpirun -np {np} projwfc.x < {tag}_projwfc.in > {tag}_projwfc.out"
        
    print(f"  [PROJWFC] Running: {cmd}")
    subprocess.run(cmd, shell=True, cwd=str(run_dir), check=True)

def main() -> None:
    parser = argparse.ArgumentParser(description="Run projwfc.x over existing QE nscf data")
    
    parser.add_argument("--prefix", default="VN_Li2S4", help="Prefix of the QE calculation files")
    parser.add_argument("--outdir", default="./tmp", help="Path to your outdir where .save/ folders are")
    parser.add_argument("--filpdos", default="VN_Li2S4_pdos", help="Prefix for output PDOS files")
    parser.add_argument("--deltae", type=float, default=0.05, help="Energy step for PDOS")
    parser.add_argument("--emin", type=float, default=-12.0, help="Minimum energy for PDOS")
    parser.add_argument("--emax", type=float, default=8.0, help="Maximum energy for PDOS")
    parser.add_argument("--np", type=int, default=4, help="Number of MPI processes for local runs")
    parser.add_argument("--run", action="store_true", help="Actually execute the projwfc.x binary")

    args = parser.parse_args()

    run_dir = Path(".").resolve()

    input_file = run_dir / f"{args.prefix}_projwfc.in"
    
    print(f"Generating input file: {input_file}")
    write_projwfc_input(
        file_path=input_file,
        tag=args.prefix,
        outdir=args.outdir,
        filpdos=args.filpdos,
        delta_e=args.deltae,
        e_min=args.emin,
        e_max=args.emax
    )

    if args.run:
        run_projwfc(run_dir, args.prefix, args.np)
        print("PDOS generation completed successfully.")
    else:
        print("\nInput file written. Run again with `--run` to execute the calculation.")

if __name__ == "__main__":
    main()
