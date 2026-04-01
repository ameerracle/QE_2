#!/usr/bin/env python3
"""
Parse QE output files to extract final gradient/energy errors and BFGS steps.
Identify structures that already meet a relaxed force convergence threshold.

NEW_THRESHOLD = 0.000778 * 2 = 0.001556 Ry/Bohr
ENERGY_THRESHOLD = 1.0E-4 Ry (default)

Usage: python check_convergence.py [--threshold VALUE] [--verbose]
"""

import re
import argparse
from pathlib import Path
from dataclasses import dataclass

@dataclass
class ConvData:
    """Convergence data extracted from output file."""
    structure: str
    file_path: Path
    gradient_error: float = None
    energy_error: float = None
    bfgs_step: int = None
    converged: bool = False
    reason: str = ""
    file_used: str = ""


def parse_output_file(filepath: Path, grad_threshold: float) -> ConvData:
    """
    Parse QE output file and extract ALL Energy error and Gradient error values.
    Check if ANY point in the file meets the threshold.
    Returns the FIRST point that meets criteria, or the FINAL state if none do.
    """
    data = ConvData(
        structure=filepath.parent.name,
        file_path=filepath,
    )
    
    if not filepath.exists():
        data.reason = "File not found"
        return data
    
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
    except Exception as e:
        data.reason = f"Could not read file: {e}"
        return data
    
    # Extract ALL Energy and Gradient error values verbatim
    energy_errors = []
    gradient_errors = []
    step_numbers = []
    
    for i, line in enumerate(lines):
        # Look for "Energy error =" line
        if "Energy error" in line and "=" in line:
            m = re.search(r'=\s*(\d+\.?\d*[Ee][+-]?\d+)', line)
            if m:
                energy_errors.append((float(m.group(1)), i))
        
        # Look for "Gradient error =" line
        if "Gradient error" in line and "=" in line:
            m = re.search(r'=\s*(\d+\.?\d*[Ee][+-]?\d+)', line)
            if m:
                gradient_errors.append((float(m.group(1)), i))
        
        # Look for "number of bfgs steps" for reference
        if "number of bfgs steps" in line:
            m = re.search(r'=\s*(\d+)', line)
            if m:
                step_numbers.append(int(m.group(1)))
    
    if not energy_errors or not gradient_errors:
        data.reason = "No Energy or Gradient error values found in file"
        return data
    
    # Check ALL pairs in order to find FIRST that meets criteria
    min_count = min(len(energy_errors), len(gradient_errors))
    
    for i in range(min_count):
        energy_error = energy_errors[i][0]
        gradient_error = gradient_errors[i][0]
        bfgs_step = step_numbers[i] if i < len(step_numbers) else i
        
        # Check if meets threshold
        if gradient_error <= grad_threshold and energy_error <= 1.0e-4:
            data.energy_error = energy_error
            data.gradient_error = gradient_error
            data.bfgs_step = bfgs_step
            data.converged = True
            data.reason = f"✓ Meets threshold at step {bfgs_step}"
            return data
    
    # If never met criteria, report FINAL state
    final_energy = energy_errors[-1][0]
    final_gradient = gradient_errors[-1][0]
    final_step = step_numbers[-1] if step_numbers else min_count - 1
    
    data.energy_error = final_energy
    data.gradient_error = final_gradient
    data.bfgs_step = final_step
    data.converged = False
    data.reason = f"Never meets criteria (final: step {final_step})"
    return data


def check_convergence(data: ConvData, grad_threshold: float, energy_threshold: float = 1.0e-4):
    """
    Check if structure meets both gradient and energy thresholds.
    Update data.converged and data.reason.
    """
    if data.gradient_error is None or data.energy_error is None:
        data.converged = False
        return
    
    grad_ok = data.gradient_error <= grad_threshold
    energy_ok = data.energy_error <= energy_threshold
    
    data.converged = grad_ok and energy_ok
    
    if not grad_ok:
        data.reason = f"Gradient {data.gradient_error:.3E} > {grad_threshold:.3E}"
    elif not energy_ok:
        data.reason = f"Energy {data.energy_error:.3E} > {energy_threshold:.3E}"
    else:
        data.reason = "CONVERGED ✓"


def main():
    parser = argparse.ArgumentParser(
        description="Check QE relaxation convergence from output files."
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.00156,
        help="Gradient threshold in Ry/Bohr (default: 0.001556 = 2x0.00078)."
    )
    parser.add_argument(
        "--energy-threshold",
        type=float,
        default=1.0e-4,
        help="Energy threshold in Ry (default: 1.0e-4)."
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print detailed info for all structures."
    )
    args = parser.parse_args()
    
    combi_dir = Path(".")
    
    # Get all combi_relax subdirectories
    relax_dirs = sorted([d for d in combi_dir.iterdir() if d.is_dir() and "_combi_relax" in d.name])
    
    if not relax_dirs:
        print("No *_combi_relax directories found in current directory.")
        return
    
    results = []
    
    print(f"Scanning {len(relax_dirs)} relaxation directories...")
    print(f"Gradient threshold: {args.threshold:.6E} Ry/Bohr")
    print(f"Energy threshold: {args.energy_threshold:.6E} Ry")
    print()
    
    for relax_dir in relax_dirs:
        structure_name = relax_dir.name.replace("_combi_relax", "")
        
        # Look for output files: first .out, then _restart.out
        out_file = relax_dir / f"{structure_name}_combi.out"
        restart_file = relax_dir / f"{structure_name}_combi_restart.out"
        
        data = None
        file_used = None
        
        if out_file.exists():
            data = parse_output_file(out_file, args.threshold)
            file_used = out_file.name
            # If main file didn't converge, try restart
            if not data.converged and not data.reason.startswith("✓"):
                if restart_file.exists():
                    data = parse_output_file(restart_file, args.threshold)
                    file_used = restart_file.name
        elif restart_file.exists():
            data = parse_output_file(restart_file, args.threshold)
            file_used = restart_file.name
        else:
            data = ConvData(
                structure=structure_name,
                file_path=relax_dir,
                reason="No .out or _restart.out file found"
            )
            file_used = "MISSING"
        
        # Store file info
        data.file_used = file_used
        results.append(data)
    
    # Separate converged and not converged
    converged = [r for r in results if r.converged]
    not_converged = [r for r in results if not r.converged]
    
    # Prepare output content
    output_lines = []
    
    output_lines.append("=" * 90)
    output_lines.append(f"CONVERGED: {len(converged)} / {len(results)}")
    output_lines.append("=" * 90)
    
    for data in converged:
        output_lines.append(f"  {data.structure:45} | File: {data.file_used:40}")
        if args.verbose:
            output_lines.append(f"    Gradient: {data.gradient_error:.3E} Ry/Bohr")
            output_lines.append(f"    Energy:   {data.energy_error:.3E} Ry")
            output_lines.append(f"    BFGS step: {data.bfgs_step}")
            output_lines.append("")
        else:
            output_lines.append(f"    Grad: {data.gradient_error:.3E} | Energy: {data.energy_error:.3E} | BFGS: {data.bfgs_step:3d}")
    
    output_lines.append("")
    output_lines.append("=" * 90)
    output_lines.append(f"NOT CONVERGED: {len(not_converged)} / {len(results)}")
    output_lines.append("=" * 90)
    
    for data in not_converged:
        output_lines.append(f"  {data.structure:45} | File: {data.file_used:40}")
        if data.gradient_error is not None:
            if args.verbose:
                output_lines.append(f"    Gradient: {data.gradient_error:.3E} Ry/Bohr")
                output_lines.append(f"    Energy:   {data.energy_error:.3E} Ry")
                output_lines.append(f"    BFGS step: {data.bfgs_step}")
            else:
                output_lines.append(f"    Grad: {data.gradient_error:.3E} | Energy: {data.energy_error:.3E} | BFGS: {data.bfgs_step:3d}")
        output_lines.append(f"    Status: {data.reason}")
    
    output_lines.append("")
    output_lines.append("=" * 90)
    output_lines.append(f"SUMMARY: {len(converged)}/{len(results)} structures ready (no re-run needed)")
    output_lines.append(f"         {len(not_converged)}/{len(results)} structures need tighter convergence or still running")
    output_lines.append("=" * 90)
    
    # Print to console
    for line in output_lines:
        print(line)
    
    # Write to file
    output_file = Path("convergence_report.txt")
    with open(output_file, "w") as f:
        f.write("\n".join(output_lines))
    print(f"\n✓ Report saved to {output_file}")


if __name__ == "__main__":
    main()
