"""Robust reader for Quantum ESPRESSO .out files with HUBBARD sections.

This parser avoids strict assumptions in ase.io.espresso and instead extracts
geometry snapshots from repeated CELL_PARAMETERS / ATOMIC_POSITIONS blocks.
It is designed for vc-relax outputs from newer QE versions where additional
HUBBARD text blocks may confuse default parsers.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Iterator, List, Optional, Sequence, Tuple, Union

import numpy as np
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import write
from ase.units import Bohr
from ase.utils import string2index

Number = Union[int, float]
IndexLike = Union[int, slice, str, None]


FLOAT_RE = re.compile(r"[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?")
ALAT_RE = re.compile(r"lattice parameter \(alat\)\s*=\s*([-+0-9.eE]+)")
ENERGY_RE = re.compile(r"!\s+total energy\s*=\s*([-+0-9.eE]+)\s+Ry")
BFGS_STEP_RE = re.compile(r"number of bfgs steps\s*=\s*(\d+)", re.IGNORECASE)


class QESnapshotList(list):
    """List of QE snapshots with convenience access to the final image."""

    def get_potential_energy(self) -> float:
        return self[-1].get_potential_energy()

    def get_forces(self) -> np.ndarray:
        return self[-1].get_forces()


def _extract_floats(line: str) -> List[float]:
    return [float(x) for x in FLOAT_RE.findall(line)]


def _normalize_unit(unit: Optional[str], default: str) -> str:
    if unit is None:
        return default
    unit = unit.strip().lower()
    aliases = {
        "angstrom": "angstrom",
        "ang": "angstrom",
        "a": "angstrom",
        "bohr": "bohr",
        "a.u.": "bohr",
        "au": "bohr",
        "alat": "alat",
        "crystal": "crystal",
    }
    return aliases.get(unit, default)


def _to_angstrom(values: np.ndarray, unit: str, alat_ang: Optional[float]) -> np.ndarray:
    if unit == "angstrom":
        return values
    if unit == "bohr":
        return values * Bohr
    if unit == "alat":
        if alat_ang is None:
            raise ValueError("Encountered alat coordinates/cell but alat was not found in output.")
        return values * alat_ang
    raise ValueError(f"Unsupported Cartesian unit: {unit}")


def _to_ev_per_angstrom(values: np.ndarray, unit: str) -> np.ndarray:
    unit = unit.strip().lower()
    if unit in {"ev/angstrom", "ev/ang", "ev/a"}:
        return values
    if unit in {"ry/au", "ry/bohr", "ry/a.u.", "ry/atomic unit"}:
        return values * 13.605693009 / Bohr
    if unit in {"ha/bohr", "hartree/bohr", "hartree/au"}:
        return values * 27.211386018 / Bohr
    raise ValueError(f"Unsupported force unit: {unit}")


def _attach_single_point_results(
    atoms: Atoms,
    *,
    energy_ry: Optional[float] = None,
    forces: Optional[np.ndarray] = None,
) -> None:
    results = {}

    if energy_ry is not None:
        energy_ev = energy_ry * 13.605693009
        atoms.info["energy_ry"] = energy_ry
        atoms.info["energy_ev"] = energy_ev
        results["energy"] = energy_ev
        results["free_energy"] = energy_ev

    if forces is not None:
        results["forces"] = forces

    if results:
        atoms.calc = SinglePointCalculator(atoms, **results)


def _parse_forces_block(lines: Sequence[str], start: int) -> Tuple[Optional[np.ndarray], str, int]:
    line = lines[start].strip()
    match = re.search(r"\(([^)]+)\)", line)
    force_unit = match.group(1).split(",")[-1].strip() if match else "Ry/au"

    forces: List[List[float]] = []
    i = start + 1
    started = False
    while i < len(lines):
        raw = lines[i].strip()
        if not raw:
            if started:
                break
            i += 1
            continue

        lowered = raw.lower()
        if lowered.startswith((
            "number of bfgs steps",
            "cell_parameters",
            "atomic_positions",
            "end final coordinates",
            "writing config-only",
            "k_points",
            "total force",
            "scf correction",
        )):
            break

        if lowered.startswith("atom"):
            vals = _extract_floats(raw)
            if len(vals) < 3:
                break
            forces.append(vals[-3:])
            started = True
            i += 1
            continue

        if started:
            break

        i += 1

    if not forces:
        return None, force_unit, i

    return np.array(forces, dtype=float), force_unit, i


def _parse_cell_block(lines: Sequence[str], start: int) -> Tuple[Optional[np.ndarray], str, int]:
    line = lines[start].strip()
    match = re.search(r"\(([^)]+)\)", line)
    cell_unit = _normalize_unit(match.group(1) if match else None, default="alat")

    mat = []
    i = start + 1
    for _ in range(3):
        if i >= len(lines):
            return None, cell_unit, i
        vals = _extract_floats(lines[i])
        if len(vals) < 3:
            return None, cell_unit, i
        mat.append(vals[:3])
        i += 1

    return np.array(mat, dtype=float), cell_unit, i


def _parse_crystal_axes_block(lines: Sequence[str], start: int) -> Tuple[Optional[np.ndarray], str, int]:
    mat = []
    i = start + 1
    for _ in range(3):
        if i >= len(lines):
            return None, "alat", i
        groups = re.findall(r"\(([^()]*)\)", lines[i])
        vals = _extract_floats(groups[-1]) if groups else []
        if len(vals) < 3:
            return None, "alat", i
        mat.append(vals[:3])
        i += 1

    return np.array(mat, dtype=float), "alat", i


def _parse_positions_block(lines: Sequence[str], start: int) -> Tuple[List[str], np.ndarray, str, int]:
    line = lines[start].strip()
    match = re.search(r"\(([^)]+)\)", line)
    pos_unit = _normalize_unit(match.group(1) if match else None, default="angstrom")

    symbols: List[str] = []
    positions: List[List[float]] = []

    i = start + 1
    while i < len(lines):
        raw = lines[i].strip()
        if not raw:
            break
        if raw.startswith(("CELL_PARAMETERS", "K_POINTS", "ATOMIC_SPECIES", "End final coordinates", "Writing config-only")):
            break

        parts = raw.split()
        if len(parts) < 4:
            break

        symbol = parts[0]
        try:
            xyz = [float(parts[1]), float(parts[2]), float(parts[3])]
        except ValueError:
            break

        symbols.append(symbol)
        positions.append(xyz)
        i += 1

    return symbols, np.array(positions, dtype=float), pos_unit, i


def _parse_qe_out_all(filename: Union[str, Path]) -> List[Atoms]:
    path = Path(filename)
    lines = path.read_text(errors="replace").splitlines()

    snapshots: List[Atoms] = []
    last_cell_raw: Optional[np.ndarray] = None
    last_cell_unit: str = "angstrom"
    alat_ang: Optional[float] = None
    last_energy_ry: Optional[float] = None
    last_forces_raw: Optional[np.ndarray] = None
    last_forces_unit: str = "Ry/au"
    last_bfgs_step: Optional[int] = None

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        m_alat = ALAT_RE.search(line)
        if m_alat:
            # QE prints alat in atomic units (Bohr).
            alat_ang = float(m_alat.group(1)) * Bohr

        m_energy = ENERGY_RE.search(line)
        if m_energy:
            last_energy_ry = float(m_energy.group(1))

        m_bfgs = BFGS_STEP_RE.search(line)
        if m_bfgs:
            last_bfgs_step = int(m_bfgs.group(1))

        if line.startswith("Forces acting on atoms"):
            forces_raw, force_unit, new_i = _parse_forces_block(lines, i)
            if forces_raw is not None:
                last_forces_raw = forces_raw
                last_forces_unit = force_unit
            i = new_i
            continue

        if line.startswith("CELL_PARAMETERS"):
            cell_raw, cell_unit, new_i = _parse_cell_block(lines, i)
            if cell_raw is not None:
                last_cell_raw = cell_raw
                last_cell_unit = cell_unit
            i = new_i
            continue

        if line.startswith("crystal axes:"):
            cell_raw, cell_unit, new_i = _parse_crystal_axes_block(lines, i)
            if cell_raw is not None:
                last_cell_raw = cell_raw
                last_cell_unit = cell_unit
            i = new_i
            continue

        if line.startswith("ATOMIC_POSITIONS"):
            symbols, pos_raw, pos_unit, new_i = _parse_positions_block(lines, i)
            i = new_i
            if len(symbols) == 0:
                continue

            if last_cell_raw is None:
                # No cell seen yet; skip this snapshot.
                continue

            try:
                cell_ang = _to_angstrom(last_cell_raw, last_cell_unit, alat_ang)
            except ValueError:
                continue

            if pos_unit == "crystal":
                atoms = Atoms(symbols=symbols, cell=cell_ang, pbc=True)
                atoms.set_scaled_positions(pos_raw)
            else:
                try:
                    pos_ang = _to_angstrom(pos_raw, pos_unit, alat_ang)
                except ValueError:
                    continue
                atoms = Atoms(symbols=symbols, positions=pos_ang, cell=cell_ang, pbc=True)

            forces_ang = None
            if last_forces_raw is not None and len(last_forces_raw) == len(symbols):
                try:
                    forces_ang = _to_ev_per_angstrom(last_forces_raw, last_forces_unit)
                except ValueError:
                    forces_ang = None

            if last_bfgs_step is not None:
                atoms.info["bfgs_step"] = last_bfgs_step
                atoms.info["force_counter"] = last_bfgs_step

            _attach_single_point_results(atoms, energy_ry=last_energy_ry, forces=forces_ang)

            snapshots.append(atoms)
            continue

        i += 1

    # Some QE runs print an updated final energy after the last geometry block.
    if snapshots and last_energy_ry is not None:
        final_atoms = snapshots[-1]
        final_forces = None
        if final_atoms.calc is not None:
            final_forces = final_atoms.calc.results.get("forces")
        _attach_single_point_results(final_atoms, energy_ry=last_energy_ry, forces=final_forces)

        if last_bfgs_step is not None:
            final_atoms.info["bfgs_step"] = last_bfgs_step
            final_atoms.info["force_counter"] = last_bfgs_step

    return snapshots


def iread_qe_hubbard_out(filename: Union[str, Path], index: IndexLike = ":") -> Iterator[Atoms]:
    snapshots = _parse_qe_out_all(filename)
    if not snapshots:
        raise ValueError(f"No geometry snapshots found in {filename}")

    if isinstance(index, str):
        index = string2index(index)

    if index is None:
        index = slice(None)

    if isinstance(index, int):
        yield snapshots[index]
        return

    for atoms in snapshots[index]:
        yield atoms


def read_qe_hubbard_out(filename: Union[str, Path], index: IndexLike = -1) -> Union[Atoms, List[Atoms]]:
    snapshots = _parse_qe_out_all(filename)
    if not snapshots:
        raise ValueError(f"No geometry snapshots found in {filename}")

    if isinstance(index, str):
        index = string2index(index)

    if index is None:
        index = -1

    if isinstance(index, slice):
        return QESnapshotList(snapshots[index])

    return snapshots[index]


# Convenience alias so you can do: from readqe import read
read = read_qe_hubbard_out


def _main() -> None:
    parser = argparse.ArgumentParser(description="Read QE .out files with HUBBARD sections.")
    parser.add_argument("input", help="Path to QE output file (e.g. vcrelax.out)")
    parser.add_argument("--index", default="-1", help="Image index, e.g. -1, 0, :, -3:")
    parser.add_argument("--write", default=None, help="Optional output structure path (extxyz/cif/etc.)")
    args = parser.parse_args()

    idx: Union[int, slice, str]
    try:
        idx = int(args.index)
    except ValueError:
        idx = args.index

    out = read_qe_hubbard_out(args.input, index=idx)

    if isinstance(out, list):
        print(f"Parsed {len(out)} snapshots from {args.input}")
        if args.write:
            write(args.write, out)
            print(f"Wrote snapshots to {args.write}")
    else:
        print(f"Parsed one snapshot from {args.input}")
        print(f"Formula: {out.get_chemical_formula()} | natoms: {len(out)}")
        print(f"Cell a,b,c (A): {out.cell.lengths()}")
        if "energy_ry" in out.info:
            print(f"Energy: {out.info['energy_ry']:.10f} Ry")
        if args.write:
            write(args.write, out)
            print(f"Wrote structure to {args.write}")


if __name__ == "__main__":
    _main()
