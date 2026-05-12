#!/usr/bin/env python3
"""Generate and optionally run Critic2 BCP/Laplacian analyses.

The driver keeps the Critic2 input minimal:

    crystal <cube>
    load <cube> id rho_tot
    reference rho_tot
    auto
    cpreport LONG
    END

It expects ASCII-safe file names and writes one .cri/.cro pair per system.
"""

from __future__ import annotations

import argparse
import shlex
import subprocess
from dataclasses import dataclass
from pathlib import Path

DEFAULT_METALS = ("TiN", "VN", "ScN", "NbN", "ZrN")
DEFAULT_ADSORBATES = ("Li2S4", "Li2S8", "S8")
DEFAULT_CUBE_SUFFIX = "_charge.cube"
DEFAULT_CUBE_FALLBACK = ".cube"


@dataclass(frozen=True)
class JobSpec:
    metal: str
    adsorbate: str
    stem: str
    cube_path: Path


def ascii_safe(text: str) -> bool:
    try:
        text.encode("ascii")
    except UnicodeEncodeError:
        return False
    return True


def system_stem(metal: str, adsorbate: str) -> str:
    stem = f"{metal}_{adsorbate}"
    if not ascii_safe(stem):
        raise ValueError(f"Non-ASCII system stem: {stem}")
    return stem


def candidate_cube_paths(base_dir: Path, stem: str, cube_suffix: str) -> list[Path]:
    candidates = [
        base_dir / stem / f"{stem}{cube_suffix}",
        base_dir / stem / f"{stem}{DEFAULT_CUBE_FALLBACK}",
        base_dir / f"{stem}{cube_suffix}",
        base_dir / f"{stem}{DEFAULT_CUBE_FALLBACK}",
    ]
    unique_candidates: list[Path] = []
    seen: set[Path] = set()
    for candidate in candidates:
        resolved = candidate.resolve()
        if resolved not in seen:
            unique_candidates.append(candidate)
            seen.add(resolved)
    return unique_candidates


def find_cube_path(base_dir: Path, stem: str, cube_suffix: str) -> Path:
    for candidate in candidate_cube_paths(base_dir, stem, cube_suffix):
        if candidate.exists():
            return candidate

    recursive_matches = sorted({path.resolve() for path in base_dir.glob(f"**/{stem}{cube_suffix}")})
    if not recursive_matches and cube_suffix != DEFAULT_CUBE_FALLBACK:
        recursive_matches = sorted({path.resolve() for path in base_dir.glob(f"**/{stem}{DEFAULT_CUBE_FALLBACK}")})

    if len(recursive_matches) == 1:
        return recursive_matches[0]
    if len(recursive_matches) > 1:
        formatted = "\n  ".join(str(path) for path in recursive_matches)
        raise RuntimeError(f"Multiple cube files match {stem}:\n  {formatted}")

    raise FileNotFoundError(f"Missing cube for {stem} under {base_dir}")


def build_jobs(base_dir: Path, metals: list[str], adsorbates: list[str], cube_suffix: str) -> list[JobSpec]:
    jobs: list[JobSpec] = []
    for metal in metals:
        for adsorbate in adsorbates:
            stem = system_stem(metal, adsorbate)
            cube_path = find_cube_path(base_dir, stem, cube_suffix)
            jobs.append(JobSpec(metal=metal, adsorbate=adsorbate, stem=stem, cube_path=cube_path))
    return jobs


def write_critic2_input(cri_path: Path, cube_path: Path) -> None:
    cube_text = cube_path.as_posix()
    lines = [
        f"crystal {cube_text}",
        f"load {cube_text} id rho_tot",
        "reference rho_tot",
        "auto",
        "cpreport LONG",
        "END",
        "",
    ]
    cri_path.write_text("\n".join(lines), encoding="utf-8")


def run_critic2(critic2_cmd: str, cri_path: Path, cro_path: Path, cwd: Path) -> None:
    cmd = shlex.split(critic2_cmd)
    with cri_path.open("r", encoding="utf-8") as cri_handle, cro_path.open("w", encoding="utf-8") as cro_handle:
        subprocess.run(
            cmd,
            cwd=str(cwd),
            stdin=cri_handle,
            stdout=cro_handle,
            stderr=subprocess.STDOUT,
            check=True,
        )


def prepare_output_dir(job: JobSpec, output_root: Path | None) -> Path:
    if output_root is None:
        return job.cube_path.parent
    return output_root / job.stem


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate/run Critic2 BCP + Laplacian reports for combined systems.")
    parser.add_argument("--base-dir", type=Path, default=Path.cwd(), help="Root directory that contains the system subdirs or cube files.")
    parser.add_argument("--output-root", type=Path, default=None, help="Optional directory for .cri/.cro output. Defaults to each cube's directory.")
    parser.add_argument("--metals", nargs="+", default=list(DEFAULT_METALS), help="Metals to process, e.g. TiN VN ScN NbN ZrN.")
    parser.add_argument("--adsorbates", nargs="+", default=list(DEFAULT_ADSORBATES), help="Adsorbates to process, e.g. Li2S4 Li2S8 S8.")
    parser.add_argument("--cube-suffix", default=DEFAULT_CUBE_SUFFIX, help="Suffix for the charge-density cube file, default: _charge.cube.")
    parser.add_argument("--critic2-cmd", default="critic2", help="Critic2 command to execute.")
    parser.add_argument("--run", action="store_true", help="Run Critic2 after writing the .cri file.")
    parser.add_argument("--write-only", action="store_true", help="Only write .cri files, never execute Critic2.")
    args = parser.parse_args()

    if args.run and args.write_only:
        raise ValueError("Choose either --run or --write-only, not both.")

    base_dir = args.base_dir.resolve()
    output_root = args.output_root.resolve() if args.output_root is not None else None

    jobs = build_jobs(base_dir, list(args.metals), list(args.adsorbates), args.cube_suffix)
    if not jobs:
        print("No systems matched the requested metals/adsorbates.")
        return

    completed: list[str] = []
    failed: list[str] = []

    for job in jobs:
        if not ascii_safe(job.cube_path.name):
            raise ValueError(f"Non-ASCII cube filename: {job.cube_path.name}")

        out_dir = prepare_output_dir(job, output_root)
        out_dir.mkdir(parents=True, exist_ok=True)

        cri_path = out_dir / f"{job.stem}.cri"
        cro_path = out_dir / f"{job.stem}.cro"

        write_critic2_input(cri_path, job.cube_path.resolve())
        print(f"Wrote {cri_path}")

        if args.run:
            try:
                run_critic2(args.critic2_cmd, cri_path, cro_path, out_dir)
            except subprocess.CalledProcessError as exc:
                failed.append(job.stem)
                print(f"[FAIL] {job.stem}: Critic2 exited with {exc.returncode}")
                continue

            completed.append(job.stem)
            print(f"[OK] {job.stem} -> {cro_path}")

    if args.run:
        print(f"Completed {len(completed)} system(s).")
        if failed:
            print(f"Failed {len(failed)} system(s): {', '.join(failed)}")
            raise SystemExit(1)


if __name__ == "__main__":
    main()