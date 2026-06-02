#!/usr/bin/env python3
"""Compare QGPU and Fortran context debug dumps.

The debug dump is structured, so this parser compares fields semantically
instead of relying on line-by-line text equality.
"""

from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


DEFAULT_CUDA = Path("/home/mcpi-02/code/Q/src/core/cuda_context_debug.txt")
DEFAULT_FORTRAN = Path(
    "/home/mcpi-02/code/qligfepv2-BenchmarkExperiments/perturbations/cdk2/"
    "2.protein/FEP_1h1q_1oiu_lab2_test/f90_input/fortran_context_debug.txt"
)

Vector = tuple[float, float, float]


@dataclass(frozen=True)
class ContextDebug:
    fresh_start: int
    n_atoms: int
    n_atoms_solute: int
    n_patoms: int
    n_qatoms: int
    n_waters: int
    n_molecules: int
    velocities: tuple[Vector, ...]
    dvelocities: tuple[Vector, ...]


class ParseError(RuntimeError):
    pass


class Cursor:
    def __init__(self, path: Path):
        self.path = path
        self.lines = path.read_text().splitlines()
        self.index = 0

    def next_line(self, label: str) -> str:
        if self.index >= len(self.lines):
            raise ParseError(f"{self.path}: expected {label}, but reached EOF")
        line = self.lines[self.index]
        self.index += 1
        return line

    @property
    def line_no(self) -> int:
        return self.index

    def remaining_nonempty(self) -> list[tuple[int, str]]:
        return [
            (i + 1, line)
            for i, line in enumerate(self.lines[self.index :], start=self.index)
            if line.strip()
        ]


def parse_int_line(cur: Cursor, label: str) -> int:
    line = cur.next_line(label)
    parts = line.split()
    if len(parts) != 1:
        raise ParseError(f"{cur.path}:{cur.line_no}: expected one integer for {label}, got {line!r}")
    try:
        return int(parts[0])
    except ValueError as exc:
        raise ParseError(f"{cur.path}:{cur.line_no}: invalid integer for {label}: {line!r}") from exc


def parse_float(token: str) -> float:
    return float(token.replace("D", "E").replace("d", "e"))


def parse_vector(line: str, path: Path, line_no: int, label: str) -> Vector:
    parts = line.split()
    if len(parts) != 3:
        raise ParseError(f"{path}:{line_no}: expected 3 floats for {label}, got {line!r}")
    try:
        return parse_float(parts[0]), parse_float(parts[1]), parse_float(parts[2])
    except ValueError as exc:
        raise ParseError(f"{path}:{line_no}: invalid vector for {label}: {line!r}") from exc


def parse_vector_block(cur: Cursor, label: str) -> tuple[Vector, ...]:
    n_vectors = parse_int_line(cur, f"number of {label}")
    vectors: list[Vector] = []
    for idx in range(n_vectors):
        line = cur.next_line(f"{label}[{idx}]")
        vectors.append(parse_vector(line, cur.path, cur.line_no, f"{label}[{idx}]"))
    return tuple(vectors)


def parse_context_debug(path: Path) -> ContextDebug:
    cur = Cursor(path)
    data = ContextDebug(
        fresh_start=parse_int_line(cur, "fresh_start"),
        n_atoms=parse_int_line(cur, "n_atoms"),
        n_atoms_solute=parse_int_line(cur, "n_atoms_solute"),
        n_patoms=parse_int_line(cur, "n_patoms"),
        n_qatoms=parse_int_line(cur, "n_qatoms"),
        n_waters=parse_int_line(cur, "n_waters"),
        n_molecules=parse_int_line(cur, "n_molecules"),
        velocities=parse_vector_block(cur, "velocities"),
        dvelocities=parse_vector_block(cur, "dvelocities"),
    )
    extra = cur.remaining_nonempty()
    if extra:
        line_no, line = extra[0]
        raise ParseError(f"{path}:{line_no}: unexpected trailing content: {line!r}")
    return data


def compare_scalars(cuda: ContextDebug, fort: ContextDebug) -> list[str]:
    issues: list[str] = []
    for field in (
        "fresh_start",
        "n_atoms",
        "n_atoms_solute",
        "n_patoms",
        "n_qatoms",
        "n_waters",
        "n_molecules",
    ):
        lhs = getattr(cuda, field)
        rhs = getattr(fort, field)
        if lhs != rhs:
            issues.append(f"{field}: cuda={lhs} fortran={rhs}")
    return issues


def compare_vectors(
    name: str,
    cuda: Sequence[Vector],
    fort: Sequence[Vector],
    atol: float,
    rtol: float,
    max_report: int,
) -> tuple[list[str], dict[str, float | int]]:
    issues: list[str] = []
    stats: dict[str, float | int] = {
        "vectors": min(len(cuda), len(fort)),
        "components": 0,
        "over_tolerance": 0,
        "max_abs": 0.0,
        "max_index": -1,
        "max_component": -1,
        "rms": 0.0,
    }

    if len(cuda) != len(fort):
        issues.append(f"{name}.count: cuda={len(cuda)} fortran={len(fort)}")

    sum_sq = 0.0
    n_components = 0
    over_tolerance = 0
    reported = 0
    max_abs = -1.0
    max_index = -1
    max_component = -1

    for idx, (c_vec, f_vec) in enumerate(zip(cuda, fort)):
        for comp_idx, (cv, fv) in enumerate(zip(c_vec, f_vec)):
            diff = abs(cv - fv)
            n_components += 1
            sum_sq += diff * diff
            if diff > max_abs:
                max_abs = diff
                max_index = idx
                max_component = comp_idx
            if diff > atol + rtol * abs(cv):
                over_tolerance += 1
                if reported < max_report:
                    comp = "xyz"[comp_idx]
                    issues.append(
                        f"{name}[{idx}].{comp}: cuda={cv:.12g} fortran={fv:.12g} "
                        f"abs_diff={diff:.6g}"
                    )
                    reported += 1

    stats["components"] = n_components
    stats["over_tolerance"] = over_tolerance
    stats["max_abs"] = max_abs if max_abs >= 0 else 0.0
    stats["max_index"] = max_index
    stats["max_component"] = max_component
    stats["rms"] = math.sqrt(sum_sq / n_components) if n_components else 0.0

    if over_tolerance > reported:
        issues.append(f"... {over_tolerance - reported} more {name} components exceed tolerance")

    return issues, stats


def print_vector_summary(name: str, stats: dict[str, float | int]) -> None:
    print(f"  {name} components:     {stats['components']}")
    print(f"  {name} over tolerance: {stats['over_tolerance']}")
    print(f"  {name} max abs diff:   {stats['max_abs']:.12g}")
    print(f"  {name} RMS diff:       {stats['rms']:.12g}")
    if stats["max_index"] != -1:
        comp = "xyz"[int(stats["max_component"])]
        print(f"  {name} max location:   index {stats['max_index']} component {comp}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cuda", nargs="?", type=Path, default=DEFAULT_CUDA)
    parser.add_argument("fortran", nargs="?", type=Path, default=DEFAULT_FORTRAN)
    parser.add_argument("--atol", type=float, default=1e-8, help="absolute tolerance for vector components")
    parser.add_argument("--rtol", type=float, default=0.0, help="relative tolerance for vector components")
    parser.add_argument("--max-report", type=int, default=20, help="maximum detailed mismatches per vector block")
    args = parser.parse_args()

    try:
        cuda = parse_context_debug(args.cuda)
        fort = parse_context_debug(args.fortran)
    except (OSError, ParseError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    scalar_issues = compare_scalars(cuda, fort)
    velocity_issues, velocity_stats = compare_vectors(
        "velocities", cuda.velocities, fort.velocities, args.atol, args.rtol, args.max_report
    )
    dvelocity_issues, dvelocity_stats = compare_vectors(
        "dvelocities", cuda.dvelocities, fort.dvelocities, args.atol, args.rtol, args.max_report
    )

    print(f"cuda:    {args.cuda}")
    print(f"fortran: {args.fortran}")
    print()
    print("Summary")
    print(f"  scalar mismatches:      {len(scalar_issues)}")
    print_vector_summary("velocity", velocity_stats)
    print_vector_summary("dvelocity", dvelocity_stats)

    sections = [
        ("Scalars", scalar_issues),
        ("Velocities", velocity_issues),
        ("DVelocities", dvelocity_issues),
    ]
    for title, issues in sections:
        if not issues:
            continue
        print()
        print(title)
        for issue in issues:
            print(f"  - {issue}")

    if scalar_issues or velocity_issues or dvelocity_issues:
        return 1

    print()
    print("OK: debug contexts match within tolerance")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
