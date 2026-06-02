#!/usr/bin/env python3
"""Compare numeric rows from two debug output files.

This is intended for quick Fortran-vs-C++ debug logs such as nonbonded pair
rows:

    i j Vela scaling coulomb_constant crg_i crg_j r

Non-numeric lines, blank lines, and separator lines are skipped.
"""

from __future__ import annotations

import argparse
import math
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class Row:
    line_no: int
    values: tuple[float, ...]
    text: str


@dataclass(frozen=True)
class Difference:
    row: int | tuple[float | int, ...]
    column: int
    lhs: float
    rhs: float
    abs_diff: float
    rel_diff: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare numeric rows from two debug files and report the largest differences."
    )
    parser.add_argument("lhs", type=Path, help="First debug file, for example C++ output.")
    parser.add_argument("rhs", type=Path, help="Second debug file, for example Fortran output.")
    parser.add_argument("--atol", type=float, default=1.0e-8, help="Absolute tolerance. Default: 1e-8.")
    parser.add_argument("--rtol", type=float, default=1.0e-6, help="Relative tolerance. Default: 1e-6.")
    parser.add_argument(
        "--key-cols",
        default="",
        help=(
            "Comma-separated zero-based columns used as row keys, for example '0,1' for atom i,j. "
            "By default rows are compared by numeric-row position."
        ),
    )
    parser.add_argument(
        "--compare-key-cols",
        action="store_true",
        help="Also compare key columns as values. By default key columns are not value-compared.",
    )
    parser.add_argument(
        "--sort-key-cols",
        action="store_true",
        help="Sort key-column values before matching, useful for unordered atom pairs i,j vs j,i.",
    )
    parser.add_argument(
        "--ignore-cols",
        default="",
        help="Comma-separated zero-based columns to ignore during value comparison.",
    )
    parser.add_argument(
        "--lhs-key-offset",
        type=int,
        default=0,
        help="Integer offset added to lhs key columns before matching. Use 1 for 0-based C++ atom ids.",
    )
    parser.add_argument(
        "--rhs-key-offset",
        type=int,
        default=0,
        help="Integer offset added to rhs key columns before matching.",
    )
    parser.add_argument(
        "--max-report",
        type=int,
        default=20,
        help="Maximum number of over-tolerance differences to print. Default: 20.",
    )
    return parser.parse_args()


def parse_column_list(value: str) -> set[int]:
    if not value.strip():
        return set()
    columns: set[int] = set()
    for part in value.split(","):
        part = part.strip()
        if not part:
            continue
        columns.add(int(part))
    return columns


def parse_float(token: str) -> float:
    return float(token.replace("D", "E").replace("d", "e"))


def parse_numeric_rows(path: Path) -> list[Row]:
    rows: list[Row] = []
    with path.open("r", encoding="utf-8") as handle:
        for line_no, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            parts = line.split()
            try:
                values = tuple(parse_float(part) for part in parts)
            except ValueError:
                continue
            rows.append(Row(line_no=line_no, values=values, text=line))
    return rows


def normalized_key(row: Row, key_cols: set[int], offset: int, sort_key_cols: bool) -> tuple[float | int, ...]:
    key: list[float | int] = []
    for col in sorted(key_cols):
        value = row.values[col]
        rounded = int(value)
        if math.isclose(value, rounded, rel_tol=0.0, abs_tol=1.0e-12):
            key.append(rounded + offset)
        else:
            key.append(value + offset)
    if sort_key_cols:
        key.sort()
    return tuple(key)


def rows_by_key(
    rows: Iterable[Row], key_cols: set[int], offset: int, sort_key_cols: bool
) -> tuple[dict[tuple[float | int, ...], Row], int]:
    counts: dict[tuple[float | int, ...], int] = defaultdict(int)
    keyed: dict[tuple[float | int, ...], Row] = {}
    skipped_short = 0
    min_width = max(key_cols) + 1
    for row in rows:
        if len(row.values) < min_width:
            skipped_short += 1
            continue
        base_key = normalized_key(row, key_cols, offset, sort_key_cols)
        occurrence = counts[base_key]
        counts[base_key] += 1
        keyed[(*base_key, occurrence)] = row
    return keyed, skipped_short


def relative_diff(lhs: float, rhs: float) -> float:
    scale = max(abs(lhs), abs(rhs), 1.0)
    return abs(lhs - rhs) / scale


def compare_row_values(
    row_id: int | tuple[float | int, ...],
    lhs: Row,
    rhs: Row,
    ignored_cols: set[int],
    atol: float,
    rtol: float,
) -> tuple[list[Difference], int, float, tuple[int, float, float] | None]:
    width = min(len(lhs.values), len(rhs.values))
    differences: list[Difference] = []
    compared = 0
    sum_sq = 0.0
    max_item: tuple[int, float, float] | None = None
    max_abs = -1.0

    for col in range(width):
        if col in ignored_cols:
            continue
        lv = lhs.values[col]
        rv = rhs.values[col]
        diff = abs(lv - rv)
        compared += 1
        sum_sq += diff * diff
        if diff > max_abs:
            max_abs = diff
            max_item = (col, lv, rv)
        if diff > atol + rtol * max(abs(lv), abs(rv)):
            differences.append(
                Difference(
                    row=row_id,
                    column=col,
                    lhs=lv,
                    rhs=rv,
                    abs_diff=diff,
                    rel_diff=relative_diff(lv, rv),
                )
            )
    return differences, compared, sum_sq, max_item


def compare_by_position(
    lhs_rows: list[Row],
    rhs_rows: list[Row],
    ignored_cols: set[int],
    atol: float,
    rtol: float,
) -> tuple[list[Difference], dict[str, float | int]]:
    differences: list[Difference] = []
    compared = 0
    sum_sq = 0.0
    max_abs = 0.0
    max_row = -1
    max_col = -1

    for row_idx, (lhs, rhs) in enumerate(zip(lhs_rows, rhs_rows)):
        row_diffs, row_compared, row_sum_sq, row_max = compare_row_values(
            row_idx, lhs, rhs, ignored_cols, atol, rtol
        )
        differences.extend(row_diffs)
        compared += row_compared
        sum_sq += row_sum_sq
        if row_max is not None:
            col, lv, rv = row_max
            diff = abs(lv - rv)
            if diff > max_abs:
                max_abs = diff
                max_row = row_idx
                max_col = col

    return differences, {
        "lhs_rows": len(lhs_rows),
        "rhs_rows": len(rhs_rows),
        "matched_rows": min(len(lhs_rows), len(rhs_rows)),
        "compared_values": compared,
        "over_tolerance": len(differences),
        "max_abs_diff": max_abs,
        "max_row": max_row,
        "max_col": max_col,
        "rms_abs_diff": math.sqrt(sum_sq / compared) if compared else 0.0,
    }


def compare_by_key(
    lhs_rows: list[Row],
    rhs_rows: list[Row],
    key_cols: set[int],
    ignored_cols: set[int],
    lhs_key_offset: int,
    rhs_key_offset: int,
    sort_key_cols: bool,
    atol: float,
    rtol: float,
) -> tuple[list[Difference], dict[str, float | int]]:
    lhs_keyed, lhs_skipped_short = rows_by_key(lhs_rows, key_cols, lhs_key_offset, sort_key_cols)
    rhs_keyed, rhs_skipped_short = rows_by_key(rhs_rows, key_cols, rhs_key_offset, sort_key_cols)
    matched_keys = sorted(set(lhs_keyed) & set(rhs_keyed))
    missing_lhs = len(set(rhs_keyed) - set(lhs_keyed))
    missing_rhs = len(set(lhs_keyed) - set(rhs_keyed))

    differences: list[Difference] = []
    compared = 0
    sum_sq = 0.0
    max_abs = 0.0
    max_key: tuple[float | int, ...] | None = None
    max_col = -1

    for key in matched_keys:
        lhs = lhs_keyed[key]
        rhs = rhs_keyed[key]
        row_diffs, row_compared, row_sum_sq, row_max = compare_row_values(
            key, lhs, rhs, ignored_cols, atol, rtol
        )
        differences.extend(row_diffs)
        compared += row_compared
        sum_sq += row_sum_sq
        if row_max is not None:
            col, lv, rv = row_max
            diff = abs(lv - rv)
            if diff > max_abs:
                max_abs = diff
                max_key = key
                max_col = col

    return differences, {
        "lhs_rows": len(lhs_rows),
        "rhs_rows": len(rhs_rows),
        "lhs_rows_skipped_short_for_key": lhs_skipped_short,
        "rhs_rows_skipped_short_for_key": rhs_skipped_short,
        "matched_rows": len(matched_keys),
        "missing_lhs": missing_lhs,
        "missing_rhs": missing_rhs,
        "compared_values": compared,
        "over_tolerance": len(differences),
        "max_abs_diff": max_abs,
        "max_row": str(max_key) if max_key is not None else "",
        "max_col": max_col,
        "rms_abs_diff": math.sqrt(sum_sq / compared) if compared else 0.0,
    }


def print_summary(stats: dict[str, float | int]) -> None:
    for key, value in stats.items():
        print(f"{key}: {value}")


def main() -> int:
    args = parse_args()
    key_cols = parse_column_list(args.key_cols)
    ignored_cols = parse_column_list(args.ignore_cols)
    if key_cols and not args.compare_key_cols:
        ignored_cols |= key_cols

    lhs_rows = parse_numeric_rows(args.lhs)
    rhs_rows = parse_numeric_rows(args.rhs)
    if not lhs_rows:
        print(f"No numeric rows found in {args.lhs}", file=sys.stderr)
        return 2
    if not rhs_rows:
        print(f"No numeric rows found in {args.rhs}", file=sys.stderr)
        return 2

    if key_cols:
        differences, stats = compare_by_key(
            lhs_rows,
            rhs_rows,
            key_cols,
            ignored_cols,
            args.lhs_key_offset,
            args.rhs_key_offset,
            args.sort_key_cols,
            args.atol,
            args.rtol,
        )
    else:
        differences, stats = compare_by_position(lhs_rows, rhs_rows, ignored_cols, args.atol, args.rtol)

    print_summary(stats)
    if differences:
        print()
        print(f"First {min(args.max_report, len(differences))} differences over tolerance:")
        for diff in sorted(differences, key=lambda item: item.abs_diff, reverse=True)[: args.max_report]:
            print(
                f"row={diff.row} col={diff.column} "
                f"lhs={diff.lhs:.12g} rhs={diff.rhs:.12g} "
                f"abs={diff.abs_diff:.12g} rel={diff.rel_diff:.12g}"
            )
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
