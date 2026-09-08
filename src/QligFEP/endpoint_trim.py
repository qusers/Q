"""Prepare and summarize an explicitly truncated two-state QFEP analysis.

Raw energies and trajectories are never edited. Lambda labels come from the
original MD inputs, not filename digits or rounded legacy output. This is an
analysis workaround, not a soft-core implementation or an overlap diagnostic.
"""
from __future__ import annotations

import argparse
from decimal import Decimal, InvalidOperation
import hashlib
import json
import math
from pathlib import Path
import shlex


def _hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _tokens(line: str) -> list[str]:
    lexer = shlex.shlex(line, posix=True)
    lexer.whitespace_split = True
    lexer.commenters = "!#"
    return list(lexer)


def _number(value: str) -> Decimal:
    try:
        result = Decimal(value.lower().replace("d", "e"))
    except InvalidOperation as error:
        raise ValueError(f"Invalid numeric value: {value}") from error
    if not result.is_finite():
        raise ValueError(f"Nonfinite value: {value}")
    return result


def _md_windows(directory: Path) -> dict[str, dict]:
    windows = {}
    for path in sorted(directory.glob("*.inp")):
        sections: dict[str, list[list[str]]] = {}
        section = None
        for line in path.read_text().splitlines():
            tokens = _tokens(line)
            if not tokens:
                continue
            if tokens[0].startswith("["):
                section = tokens[0].lower()
            elif section:
                sections.setdefault(section, []).append(tokens)
        energy = [row[1:] for row in sections.get("[files]", []) if row[0].lower() == "energy"]
        if not energy or "[lambdas]" not in sections:
            continue
        if len(energy) != 1 or len(energy[0]) != 1:
            raise ValueError(f"Ambiguous energy declaration: {path}")
        weights = [_number(v) for row in sections["[lambdas]"] for v in row]
        if len(weights) != 2 or any(v < 0 or v > 1 for v in weights) or sum(weights) != 1:
            raise ValueError(f"Expected normalized two-state mapping weights in {path}")
        # QFEP prints lambda with six decimals; finer labels cannot be checked
        # reliably against the resulting text summary by this temporary utility.
        if any(v != v.quantize(Decimal("0.000001")) for v in weights):
            raise ValueError(f"Mapping weights exceed QFEP's six-decimal output precision: {path}")
        name = Path(energy[0][0]).name
        if name in windows:
            raise ValueError(f"Multiple MD inputs declare energy basename {name}")
        windows[name] = {"lambdas": [str(v) for v in weights],
                         "md_input": str(path.resolve()), "md_input_sha256": _hash(path)}
    return windows


def prepare(qfep_input: Path, md_input_dir: Path, output_dir: Path,
            lambda_min: str | None = None, lambda_max: str | None = None,
            energy_dir: Path | None = None) -> dict:
    """Exclude exact endpoint files, optionally selecting a fixed interior interval."""
    qfep_input = qfep_input.resolve()
    energy_base = energy_dir.resolve() if energy_dir is not None else qfep_input.parent
    output_dir = output_dir.resolve()
    if output_dir.exists():
        raise FileExistsError(f"Output directory already exists: {output_dir}")
    lines = [line for line in qfep_input.read_text().splitlines() if _tokens(line)]
    if not lines:
        raise ValueError("Empty QFEP input")
    count = int(_tokens(lines[0])[0])
    if count < 2 or len(lines)-count < 8:
        raise ValueError("Invalid QFEP energy-file count/header")
    header, filenames = lines[:-count], lines[-count:]
    state_control = _tokens(header[1])
    if state_control != ["2", "0"] or _tokens(header[6]) != ["0"]:
        raise ValueError("Only two-state QFEP without off-diagonal coupling is supported")
    analysis_control = _tokens(header[2])
    if len(analysis_control) < 2 or (len(analysis_control) > 2 and int(analysis_control[2]) != 0):
        raise ValueError("Full-Hamiltonian QFEP (gas=0) is required")
    if _number(analysis_control[0]) <= 0 or int(analysis_control[1]) < 0:
        raise ValueError("QFEP requires positive kT and nonnegative equilibration discard")
    if (lambda_min is None) != (lambda_max is None):
        raise ValueError("Supply both --lambda-min and --lambda-max, or neither")
    lower = _number(lambda_min) if lambda_min is not None else Decimal(0)
    upper = _number(lambda_max) if lambda_max is not None else Decimal(1)
    if not 0 <= lower < upper <= 1:
        raise ValueError("Lambda bounds must satisfy 0 <= min < max <= 1")
    if lambda_min is not None and (lower == 0 or upper == 1):
        raise ValueError("Explicit interval bounds must be strictly interior")
    known = _md_windows(md_input_dir)
    windows = []
    for line in filenames:
        tokens = _tokens(line)
        if len(tokens) != 1:
            raise ValueError(f"Expected a single energy filename: {line}")
        source = (energy_base / tokens[0]).resolve()
        if source.name not in known:
            raise ValueError(f"No unique original MD mapping for {source.name}")
        window = {**known[source.name], "source_energy": str(source)}
        weight = _number(window["lambdas"][0])
        window["excluded_reason"] = ("exact_endpoint" if weight in (0, 1)
            else "outside_explicit_interval" if weight < lower or weight > upper else None)
        windows.append(window)
    labels = [_number(w["lambdas"][0]) for w in windows]
    differences = [b-a for a,b in zip(labels,labels[1:])]
    if not (all(d > 0 for d in differences) or all(d < 0 for d in differences)):
        raise ValueError("Original QFEP windows must be unique and strictly monotonic")
    retained = [w for w in windows if w["excluded_reason"] is None]
    if len(retained) < 2:
        raise ValueError("At least two sampled interior windows are required")
    retained_labels = [_number(w["lambdas"][0]) for w in retained]
    if lambda_min is not None and (min(retained_labels) != lower or max(retained_labels) != upper):
        raise ValueError("Explicit interval endpoints must be present in the sampled ladder")
    for window in retained:
        source = Path(window["source_energy"])
        if not source.is_file() or source.stat().st_size == 0:
            raise ValueError(f"Missing or empty retained energy file: {source}")
        window["source_energy_sha256"] = _hash(source)
    start, end = retained[0]["lambdas"], retained[-1]["lambdas"]
    scale = _number(end[1])-_number(start[1])
    manifest = {
        "schema_version": 1, "analysis": "truncated_two_state_BAR",
        "full_endpoint_free_energy": False,
        "source_qfep_input": str(qfep_input), "source_qfep_input_sha256": _hash(qfep_input),
        "energy_base": str(energy_base),
        "kT_kcal_mol": str(_number(analysis_control[0])), "gas": 0,
        "start_lambdas": start, "end_lambdas": end,
        "state_constant_gap_multiplier": str(scale),
        "state_constant_rule": "For a separately added state constant g, add multiplier*(g2-g1); never add it twice if already in the energy files.",
        "retained": retained, "excluded": [w for w in windows if w["excluded_reason"]],
        "qualification": "Overlap, ESS, stationarity and cross-leg/state-identity checks remain required. No high-energy frames or interior pairs were removed.",
    }
    output_dir.mkdir(parents=True)
    names = []
    for i, window in enumerate(retained):
        name = f"window_{i:04d}.en"  # Within QFEP's 80-character filename limit.
        (output_dir / name).symlink_to(window["source_energy"])
        names.append(name)
    text = "\n".join([str(len(retained)), *header[1:], *names])+"\n"
    (output_dir / "qfep.inp").write_text(text)
    manifest["prepared_qfep_input_sha256"] = _hash(output_dir / "qfep.inp")
    (output_dir / "endpoint-trim.json").write_text(json.dumps(manifest, indent=2)+"\n")
    return manifest


def summarize(directory: Path) -> dict:
    """Read BAR's final cumulative value from a separately rerun trimmed input."""
    manifest = json.loads((directory / "endpoint-trim.json").read_text())
    if _hash(directory / "qfep.inp") != manifest["prepared_qfep_input_sha256"]:
        raise ValueError("Prepared QFEP input changed after interval selection")
    for i, window in enumerate(manifest["retained"]):
        if _hash(directory / f"window_{i:04d}.en") != window["source_energy_sha256"]:
            raise ValueError("Retained energy file changed after preparation")
    rows = []
    active = False
    for line in (directory / "qfep.out").read_text().splitlines():
        stripped = line.strip()
        if stripped.startswith("# Part 6: BAR"):
            if active or rows:
                raise ValueError("Multiple BAR sections in QFEP output")
            active = True
        elif stripped.startswith("# Part"):
            active = False
        elif active and stripped and not stripped.startswith("#"):
            fields = stripped.split()
            if len(fields) != 3:
                raise ValueError(f"Malformed BAR row: {stripped}")
            numbers = [float(_number(value)) for value in fields]
            if not all(math.isfinite(v) for v in numbers):
                raise ValueError("Nonfinite BAR result")
            rows.append(numbers)
    if len(rows) != len(manifest["retained"]):
        raise ValueError("Incomplete BAR output or wrong number of retained windows")
    for row, window in zip(rows, manifest["retained"]):
        if abs(row[0]-float(window["lambdas"][0])) > 5.1e-7:
            raise ValueError("QFEP energy-file mapping disagrees with original MD input")
    if abs(rows[0][1])+abs(rows[0][2]) > 1e-6:
        raise ValueError("BAR cumulative origin must be zero")
    for previous, current in zip(rows, rows[1:]):
        if abs(current[2]-previous[2]-current[1]) > .0016:
            raise ValueError("Inconsistent BAR increment/cumulative columns")
    return {"analysis": manifest["analysis"], "full_endpoint_free_energy": False,
            "start_lambdas": manifest["start_lambdas"], "end_lambdas": manifest["end_lambdas"],
            "kT_kcal_mol": manifest["kT_kcal_mol"], "bar_kcal_mol": rows[-1][2],
            "retained_pairs": len(rows)-1, "reported_precision_kcal_mol": .001,
            "state_constant_gap_multiplier": manifest["state_constant_gap_multiplier"],
            "qfep_output_sha256": _hash(directory / "qfep.out"),
            "qualification": manifest["qualification"]}


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    prep = commands.add_parser("prepare", help="Create a new analysis directory with interior energy-file links")
    prep.add_argument("qfep_input", type=Path)
    prep.add_argument("--md-input-dir", type=Path, required=True)
    prep.add_argument("--output-dir", type=Path, required=True)
    prep.add_argument("--energy-dir", type=Path,
                      help="Resolve relative energy filenames here (default: QFEP input directory)")
    prep.add_argument("--lambda-min")
    prep.add_argument("--lambda-max")
    read = commands.add_parser("summarize", help="Print a labelled result after rerunning QFEP in that directory")
    read.add_argument("directory", type=Path)
    args = parser.parse_args(argv)
    try:
        if args.command == "prepare":
            result = prepare(args.qfep_input, args.md_input_dir, args.output_dir,
                             args.lambda_min, args.lambda_max, args.energy_dir)
        else:
            result = summarize(args.directory)
    except (ValueError, OSError) as error:
        parser.exit(2, f"{error}\n")
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
