#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable


@dataclass
class TotalEnergyFrame:
    interval: int
    temperature: float | None
    ukin: float
    upot: float
    utot: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build and run QDyn in default and QDYN_SPFP modes, collect logs, "
            "parse per-interval energies.csv output, and compare energy differences."
        )
    )
    parser.add_argument(
        "--test",
        required=True,
        help=(
            "Path to a generated QDyn case directory containing md.csv/topo.csv/etc, "
            "or an outer test directory containing TEST/<case>."
        ),
    )
    parser.add_argument(
        "--backend",
        choices=("gpu", "cpu"),
        default="gpu",
        help="Run the GPU or CPU backend. Default: gpu.",
    )
    parser.add_argument(
        "--output-dir",
        help=(
            "Directory for staged runs, logs, parsed outputs, and comparison report. "
            "Default: <src/core>/precision_compare_runs/<case>_<timestamp>."
        ),
    )
    parser.add_argument(
        "--build-root",
        help=(
            "Root directory for CMake build trees. "
            "Default: <src/core>/precision_compare_builds."
        ),
    )
    parser.add_argument(
        "--double-binary",
        help="Path to an already-built default-precision qdyn binary. Skips default build if provided.",
    )
    parser.add_argument(
        "--spfp-binary",
        help="Path to an already-built QDYN_SPFP qdyn binary. Skips SPFP build if provided.",
    )
    parser.add_argument(
        "--skip-build",
        action="store_true",
        help="Require --double-binary and --spfp-binary instead of building binaries.",
    )
    parser.add_argument(
        "--keep-work",
        action="store_true",
        help="Keep staged copied case directories even if a run fails.",
    )
    return parser.parse_args()


def run_command(cmd: list[str], cwd: Path, log_path: Path | None = None) -> None:
    if log_path is None:
        subprocess.run(cmd, cwd=cwd, check=True)
        return

    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as log_file:
        process = subprocess.run(
            cmd,
            cwd=cwd,
            stdout=log_file,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
    if process.returncode != 0:
        raise subprocess.CalledProcessError(process.returncode, cmd)


def find_case_dir(path: Path) -> Path:
    path = path.resolve()
    if (path / "md.csv").is_file() and (path / "topo.csv").is_file():
        return path

    test_root = path / "TEST"
    if test_root.is_dir():
        candidates = sorted(
            candidate
            for candidate in test_root.iterdir()
            if candidate.is_dir() and (candidate / "md.csv").is_file() and (candidate / "topo.csv").is_file()
        )
        if len(candidates) == 1:
            return candidates[0]
        if not candidates:
            raise FileNotFoundError(f"No generated test case with md.csv/topo.csv found under {test_root}")
        raise RuntimeError(
            f"Multiple generated test cases found under {test_root}. "
            f"Pass the exact case directory instead."
        )

    raise FileNotFoundError(
        f"Could not find a QDyn case directory from {path}. "
        "Expected md.csv/topo.csv directly or under TEST/<case>."
    )


def core_dir_from_script() -> Path:
    return Path(__file__).resolve().parents[1]


def default_output_dir(core_dir: Path, case_name: str) -> Path:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return core_dir / "precision_compare_runs" / f"{case_name}_{timestamp}"


def default_build_root(core_dir: Path) -> Path:
    return core_dir / "precision_compare_builds"


def ensure_clean_output_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def configure_and_build(core_dir: Path, build_dir: Path, spfp: bool, log_path: Path) -> Path:
    build_dir.mkdir(parents=True, exist_ok=True)
    configure_cmd = [
        "cmake",
        "-S",
        str(core_dir),
        "-B",
        str(build_dir),
        "-DCMAKE_BUILD_TYPE=Release",
    ]
    if spfp:
        configure_cmd.append("-DQDYN_SPFP=ON")

    build_cmd = ["cmake", "--build", str(build_dir), "-j", "4"]
    run_command(configure_cmd, cwd=core_dir, log_path=log_path)
    run_command(build_cmd, cwd=core_dir, log_path=log_path)

    binary = build_dir / "qdyn"
    if not binary.is_file():
        raise FileNotFoundError(f"Expected built binary at {binary}")
    return binary


def stage_case(src_case_dir: Path, dst_case_dir: Path) -> Path:
    if dst_case_dir.exists():
        shutil.rmtree(dst_case_dir)
    shutil.copytree(src_case_dir, dst_case_dir)

    output_dir = dst_case_dir / "output"
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    return dst_case_dir


def run_case(binary: Path, backend: str, case_dir: Path, log_path: Path) -> Path:
    cmd = [str(binary)]
    if backend == "gpu":
        cmd.append("--gpu")
    cmd.append(str(case_dir))
    run_command(cmd, cwd=case_dir, log_path=log_path)

    energies_path = case_dir / "output" / "energies.csv"
    if not energies_path.is_file():
        raise FileNotFoundError(f"Expected energies file at {energies_path}")
    return energies_path


def parse_energies_csv(path: Path) -> list[TotalEnergyFrame]:
    frames: list[TotalEnergyFrame] = []
    current_interval: int | None = None
    current_temperature: float | None = None
    current_ukin: float | None = None
    current_upot: float | None = None
    current_utot: float | None = None
    block: str | None = None

    def append_current_frame() -> None:
        if current_interval is None:
            return
        if current_ukin is None or current_upot is None or current_utot is None:
            raise ValueError(f"Incomplete [total] block for interval {current_interval} in {path}")
        frames.append(
            TotalEnergyFrame(
                interval=current_interval,
                temperature=current_temperature,
                ukin=current_ukin,
                upot=current_upot,
                utot=current_utot,
            )
        )

    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line == "lambdas":
                block = None
                continue
            if line.startswith("interval "):
                append_current_frame()
                current_interval = int(line.split()[1])
                current_temperature = None
                current_ukin = None
                current_upot = None
                current_utot = None
                block = None
                continue
            if line.startswith("[") and line.endswith("]"):
                block = line[1:-1].lower()
                continue

            if current_interval is None:
                continue

            if block not in {"temperature", "total"}:
                continue

            parts = line.split()
            if len(parts) < 2:
                continue
            key = parts[0]
            if block == "temperature" and key != "Temp":
                continue
            if block == "total" and key not in {"Ukin", "Upot", "Utot"}:
                continue

            try:
                value = float(parts[1])
            except ValueError:
                continue

            if block == "temperature" and key == "Temp":
                current_temperature = value
            elif block == "total":
                if key == "Ukin":
                    current_ukin = value
                elif key == "Upot":
                    current_upot = value
                elif key == "Utot":
                    current_utot = value

    append_current_frame()

    if not frames:
        raise ValueError(f"No intervals parsed from {path}")
    return frames


def write_parsed_energy_csv(path: Path, frames: Iterable[TotalEnergyFrame]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["interval", "temperature", "Ukin", "Upot", "Utot"])
        for frame in frames:
            writer.writerow([frame.interval, frame.temperature, frame.ukin, frame.upot, frame.utot])


def compare_frames(double_frames: list[TotalEnergyFrame], spfp_frames: list[TotalEnergyFrame]) -> dict:
    if len(double_frames) != len(spfp_frames):
        raise ValueError(
            f"Frame count mismatch: double={len(double_frames)} spfp={len(spfp_frames)}"
        )

    rows = []
    max_abs = {"Ukin": 0.0, "Upot": 0.0, "Utot": 0.0}
    sum_sq = {"Ukin": 0.0, "Upot": 0.0, "Utot": 0.0}

    for double_frame, spfp_frame in zip(double_frames, spfp_frames):
        if double_frame.interval != spfp_frame.interval:
            raise ValueError(
                f"Interval mismatch: double={double_frame.interval} spfp={spfp_frame.interval}"
            )

        diff_ukin = spfp_frame.ukin - double_frame.ukin
        diff_upot = spfp_frame.upot - double_frame.upot
        diff_utot = spfp_frame.utot - double_frame.utot

        row = {
            "interval": double_frame.interval,
            "double_temperature": double_frame.temperature,
            "spfp_temperature": spfp_frame.temperature,
            "double_Ukin": double_frame.ukin,
            "spfp_Ukin": spfp_frame.ukin,
            "diff_Ukin": diff_ukin,
            "abs_diff_Ukin": abs(diff_ukin),
            "double_Upot": double_frame.upot,
            "spfp_Upot": spfp_frame.upot,
            "diff_Upot": diff_upot,
            "abs_diff_Upot": abs(diff_upot),
            "double_Utot": double_frame.utot,
            "spfp_Utot": spfp_frame.utot,
            "diff_Utot": diff_utot,
            "abs_diff_Utot": abs(diff_utot),
        }
        rows.append(row)

        for key, value in (("Ukin", diff_ukin), ("Upot", diff_upot), ("Utot", diff_utot)):
            max_abs[key] = max(max_abs[key], abs(value))
            sum_sq[key] += value * value

    count = len(rows)
    rms = {key: math.sqrt(sum_sq[key] / count) for key in sum_sq}
    max_abs_frames: dict[str, list[dict]] = {}
    diff_field_map = {"Ukin": "diff_Ukin", "Upot": "diff_Upot", "Utot": "diff_Utot"}
    abs_diff_field_map = {"Ukin": "abs_diff_Ukin", "Upot": "abs_diff_Upot", "Utot": "abs_diff_Utot"}

    for metric in ("Ukin", "Upot", "Utot"):
        abs_field = abs_diff_field_map[metric]
        diff_field = diff_field_map[metric]
        metric_hits = []
        for row in rows:
            if math.isclose(row[abs_field], max_abs[metric], rel_tol=0.0, abs_tol=1e-15):
                metric_hits.append(
                    {
                        "interval": row["interval"],
                        "double": row[f"double_{metric}"],
                        "spfp": row[f"spfp_{metric}"],
                        "diff": row[diff_field],
                        "abs_diff": row[abs_field],
                    }
                )
        max_abs_frames[metric] = metric_hits

    return {
        "frame_count": count,
        "max_abs_diff": max_abs,
        "max_abs_diff_frames": max_abs_frames,
        "rms_diff": rms,
        "rows": rows,
    }


def write_comparison_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError("No comparison rows to write")
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    core_dir = core_dir_from_script()
    case_dir = find_case_dir(Path(args.test))
    case_name = case_dir.name

    output_dir = Path(args.output_dir).resolve() if args.output_dir else default_output_dir(core_dir, case_name)
    build_root = Path(args.build_root).resolve() if args.build_root else default_build_root(core_dir)
    ensure_clean_output_dir(output_dir)

    double_build_log = output_dir / "double_build.log"
    spfp_build_log = output_dir / "spfp_build.log"

    try:
        if args.skip_build:
            if not args.double_binary or not args.spfp_binary:
                raise ValueError("--skip-build requires --double-binary and --spfp-binary")

        if args.double_binary:
            double_binary = Path(args.double_binary).resolve()
        else:
            double_binary = configure_and_build(
                core_dir=core_dir,
                build_dir=build_root / "double",
                spfp=False,
                log_path=double_build_log,
            )

        if args.spfp_binary:
            spfp_binary = Path(args.spfp_binary).resolve()
        else:
            spfp_binary = configure_and_build(
                core_dir=core_dir,
                build_dir=build_root / "spfp",
                spfp=True,
                log_path=spfp_build_log,
            )

        staged_root = output_dir / "runs"
        double_case = stage_case(case_dir, staged_root / "double" / case_name)
        spfp_case = stage_case(case_dir, staged_root / "spfp" / case_name)

        double_log = output_dir / "double_run.log"
        spfp_log = output_dir / "spfp_run.log"

        double_energy_path = run_case(double_binary, args.backend, double_case, double_log)
        spfp_energy_path = run_case(spfp_binary, args.backend, spfp_case, spfp_log)

        double_frames = parse_energies_csv(double_energy_path)
        spfp_frames = parse_energies_csv(spfp_energy_path)

        write_parsed_energy_csv(output_dir / "double_total_energies.csv", double_frames)
        write_parsed_energy_csv(output_dir / "spfp_total_energies.csv", spfp_frames)

        comparison = compare_frames(double_frames, spfp_frames)
        write_comparison_csv(output_dir / "energy_diff.csv", comparison["rows"])

        summary = {
            "case_dir": str(case_dir),
            "backend": args.backend,
            "double_binary": str(double_binary),
            "spfp_binary": str(spfp_binary),
            "double_run_log": str(double_log),
            "spfp_run_log": str(spfp_log),
            "double_energy_file": str(double_energy_path),
            "spfp_energy_file": str(spfp_energy_path),
            "comparison_csv": str(output_dir / "energy_diff.csv"),
            "double_total_csv": str(output_dir / "double_total_energies.csv"),
            "spfp_total_csv": str(output_dir / "spfp_total_energies.csv"),
            "frame_count": comparison["frame_count"],
            "max_abs_diff": comparison["max_abs_diff"],
            "max_abs_diff_frames": comparison["max_abs_diff_frames"],
            "rms_diff": comparison["rms_diff"],
        }

        with (output_dir / "summary.json").open("w", encoding="utf-8") as handle:
            json.dump(summary, handle, indent=2)

        print(f"Case: {case_dir}")
        print(f"Backend: {args.backend}")
        print(f"Output: {output_dir}")
        print(f"Frames compared: {comparison['frame_count']}")
        print("Max abs diff:")
        print(
            f"  Ukin={comparison['max_abs_diff']['Ukin']:.10f} "
            f"Upot={comparison['max_abs_diff']['Upot']:.10f} "
            f"Utot={comparison['max_abs_diff']['Utot']:.10f}"
        )
        print("Max abs diff frames:")
        for metric in ("Ukin", "Upot", "Utot"):
            hits = comparison["max_abs_diff_frames"][metric]
            hit_text = ", ".join(
                f"interval {hit['interval']} (double={hit['double']:.6f}, spfp={hit['spfp']:.6f}, diff={hit['diff']:+.6f})"
                for hit in hits
            )
            print(f"  {metric}: {hit_text}")
        print("RMS diff:")
        print(
            f"  Ukin={comparison['rms_diff']['Ukin']:.10f} "
            f"Upot={comparison['rms_diff']['Upot']:.10f} "
            f"Utot={comparison['rms_diff']['Utot']:.10f}"
        )
        print(f"Per-step diff CSV: {output_dir / 'energy_diff.csv'}")
        return 0

    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        if not args.keep_work:
            staged_root = output_dir / "runs"
            if staged_root.exists():
                shutil.rmtree(staged_root)
        return 1


if __name__ == "__main__":
    sys.exit(main())
