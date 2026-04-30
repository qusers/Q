#!/usr/bin/env python3

import argparse
import csv
import io
import json
import math
import os
import shutil
import sys
from contextlib import redirect_stdout
from datetime import datetime
from pathlib import Path
from statistics import mean

os.environ.setdefault("MPLCONFIGDIR", "/tmp/qgpu-benchmark-matplotlib")

import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt

from benchmark_test import (
    ROOT,
    command_text,
    prepare_qgpu_input,
    prepare_restart_with_qdyn_test,
    resolve_fortran_bin,
    resolve_qgpu_bin,
    resolve_test_data,
    run_timed,
    write_md_input,
)

sys.path.insert(0, str(ROOT / "src" / "Qgpu"))

import compare  # noqa: E402
import energy as ENERGY  # noqa: E402


def default_collect_out(test_name):
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return ROOT / "benchmark-qgpu" / "results" / f"{stamp}_{test_name}_correctness"


def run_qgpu_once(qgpu_bin, prepared_data_dir, run_dir):
    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True)
    data_dir = run_dir / prepared_data_dir.name
    shutil.copytree(prepared_data_dir, data_dir)

    stdout_path = run_dir / "qgpu.log"
    stderr_path = run_dir / "qgpu.err"
    args = [str(qgpu_bin), "--gpu", str(data_dir)]
    return_code, wall_seconds = run_timed(args, ROOT, stdout_path, stderr_path)
    if return_code != 0:
        raise RuntimeError(
            "QGPU correctness run failed. "
            f"Command: {command_text(args)} Logs: stdout={stdout_path} stderr={stderr_path}"
        )
    return data_dir, {
        "command": command_text(args),
        "return_code": return_code,
        "wall_seconds": wall_seconds,
        "stdout": str(stdout_path),
        "stderr": str(stderr_path),
    }


def load_qgpu_energy(qgpu_data_dir):
    energy_path = Path(qgpu_data_dir) / "output" / "energies.csv"
    if not energy_path.exists():
        raise FileNotFoundError(f"QGPU energy file not found: {energy_path}")
    return ENERGY.Read_Energy(str(energy_path), 0).QDYN(), energy_path


def load_fortran_energy(fortran_dir):
    q_data_path = Path(fortran_dir) / "Q_data.json"
    if not q_data_path.exists():
        raise FileNotFoundError(f"Fortran energy JSON not found: {q_data_path}")
    with open(q_data_path, encoding="utf-8") as json_f:
        return json.load(json_f), q_data_path


def find_prepared_qgpu_dir(reference_dir):
    prepare_root = Path(reference_dir) / "qgpu_prepare" / "TEST"
    if not prepare_root.is_dir():
        raise FileNotFoundError(f"Prepared QGPU TEST directory not found: {prepare_root}")
    candidates = sorted(path for path in prepare_root.iterdir() if path.is_dir() and (path / "md.csv").exists())
    if len(candidates) != 1:
        shown = ", ".join(str(path) for path in candidates)
        raise RuntimeError(f"Expected exactly one prepared QGPU directory under {prepare_root}; found: {shown}")
    return candidates[0]


def copy_reference_inputs(reference_dir, out_dir):
    reference_dir = Path(reference_dir).expanduser().resolve()
    source_fortran_dir = reference_dir / "fortran_reference"
    if not (source_fortran_dir / "Q_data.json").exists():
        raise FileNotFoundError(f"Fortran reference Q_data.json not found: {source_fortran_dir / 'Q_data.json'}")

    source_prepared_dir = find_prepared_qgpu_dir(reference_dir)
    fortran_dir = out_dir / "fortran_reference"
    prep_dir = out_dir / "qgpu_prepare"
    prepared_data_dir = prep_dir / "TEST" / source_prepared_dir.name

    if fortran_dir.exists():
        shutil.rmtree(fortran_dir)
    if prep_dir.exists():
        shutil.rmtree(prep_dir)

    shutil.copytree(source_fortran_dir, fortran_dir)
    shutil.copytree(source_prepared_dir, prepared_data_dir)
    return fortran_dir, prep_dir, prepared_data_dir, reference_dir


def build_correctness_rows(fortran_data, qgpu_data, tolerance):
    compare.ENERGY_TOLERANCE = tolerance
    rows = []
    frames = sorted(int(key) for key in fortran_data.keys() if key.isdigit())
    for frame in frames:
        if frame >= len(qgpu_data):
            continue
        with redirect_stdout(io.StringIO()):
            passed, fortran_values, qgpu_values = compare.compare_energies(
                fortran_data[str(frame)],
                qgpu_data[frame],
            )
        for term, fortran_value, qgpu_value in zip(compare.header, fortran_values, qgpu_values):
            if math.isnan(fortran_value) or math.isnan(qgpu_value):
                continue
            abs_error = abs(fortran_value - qgpu_value)
            rel_error = abs_error / abs(fortran_value) if fortran_value != 0 else ""
            rows.append(
                {
                    "frame": frame,
                    "term": term,
                    "fortran": fortran_value,
                    "qgpu": qgpu_value,
                    "abs_error": abs_error,
                    "rel_error": rel_error,
                    "passed_tolerance": abs_error <= tolerance,
                    "frame_passed": passed,
                }
            )
    if not rows:
        raise RuntimeError("No comparable energy rows were produced.")
    return rows


def summarize_rows(rows, tolerance):
    abs_errors = [float(row["abs_error"]) for row in rows]
    by_term = {}
    for row in rows:
        by_term.setdefault(row["term"], []).append(float(row["abs_error"]))

    term_summary = []
    for term, values in sorted(by_term.items()):
        term_summary.append(
            {
                "term": term,
                "max_abs_error": max(values),
                "mean_abs_error": mean(values),
                "rmse": math.sqrt(mean([value * value for value in values])),
            }
        )

    return {
        "tolerance": tolerance,
        "frames": sorted({int(row["frame"]) for row in rows}),
        "terms": len(by_term),
        "rows": len(rows),
        "max_abs_error": max(abs_errors),
        "mean_abs_error": mean(abs_errors),
        "rmse": math.sqrt(mean([value * value for value in abs_errors])),
        "passed": all(float(row["abs_error"]) <= tolerance for row in rows),
        "term_summary": term_summary,
    }


def write_outputs(rows, summary, out_dir, metadata):
    terms_csv = out_dir / "correctness_terms.csv"
    summary_json = out_dir / "correctness_summary.json"

    with open(terms_csv, "w", newline="", encoding="utf-8") as csv_f:
        fieldnames = [
            "frame",
            "term",
            "fortran",
            "qgpu",
            "abs_error",
            "rel_error",
            "passed_tolerance",
            "frame_passed",
        ]
        writer = csv.DictWriter(csv_f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    payload = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "metadata": metadata,
        "summary": summary,
    }
    with open(summary_json, "w", encoding="utf-8") as json_f:
        json.dump(payload, json_f, indent=2)

    return terms_csv, summary_json


def collect(args):
    out_dir = Path(args.out).expanduser().resolve() if args.out else default_collect_out(args.test)
    out_dir.mkdir(parents=True, exist_ok=True)

    qgpu_bin = resolve_qgpu_bin(args.qgpu_bin)
    qgpu_run_dir = out_dir / "qgpu_run"

    if args.reference_dir:
        print(f"Reusing Fortran/QGPU prepared reference from {args.reference_dir}")
        fortran_dir, prep_dir, prepared_data_dir, reference_dir = copy_reference_inputs(args.reference_dir, out_dir)
        prep_fortran_bin = None
    else:
        default_prep_fortran_bin = (
            ROOT / "src" / "q6" / "bin" / "q6" / "qdynp"
            if args.prep_fortran_mpi_procs is not None
            else ROOT / "src" / "q6" / "bin" / "q6" / "qdyn_test"
        )
        prep_fortran_bin = resolve_fortran_bin(args.prep_fortran_bin or default_prep_fortran_bin)
        data = resolve_test_data(args.test, args.steps, args.lambda_name, args.shake)

        fortran_dir = out_dir / "fortran_reference"
        prep_dir = out_dir / "qgpu_prepare"
        fortran_dir.mkdir(parents=True, exist_ok=True)

        if args.prep_fortran_mpi_procs is None:
            print(f"Preparing Fortran reference for {args.test}")
        else:
            print(
                f"Preparing Fortran reference for {args.test} "
                f"with {args.prep_fortran_mpi_procs} MPI rank(s)"
            )
        write_md_input(data, fortran_dir)
        prepare_restart_with_qdyn_test(
            data,
            prep_fortran_bin,
            fortran_dir,
            mpi_procs=args.prep_fortran_mpi_procs,
            mpirun_bin=args.mpirun_bin,
            mpirun_args=args.mpirun_args,
        )

        print("Preparing QGPU input")
        prepared_data_dir = prepare_qgpu_input(data, fortran_dir, prep_dir)
        reference_dir = None

    print("Running QGPU correctness simulation")
    qgpu_data_dir, qgpu_run = run_qgpu_once(qgpu_bin, prepared_data_dir, qgpu_run_dir)

    fortran_data, fortran_energy_path = load_fortran_energy(fortran_dir)
    qgpu_data, qgpu_energy_path = load_qgpu_energy(qgpu_data_dir)
    rows = build_correctness_rows(fortran_data, qgpu_data, args.tolerance)
    summary = summarize_rows(rows, args.tolerance)

    terms_csv, summary_json = write_outputs(
        rows,
        summary,
        out_dir,
        {
            "test": args.test,
            "steps": args.steps,
            "lambda": args.lambda_name,
            "shake": args.shake,
            "qgpu_bin": str(qgpu_bin),
            "prep_fortran_bin": str(prep_fortran_bin) if prep_fortran_bin is not None else None,
            "prep_fortran_mpi_procs": args.prep_fortran_mpi_procs,
            "mpirun_bin": args.mpirun_bin,
            "mpirun_args": args.mpirun_args,
            "reference_dir": str(reference_dir) if reference_dir is not None else None,
            "prepared_qgpu_input": str(prepared_data_dir),
            "fortran_energy": str(fortran_energy_path),
            "qgpu_energy": str(qgpu_energy_path),
            "qgpu_run": qgpu_run,
        },
    )

    print(f"Terms CSV: {terms_csv}")
    print(f"Summary JSON: {summary_json}")
    print(
        f"max |delta E| = {summary['max_abs_error']:.6g} kcal/mol; "
        f"RMSE = {summary['rmse']:.6g}; passed = {summary['passed']}"
    )
    return 0


def load_rows(csv_path):
    rows = []
    with open(csv_path, newline="", encoding="utf-8") as csv_f:
        reader = csv.DictReader(csv_f)
        for row in reader:
            row["frame"] = int(row["frame"])
            row["fortran"] = float(row["fortran"])
            row["qgpu"] = float(row["qgpu"])
            row["abs_error"] = float(row["abs_error"])
            rows.append(row)
    if not rows:
        raise RuntimeError(f"No rows found in {csv_path}")
    return rows


def select_term_rows(rows, term):
    selected = [row for row in rows if row["term"] == term]
    if not selected:
        terms = ", ".join(sorted({row["term"] for row in rows}))
        raise ValueError(f"Term '{term}' not found. Available terms: {terms}")
    return sorted(selected, key=lambda row: row["frame"])


def plot(args):
    rows = load_rows(Path(args.csv).expanduser().resolve())
    selected = select_term_rows(rows, args.term)

    frames = [row["frame"] for row in selected]
    fortran_values = [row["fortran"] for row in selected]
    qgpu_values = [row["qgpu"] for row in selected]
    abs_errors = [row["abs_error"] for row in selected]
    rel_errors_pct = [
        (row["abs_error"] / abs(row["fortran"]) * 100.0) if row["fortran"] != 0 else 0.0
        for row in selected
    ]
    max_abs_error = max(abs_errors)
    mean_abs_error = mean(abs_errors)
    rmse = math.sqrt(mean([value * value for value in abs_errors]))
    max_rel_error = max(rel_errors_pct)
    mean_rel_error = mean(rel_errors_pct)

    if args.error_mode == "relative":
        plotted_errors = rel_errors_pct
        error_ylabel = "Relative error (%)"
        tolerance = args.tolerance
        tolerance_label = "rel. tolerance"
    else:
        plotted_errors = abs_errors
        error_ylabel = "|delta E|"
        tolerance = args.tolerance
        tolerance_label = "tolerance"

    out_path = Path(args.out).expanduser().resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(9.8, 4.2))
    grid = fig.add_gridspec(2, 2, width_ratios=[4.2, 1.45], height_ratios=[2.3, 1.3])
    ax_energy = fig.add_subplot(grid[0, 0])
    ax_error = fig.add_subplot(grid[1, 0], sharex=ax_energy)
    ax_panel = fig.add_subplot(grid[:, 1])

    ax_energy.plot(frames, fortran_values, color="#4a4a4a", linewidth=1.8, label="Fortran")
    ax_energy.plot(frames, qgpu_values, color="#0b71c8", linestyle="--", linewidth=1.6, label="QGPU")
    ax_energy.set_title(args.title, loc="left", fontsize=13, weight="bold", color="#113b5f")
    ax_energy.set_ylabel(f"{args.term} (kcal/mol)")
    ax_energy.grid(axis="y", color="#e5e8ee", linewidth=0.8)
    ax_energy.legend(frameon=False, loc="best", fontsize=8)
    ax_energy.spines["top"].set_visible(False)
    ax_energy.spines["right"].set_visible(False)

    ax_error.plot(frames, plotted_errors, color="#d62728", linewidth=1.6)
    ax_error.fill_between(frames, plotted_errors, color="#d62728", alpha=0.13)
    if tolerance is not None:
        ax_error.axhline(tolerance, color="#777777", linestyle=":", linewidth=1.0, label=tolerance_label)
        ax_error.legend(frameon=False, loc="best", fontsize=8)
    ax_error.set_xlabel("MD step")
    ax_error.set_ylabel(error_ylabel)
    ax_error.grid(axis="y", color="#e5e8ee", linewidth=0.8)
    ax_error.spines["top"].set_visible(False)
    ax_error.spines["right"].set_visible(False)

    ax_panel.set_facecolor("#eef5fd")
    for spine in ax_panel.spines.values():
        spine.set_color("#8ab9ef")
    ax_panel.set_xticks([])
    ax_panel.set_yticks([])
    if args.error_mode == "relative":
        ax_panel.text(0.5, 0.84, "Consistency", ha="center", va="center", fontsize=13, weight="bold", color="#0b3970")
        ax_panel.text(0.5, 0.64, f"{max_rel_error:.3f}%", ha="center", va="center", fontsize=24, weight="bold", color="#003c7f")
        ax_panel.text(0.5, 0.50, "max rel. error", ha="center", va="center", fontsize=10, color="#0b3970")
        ax_panel.axhline(0.36, xmin=0.15, xmax=0.85, color="#8ab9ef", linewidth=0.8)
        ax_panel.text(0.5, 0.25, f"mean {mean_rel_error:.3f}%", ha="center", va="center", fontsize=11, weight="bold", color="#0b3970")
        ax_panel.text(0.5, 0.13, f"abs RMSE {rmse:.2e}", ha="center", va="center", fontsize=9, color="#0b3970")
    else:
        ax_panel.text(0.5, 0.82, "Agreement", ha="center", va="center", fontsize=13, weight="bold", color="#0b3970")
        ax_panel.text(0.5, 0.62, f"{max_abs_error:.2e}", ha="center", va="center", fontsize=22, weight="bold", color="#003c7f")
        ax_panel.text(0.5, 0.48, "max |delta E|", ha="center", va="center", fontsize=10, color="#0b3970")
        ax_panel.axhline(0.34, xmin=0.15, xmax=0.85, color="#8ab9ef", linewidth=0.8)
        ax_panel.text(0.5, 0.23, f"RMSE {rmse:.2e}", ha="center", va="center", fontsize=11, weight="bold", color="#0b3970")
        ax_panel.text(0.5, 0.12, f"mean {mean_abs_error:.2e}", ha="center", va="center", fontsize=10, color="#0b3970")

    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)
    print(f"Plot written to: {out_path}")
    return 0


def positive_int(value):
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be >= 1")
    return parsed


def nonnegative_float(value):
    parsed = float(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("must be >= 0")
    return parsed


def parse_args():
    parser = argparse.ArgumentParser(description="Collect and plot Fortran vs QGPU energy correctness.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    collect_parser = subparsers.add_parser("collect", help="Run a correctness benchmark and write CSV data.")
    collect_parser.add_argument("--test", required=True, help="runTEST.py test name.")
    collect_parser.add_argument("--steps", type=positive_int, required=True, help="MD steps.")
    collect_parser.add_argument("--lambda", dest="lambda_name", default=None, help="Perturbation lambda suffix, e.g. eq5.")
    collect_parser.add_argument("--shake", action="store_true", help="Enable shake.")
    collect_parser.add_argument("--out", help="Output directory.")
    collect_parser.add_argument("--qgpu-bin", help="Path to QGPU qdyn binary.")
    collect_parser.add_argument(
        "--reference-dir",
        help="Existing correctness result directory containing fortran_reference/ and qgpu_prepare/ to reuse.",
    )
    collect_parser.add_argument(
        "--prep-fortran-bin",
        default=None,
        help="Path to Fortran binary used to generate reference data. Defaults to qdynp with MPI, otherwise qdyn_test.",
    )
    collect_parser.add_argument(
        "--prep-fortran-mpi-procs",
        type=positive_int,
        default=None,
        help="Run the Fortran reference preparation through mpirun with this many MPI ranks.",
    )
    collect_parser.add_argument(
        "--mpirun-bin",
        default="mpirun",
        help="MPI launcher to use with --prep-fortran-mpi-procs. Defaults to mpirun.",
    )
    collect_parser.add_argument(
        "--mpirun-args",
        default=None,
        help='Extra MPI launcher arguments, quoted as one string, e.g. "--bind-to core".',
    )
    collect_parser.add_argument(
        "--tolerance",
        type=nonnegative_float,
        default=1e-3,
        help="Absolute energy tolerance in kcal/mol for pass/fail summary.",
    )

    plot_parser = subparsers.add_parser("plot", help="Plot correctness from correctness_terms.csv.")
    plot_parser.add_argument("csv", help="correctness_terms.csv from collect.")
    plot_parser.add_argument("--out", required=True, help="Output PNG path.")
    plot_parser.add_argument("--term", default="total-Utot", help="Energy term to plot.")
    plot_parser.add_argument(
        "--title",
        default="Long-Run Energy Consistency",
        help="Plot title.",
    )
    plot_parser.add_argument(
        "--error-mode",
        choices=["absolute", "relative"],
        default="absolute",
        help="Plot absolute kcal/mol error or relative percent error.",
    )
    plot_parser.add_argument(
        "--tolerance",
        type=nonnegative_float,
        default=None,
        help="Optional horizontal tolerance line on the error panel. Units follow --error-mode.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    if args.command == "collect":
        return collect(args)
    if args.command == "plot":
        return plot(args)
    raise SystemExit(f"Unknown command: {args.command}")


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
