#!/usr/bin/env python3

import argparse
import csv
import json
import math
import os
import sys
from datetime import datetime
from pathlib import Path
from statistics import median

os.environ.setdefault("MPLCONFIGDIR", "/tmp/qgpu-benchmark-matplotlib")

import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt

from benchmark_test import (
    ROOT,
    ns_per_day,
    prepare_qgpu_input,
    prepare_restart_with_qdyn_test,
    resolve_fortran_bin,
    resolve_qgpu_bin,
    resolve_test_data,
    run_fortran_repeats,
    run_qgpu_repeats,
    write_md_input,
)
from benchmark_nsday import run_concurrency_batch


RESTART_INIT_STEPS = 1


def default_collect_out():
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return ROOT / "benchmark-qgpu" / "results" / f"{stamp}_system_scaling"


def count_atoms(prepared_data_dir):
    coords_path = Path(prepared_data_dir) / "coords.csv"
    if not coords_path.exists():
        raise FileNotFoundError(f"coords.csv not found: {coords_path}")
    with open(coords_path, encoding="utf-8") as coords_f:
        return int(coords_f.readline().strip())


def successful_times(records):
    return [float(record["wall_seconds"]) for record in records if int(record["return_code"]) == 0]


def parse_optional_float(value):
    if value in (None, ""):
        return float("nan")
    return float(value)


def write_raw_records(records, out_dir):
    path = out_dir / "system_scaling_raw.csv"
    fieldnames = [
        "test",
        "runner",
        "repeat",
        "command",
        "return_code",
        "wall_seconds",
        "steps",
        "ns_per_day",
        "stdout",
        "stderr",
    ]
    with open(path, "w", newline="", encoding="utf-8") as csv_f:
        writer = csv.DictWriter(csv_f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)
    return path


def write_qgpu_concurrency_records(records, out_dir):
    path = out_dir / "system_scaling_qgpu_concurrency.csv"
    fieldnames = [
        "test",
        "label",
        "concurrency",
        "repeat",
        "steps",
        "batch_wall_seconds",
        "total_ns_per_day",
        "mean_process_ns_per_day",
        "failed_processes",
        "command",
    ]
    with open(path, "w", newline="", encoding="utf-8") as csv_f:
        writer = csv.DictWriter(csv_f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)
    return path


def write_summary(rows, out_dir, metadata):
    summary_csv = out_dir / "system_scaling.csv"
    meta_json = out_dir / "system_scaling_meta.json"

    fieldnames = [
        "test",
        "atoms",
        "steps",
        "fortran_wall_median_s",
        "qgpu_wall_median_s",
        "fortran_ns_per_day",
        "qgpu_ns_per_day",
        "qgpu_best_concurrency",
        "speedup_x",
        "fortran_repeats",
        "qgpu_repeats",
    ]
    with open(summary_csv, "w", newline="", encoding="utf-8") as csv_f:
        writer = csv.DictWriter(csv_f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    with open(meta_json, "w", encoding="utf-8") as json_f:
        json.dump(metadata, json_f, indent=2)

    return summary_csv, meta_json


def run_qgpu_concurrency_sweep(args, test_name, qgpu_bin, prepared_data_dir, qgpu_runs_dir):
    batch_rows = []
    process_rows = []
    for concurrency in args.concurrency:
        for repeat in range(1, args.repeat + 1):
            run_dir = qgpu_runs_dir / f"c{concurrency:03d}" / f"repeat_{repeat:03d}"
            print(f"Running QGPU for {test_name}: concurrency={concurrency}, repeat={repeat}")
            batch_row, rows = run_concurrency_batch(
                qgpu_bin=qgpu_bin,
                prepared_data_dir=prepared_data_dir,
                run_dir=run_dir,
                concurrency=concurrency,
                steps=args.steps,
                label=test_name,
                repeat=repeat,
            )
            batch_row["test"] = test_name
            batch_rows.append(batch_row)
            process_rows.extend(rows)
            if batch_row["failed_processes"]:
                raise RuntimeError(
                    f"{batch_row['failed_processes']} QGPU process(es) failed for {test_name} "
                    f"at concurrency {concurrency}, repeat {repeat}. Logs are under {run_dir}"
                )
    return batch_rows, process_rows


def collect_one_test(args, test_name, out_dir, fortran_bin, prep_fortran_bin, qgpu_bin):
    test_dir = out_dir / test_name
    fortran_dir = test_dir / "fortran"
    prep_dir = test_dir / "qgpu_prepare"
    qgpu_runs_dir = test_dir / "qgpu_runs"
    fortran_dir.mkdir(parents=True, exist_ok=True)

    data = resolve_test_data(test_name, args.steps, args.lambda_name, args.shake)
    fortran_records = []
    fortran_times = []

    if args.gpu_only:
        init_data = resolve_test_data(test_name, RESTART_INIT_STEPS, args.lambda_name, args.shake)
        print(f"Preparing QGPU restart for {test_name} with {RESTART_INIT_STEPS} MD step(s)")
        write_md_input(init_data, fortran_dir)
        prepare_restart_with_qdyn_test(init_data, prep_fortran_bin, fortran_dir)

        print(f"Writing QGPU benchmark input for {test_name} with {args.steps} MD step(s)")
        write_md_input(data, fortran_dir)
    else:
        print(f"Preparing {test_name}")
        write_md_input(data, fortran_dir)

        print(f"Running Fortran qdyn for {test_name} ({args.repeat} repeat(s))")
        fortran_records, fortran_ok = run_fortran_repeats(data, fortran_bin, fortran_dir, args.repeat, args.steps)
        if not fortran_ok:
            return None, fortran_records, []
        fortran_times = successful_times(fortran_records)

        print(f"Preparing QGPU input for {test_name}")
        prepare_restart_with_qdyn_test(data, prep_fortran_bin, fortran_dir)

    prepared_data_dir = prepare_qgpu_input(data, fortran_dir, prep_dir)
    atoms = count_atoms(prepared_data_dir)

    qgpu_concurrency_rows = []
    if args.concurrency:
        qgpu_records = []
        qgpu_concurrency_rows, _ = run_qgpu_concurrency_sweep(
            args, test_name, qgpu_bin, prepared_data_dir, qgpu_runs_dir
        )
    else:
        print(f"Running QGPU for {test_name} ({args.repeat} repeat(s))")
        qgpu_records = run_qgpu_repeats(data, qgpu_bin, prepared_data_dir, qgpu_runs_dir, args.repeat, args.steps)

    qgpu_times = successful_times(qgpu_records)
    if args.concurrency:
        successful_batches = [row for row in qgpu_concurrency_rows if int(row["failed_processes"]) == 0]
        if not successful_batches:
            return None, [*fortran_records, *qgpu_records], qgpu_concurrency_rows
        best_qgpu = max(successful_batches, key=lambda row: float(row["total_ns_per_day"]))
        qgpu_wall = float(best_qgpu["batch_wall_seconds"])
        qgpu_ns_day = float(best_qgpu["total_ns_per_day"])
        qgpu_best_concurrency = int(best_qgpu["concurrency"])
        qgpu_repeat_count = len(qgpu_concurrency_rows)
    else:
        if not qgpu_times:
            return None, [*fortran_records, *qgpu_records], qgpu_concurrency_rows
        qgpu_wall = median(qgpu_times)
        qgpu_ns_day = ns_per_day(args.steps, qgpu_wall)
        qgpu_best_concurrency = 1
        qgpu_repeat_count = len(qgpu_records)

    if not args.gpu_only and not fortran_times:
        return None, [*fortran_records, *qgpu_records], qgpu_concurrency_rows

    fortran_wall = median(fortran_times) if fortran_times else None
    fortran_ns_day = ns_per_day(args.steps, fortran_wall) if fortran_wall is not None else None
    row = {
        "test": test_name,
        "atoms": atoms,
        "steps": args.steps,
        "fortran_wall_median_s": fortran_wall if fortran_wall is not None else "",
        "qgpu_wall_median_s": qgpu_wall,
        "fortran_ns_per_day": fortran_ns_day if fortran_ns_day is not None else "",
        "qgpu_ns_per_day": qgpu_ns_day,
        "qgpu_best_concurrency": qgpu_best_concurrency,
        "speedup_x": qgpu_ns_day / fortran_ns_day if fortran_ns_day is not None and fortran_ns_day > 0 else "",
        "fortran_repeats": len(fortran_records),
        "qgpu_repeats": qgpu_repeat_count,
    }
    return row, [*fortran_records, *qgpu_records], qgpu_concurrency_rows


def collect(args):
    out_dir = Path(args.out).expanduser().resolve() if args.out else default_collect_out()
    out_dir.mkdir(parents=True, exist_ok=True)
    fortran_bin = None if args.gpu_only else resolve_fortran_bin(args.fortran_bin)
    prep_fortran_bin = resolve_fortran_bin(args.prep_fortran_bin)
    qgpu_bin = resolve_qgpu_bin(args.qgpu_bin)

    rows = []
    raw_records = []
    qgpu_concurrency_records = []
    try:
        for test_name in args.test:
            row, records, concurrency_records = collect_one_test(
                args, test_name, out_dir, fortran_bin, prep_fortran_bin, qgpu_bin
            )
            raw_records.extend(records)
            qgpu_concurrency_records.extend(concurrency_records)
            write_raw_records(raw_records, out_dir)
            if args.concurrency:
                write_qgpu_concurrency_records(qgpu_concurrency_records, out_dir)
            if row is not None:
                rows.append(row)
                write_summary(
                    rows,
                    out_dir,
                    {
                        "created_at": datetime.now().isoformat(timespec="seconds"),
                        "tests": args.test,
                        "steps": args.steps,
                        "repeat": args.repeat,
                        "gpu_only": args.gpu_only,
                        "concurrency": args.concurrency,
                        "fortran_bin": str(fortran_bin) if fortran_bin is not None else None,
                        "prep_fortran_bin": str(prep_fortran_bin),
                        "qgpu_bin": str(qgpu_bin),
                    },
                )
    finally:
        raw_path = write_raw_records(raw_records, out_dir)
        concurrency_path = (
            write_qgpu_concurrency_records(qgpu_concurrency_records, out_dir) if args.concurrency else None
        )

    failures = [record for record in raw_records if int(record["return_code"]) != 0]
    if failures:
        first = failures[0]
        raise RuntimeError(
            f"{first['runner']} failed for {first['test']} repeat {first['repeat']}. "
            f"Logs: stdout={first['stdout']} stderr={first['stderr']}; raw CSV: {raw_path}"
        )

    summary_csv, meta_json = write_summary(
        rows,
        out_dir,
        {
            "created_at": datetime.now().isoformat(timespec="seconds"),
            "tests": args.test,
            "steps": args.steps,
            "repeat": args.repeat,
            "gpu_only": args.gpu_only,
            "concurrency": args.concurrency,
            "fortran_bin": str(fortran_bin) if fortran_bin is not None else None,
            "prep_fortran_bin": str(prep_fortran_bin),
            "qgpu_bin": str(qgpu_bin),
        },
    )
    print(f"Summary CSV: {summary_csv}")
    print(f"Raw CSV: {raw_path}")
    if concurrency_path is not None:
        print(f"QGPU concurrency CSV: {concurrency_path}")
    print(f"Metadata JSON: {meta_json}")
    return 0


def load_rows(csv_path):
    rows = []
    with open(csv_path, newline="", encoding="utf-8") as csv_f:
        reader = csv.DictReader(csv_f)
        for row in reader:
            parsed = dict(row)
            for key in [
                "atoms",
                "steps",
                "fortran_wall_median_s",
                "qgpu_wall_median_s",
                "fortran_ns_per_day",
                "qgpu_ns_per_day",
                "qgpu_best_concurrency",
                "speedup_x",
            ]:
                parsed[key] = parse_optional_float(parsed.get(key))
            rows.append(parsed)
    if not rows:
        raise RuntimeError(f"No rows found in {csv_path}")
    return rows


def fmt_atoms(atoms):
    atoms = int(atoms)
    if atoms >= 1000:
        return f"{atoms / 1000:.1f}k atoms"
    return f"{atoms} atoms"


def annotate_bars(ax, bars, formatter):
    for bar in bars:
        height = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            height,
            formatter(height),
            ha="center",
            va="bottom",
            fontsize=8,
            weight="bold",
        )


def plot_speedup(rows, out_path, title):
    if not any(math.isfinite(row["speedup_x"]) for row in rows):
        raise RuntimeError("speedup plot requires Fortran data. Use --metric nsday for --gpu-only results.")
    labels = [row["test"] for row in rows]
    speedups = [row["speedup_x"] for row in rows]
    atoms = [row["atoms"] for row in rows]

    fig, (ax, panel) = plt.subplots(
        1,
        2,
        figsize=(9.2, 3.3),
        gridspec_kw={"width_ratios": [4.3, 1.55]},
    )
    x = range(len(rows))
    bars = ax.bar(x, speedups, color="#0b71c8", width=0.62)
    annotate_bars(ax, bars, lambda value: f"{value:.1f}x")
    ax.set_title(title, loc="left", fontsize=13, weight="bold", color="#113b5f")
    ax.set_ylabel("Speedup vs Fortran (x)")
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels)
    ax.grid(axis="y", color="#e5e8ee", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for xpos, atom_count in zip(x, atoms):
        ax.text(xpos, -0.08, fmt_atoms(atom_count), transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=8)

    best = max(rows, key=lambda row: row["speedup_x"])
    panel.set_facecolor("#eef5fd")
    for spine in panel.spines.values():
        spine.set_color("#8ab9ef")
    panel.set_xticks([])
    panel.set_yticks([])
    panel.text(0.5, 0.80, "Best", ha="center", va="center", fontsize=12, weight="bold", color="#0b3970")
    panel.text(0.5, 0.55, f"{best['speedup_x']:.1f}x", ha="center", va="center", fontsize=30, weight="bold", color="#003c7f")
    panel.text(0.5, 0.35, "speedup", ha="center", va="center", fontsize=13, weight="bold", color="#0b3970")
    panel.text(0.5, 0.18, best["test"], ha="center", va="center", fontsize=10, color="#0b3970")
    panel.text(0.5, 0.08, fmt_atoms(best["atoms"]), ha="center", va="center", fontsize=9, color="#0b3970")

    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def plot_nsday(rows, out_path, title):
    labels = [row["test"] for row in rows]
    x = list(range(len(rows)))
    width = 0.34

    fig, ax = plt.subplots(figsize=(8.6, 3.5))
    fortran = [row["fortran_ns_per_day"] for row in rows]
    qgpu = [row["qgpu_ns_per_day"] for row in rows]
    has_fortran = any(math.isfinite(value) for value in fortran)
    has_concurrency = any(math.isfinite(row["qgpu_best_concurrency"]) and row["qgpu_best_concurrency"] > 1 for row in rows)
    qgpu_label = "QGPU best total" if has_concurrency else "QGPU"
    if has_fortran:
        bars_cpu = ax.bar([i - width / 2 for i in x], fortran, width, label="Fortran CPU", color="#9b9b9b")
        bars_gpu = ax.bar([i + width / 2 for i in x], qgpu, width, label=qgpu_label, color="#0b71c8")
        annotate_bars(ax, bars_cpu, lambda value: f"{value:.1f}")
    else:
        bars_gpu = ax.bar(x, qgpu, width * 1.55, label=qgpu_label, color="#0b71c8")
    annotate_bars(ax, bars_gpu, lambda value: f"{value:.1f}")
    ax.set_title(title, loc="left", fontsize=13, weight="bold", color="#113b5f")
    ax.set_ylabel("Best total ns/day" if has_concurrency else "ns/day")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.grid(axis="y", color="#e5e8ee", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, loc="best")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for xpos, row in zip(x, rows):
        ax.text(xpos, -0.08, fmt_atoms(row["atoms"]), transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=8)

    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def plot_atoms(rows, out_path, title):
    fig, ax = plt.subplots(figsize=(6.5, 3.8))
    xs = [row["atoms"] for row in rows]
    has_speedup = any(math.isfinite(row["speedup_x"]) for row in rows)
    value_key = "speedup_x" if has_speedup else "qgpu_ns_per_day"
    ys = [row[value_key] for row in rows]
    ax.plot(xs, ys, color="#0b71c8", marker="o", linewidth=1.8)
    for row in rows:
        suffix = f"{row['speedup_x']:.1f}x" if has_speedup else f"{row['qgpu_ns_per_day']:.1f} ns/day"
        ax.text(row["atoms"], row[value_key], f" {row['test']} ({suffix})", va="center", fontsize=8)
    ax.set_title(title, loc="left", fontsize=13, weight="bold", color="#113b5f")
    ax.set_xlabel("Atoms")
    ax.set_ylabel("Speedup vs Fortran (x)" if has_speedup else "QGPU ns/day")
    ax.grid(True, color="#e5e8ee", linewidth=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def plot(args):
    rows = load_rows(Path(args.csv).expanduser().resolve())
    rows.sort(key=lambda row: row["atoms"] if args.sort == "atoms" else row["test"])
    out_path = Path(args.out).expanduser().resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if args.metric == "speedup":
        plot_speedup(rows, out_path, args.title)
    elif args.metric == "nsday":
        plot_nsday(rows, out_path, args.title)
    elif args.metric == "atoms":
        plot_atoms(rows, out_path, args.title)
    else:
        raise SystemExit(f"Unknown metric: {args.metric}")

    print(f"Plot written to: {out_path}")
    return 0


def positive_int(value):
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be >= 1")
    return parsed


def parse_args():
    parser = argparse.ArgumentParser(description="Collect and plot QGPU scaling across molecular systems.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    collect_parser = subparsers.add_parser("collect", help="Run Fortran/QGPU benchmark for multiple tests.")
    collect_parser.add_argument("--test", nargs="+", required=True, help="runTEST.py test names.")
    collect_parser.add_argument("--steps", type=positive_int, required=True, help="MD steps.")
    collect_parser.add_argument("--lambda", dest="lambda_name", default=None, help="Perturbation lambda suffix, e.g. eq5.")
    collect_parser.add_argument("--shake", action="store_true", help="Enable shake.")
    collect_parser.add_argument("--repeat", type=positive_int, default=1, help="Repeats per runner per system.")
    collect_parser.add_argument(
        "--concurrency",
        type=positive_int,
        nargs="+",
        help="Concurrent QGPU instance counts to sweep; summary uses the maximum total ns/day.",
    )
    collect_parser.add_argument(
        "--gpu-only",
        action="store_true",
        help="Skip timed Fortran qdyn runs and collect only QGPU performance.",
    )
    collect_parser.add_argument("--out", help="Output directory.")
    collect_parser.add_argument(
        "--fortran-bin",
        default=str(ROOT / "src" / "q6" / "bin" / "q6" / "qdyn"),
        help="Path to production Fortran qdyn binary.",
    )
    collect_parser.add_argument(
        "--prep-fortran-bin",
        default=str(ROOT / "src" / "q6" / "bin" / "q6" / "qdyn_test"),
        help="Path to qdyn_test used only to prepare QGPU restart CSVs.",
    )
    collect_parser.add_argument("--qgpu-bin", help="Path to QGPU qdyn binary.")

    plot_parser = subparsers.add_parser("plot", help="Plot system scaling from system_scaling.csv.")
    plot_parser.add_argument("csv", help="system_scaling.csv from collect.")
    plot_parser.add_argument("--out", required=True, help="Output PNG path.")
    plot_parser.add_argument(
        "--metric",
        choices=["speedup", "nsday", "atoms"],
        default="speedup",
        help="Plot style.",
    )
    plot_parser.add_argument("--sort", choices=["atoms", "test"], default="atoms", help="System order.")
    plot_parser.add_argument("--title", default="Performance Across Molecular Systems", help="Plot title.")
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
