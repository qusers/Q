#!/usr/bin/env python3

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
from statistics import median

os.environ.setdefault("MPLCONFIGDIR", "/tmp/qgpu-benchmark-matplotlib")

import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt

from benchmark_test import (
    ROOT,
    command_text,
    resolve_qgpu_bin,
    resolve_test_data,
    write_md_input,
)
from benchmark_input import ns_per_day, read_inp_settings, stage_inp_input


def default_collect_out(label):
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    safe_label = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in label)
    return ROOT / "benchmark-qgpu" / "results" / f"{stamp}_{safe_label}_nsday"


def prepare_from_test(args, out_dir):
    benchmark_data = resolve_test_data(args.test, args.steps, args.lambda_name, args.shake)
    input_dir = out_dir / "prepare" / args.test
    input_dir.mkdir(parents=True, exist_ok=True)

    print(f"Writing native QGPU input for {args.test} with {args.steps} MD step(s)")
    write_md_input(benchmark_data, input_dir)
    return input_dir / "eq1.inp"


def resolve_collect_input(args, out_dir):
    if args.input:
        settings = read_inp_settings(args.input)
        if args.steps is not None and args.steps != settings["steps"]:
            raise ValueError(
                f"--steps {args.steps} does not match {settings['steps']} in {settings['path']}"
            )
        return settings

    if not args.test:
        raise SystemExit("collect requires --test or --input.")
    if args.steps is None:
        raise SystemExit("collect with --test requires --steps.")
    return read_inp_settings(prepare_from_test(args, out_dir))


def cleanup_successful_task_data(processes):
    removed = 0
    for item in processes:
        if item["return_code"] != 0:
            continue
        input_dir = item["input_dir"]
        if input_dir.exists():
            shutil.rmtree(input_dir)
            removed += 1
    return removed


def run_concurrency_batch(
    qgpu_bin,
    input_file,
    run_dir,
    concurrency,
    steps,
    stepsize_fs,
    label,
    repeat,
    cleanup_run_data=False,
):
    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True)

    command_template = None
    launch_specs = []
    input_file = read_inp_settings(input_file)["path"]
    for index in range(1, concurrency + 1):
        proc_dir = run_dir / f"proc_{index:03d}"
        proc_dir.mkdir(parents=True)
        staged_input = stage_inp_input(input_file, proc_dir)

        stdout_path = proc_dir / "qgpu.log"
        stderr_path = proc_dir / "qgpu.err"
        args = [str(qgpu_bin), "--gpu", str(staged_input)]
        command_template = command_text([str(qgpu_bin), "--gpu", "<input.inp>"])
        launch_specs.append(
            {
                "index": index,
                "args": args,
                "stdout": stdout_path,
                "stderr": stderr_path,
                "input_dir": staged_input.parent,
                "command": command_text(args),
            }
        )

    processes = []
    process_rows = []
    batch_start = time.perf_counter()
    for spec in launch_specs:
        stdout_f = open(spec["stdout"], "w", encoding="utf-8")
        stderr_f = open(spec["stderr"], "w", encoding="utf-8")
        proc_start = time.perf_counter()
        process = subprocess.Popen(spec["args"], cwd=input_file.parent, stdout=stdout_f, stderr=stderr_f)
        processes.append(
            {
                "index": spec["index"],
                "process": process,
                "stdout_file": stdout_f,
                "stderr_file": stderr_f,
                "stdout": spec["stdout"],
                "stderr": spec["stderr"],
                "input_dir": spec["input_dir"],
                "start": proc_start,
                "command": spec["command"],
            }
        )

    remaining = set(range(len(processes)))
    while remaining:
        for item_index in list(remaining):
            item = processes[item_index]
            return_code = item["process"].poll()
            if return_code is None:
                continue
            item["return_code"] = return_code
            item["end"] = time.perf_counter()
            item["stdout_file"].close()
            item["stderr_file"].close()
            remaining.remove(item_index)
        if remaining:
            time.sleep(0.01)

    for item in processes:
        wall_seconds = item["end"] - item["start"]
        process_rows.append(
            {
                "label": label,
                "concurrency": concurrency,
                "repeat": repeat,
                "process_index": item["index"],
                "return_code": item["return_code"],
                "process_wall_seconds": wall_seconds,
                "process_ns_per_day": ns_per_day(steps, stepsize_fs, wall_seconds) or "",
                "stdout": str(item["stdout"]),
                "stderr": str(item["stderr"]),
                "command": item["command"],
            }
        )

    batch_wall_seconds = time.perf_counter() - batch_start
    failed = sum(1 for row in process_rows if row["return_code"] != 0)

    if cleanup_run_data:
        removed_task_data = cleanup_successful_task_data(processes)
        if removed_task_data:
            print(f"Removed copied QGPU data for {removed_task_data} successful task(s) under {run_dir}")

    total_ns_per_day = concurrency * ns_per_day(steps, stepsize_fs, batch_wall_seconds)
    mean_process_ns_per_day = (
        sum(float(row["process_ns_per_day"]) for row in process_rows if row["process_ns_per_day"] != "")
        / len(process_rows)
    )
    return {
        "label": label,
        "concurrency": concurrency,
        "repeat": repeat,
        "steps": steps,
        "batch_wall_seconds": batch_wall_seconds,
        "total_ns_per_day": total_ns_per_day,
        "mean_process_ns_per_day": mean_process_ns_per_day,
        "failed_processes": failed,
        "command": command_template,
    }, process_rows


def write_collect_outputs(batch_rows, process_rows, out_dir, meta):
    summary_csv = out_dir / "nsday_summary.csv"
    process_csv = out_dir / "nsday_processes.csv"
    meta_json = out_dir / "nsday_meta.json"

    with open(summary_csv, "w", newline="", encoding="utf-8") as csv_f:
        fieldnames = [
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
        writer = csv.DictWriter(csv_f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(batch_rows)

    with open(process_csv, "w", newline="", encoding="utf-8") as csv_f:
        fieldnames = [
            "label",
            "concurrency",
            "repeat",
            "process_index",
            "return_code",
            "process_wall_seconds",
            "process_ns_per_day",
            "stdout",
            "stderr",
            "command",
        ]
        writer = csv.DictWriter(csv_f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(process_rows)

    with open(meta_json, "w", encoding="utf-8") as json_f:
        json.dump(meta, json_f, indent=2)

    return summary_csv, process_csv, meta_json


def collect(args):
    label = args.label or args.test or Path(args.input).stem
    out_dir = Path(args.out).expanduser().resolve() if args.out else default_collect_out(label)
    out_dir.mkdir(parents=True, exist_ok=True)
    qgpu_bin = resolve_qgpu_bin(args.qgpu_bin)
    settings = resolve_collect_input(args, out_dir)
    input_file = settings["path"]
    steps = settings["steps"]
    stepsize_fs = settings["stepsize_fs"]

    batch_rows = []
    process_rows = []
    for concurrency in args.concurrency:
        for repeat in range(1, args.repeat + 1):
            run_dir = out_dir / "runs" / f"c{concurrency:03d}" / f"repeat_{repeat:03d}"
            print(f"Running {label}: concurrency={concurrency}, repeat={repeat}")
            batch_row, rows = run_concurrency_batch(
                qgpu_bin=qgpu_bin,
                input_file=input_file,
                run_dir=run_dir,
                concurrency=concurrency,
                steps=steps,
                stepsize_fs=stepsize_fs,
                label=label,
                repeat=repeat,
                cleanup_run_data=not args.keep_run_data,
            )
            batch_rows.append(batch_row)
            process_rows.extend(rows)
            if batch_row["failed_processes"]:
                summary_csv, process_csv, meta_json = write_collect_outputs(
                    batch_rows,
                    process_rows,
                    out_dir,
                    {
                        "created_at": datetime.now().isoformat(timespec="seconds"),
                        "label": label,
                        "qgpu_bin": str(qgpu_bin),
                        "input": str(input_file),
                        "steps": steps,
                        "stepsize_fs": stepsize_fs,
                        "keep_run_data": args.keep_run_data,
                    },
                )
                raise RuntimeError(
                    f"{batch_row['failed_processes']} process(es) failed at concurrency "
                    f"{concurrency}, repeat {repeat}. Summary: {summary_csv}; processes: {process_csv}; meta: {meta_json}"
                )
            if args.pause_seconds > 0:
                time.sleep(args.pause_seconds)

    summary_csv, process_csv, meta_json = write_collect_outputs(
        batch_rows,
        process_rows,
        out_dir,
        {
            "created_at": datetime.now().isoformat(timespec="seconds"),
            "label": label,
            "test": args.test,
            "input": str(input_file),
            "qgpu_bin": str(qgpu_bin),
            "steps": steps,
            "stepsize_fs": stepsize_fs,
            "concurrency": args.concurrency,
            "repeat": args.repeat,
            "keep_run_data": args.keep_run_data,
        },
    )
    print(f"Summary CSV: {summary_csv}")
    print(f"Process CSV: {process_csv}")
    print(f"Metadata JSON: {meta_json}")
    return 0


def load_plot_series(csv_paths, metric):
    series = {}
    for csv_path in csv_paths:
        with open(csv_path, newline="", encoding="utf-8") as csv_f:
            reader = csv.DictReader(csv_f)
            for row in reader:
                if int(row.get("failed_processes") or 0) != 0:
                    continue
                label = row["label"]
                concurrency = int(row["concurrency"])
                value = float(row[metric])
                series.setdefault(label, {}).setdefault(concurrency, []).append(value)

    plotted = []
    for label, by_concurrency in sorted(series.items()):
        xs = sorted(by_concurrency)
        ys = [median(by_concurrency[x]) for x in xs]
        plotted.append({"label": label, "xs": xs, "ys": ys})
    if not plotted:
        raise RuntimeError("No successful rows found in input CSV file(s).")
    return plotted


def plot(args):
    metric = args.metric
    metric_labels = {
        "total_ns_per_day": "Total Throughput (ns/day)",
        "mean_process_ns_per_day": "Mean Per-Process Throughput (ns/day)",
    }
    series = load_plot_series([Path(path).expanduser().resolve() for path in args.csv], metric)
    out_path = Path(args.out).expanduser().resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fig, (ax, panel) = plt.subplots(
        1,
        2,
        figsize=(9.5, 3.2),
        gridspec_kw={"width_ratios": [4.4, 1.55]},
    )
    palette = ["#1f77b4", "#43a047", "#f57c00", "#7b1fa2", "#00838f"]
    all_points = []
    all_xs = sorted({x for item in series for x in item["xs"]})
    for index, item in enumerate(series):
        color = palette[index % len(palette)]
        ax.plot(item["xs"], item["ys"], marker="o", linewidth=1.8, markersize=4.5, color=color, label=item["label"])
        for x, y in zip(item["xs"], item["ys"]):
            all_points.append((y, item["label"], x))
            if len(item["xs"]) == 1:
                x_offset = 0
                ha = "center"
            elif x == item["xs"][0]:
                x_offset = 6
                ha = "left"
            elif x == item["xs"][-1]:
                x_offset = -6
                ha = "right"
            else:
                x_offset = 0
                ha = "center"
            ax.annotate(
                f"{y:.1f}",
                xy=(x, y),
                xytext=(x_offset, 6),
                textcoords="offset points",
                ha=ha,
                va="bottom",
                fontsize=8,
                weight="bold",
                color="#253142",
            )

    y_values = [point[0] for point in all_points]
    y_min = min(y_values)
    y_max = max(y_values)
    y_span = y_max - y_min
    y_pad = max(y_span * 0.22, y_max * 0.035, 0.5)
    ax.set_ylim(max(0, y_min - y_pad * 0.35), y_max + y_pad)
    ax.set_xticks(all_xs)
    if len(all_xs) == 1:
        ax.set_xlim(all_xs[0] - 0.5, all_xs[0] + 0.5)
    else:
        x_pad = max((all_xs[-1] - all_xs[0]) * 0.06, 0.25)
        ax.set_xlim(all_xs[0] - x_pad, all_xs[-1] + x_pad)

    ax.text(0.0, 1.14, args.title, transform=ax.transAxes, fontsize=13, weight="bold", color="#0f5f18")
    ax.text(0.0, 1.07, args.subtitle, transform=ax.transAxes, fontsize=9, style="italic", color="#253142")
    ax.set_xlabel("Concurrent Simulations")
    ax.set_ylabel(metric_labels[metric])
    ax.grid(axis="y", color="#e3e7ed", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if len(series) > 1:
        ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.0, 1.01), ncols=2, fontsize=8)

    best_points = sorted(all_points, reverse=True)
    best = best_points[0]
    second = None
    seen_labels = {best[1]}
    for point in best_points[1:]:
        if point[1] not in seen_labels:
            second = point
            break

    panel.set_facecolor("#edf7eb")
    for spine in panel.spines.values():
        spine.set_color("#a3d39b")
    panel.set_xticks([])
    panel.set_yticks([])
    panel.text(0.5, 0.80, "Up to", ha="center", va="center", fontsize=11, weight="bold", color="#14751c")
    panel.text(0.5, 0.55, f"{best[0]:.1f}", ha="center", va="center", fontsize=30, weight="bold", color="#14751c")
    panel.text(0.5, 0.35, "ns/day", ha="center", va="center", fontsize=13, weight="bold", color="#14751c")
    panel.text(0.5, 0.20, f"{best[1]}", ha="center", va="center", fontsize=9, color="#253142")
    if second is not None:
        panel.axhline(0.12, xmin=0.12, xmax=0.88, color="#7fbf79", linewidth=0.8)
        panel.text(0.5, 0.05, f"{second[0]:.1f} ns/day", ha="center", va="bottom", fontsize=10, weight="bold", color="#14751c")

    fig.tight_layout(rect=(0, 0, 1, 0.9))
    fig.savefig(out_path, dpi=220)
    plt.close(fig)
    print(f"Plot written to: {out_path}")
    return 0


def positive_int(value):
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be >= 1")
    return parsed


def parse_args():
    parser = argparse.ArgumentParser(description="Collect and plot QGPU concurrency throughput in ns/day.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    collect_parser = subparsers.add_parser("collect", help="Run QGPU concurrency benchmark and write CSV data.")
    input_group = collect_parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--test", help="runTEST.py test name to prepare and benchmark.")
    input_group.add_argument("--input", help="Existing native Q .inp file.")
    collect_parser.add_argument("--steps", type=positive_int, help="MD steps. Required with --test; validated with --input.")
    collect_parser.add_argument("--lambda", dest="lambda_name", default=None, help="Perturbation lambda suffix, e.g. eq5.")
    collect_parser.add_argument("--shake", action="store_true", help="Enable shake when preparing from --test.")
    collect_parser.add_argument(
        "--concurrency",
        type=positive_int,
        nargs="+",
        default=[1, 2, 4, 8],
        help="Concurrent QGPU simulations to run.",
    )
    collect_parser.add_argument("--repeat", type=positive_int, default=1, help="Repeats per concurrency level.")
    collect_parser.add_argument("--label", help="Series label written into the CSV, e.g. 'A100 (thrombin)'.")
    collect_parser.add_argument("--out", help="Output directory.")
    collect_parser.add_argument("--qgpu-bin", help="Path to QGPU qdyn binary.")
    collect_parser.add_argument("--pause-seconds", type=float, default=0.0, help="Pause between batches.")
    collect_parser.add_argument(
        "--keep-run-data",
        action="store_true",
        help=(
            "Keep per-task copied QGPU input/output directories. By default successful task data is "
            "deleted after logs and timing are recorded."
        ),
    )

    plot_parser = subparsers.add_parser("plot", help="Plot ns/day vs concurrency from one or more CSV files.")
    plot_parser.add_argument("csv", nargs="+", help="One or more nsday_summary.csv files from collect.")
    plot_parser.add_argument("--out", required=True, help="Output PNG path.")
    plot_parser.add_argument(
        "--metric",
        choices=["total_ns_per_day", "mean_process_ns_per_day"],
        default="total_ns_per_day",
        help="Y-axis metric.",
    )
    plot_parser.add_argument("--title", default="Multi-System Concurrency (MPS)", help="Plot title.")
    plot_parser.add_argument(
        "--subtitle",
        default="Total simulation throughput at different concurrency levels",
        help="Plot subtitle.",
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
