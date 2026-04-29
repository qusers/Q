#!/usr/bin/env python3

import argparse
import csv
import json
import os
import shlex
import shutil
import subprocess
import sys
import time
from contextlib import contextmanager
from datetime import datetime
from pathlib import Path
from statistics import median

os.environ.setdefault("MPLCONFIGDIR", "/tmp/qgpu-benchmark-matplotlib")

import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
TIME_STEP_NS = 2e-6

sys.path.insert(0, str(ROOT / "test"))
sys.path.insert(0, str(ROOT / "src" / "qligfep-newbin-unfinished"))

import runTEST  # noqa: E402
import qdyn as qdyn_prepare  # noqa: E402


@contextmanager
def pushd(path):
    previous = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


def abs_path(path):
    if path is None:
        return None
    return str(Path(path).expanduser().resolve())


def command_text(args):
    return " ".join(shlex.quote(str(arg)) for arg in args)


def split_mpirun_args(args):
    if args is None:
        return []
    if isinstance(args, str):
        return shlex.split(args)
    return [str(arg) for arg in args]


def build_fortran_command(fortran_bin, input_file, mpi_procs=None, mpirun_bin="mpirun", mpirun_args=None):
    command = [str(fortran_bin), input_file]
    if mpi_procs is None:
        return command
    return [str(mpirun_bin), "-np", str(mpi_procs), *split_mpirun_args(mpirun_args), *command]


def resolve_qgpu_bin(path):
    if path:
        candidate = Path(path).expanduser()
        if not candidate.is_absolute():
            candidate = (Path.cwd() / candidate).resolve()
        if not candidate.exists():
            raise FileNotFoundError(f"QGPU binary not found: {candidate}")
        return candidate

    for candidate in (ROOT / "bin" / "qdyn", ROOT / "src" / "core" / "qdyn"):
        if candidate.exists():
            return candidate
    raise FileNotFoundError(
        "QGPU binary not found. Expected bin/qdyn or src/core/qdyn, "
        "or pass --qgpu-bin."
    )


def resolve_fortran_bin(path):
    candidate = Path(path).expanduser()
    if not candidate.is_absolute():
        candidate = (Path.cwd() / candidate).resolve()
    if not candidate.exists():
        raise FileNotFoundError(f"Fortran binary not found: {candidate}")
    return candidate


def resolve_test_data(test_name, steps, lambda_name, shake):
    testinfo = runTEST.get_default_testinfo()
    if test_name not in testinfo:
        available = ", ".join(sorted(testinfo))
        raise ValueError(f"Unknown test '{test_name}'. Available tests: {available}")

    topdir = ROOT / "test" / "data" / "topology"
    inputdir = ROOT / "test" / "data" / "inputs"
    info = testinfo[test_name]
    topfile = info[0]
    if len(info) >= 3 and lambda_name is not None:
        stem, suffix = topfile.rsplit(".", 1)
        topfile = f"{stem}_{lambda_name}.{suffix}"

    data = {
        "avg": False,
        "curtest": None,
        "fep_path": None,
        "inputdir": str(inputdir),
        "lambda": lambda_name,
        "plot": False,
        "restraints_path": None,
        "shake": shake,
        "test": test_name,
        "testinfo": testinfo,
        "timestep": str(steps),
        "topdir": str(topdir),
        "topfile": topfile,
        "topology_path": str(topdir / topfile),
        "verbose": False,
    }
    if len(info) >= 3:
        data["fep_path"] = str(inputdir / info[2])
    if len(info) >= 4:
        data["restraints_path"] = str(inputdir / info[3])

    required = [Path(data["topology_path"])]
    if data["fep_path"] is not None:
        required.append(Path(data["fep_path"]))
    if data["restraints_path"] is not None:
        required.append(Path(data["restraints_path"]))
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("Required input file(s) not found: " + ", ".join(missing))

    return data


def run_timed(args, cwd, stdout_path, stderr_path):
    start = time.perf_counter()
    with open(stdout_path, "w", encoding="utf-8") as stdout_f, open(
        stderr_path, "w", encoding="utf-8"
    ) as stderr_f:
        completed = subprocess.run(args, cwd=cwd, stdout=stdout_f, stderr=stderr_f)
    wall_seconds = time.perf_counter() - start
    return completed.returncode, wall_seconds


def ns_per_day(steps, wall_seconds):
    if wall_seconds <= 0:
        return None
    return steps * TIME_STEP_NS * 86400 / wall_seconds


def write_md_input(data, fortran_dir):
    data["curtest"] = str(fortran_dir)
    with pushd(fortran_dir):
        runTEST.create_MD_input(data)


def run_fortran_repeats(
    data,
    fortran_bin,
    fortran_dir,
    repeat,
    steps,
    mpi_procs=None,
    mpirun_bin="mpirun",
    mpirun_args=None,
):
    records = []
    saw_success = False

    for index in range(1, repeat + 1):
        stdout_name = "fortran.log" if repeat == 1 else f"fortran_{index}.log"
        stderr_name = "fortran.err" if repeat == 1 else f"fortran_{index}.err"
        stdout_path = fortran_dir / stdout_name
        stderr_path = fortran_dir / stderr_name
        args = build_fortran_command(
            fortran_bin,
            "eq1.inp",
            mpi_procs=mpi_procs,
            mpirun_bin=mpirun_bin,
            mpirun_args=mpirun_args,
        )
        return_code, wall_seconds = run_timed(args, fortran_dir, stdout_path, stderr_path)
        if return_code == 0:
            saw_success = True
        records.append(
            {
                "test": data["test"],
                "runner": "fortran",
                "repeat": index,
                "command": command_text(args),
                "return_code": return_code,
                "wall_seconds": wall_seconds,
                "steps": steps,
                "ns_per_day": ns_per_day(steps, wall_seconds),
                "stdout": str(stdout_path),
                "stderr": str(stderr_path),
            }
        )

    return records, saw_success


def prepare_restart_with_qdyn_test(data, prep_fortran_bin, fortran_dir, prep_steps=None):
    input_path = fortran_dir / "eq1.inp"
    original_input = input_path.read_text(encoding="utf-8")
    parse_data = data
    if prep_steps is not None:
        prep_data = dict(data)
        prep_data["timestep"] = str(prep_steps)
        write_md_input(prep_data, fortran_dir)
        parse_data = prep_data

    stdout_path = fortran_dir / "restart_prep_qdyn_test.log"
    stderr_path = fortran_dir / "restart_prep_qdyn_test.err"
    args = [str(prep_fortran_bin), "eq1.inp"]
    try:
        return_code, _ = run_timed(args, fortran_dir, stdout_path, stderr_path)
        if return_code != 0:
            raise RuntimeError(
                "QGPU restart preparation failed. "
                f"Command: {command_text(args)} Logs: stdout={stdout_path} stderr={stderr_path}"
            )

        shutil.copyfile(stdout_path, fortran_dir / "eq1.log")
        with pushd(fortran_dir):
            runTEST.Parse_Q6_data(parse_data)
    finally:
        if prep_steps is not None:
            input_path.write_text(original_input, encoding="utf-8")


def prepare_qgpu_input(data, fortran_dir, prep_dir):
    prep_dir.mkdir(parents=True, exist_ok=True)
    restart_dir = prep_dir / "restart"
    restart_dir.mkdir(exist_ok=True)
    shutil.copyfile(fortran_dir / "coords.csv", restart_dir / "coords.csv")
    shutil.copyfile(fortran_dir / "velocities.csv", restart_dir / "velocities.csv")

    top_stem = Path(data["topfile"]).stem
    wd_rel = f"TEST/{top_stem}"
    with pushd(prep_dir):
        qdyn_prepare.Create_Environment(top=data["topology_path"], wd=wd_rel)
        qdyn_prepare.Prepare_Topology(top=data["topology_path"], wd=wd_rel)
        qdyn_prepare.Prepare_MD(top=data["topology_path"], md=str(fortran_dir / "eq1.inp"), wd=wd_rel)
        qdyn_prepare.Prepare_FEP(
            fepfile=data["fep_path"],
            wd=wd_rel,
            top=data["topology_path"],
        )
        qdyn_prepare.Read_Restart(restart=str(restart_dir), wd=wd_rel, top=data["topology_path"])

    prepared_data_dir = prep_dir / wd_rel
    if not (prepared_data_dir / "md.csv").exists():
        raise RuntimeError(f"Prepared QGPU data is missing md.csv: {prepared_data_dir}")
    return prepared_data_dir


def run_qgpu_repeats(data, qgpu_bin, prepared_data_dir, qgpu_runs_dir, repeat, steps):
    qgpu_runs_dir.mkdir(parents=True, exist_ok=True)
    records = []

    for index in range(1, repeat + 1):
        run_dir = qgpu_runs_dir / f"repeat_{index:03d}"
        data_dir = run_dir / prepared_data_dir.name
        if run_dir.exists():
            shutil.rmtree(run_dir)
        run_dir.mkdir(parents=True)
        shutil.copytree(prepared_data_dir, data_dir)

        stdout_path = run_dir / "qgpu.log"
        stderr_path = run_dir / "qgpu.err"
        args = [str(qgpu_bin), "--gpu", str(data_dir)]
        return_code, wall_seconds = run_timed(args, ROOT, stdout_path, stderr_path)
        records.append(
            {
                "test": data["test"],
                "runner": "qgpu",
                "repeat": index,
                "command": command_text(args),
                "return_code": return_code,
                "wall_seconds": wall_seconds,
                "steps": steps,
                "ns_per_day": ns_per_day(steps, wall_seconds),
                "stdout": str(stdout_path),
                "stderr": str(stderr_path),
            }
        )

    return records


def write_summary_csv(records, out_dir):
    csv_path = out_dir / "summary.csv"
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
    with open(csv_path, "w", newline="", encoding="utf-8") as csv_f:
        writer = csv.DictWriter(csv_f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)
    return csv_path


def read_summary_csv(csv_path):
    records = []
    with open(csv_path, newline="", encoding="utf-8") as csv_f:
        reader = csv.DictReader(csv_f)
        for row in reader:
            parsed = dict(row)
            parsed["repeat"] = int(parsed["repeat"])
            parsed["return_code"] = int(parsed["return_code"])
            parsed["wall_seconds"] = float(parsed["wall_seconds"])
            parsed["steps"] = int(parsed["steps"])
            parsed["ns_per_day"] = float(parsed["ns_per_day"]) if parsed.get("ns_per_day") else None
            records.append(parsed)
    if not records:
        raise RuntimeError(f"No records found in {csv_path}")
    return records


def summarize(records, args, qgpu_bin, fortran_bin, prep_fortran_bin):
    by_test = {}
    for record in records:
        by_test.setdefault(record["test"], {}).setdefault(record["runner"], []).append(record)

    tests = []
    for test_name in sorted(by_test):
        fortran_records = by_test[test_name].get("fortran", [])
        qgpu_records = by_test[test_name].get("qgpu", [])
        fortran_ok = [float(r["wall_seconds"]) for r in fortran_records if int(r["return_code"]) == 0]
        qgpu_ok = [float(r["wall_seconds"]) for r in qgpu_records if int(r["return_code"]) == 0]
        if not fortran_ok or not qgpu_ok:
            continue
        fortran_median = median(fortran_ok)
        qgpu_median = median(qgpu_ok)
        speedup = fortran_median / qgpu_median if qgpu_median > 0 else None
        tests.append(
            {
                "test": test_name,
                "fortran_median_seconds": fortran_median,
                "qgpu_median_seconds": qgpu_median,
                "speedup_x": speedup,
                "improvement_pct": (speedup - 1) * 100 if speedup is not None else None,
                "fortran_repeats": len(fortran_records),
                "qgpu_repeats": len(qgpu_records),
            }
        )

    return {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "args": {
            "test": args.test,
            "steps": args.steps,
            "lambda": args.lambda_name,
            "shake": args.shake,
            "repeat": args.repeat,
            "restart_prep_steps": getattr(args, "restart_prep_steps", None),
            "fortran_mpi_procs": getattr(args, "fortran_mpi_procs", None),
            "mpirun_bin": getattr(args, "mpirun_bin", None),
            "mpirun_args": getattr(args, "mpirun_args", None),
        },
        "binaries": {
            "fortran": str(fortran_bin),
            "restart_prep_fortran": str(prep_fortran_bin),
            "qgpu": str(qgpu_bin),
        },
        "tests": tests,
    }


def summarize_for_plot(records):
    args = argparse.Namespace(
        test=sorted({record["test"] for record in records}),
        steps=sorted({int(record["steps"]) for record in records}),
        lambda_name=None,
        shake=None,
        repeat=None,
    )
    return summarize(records, args, qgpu_bin="<from summary.csv>", fortran_bin="<from summary.csv>", prep_fortran_bin="")


def write_summary_json(summary, out_dir):
    json_path = out_dir / "summary.json"
    with open(json_path, "w", encoding="utf-8") as json_f:
        json.dump(summary, json_f, indent=2)
    return json_path


def plot_speedup(summary, out_dir):
    tests = summary["tests"]
    if not tests:
        return None

    fig_width = max(8.0, 2.0 + len(tests) * 1.2)
    fig, (ax, panel) = plt.subplots(
        1,
        2,
        figsize=(fig_width, 3.0),
        gridspec_kw={"width_ratios": [3.6, 1.8]},
    )

    x_positions = list(range(len(tests)))
    width = 0.34
    fortran_times = [item["fortran_median_seconds"] for item in tests]
    qgpu_times = [item["qgpu_median_seconds"] for item in tests]
    labels = [item["test"] for item in tests]

    ax.bar([x - width / 2 for x in x_positions], fortran_times, width, label="Fortran", color="#9b9b9b")
    ax.bar([x + width / 2 for x in x_positions], qgpu_times, width, label="QGPU", color="#0b71c8")
    ax.set_title("Execution Time (s)", fontsize=11, weight="bold")
    ax.set_ylabel("Time (s)")
    ax.set_xticks(x_positions)
    ax.set_xticklabels(labels, rotation=0 if len(labels) <= 3 else 30, ha="center")
    ax.legend(frameon=False, loc="upper right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#e7e7e7", linewidth=0.8)
    ax.set_axisbelow(True)

    for x, value in zip([x - width / 2 for x in x_positions], fortran_times):
        ax.text(x, value, f"{value:.1f}", ha="center", va="bottom", fontsize=8, weight="bold")
    for x, value in zip([x + width / 2 for x in x_positions], qgpu_times):
        ax.text(x, value, f"{value:.1f}", ha="center", va="bottom", fontsize=8, weight="bold")

    if len(tests) == 1:
        ymax = max(fortran_times[0], qgpu_times[0])
        ax.annotate(
            "",
            xy=(x_positions[0] + width / 2, qgpu_times[0] + ymax * 0.15),
            xytext=(x_positions[0] - width / 2, fortran_times[0] * 0.85),
            arrowprops={"arrowstyle": "->", "linestyle": "--", "color": "#0b71c8", "lw": 1.2},
        )

    best = max(tests, key=lambda item: item["speedup_x"] or 0)
    panel.set_facecolor("#eef5fd")
    for spine in panel.spines.values():
        spine.set_color("#8ab9ef")
    panel.set_xticks([])
    panel.set_yticks([])
    panel.text(0.5, 0.82, "Up to", ha="center", va="center", fontsize=13, weight="bold", color="#0b3970")
    panel.text(
        0.5,
        0.52,
        f"{best['speedup_x']:.1f}x",
        ha="center",
        va="center",
        fontsize=32,
        weight="bold",
        color="#003c7f",
    )
    panel.text(0.5, 0.28, "speedup", ha="center", va="center", fontsize=14, weight="bold", color="#0b3970")
    panel.text(0.5, 0.12, "(vs. Fortran)", ha="center", va="center", fontsize=10, color="#0b3970")

    fig.tight_layout()
    png_path = out_dir / "speedup.png"
    fig.savefig(png_path, dpi=200)
    plt.close(fig)
    return png_path


def plot_summary_csv(args):
    csv_path = Path(args.csv).expanduser().resolve()
    records = read_summary_csv(csv_path)
    summary = summarize_for_plot(records)
    if args.out:
        out_path = Path(args.out).expanduser().resolve()
        out_dir = out_path.parent
        out_dir.mkdir(parents=True, exist_ok=True)
        png_path = plot_speedup(summary, out_dir)
        if png_path != out_path:
            if out_path.exists():
                out_path.unlink()
            png_path.rename(out_path)
            png_path = out_path
    else:
        png_path = plot_speedup(summary, csv_path.parent)
    if png_path is None:
        raise RuntimeError("No successful Fortran/QGPU pairs found to plot.")
    print(f"Speedup plot: {png_path}")
    return 0


def default_out_dir(test_names):
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    label = test_names[0] if len(test_names) == 1 else "multi"
    return ROOT / "benchmark-qgpu" / "results" / f"{stamp}_{label}"


def parse_args():
    if len(sys.argv) > 1 and sys.argv[1] == "plot":
        parser = argparse.ArgumentParser(description="Plot benchmark_test.py speedup from an existing summary.csv.")
        parser.add_argument("command", choices=["plot"])
        parser.add_argument("csv", help="summary.csv written by benchmark_test.py.")
        parser.add_argument("--out", help="Output PNG path. Defaults to speedup.png next to the CSV.")
        return parser.parse_args()

    parser = argparse.ArgumentParser(description="Benchmark Fortran vs QGPU for runTEST.py test cases.")
    parser.add_argument("--test", nargs="+", help="Test name(s) from test/runTEST.py.")
    parser.add_argument("--list-tests", action="store_true", help="List available tests and exit.")
    parser.add_argument("--steps", type=int, help="MD steps to write into eq1.inp.")
    parser.add_argument("--lambda", dest="lambda_name", default=None, help="Perturbation lambda suffix, e.g. eq5.")
    parser.add_argument("--shake", action="store_true", help="Enable shake in generated MD input.")
    parser.add_argument("--repeat", type=int, default=1, help="Number of repeats for each runner.")
    parser.add_argument("--out", default=None, help="Output directory.")
    parser.add_argument(
        "--restart-prep-steps",
        type=int,
        default=1,
        help="MD steps used only for qdyn_test restart preparation. Defaults to 1.",
    )
    parser.add_argument(
        "--fortran-bin",
        default=str(ROOT / "src" / "q6" / "bin" / "q6" / "qdyn"),
        help="Path to production Fortran qdyn/qdynp binary used for timed Fortran runs.",
    )
    parser.add_argument(
        "--fortran-mpi-procs",
        type=int,
        default=None,
        help="Run the timed Fortran binary through mpirun with this many MPI ranks.",
    )
    parser.add_argument(
        "--mpirun-bin",
        default="mpirun",
        help="MPI launcher to use with --fortran-mpi-procs. Defaults to mpirun.",
    )
    parser.add_argument(
        "--mpirun-args",
        default=None,
        help='Extra MPI launcher arguments, quoted as one string, e.g. "--bind-to core".',
    )
    parser.add_argument(
        "--prep-fortran-bin",
        default=str(ROOT / "src" / "q6" / "bin" / "q6" / "qdyn_test"),
        help="Path to qdyn_test binary used only to prepare QGPU restart CSVs.",
    )
    parser.add_argument("--qgpu-bin", default=None, help="Path to QGPU qdyn binary.")
    return parser.parse_args()


def validate_args(args):
    if getattr(args, "command", None) == "plot":
        return
    if args.list_tests:
        return
    if not args.test:
        raise SystemExit("--test is required unless --list-tests is used.")
    if args.steps is None:
        raise SystemExit("--steps is required unless --list-tests is used.")
    if args.steps < 1:
        raise SystemExit("--steps must be >= 1.")
    if args.repeat < 1:
        raise SystemExit("--repeat must be >= 1.")
    if args.restart_prep_steps < 1:
        raise SystemExit("--restart-prep-steps must be >= 1.")
    if args.fortran_mpi_procs is not None and args.fortran_mpi_procs < 1:
        raise SystemExit("--fortran-mpi-procs must be >= 1.")


def main():
    args = parse_args()
    validate_args(args)

    if getattr(args, "command", None) == "plot":
        return plot_summary_csv(args)

    testinfo = runTEST.get_default_testinfo()
    if args.list_tests:
        for test_name in sorted(testinfo):
            print(test_name)
        return 0

    qgpu_bin = resolve_qgpu_bin(args.qgpu_bin)
    fortran_bin = resolve_fortran_bin(args.fortran_bin)
    prep_fortran_bin = resolve_fortran_bin(args.prep_fortran_bin)
    out_dir = Path(args.out).expanduser().resolve() if args.out else default_out_dir(args.test)
    out_dir.mkdir(parents=True, exist_ok=True)

    all_records = []
    try:
        for test_name in args.test:
            test_dir = out_dir / test_name
            fortran_dir = test_dir / "fortran"
            prep_dir = test_dir / "qgpu_prepare"
            qgpu_runs_dir = test_dir / "qgpu_runs"
            fortran_dir.mkdir(parents=True, exist_ok=True)

            data = resolve_test_data(test_name, args.steps, args.lambda_name, args.shake)
            print(f"Preparing Fortran input for {test_name} in {fortran_dir}")
            write_md_input(data, fortran_dir)

            if args.fortran_mpi_procs is None:
                print(f"Running Fortran for {test_name} ({args.repeat} repeat(s))")
            else:
                print(
                    f"Running Fortran for {test_name} with {args.fortran_mpi_procs} MPI rank(s) "
                    f"({args.repeat} repeat(s))"
                )
            fortran_records, fortran_ok = run_fortran_repeats(
                data,
                fortran_bin,
                fortran_dir,
                args.repeat,
                args.steps,
                mpi_procs=args.fortran_mpi_procs,
                mpirun_bin=args.mpirun_bin,
                mpirun_args=args.mpirun_args,
            )
            all_records.extend(fortran_records)
            if not fortran_ok:
                continue

            print(f"Preparing QGPU restart with qdyn_test for {test_name} ({args.restart_prep_steps} step(s))")
            prepare_restart_with_qdyn_test(
                data,
                prep_fortran_bin,
                fortran_dir,
                prep_steps=args.restart_prep_steps,
            )

            print(f"Preparing QGPU CSV input for {test_name}")
            prepared_data_dir = prepare_qgpu_input(data, fortran_dir, prep_dir)

            print(f"Running QGPU for {test_name} ({args.repeat} repeat(s))")
            all_records.extend(
                run_qgpu_repeats(data, qgpu_bin, prepared_data_dir, qgpu_runs_dir, args.repeat, args.steps)
            )
    finally:
        write_summary_csv(all_records, out_dir)

    failures = [record for record in all_records if record["return_code"] != 0]
    if failures:
        first = failures[0]
        raise RuntimeError(
            f"{first['runner']} failed for {first['test']} repeat {first['repeat']}. "
            f"Logs: stdout={first['stdout']} stderr={first['stderr']}"
        )

    summary = summarize(all_records, args, qgpu_bin, fortran_bin, prep_fortran_bin)
    csv_path = write_summary_csv(all_records, out_dir)
    json_path = write_summary_json(summary, out_dir)
    png_path = plot_speedup(summary, out_dir)

    print(f"Summary CSV: {csv_path}")
    print(f"Summary JSON: {json_path}")
    if png_path is not None:
        print(f"Speedup plot: {png_path}")
    for item in summary["tests"]:
        print(
            f"{item['test']}: Fortran {item['fortran_median_seconds']:.3f}s, "
            f"QGPU {item['qgpu_median_seconds']:.3f}s, speedup {item['speedup_x']:.2f}x"
        )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, RuntimeError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
