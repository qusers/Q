import os
import time
import subprocess
import multiprocessing
import threading
import shutil
import json
import sys
from datetime import datetime
import glob
import csv
import psutil
from benchmark_report import make_html_report

TIME_STEP = 2 * 1e-6  # 2fs per step = 0.000002 ns per step

 
def _collect_gpu_mem_mb_for_pids(pids):
    """
    Returns the summed GPU memory (MB) used by the given PIDs, if NVIDIA tools exist.
    Otherwise returns None.
    """
    if shutil.which("nvidia-smi") is None or not pids:
        return None
    try:
        # Query all compute app processes once, map pid->MB
        out = subprocess.check_output(
            ["nvidia-smi", "--query-compute-apps=pid,used_memory",
             "--format=csv,noheader,nounits"],
            stderr=subprocess.DEVNULL,
            text=True
        )
        pid2mb = {}
        for line in out.strip().splitlines():
            parts = [x.strip() for x in line.split(",")]
            if len(parts) != 2:
                continue
            pid_str, mem_mb_str = parts
            try:
                pid2mb[int(pid_str)] = float(mem_mb_str)
            except ValueError:
                pass
        total = 0.0
        found_any = False
        for pid in pids:
            if pid in pid2mb:
                total += pid2mb[pid]
                found_any = True
        return total if found_any else 0.0
    except Exception:
        return None



def _sample_gpu_util_pct_nvsmi():
    """Fallback: query overall GPU util (%) via nvidia-smi; returns average across GPUs."""
    if shutil.which("nvidia-smi") is None:
        return None
    try:
        out = subprocess.check_output(
            ["nvidia-smi", "--query-gpu=utilization.gpu", "--format=csv,noheader,nounits"],
            stderr=subprocess.DEVNULL, text=True
        )
        vals = []
        for line in out.strip().splitlines():
            try:
                vals.append(float(line.strip()))
            except Exception:
                pass
        if not vals:
            return 0.0
        return sum(vals) / len(vals)
    except Exception:
        return None

def _sample_gpu_mem_util_pct_nvsmi():
    if shutil.which("nvidia-smi") is None:
        return None
    try:
        out = subprocess.check_output(
            ["nvidia-smi", "--query-gpu=memory.used,memory.total",
             "--format=csv,noheader,nounits"],
            stderr=subprocess.DEVNULL, text=True
        )
        vals = []
        for line in out.strip().splitlines():
            parts = [p.strip() for p in line.split(",")]
            if len(parts) != 2: 
                continue
            used, total = parts
            try:
                used = float(used); total = float(total)
                if total > 0:
                    vals.append(100.0 * used/total)
            except Exception:
                pass
        return sum(vals)/len(vals) if vals else 0.0
    except Exception:
        return None


def _descendant_pids(root_pid):
    """Return the set of descendant PIDs (including root) using psutil if available; else just root."""
    if psutil is None:
        return {root_pid}
    try:
        root = psutil.Process(root_pid)
        kids = root.children(recursive=True)
        return {root_pid, *(c.pid for c in kids)}
    except Exception:
        return {root_pid}

def _monitor_peaks(root_pid, stop_event, sample_interval=0.2):
    """
    Periodically sample CPU RAM and GPU VRAM for all descendants of root_pid.
    Returns a dict-like object to be filled with peaks.
    """
    peaks = {"peak_rss_bytes": 0, "peak_gpu_mem_mb": 0.0, 
             "gpu_util_peak_pct": 0.0, "gpu_util_sum_pct": 0.0, "gpu_util_samples": 0, 
             "gpu_mem_util_peak_pct": 0.0, "gpu_mem_util_sum_pct": 0.0, "gpu_mem_util_samples": 0}
    while not stop_event.is_set():
        pids = _descendant_pids(root_pid)
        # --- CPU RAM (RSS) ---
        if psutil is not None:
            rss_sum = 0
            for pid in list(pids):
                try:
                    rss_sum += psutil.Process(pid).memory_info().rss
                except Exception:
                    pass
            if rss_sum > peaks["peak_rss_bytes"]:
                peaks["peak_rss_bytes"] = rss_sum

        # --- GPU VRAM (MB) ---
        gpu_mb = _collect_gpu_mem_mb_for_pids(pids)
        if gpu_mb is not None and gpu_mb > peaks["peak_gpu_mem_mb"]:
            peaks["peak_gpu_mem_mb"] = gpu_mb
        
        # --- GPU utilization (%) --- 
        util = _sample_gpu_util_pct_nvsmi()
        if util is not None:
            if util > peaks["gpu_util_peak_pct"]:
                peaks["gpu_util_peak_pct"] = util
            peaks["gpu_util_sum_pct"] += util
            peaks["gpu_util_samples"] += 1
        
        
        # --- GPU memory utilization (%) ---
        mem_util = _sample_gpu_mem_util_pct_nvsmi()
        if mem_util is not None:
            if mem_util > peaks["gpu_mem_util_peak_pct"]:
                peaks["gpu_mem_util_peak_pct"] = mem_util
            peaks["gpu_mem_util_sum_pct"] += mem_util
            peaks["gpu_mem_util_samples"] += 1
        
        time.sleep(sample_interval)
    return peaks


def run_program(program_index, program_command, output_file, steps):
    """
    Run a program and write its output to a specified file, recording execution time and memory peaks.
    """
    program_command = os.path.expanduser(program_command)
    output_file = os.path.expanduser(output_file)
    err_file = os.path.splitext(output_file)[0] + ".err"
    metrics_file = os.path.splitext(output_file)[0] + ".metrics.json"

    # Open output and error files
    out_f = open(output_file, "w")
    err_f = open(err_file, "w")

    # Start process
    # Note: preexec_fn=os.setsid makes a new process group on POSIX so children are easier to find.
    # Safe on Linux/macOS. If you need Windows support, drop preexec_fn and accept per-PID monitoring.
    start_time = time.time()
    try:
        process = subprocess.Popen(
            program_command,
            stdout=out_f,
            stderr=err_f,
            shell=True,
            preexec_fn=os.setsid if hasattr(os, "setsid") else None,
        )

        # Start resource monitor
        stop_event = threading.Event()
        peaks_holder = {}

        def monitor():
            peaks = _monitor_peaks(process.pid, stop_event)
            peaks_holder.update(peaks)

        t = threading.Thread(target=monitor, daemon=True)
        t.start()

        # Wait for completion
        retcode = process.wait()
        end_time = time.time()
        stop_event.set()
        t.join(timeout=2.0)

        # Close files so they flush
        out_f.close()
        err_f.close()

        execution_time = end_time - start_time

        
        gpu_util_mean = None
        gpu_util_peak = None
        s = peaks_holder.get("gpu_util_samples", 0)
        if s > 0:
            gpu_util_mean = peaks_holder["gpu_util_sum_pct"] / s
        gpu_util_peak = peaks_holder.get("gpu_util_peak_pct", None)

        gpu_mem_util_mean = (peaks_holder["gpu_mem_util_sum_pct"] / peaks_holder["gpu_mem_util_samples"]
                     if peaks_holder.get("gpu_mem_util_samples") else None)
        gpu_mem_util_peak = peaks_holder.get("gpu_mem_util_peak_pct")

        # Prepare metrics
        metrics = {
            "timestamp": datetime.now().isoformat(timespec="seconds"),
            "program_index": program_index + 1,
            "command": program_command,
            "output": os.path.abspath(output_file),
            "stderr": os.path.abspath(err_file),
            "return_code": retcode,
            "steps": steps,
            "time": {
                "wall_seconds": execution_time
            },
            "gpu": {
                "util_mean_pct": gpu_util_mean,
                "util_peak_pct": gpu_util_peak,
                "mem_util_mean_pct": gpu_mem_util_mean,
                "mem_util_peak_pct": gpu_mem_util_peak
            },
            "memory": {
                "peak_rss_bytes": peaks_holder.get("peak_rss_bytes", None),
                "peak_rss_mb": (peaks_holder.get("peak_rss_bytes", 0) or 0) / (1024 * 1024),
                "peak_gpu_mem_mb": peaks_holder.get("peak_gpu_mem_mb", None)
            },
            "env": {
                "hostname": os.uname().nodename if hasattr(os, "uname") else None
            }
        }

        with open(metrics_file, "w") as mf:
            json.dump(metrics, mf, indent=2)

        # Human-friendly print
        rss_mb = metrics["memory"]["peak_rss_mb"]
        gpu_mb = metrics["memory"]["peak_gpu_mem_mb"]
        gpu_str = f", GPU peak: {gpu_mb:.1f} MB" if gpu_mb is not None else ""
        print(
            f"Process {program_index + 1} finished in {execution_time:.2f}s "
            f"(RSS peak: {rss_mb:.1f} MB{gpu_str}). Metrics: {metrics_file}"
        )

    except Exception as e:
        try:
            out_f.close(); err_f.close()
        except Exception:
            pass
        print(f"Process {program_index + 1} failed with error: {e}", file=sys.stderr)


def work(max_procs, logs_dir, command, steps):
    tasks = []
    for i in range(max_procs):
        output_file = os.path.join(logs_dir, f"{i+1:03d}.log")
        tasks.append((i, command, output_file, steps))

    # Execute with a process pool
    with multiprocessing.Pool(processes=max_procs) as pool:
        results = [pool.apply_async(run_program, t) for t in tasks]
        for r in results:
            r.get()

    # Collate metrics -> JSONL and CSV
    metrics_paths = sorted(glob.glob(os.path.join(logs_dir, "*.metrics.json")))
    summaries = []
    for p in metrics_paths:
        try:
            with open(p, "r") as f:
                m = json.load(f)
            # use filename index as run id
            base = os.path.basename(p)
            run_id = os.path.splitext(base)[0]  # e.g., "001"
            m["run_id"] = run_id
            summaries.append(m)
        except Exception as e:
            print(f"Warning: could not read metrics {p}: {e}")

    summary_jsonl = os.path.join(logs_dir, "summary.jsonl")
    with open(summary_jsonl, "w") as f:
        for m in summaries:
            f.write(json.dumps(m) + "\n")

    def _get(d, dotted, default=None):
        cur = d
        for k in dotted.split("."):
            if isinstance(cur, dict) and k in cur:
                cur = cur[k]
            else:
                return default
        return cur

    summary_csv = os.path.join(logs_dir, "summary.csv")
    fieldnames = [
        "run_id", "program_index", "return_code",
        "time.wall_seconds", "memory.peak_rss_mb", "memory.peak_gpu_mem_mb", 
        "gpu.util_mean_pct", "gpu.util_peak_pct",
        "gpu.mem_util_mean_pct", "gpu.mem_util_peak_pct",
        "output", "stderr", "command", "timestamp", "steps", "ns_per_day"
    ]

    
    
    with open(summary_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for m in summaries:
            w.writerow({
                "run_id": m.get("run_id"),
                "program_index": m.get("program_index"),
                "return_code": m.get("return_code"),
                "time.wall_seconds": _get(m, "time.wall_seconds"),
                "memory.peak_rss_mb": _get(m, "memory.peak_rss_mb"),
                "memory.peak_gpu_mem_mb": _get(m, "memory.peak_gpu_mem_mb"),
                "gpu.util_mean_pct": _get(m, "gpu.util_mean_pct"),
                "gpu.util_peak_pct": _get(m, "gpu.util_peak_pct"),
                "gpu.mem_util_mean_pct": _get(m, "gpu.mem_util_mean_pct"),
                "gpu.mem_util_peak_pct": _get(m, "gpu.mem_util_peak_pct"),
                "output": m.get("output"),
                "stderr": m.get("stderr"),
                "command": m.get("command"),
                "timestamp": m.get("timestamp"),
                "steps": m.get("steps"),
                "ns_per_day": (m.get("steps") * TIME_STEP * 86400) / _get(m, "time.wall_seconds") if _get(m, "time.wall_seconds") else None, 
            })

    print("All processes finished.")
    print(f"Logs dir: {logs_dir}")
    print(f"Summary files: {summary_csv} , {summary_jsonl}")



def run(args):
    data_dir = os.path.expanduser(args.data_dir)   # e.g., TEST/water
    bin_path = os.path.expanduser(args.bin)        # e.g., /path/to/qdyn
    max_procs = int(args.max_processes)

    if not os.path.isdir(data_dir):
        raise FileNotFoundError(f"data_dir not found: {data_dir}")
    if not os.path.exists(bin_path):
        raise FileNotFoundError(f"bin not found: {bin_path}")

    # clear previous logs
    logs_root = os.path.join(os.getcwd(), "benchmark_logs")
    if os.path.exists(logs_root):
        print(f"Removing previous logs at: {logs_root}")
        shutil.rmtree(logs_root)
    
    # data_dir/md.csv should exist
    steps = None
    if not os.path.exists(os.path.join(data_dir, "md.csv")):
        raise FileNotFoundError(f"data_dir/md.csv not found: {data_dir}/md.csv")
    with open(os.path.join(data_dir, "md.csv"), "r") as f:
        # get steps from md.csv
        for line in f:
            if line.startswith("steps;"):
                parts = line.strip().split(";")
                if len(parts) >= 2:
                    steps = int(parts[1])
                    print(f"MD steps found in md.csv: {steps}")
                    break
    if steps is None:
        raise RuntimeError(f"Could not find 'steps' in {data_dir}/md.csv")
    
    current_dir = os.getcwd()
    # Run cpu base line
    print("Running CPU baseline benchmark (1 process)...")
    logs_dir = os.path.join(current_dir, f"benchmark_logs/cpu_baseline")
    os.makedirs(logs_dir, exist_ok=True)
    work(1, logs_dir, f'"{bin_path}" "{data_dir}"', steps)
    
    for process_num in range(1, max_procs + 1):
        print(f"Will run {process_num} processes in parallel:")
        logs_dir = os.path.join(current_dir, f"benchmark_logs/{process_num:02d}_procs")
        os.makedirs(logs_dir, exist_ok=True)
        # One command, repeated N times in parallel
        command = f'"{bin_path}" --gpu "{data_dir}"'
        print(f"Command: {command}")
        print(f"Running {process_num} parallel processes...")
        work(process_num, logs_dir, command, steps)
        
        print(f"Completed benchmark for {process_num} processes.\n")
        time.sleep(2)  # brief pause between runs
    
    
    
            
    # generate report
    out_html = os.path.join(current_dir, "benchmark_report.html")
    make_html_report(logs_root, out_html)

