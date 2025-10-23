import glob
import csv
import re
import statistics
import math
import os
from datetime import datetime
from jinja2 import Environment, FileSystemLoader


def _pct(values, p):
    if not values: return float("nan")
    xs = sorted(values)
    k = max(1, int(round(p/100.0 * len(xs))))
    return xs[k-1]

def make_html_report(logs_root: str, out_html: str):
    cpu_baseline = None
    cpu_csv = os.path.join(logs_root, "cpu_baseline", "summary.csv")
    if os.path.exists(cpu_csv):
        with open(cpu_csv, newline="") as f:
            reader = csv.DictReader(f)
            times = [float(row.get("time.wall_seconds") or float("nan")) for row in reader]
        if times:
            cpu_baseline = max(times)
    print(f"CPU baseline time: {cpu_baseline:.2f} seconds" if cpu_baseline else "No CPU baseline found")
    
    
    # collect rows from each *_procs/summary.csv
    runs = []
    for folder in sorted(glob.glob(os.path.join(logs_root, "*_procs"))):
        m = re.search(r"(\d+)_procs$", os.path.basename(folder))
        if not m: continue
        procs = int(m.group(1))
        csv_path = os.path.join(folder, "summary.csv")
        if not os.path.exists(csv_path): continue
        with open(csv_path, newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                def fget(k):
                    v = row.get(k, "")
                    try: return float(v)
                    except: return float("nan")
                runs.append({
                    "procs": procs,
                    "wall": fget("time.wall_seconds"),
                    "rss_mb": fget("memory.peak_rss_mb"),
                    "gpu_mb": fget("memory.peak_gpu_mem_mb"),
                    "util_mean": fget("gpu.util_mean_pct"),
                    "util_peak": fget("gpu.util_peak_pct"),
                    "mem_util_mean": fget("gpu.mem_util_mean_pct"),
                    "mem_util_peak": fget("gpu.mem_util_peak_pct"),
                    "rc": int(row.get("return_code") or 0),
                    "ns_per_day": fget("ns_per_day"),

                })
    if not runs:
        raise RuntimeError(f"No runs found under {logs_root}")

    # aggregate by procs
    by = {}
    for r in runs:
        by.setdefault(r["procs"], []).append(r)

    procs_list = sorted(by.keys())
    wall_mean = []
    wall_p95  = []
    rss_mean  = []
    gpu_mean  = []
    table_rows = []
    gpu_util_mean = []
    gpu_util_peak = []
    util_mem_mean = []
    util_mem_peak = []
    speedups = []
    ns_per_day = []

    for p in procs_list:
        vals = by[p]
        walls = [v["wall"] for v in vals if math.isfinite(v["wall"])]
        rss   = [v["rss_mb"] for v in vals if math.isfinite(v["rss_mb"])]
        gpu   = [v["gpu_mb"] for v in vals if math.isfinite(v["gpu_mb"])]
        util_mean = [v["util_mean"] for v in vals if math.isfinite(v["util_mean"])]
        util_peak = [v["util_peak"] for v in vals if math.isfinite(v["util_peak"])]
        mem_means = [v["mem_util_mean"] for v in vals if math.isfinite(v["mem_util_mean"])]
        mem_peaks = [v["mem_util_peak"] for v in vals if math.isfinite(v["mem_util_peak"])]
        tmp_ns_per_day = [v["ns_per_day"] for v in vals if math.isfinite(v["ns_per_day"])]

        
        rc_bad = sum(1 for v in vals if v["rc"] != 0)

        w_mean = statistics.mean(walls) if walls else float("nan")
        wall_mean.append(w_mean)
        wall_p95.append(_pct(walls, 95) if walls else float("nan"))
        rss_mean.append(statistics.mean(rss) if rss else float("nan"))
        gpu_mean.append(statistics.mean(gpu) if gpu else float("nan"))

        gpu_util_mean.append(statistics.mean(util_mean) if util_mean else float("nan"))
        gpu_util_peak.append(statistics.mean(util_peak) if util_peak else float("nan"))
        util_mem_mean.append(statistics.mean(mem_means) if mem_means else float("nan"))
        util_mem_peak.append(statistics.mean(mem_peaks) if mem_peaks else float("nan"))
        ns_per_day.append(statistics.mean(tmp_ns_per_day) if tmp_ns_per_day else float("nan"))
        
        
        Tn = max(walls) if walls else float("nan")
        speedup = None
        if cpu_baseline and math.isfinite(Tn) and Tn > 0:
            speedup = cpu_baseline / (Tn / p)
        
        speedups.append(speedup if speedup is not None else float("nan"))

        table_rows.append({
            "procs": p,
            "runs": len(vals),
            "rc_bad": rc_bad,
            "wall_mean": w_mean if walls else float("nan"),
            "wall_p95": _pct(walls, 95) if walls else float("nan"),
            "wall_min": min(walls) if walls else float("nan"),
            "wall_max": max(walls) if walls else float("nan"),
            "rss_mean": statistics.mean(rss) if rss else float("nan"),
            "gpu_mean": statistics.mean(gpu) if gpu else float("nan"),

            "gpu_util_mean": statistics.mean(util_mean) if util_mean else float("nan"),
            "gpu_util_peak": statistics.mean(util_peak) if util_peak else float("nan"),
            "vram_util_mean": statistics.mean(mem_means) if mem_means else float("nan"),
            "vram_util_peak": statistics.mean(mem_peaks) if mem_peaks else float("nan"),
            "Tn": Tn,
            "speedup": speedup,
            "ns_per_day": statistics.mean(tmp_ns_per_day) if tmp_ns_per_day else float("nan"),
        })


    # render template (only data passed)
    template_path = os.path.join(os.path.dirname(__file__), "benchmark_report.html.j2")
    env = Environment(loader=FileSystemLoader(os.path.dirname(template_path) or "."), autoescape=True)
    tpl = env.get_template(os.path.basename(template_path))
    html_out = tpl.render(
        title="Benchmark Report",
        logs_root=os.path.abspath(logs_root),
        generated_at=datetime.now().isoformat(timespec="seconds"),
        payload={
            "procs": procs_list,
            "wall_mean": wall_mean,
            "wall_p95": wall_p95,
            "rss_mean": rss_mean,
            "gpu_mean": gpu_mean,
            "gpu_util_mean": gpu_util_mean,
            "gpu_util_peak": gpu_util_peak,
            "gpu_mem_util_mean": util_mem_mean,
            "gpu_mem_util_peak": util_mem_peak,
            "speedup": speedups,
            "ns_per_day": ns_per_day,
            "table_rows": table_rows,
        },
    )
    with open(out_html, "w", encoding="utf-8") as f:
        f.write(html_out)
    print(f"Report written to: {out_html}")

