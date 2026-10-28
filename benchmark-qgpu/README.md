# QGPU Performance Benchmark Guide: CDK2 Example

This guide explains how to benchmark QGPU with the CDK2 `eq5.inp` input using `benchmark-qgpu/main.py`. It measures runtime, throughput, and resource usage at different concurrency levels and generates an HTML report. Run the commands in Bash on Linux or WSL.

## Environment requirements

Prepare the following tools and make sure their commands are available in your terminal:

- **Python 3**, with `psutil` and `jinja2` installed.
- **CUDA compiler `nvcc`**, provided by the CUDA Toolkit, and a compatible C++ compiler.
- **An NVIDIA GPU and driver**, with `nvidia-smi` available for resource monitoring.
- **Git and GNU Make**, for obtaining and compiling the code.

This workflow uses QGPU's `bin/qdyn` for both the CPU baseline and GPU runs. It does not require a Fortran executable or Matplotlib.

## Getting the code and test data

If you already followed the [CDK2 energy comparison guide](../test/QGPU_TESTING.md), reuse those repositories. Otherwise, clone them under the same parent directory:

```bash
git clone https://github.com/qusers/Q.git
git clone https://github.com/goodstudyqaq/qligfepv2-BenchmarkExperiments.git
cd Q
git switch feature/qgpu
```

Run the remaining commands from the Q repository root.

## 1. What the benchmark runs

The script first runs one simulation using QGPU's CPU mode as a baseline. It then runs independent GPU simulations at each requested concurrency level.

For example, `--concurrency 1 2 4 8` runs four GPU batches: one simulation, two simultaneous simulations, four simultaneous simulations, and eight simultaneous simulations. These are separate copies of the same simulation; they do not divide one simulation across processes. The script does not automatically assign processes to different GPUs or enable CUDA MPS.

Each concurrency level is run once. The script pauses for 30 seconds after each GPU batch, then generates the report after all batches finish. It measures performance; it does not compare energy agreement with QFortran.

## 2. Dependencies and compilation

### 2.1 Python dependencies

```bash
python -m pip install psutil jinja2
```

`psutil` monitors process memory, and `jinja2` generates the HTML report.

### 2.2 Compiling QGPU

From the Q repository root, run:

```bash
cd src/core
make double
cd ../..
```

This builds `bin/qdyn` with `SPFP=0`. The current Makefile uses `-arch=sm_89` by default; adjust it for your target GPU. Use `make spfp` instead if you want to benchmark an SPFP build.

## 3. CDK2 input files

The input files are in the `qligfepv2-BenchmarkExperiments` repository. Relative to the Q repository root, the directory is:

```text
../qligfepv2-BenchmarkExperiments/perturbations/experiment/cdk2_test/
```

| File | Purpose |
| --- | --- |
| `eq5.inp` | MD settings and references to the other files |
| `dualtop.top` | System topology |
| `eq4.re` | Input restart |
| `FEP1.fep` | FEP definitions |

Keep these files together. The benchmark accepts `eq5.inp` directly and reads `steps` and `stepsize` from its `[MD]` section. No CSV conversion is needed.

Before running, choose the number of steps in `eq5.inp`. The original input specifies 50,000 steps, which also applies to the CPU baseline. A short run can check that the workflow starts correctly, but a 10-step test is too short for a useful performance measurement: process startup and initialization can dominate the timing. For performance measurements, increase the step count until the measured throughput is reasonably stable across repeated runs, keeping the same settings when comparing concurrency levels.

The throughput calculation uses the actual `stepsize` in the input. For the original CDK2 input, this is 2 fs.

## 4. Running the benchmark

### 4.1 Selecting concurrency levels

From the Q repository root, run:

```bash
python benchmark-qgpu/main.py \
  --input ../qligfepv2-BenchmarkExperiments/perturbations/experiment/cdk2_test/eq5.inp \
  --bin ./bin/qdyn \
  --concurrency 1 2 4 8
```

The script runs the CPU baseline, then the selected GPU batches. Each process receives its own copy of the `.inp` file. For this CDK2 input, trajectory and final restart outputs are written under that process's directory, while the topology, FEP, and input restart are read from the source directory.

### 4.2 Testing every concurrency level up to a limit

To test 1, 2, 3, and 4 simultaneous GPU simulations, use:

```bash
python benchmark-qgpu/main.py \
  --input ../qligfepv2-BenchmarkExperiments/perturbations/experiment/cdk2_test/eq5.inp \
  --bin ./bin/qdyn \
  --max-processes 4
```

Use either `--concurrency` or `--max-processes`; they cannot be combined. For an initial run, `--concurrency 1` is enough to check the CPU baseline and a single GPU simulation.

## 5. Output files

The script writes its results under the directory from which it was launched:

| Output | Contents |
| --- | --- |
| `benchmark_report.html` | Performance charts and summary table; open it in a browser |
| `benchmark_logs/cpu_baseline/` | CPU baseline results |
| `benchmark_logs/01_procs/`, `02_procs/`, etc. | Results for each selected GPU concurrency level |

Within each batch directory, individual processes have numbered folders such as `001/`. Each contains `qdyn.log`, `qdyn.err`, `qdyn.metrics.json`, and an `input/` directory with the staged input and its simulation outputs. Each batch also has `summary.csv` and `summary.jsonl` files.

The script keeps these files automatically; no `-k` option is needed. **Running the benchmark again deletes the existing `benchmark_logs/` directory, and a successful run replaces `benchmark_report.html`.** Save previous results before rerunning if you want to compare runs. If a Qdyn process returns a nonzero exit code, the script stops with an error after saving that batch's summaries; inspect its `qdyn.log` and `qdyn.err`. An HTML report left from an earlier run is not a report for the failed run.

## 6. Reading the report

| Metric | Meaning |
| --- | --- |
| Wall mean / p95 / min / max | Runtime statistics for the processes at that concurrency level, in seconds |
| RSS mean (MB) | Mean of the measured per-process peak host-memory usage |
| GPU mean (MB) | Mean of the measured per-process peak GPU-memory usage |
| GPU util (%) | Sampled GPU compute utilization |
| VRAM util (%) | Sampled fraction of GPU memory occupied |
| Speedup | Aggregate throughput speedup relative to the single CPU baseline |
| Total ns/day | Combined simulated nanoseconds per day of wall-clock time |
| Non-zero RC | Number of processes that returned a nonzero exit code |

For each process, throughput is calculated as:

```text
ns/day = steps × stepsize_fs × 10^-6 × 86400 / wall_seconds
```

The report's total ns/day is the mean per-process ns/day multiplied by the number of concurrent simulations. Its speedup is calculated as:

```text
speedup = CPU baseline time × concurrency / longest GPU process time
```

At concurrency 1, this is the CPU-to-GPU runtime ratio. At higher concurrency levels, it describes aggregate throughput rather than the speed of an individual simulation. Wall time includes process startup and initialization.

Use total ns/day to identify which tested concurrency level provides the highest combined throughput. Increasing concurrency may improve total throughput while making each individual simulation slower. GPU and VRAM utilization are device-wide measurements, averaged across the GPUs reported by `nvidia-smi`, so other GPU workloads can affect them.

## Appendix: Command-line options

| Option | Meaning |
| --- | --- |
| `--input PATH` | Required; path to a native Q `.inp` file |
| `--bin PATH` | Required; path to the QGPU executable |
| `--concurrency N [N ...]` | Specific concurrency levels, such as `1 2 4 8`; duplicate values are removed and levels run in ascending order |
| `--max-processes N` | Tests every concurrency level from 1 through N; `--max_processes` is also accepted |
| `--help` | Displays command-line help |

Either `--concurrency` or `--max-processes` is required. The step count and step size are taken from the input file.
