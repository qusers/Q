# QGPU core

The CUDA/C++ molecular dynamics engine of Q (QGPU), with two front ends:

- `qdyn`, the command-line program. It runs one Q `.inp` file, like Fortran Q's `qdyn`.
- `qgpu`, a Python package. It loads a system once and runs many stages in one process, without `.inp` files.

Input settings and their defaults follow Fortran Q6 (`src/q6/md.f90`).

## Layout

| Path | Contents |
|---|---|
| `common/` | Parser (`inp_parser`), simulation context, `Handler` (system setup and per-stage configuration), forces, output |
| `cuda/` | GPU kernels and `CudaHandler` |
| `cpu/` | CPU implementations and `CpuHandler` |
| `main.cu` | `qdyn` entry point |
| `python/` | `qgpu` package: `bindings.cpp` (pybind11), `src/qgpu/`, `setup.py`, tests in `test/` |
| `test/shake/` | C++ regression tests (CTest) |
| `CMakeLists.txt` | The build for everything above |

## Requirements

- CUDA Toolkit (`nvcc`) and an NVIDIA GPU
- A C++17 compiler
- CMake 3.24 or newer, and GNU Make
- For the Python package: Python 3.10 or newer. pip installs the remaining build tools.

## Building `qdyn`

Run the commands from the repository root.

```bash
cmake -S src/core -B src/core/build -DCMAKE_CUDA_ARCHITECTURES=89 -DQDYN_SPFP=OFF -DQDYN_BUILD_TESTS=OFF
```

```bash
cmake --build src/core/build --target move --parallel
```

This builds `qdyn` and copies it to `bin/qdyn`.

- The first command (configure) is only needed once, or after changing an option.
- After editing code, rerun only the second command. It recompiles just what changed.

For SPFP, use a separate build folder and run the binary from there. The `move` target always writes `bin/qdyn`, so building it here would replace the double-precision binary.

```bash
cmake -S src/core -B src/core/build_spfp -DCMAKE_CUDA_ARCHITECTURES=89 -DQDYN_SPFP=ON -DQDYN_BUILD_TESTS=OFF
```

```bash
cmake --build src/core/build_spfp --target qdyn --parallel
```

### Options

| Option | Default | Meaning |
|---|---|---|
| `CMAKE_CUDA_ARCHITECTURES` | `89` | GPU architecture: `89` RTX 40-series, `86` RTX 30-series, `80` A100, `90` H100 |
| `QDYN_SPFP` | `OFF` | `ON`: single-precision forces with fixed-point accumulation; `OFF`: double precision |
| `QDYN_BUILD_TESTS` | `ON` | Build the CTest targets in `test/` |
| `QGPU_BUILD_PYTHON` | `OFF` | Build the Python modules (`setup.py` sets this; there is normally no need to set it by hand) |
| `CMAKE_BUILD_TYPE` | `Release` | `Debug` adds `-G -g` for CUDA |

CUDA code is compiled with `--fmad=false`, in both `qdyn` and the Python modules. Results can therefore differ in the last digits from a build that allows fused multiply-add.

## Running `qdyn`

```bash
bin/qdyn --gpu eq1.inp > eq1.log
```

- Without `--gpu`, the CPU implementation is used.
- The energy summaries go to stdout. The energy, trajectory and restart files are named in the `[files]` section of the `.inp` file.

## The Python package

### Install

```bash
cd src/core/python && CMAKE_BUILD_PARALLEL_LEVEL=16 python -m pip install .
```

| Variable | Default | Meaning |
|---|---|---|
| `QGPU_CUDA_ARCH` | `89` | GPU architecture, as `CMAKE_CUDA_ARCHITECTURES` above |
| `QGPU_PRECISIONS` | `double,spfp` | Which modules to build. `double` alone halves the build time. |
| `CMAKE_BUILD_PARALLEL_LEVEL` | `4` | Parallel compile jobs |

- The install is a copy of the package, so reinstall after every change, to C++ or Python files alike.
- The CMake build in `python/build/` is kept between installs, so only changed files are recompiled.
- Build on the machine type you run on. A build for one architecture does not run on another GPU generation.

### Use

```python
from qgpu import Config, Engine, DistRestraint

config = Config(
    topology="dualtop.top", fep="FEP1.fep",
    seed=1, initial_temperature=1,
    shell_radius=25, lrf=True,
    steps=10000, stepsize=0.1, temperature=1, bath_coupling=0.2,
    lambdas=(0.5, 0.5),
    dist_restraints=(DistRestraint(9033, 9078, 0.0, 0.1, 1.5),),
)
engine = Engine(config, device=0, precision="double")   # or precision="spfp"
engine.run()                                            # gentle start; the 2 fs stage below needs it

engine.update_config(steps=5000, stepsize=2.0, temperature=50, bath_coupling=2.0, constrain_hydrogens=True)
engine.run()                                            # continues from the current state

engine.update_config(lambdas=(0.6, 0.4), energy=10, energy_file="md_0600_0400.en")
engine.run()
```

- **`Config`** has the same settings, names and meanings as the `.inp` keys.
  - A setting left out is filled in by the C++ parser with the same default an `.inp` file would get. Python keeps no default values of its own.
  - `steps`, `stepsize` and `temperature` are required.
  - A fresh start (no `restart`) needs `initial_temperature` and a `seed`.
- **`Config.from_inp("eq5.inp", steps=100)`** reads an existing `.inp` file, with overrides.
- **`engine.update_config(...)`** changes the settings that may change between runs. Its parameters are listed explicitly, so editors can autocomplete and check them. Settings that define the system (`shell_radius`, `lrf`, `restart`, cut-offs, ...) cannot change: create a new `Engine` for those.
- **`engine.get_state()` / `engine.set_state(state)`** save and restore coordinates, velocities and the water-polarisation state, e.g. to start two FEP branches from the same midpoint.
- **Errors** are raised as `qgpu.InputError` (bad input), `qgpu.ConstraintError` (SHAKE/LINCS failure), `qgpu.CudaError`, or `qgpu.QgpuError`.
  - A rejected `update_config()` leaves the engine usable.
  - After an error during `run()`, call `set_state()` or create a new `Engine`.

`Engine` runs inside the calling Python process, and its native output goes to that process's stdout.

- Run one `Engine` per process.
- For parallel replicates, use separate processes. For example, use `concurrent.futures.ProcessPoolExecutor` with the `spawn` start method, as in `python/test/test_qligfep_workflow.py`.
- Do not run several engines in threads of one process.

### Examples

- `python/test/test_qligfep_one_replicate.py`: one replicate of the qligfep FEP workflow (eq1 to eq5, then 101 lambda windows). It is written entirely in Python and writes `md_XXXX_YYYY.en` files for qfep.
- `python/test/test_qligfep_workflow.py`: all replicates in parallel processes.
- `python/test/test_short_simulation.py`: run an existing `.inp` file.

## Tests

Configure with `QDYN_BUILD_TESTS=ON` (the default), build, then run:

```bash
ctest --test-dir src/core/build --output-on-failure
```
