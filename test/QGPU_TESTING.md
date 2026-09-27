# QGPU and QFortran Testing Guide: CDK2 Example

This guide explains how to run QGPU and QFortran with the CDK2 `eq5.inp` input and compare their energies using `test/runTEST.py`. Run the commands in Bash on Linux or WSL.

The guide is based on the current `test/runTEST.py`, `src/Qgpu/compare.py`, both programs' makefiles, and the CDK2 input files. The commands have been checked statically. No compilation or molecular dynamics tests were performed when preparing this guide, so it does not establish that this test case passes.

## Environment requirements

Before starting, prepare the following tools and make sure their commands are available in your terminal:

- **Python 3**: required to run the test script, with NumPy and Matplotlib installed (see Section 2.1).
- **Fortran compiler `gfortran`**: required to compile the QFortran test executable.
- **CUDA compiler `nvcc`**: provided by the CUDA Toolkit and required to compile QGPU. A C++ compiler compatible with the CUDA Toolkit is also required.
- **Git and GNU Make**: required to obtain the code and run the build commands.

Running tests in GPU mode also requires an NVIDIA GPU and a suitable driver.

## Getting the code and test data

Clone both repositories into your chosen project directory, then enter the Q repository and switch to the `feature/qgpu` branch:

```bash
git clone https://github.com/qusers/Q.git
git clone https://github.com/goodstudyqaq/qligfepv2-BenchmarkExperiments.git
cd Q
git switch feature/qgpu
```

Keep the two repositories under the same parent directory. Run subsequent commands from the Q repository root; the build instructions also return to this directory when finished. Test data is accessed through `../qligfepv2-BenchmarkExperiments/`.

## 1. Purpose and scope

The test uses the same topology, FEP, restart, and MD settings to run the following sequence:

1. QFortran: the Q6 test executable, `src/q6/bin/q6/qdyn_test`.
2. QGPU: `bin/qdyn`, running on the CPU or GPU according to the selected option.
3. Energy comparison: the script parses both logs, compares energy components for each frame, and reports `TRUE` or `FALSE`.

This workflow checks energy agreement for the given input. It does not directly compare atomic coordinates, forces, trajectory files, or final free energies, and it does not measure execution speed. A pass only means that the compared frames and energy components satisfy the specified tolerance.

## 2. Dependencies and compilation

### 2.1 Python dependencies

Install the two Python packages required by the test script, NumPy and Matplotlib:

```bash
python -m pip install numpy matplotlib
```

### 2.2 Compiling the QFortran test executable

From the Q repository root, run:

```bash
cd src/q6
make test
cd ../..
```

`make test` uses `gfortran` by default and builds and installs `src/q6/bin/q6/qdyn_test`. The last command returns to the repository root for the following steps.

### 2.3 Compiling QGPU

From the Q repository root, run:

```bash
cd src/core
make double
cd ../..
```

`make double` cleans the build and rebuilds with `SPFP=0`, producing `bin/qdyn`. The last command returns to the repository root. The current Makefile includes `-arch=sm_89` by default; adjust the build settings for your target GPU. To test an SPFP build, use `make spfp` and note this in your test records.

To view the test script's options, run from the repository root:

```bash
python test/runTEST.py --help
```

## 3. CDK2 test inputs

The test files are in `perturbations/experiment/cdk2_test/` inside the `qligfepv2-BenchmarkExperiments` repository. Relative to the Q repository root, the path is:

```text
../qligfepv2-BenchmarkExperiments/perturbations/experiment/cdk2_test/
```

This guide uses `eq5.inp` in that directory as the test input. The directory contains four files:

| File | Purpose |
| --- | --- |
| `eq5.inp` | MD settings and references to the required files |
| `dualtop.top` | System topology |
| `eq4.re` | Input restart specified by `eq5.inp`, used to continue from an existing state |
| `FEP1.fep` | FEP definitions corresponding to the atom numbering and states in the topology |

Treat these four files as one input set. When replacing the topology, FEP, or restart file, ensure that they describe a matching system.

The key settings in the current `eq5.inp` are:

| Setting | Current value | Meaning for the test |
| --- | --- | --- |
| `steps` | 50000 | The full input runs for 50,000 steps |
| `stepsize` | 2.0 fs | Gives a nominal simulation length of 100 ps |
| `temperature` | 298 K | Target temperature |
| `bath_coupling` | 10.0 | Temperature coupling parameter |
| `shake_hydrogens / shake_solute / shake_solvent` | All `on` | Enables the corresponding SHAKE constraints |
| `lrf` | `off` | Disables LRF |
| `separate_scaling` | `on` | Enables separate scaling |
| All cut-offs | 99 | Retains the input's cut-off settings |
| `shell_radius / shell_force` | 25 / 10.0 | Spherical boundary settings |
| `radial_force / polarisation_force` | 60.0 / 20.0 | Solvent restraint parameters, with `polarisation on` |
| `output` | 1 | Outputs every step for frame-by-frame comparison |
| `trajectory` | 100 | Trajectory output interval |
| `non_bond` | 25 | Value in the input file; the current `feature/qgpu` branch does not implement nonbonded updates at this interval |
| `[lambdas]` | `0.500 0.500` | Weights of the two states |
| `[distance_restraints]` | 24 entries | Retains the atom pairs and restraint parameters in the input |

The `[files]` section references `dualtop.top`, `eq4.re`, and `FEP1.fep`, and names the trajectory and final restart outputs `eq5.dcd` and `eq5.re`.

## 4. Running the test

### 4.1 Comparing GPU mode with QFortran

Run directly from the Q repository root:

```bash
python test/runTEST.py \
  --inp ../qligfepv2-BenchmarkExperiments/perturbations/experiment/cdk2_test/eq5.inp \
  -a gpu -k All
```

The script automatically copies the inputs, runs QFortran and QGPU, and compares their energies. `-a gpu` selects GPU mode, and `-k All` keeps the test files. Results are stored in `1_eq5/` under the current directory, while comparison outcomes appear in the terminal. Running the command again overwrites this results directory.

Before running, reduce `steps` in `eq5.inp` as described in Section 5. The original input specifies 50,000 steps, and QFortran runs slowly.

### 4.2 Using CPU mode (optional)

To compare QGPU's CPU mode with QFortran, change `-a gpu` to `-a cpu`:

```bash
python test/runTEST.py \
  --inp ../qligfepv2-BenchmarkExperiments/perturbations/experiment/cdk2_test/eq5.inp \
  -a cpu -k All
```

This command also reruns QFortran and overwrites the same `1_eq5/` directory. To retain the GPU results, save that directory first or use `-w` to specify another existing working directory using an absolute path.

## 5. Reducing the number of steps before running

**Before executing the commands in Section 4, edit `steps` under `[MD]` in `eq5.inp`.** The original input specifies 50,000 steps. Since QFortran runs slowly, start with 10 steps for routine checks:

```text
[MD]
steps                     10
```

Change only the `steps` line and keep the other settings. The script uses the step count from the input file; `-t` does not override `steps` when using `--inp`.

Ten steps can be used to check input handling and initial energies. Restore 50,000 steps only when a full test is needed, and allow for a longer runtime.

Keep `output 1` in `[intervals]` to support frame-by-frame energy comparison.

## 6. Interpreting the results

### 6.1 Energy components compared

The script calls `src/Qgpu/compare.py` to check the following categories in each compared frame:

| Category | Main components |
| --- | --- |
| Nonbonded interactions | Electrostatic and van der Waals energies for solute–solute, solute–solvent, solvent–solvent, and Q-atom interactions |
| Bonded interactions | Bond, angle, torsion, and improper energies for solute, solvent, and Q-atoms |
| Restraint energies | Total, Ufix, Uradx, Upolx, Ushell, Upres |
| Total energies | Utot, Upot, Ukin |

If the entire `Q-atom` or `solvent` group is absent from the QFortran data, the corresponding optional comparisons are skipped. Per-state details in QGPU's `[q-energies]` section are parsed but are not independently used to determine whether the test passes. Temperature is also excluded from the pass/fail decision.

### 6.2 Meaning of the tolerance

Before comparison, the script formats each corresponding QGPU energy to two decimal places and compares it with the value in the QFortran log:

```text
abs(QFortran log value - QGPU value rounded to two decimal places) <= tolerance
```

`--tolerance` is an absolute energy difference in kcal/mol. The default, `0.0`, requires these log values to match; it does not require bitwise equality of the internal floating-point results. For example, `--tolerance 0.01` allows an absolute difference of up to 0.01 kcal/mol. This illustrates the option's meaning and is not a recommended acceptance threshold.

### 6.3 Pass and fail output

When all frames that were compared pass, the script prints:

```text
Passed test? TRUE
```

When a difference is found, it prints the energy component, both values, and the following message for the affected frame:

```text
Compared energies for frame ...
Passed test? FALSE
```

The labels `Q5` and `Q7` in mismatch messages are historical names. In the current workflow, they refer to QFortran/Q6 and QGPU, respectively. If a QGPU frame is missing, the script prints `Missing QGPU frame for Q6 frame ...`.

**Confirm that execution completed normally, a pass message is present, and no failure or exception messages were reported.** An energy comparison failure does not currently cause the script to return a nonzero exit status, so `$?` or a job scheduler's success status alone cannot establish that the test passed.

## Appendix A: Common options

| Option | Purpose and behavior in direct input mode |
| --- | --- |
| `--inp PATH [PATH ...]` | Runs one or more input files directly, preparing a separate directory for each; does not automatically chain the output of one input into the next |
| `-a gpu` / `-a cpu` | Required; selects QGPU's execution mode. QFortran runs in either case |
| `--qgpu-input inp` | Input parsing mode required for direct input tests; `inp` is the default |
| `-w PATH` | Existing parent working directory; an absolute path is recommended |
| `-k All` | Keeps the test files for inspection |
| `--tolerance VALUE` | Absolute energy comparison tolerance; default 0.0 |
| `--verbose` | Also prints the QGPU log to the script's output |
| `--avg` | Prints energy means and standard deviations |
| `--plot` | Saves `$WORKDIR/Utot.png` |
| `--seed / --temperature / --fep-file` | Replaces the corresponding placeholders only; does not override literal values in the current CDK2 input |

## Appendix B: Energy field mapping

| QFortran log group | QGPU field | Corresponding components |
| --- | --- | --- |
| First two `solute` entries | `nonbonded.pp` | Electrostatic, vdW |
| Last four `solute` entries | `bonded.p` | Bond, angle, torsion, improper |
| First two `solvent` entries | `nonbonded.ww` | Electrostatic, vdW |
| Last four `solvent` entries | `bonded.w` | Bond, angle, torsion, improper |
| `solute-solvent` | `nonbonded.pw` | Electrostatic, vdW |
| First two `Q-atom` entries | `nonbonded.qx` | Electrostatic, vdW |
| Last four `Q-atom` entries | `bonded.qp` | Bond, angle, torsion, improper |
| `restraints` | `restraint` | Total, Ufix, Uradx, Upolx, Ushell, Upres, in order |
| `SUM` | `total` | Utot, Upot, Ukin, in order |
