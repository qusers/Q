# QligFEP

**Q** is a set of Molecular Dynamics (MD) tools tailored to the following specific kinds of free energy calculations:

1. Free Energy Perturbation (FEP)
2. Empirical Valence Bond (EVB)
3. Linear Interaction Energies (LIE)

This repository is devoted to **QligFEP**, an automated workflow for small molecule free energy calculations in Q.

## Table of Contents

- [⚙️ Installation](#️-installation)
  - [Linux](#linux)
  - [MacOS](#macos)
  - [Compiling Q for HPC (MPI support)](#compiling-q-for-hpc-mpi-support)
  - [Compiling Q for local use (non-MPI)](#compiling-q-for-local-use-non-mpi)
  - [Setting up HPC configurations](#setting-up-hpc-configurations)
- [⌨️ Command line interface (CLI)](#️-command-line-interface-cli)
- [📊 Benchmarking](#-benchmarking)
- [📚 Citations](#-citations)
- [⏩ Q-GPU](#-q-gpu)
  - [Note to the current version](#note-to-the-current-version)
  - [Testing](#testing)

## ⚙️ Installation

We recommend that you use `mamba` or, preferably, its lightweight version `micromamba`. Please check this link on how to [install it](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html).

Once you have `micromamba` installed and have already cloned this repo, you can create the environment with:

### Linux
```bash
micromamba create -n qligfep_new python=3.11
micromamba activate qligfep_new
micromamba install gfortran=11.3.0 openff-toolkit=0.16.4 "openff-utilities>=0.1.12" openff-forcefields=2024.09.0 openmm=8.1.1 "openff-nagl>=0.3.8" lomap2 kartograf michellab::fkcombu -c conda-forge --yes
```

Now that you have the environment ready and activated, [clone the repository](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository), enter the `Q` directory with `cd Q/`, and install qligfep:
```bash
python -m pip install -e .
```

The `qprep` Fortran binary will be automatically compiled during installation.

<details>
<summary>To install everything in one line...</summary>

```bash
micromamba create -n qligfep_new python=3.11 gfortran=11.3.0 openff-toolkit=0.16.4 "openff-utilities>=0.1.12" openff-forcefields=2024.09.0 openmm=8.1.1 "openff-nagl>=0.3.8" lomap2 kartograf michellab::fkcombu -c conda-forge --yes && micromamba activate qligfep_new && python -m pip install -e .
```
</details>

### MacOS

Similar to Linux, [clone the repository](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository), enter the `Q` directory with `cd Q/`, create the environment and install:

``` bash
micromamba create -n qligfep_new python=3.11 gfortran=11.3.0 openff-toolkit=0.16.4 "openff-utilities>=0.1.12" openff-forcefields=2024.09.0 openmm=8.1.1 "openff-nagl>=0.3.8" lomap2 kartograf davidararipe::kcombu_bss -c conda-forge --yes
micromamba activate qligfep_new
python -m pip install joblib scipy tqdm
python -m pip install -e .
```

The `qprep` Fortran binary will be automatically compiled during installation.

<details>
<summary>To install everything in one line...</summary>

```bash
micromamba create -n qligfep_new python=3.11 gfortran=11.3.0 openff-toolkit=0.16.4 "openff-utilities>=0.1.12" openff-forcefields=2024.09.0 openmm=8.1.1 lomap2 kartograf davidararipe::kcombu_bss -c conda-forge --yes && micromamba activate qligfep_new && python -m pip install joblib scipy tqdm && python -m pip install -e .
```
</details>

### Compiling Q for HPC (MPI support)

> [!IMPORTANT]
> The current Q implementation relies on `slurm` for job management and submission. The basic `qprep` tool for topology creation is automatically compiled during pip installation and is sufficient for preparing inputs. When submitting jobs, QligFEP uses the MPI-enabled `qdynp` program (_p for parallel_) to run the molecular dynamics simulations. To actually run these simulations, you need to compile Q as described below:

On your HPC system, load the appropriate modules (system-dependent). We recommend using the GCC compiler suite and OpenMPI, as those are commonly available and compatible with `qdynp`. To check for module availability, use the command `module spider openmpi` or `module avail openmpi`.

In the output, look for a version compiled with GCC (e.g., `OpenMPI/4.1.4-GCC-11.3.0`) and load it using the `module load` command. After loading the module, navigate to the `src/q6` folder in the Q repository and compile both the serial and MPI versions of Q with the commands `make all` and `make mpi`. In the example, we show how to do this on the Snellius, the Dutch national supercomputer:

_Example_:
```bash
module load 2021
module load gompi/2021a
```

Then compile the serial and the MPI-enabled versions of Q:
```bash
cd src/q6
make all COMP=gcc
make mpi COMP=gcc
```

> [!TIP]
> Module names and versions are system-dependent. When in doubt, reach out to your system administrator. In general, we recommend finding an OpenMPI module compiled with GCC version 11.3.0. Users can also refer to the `settings.py` file in this repository, which outlines the modules we used on other HPC systems, as [described below](#setting-up-hpc-configurations).

### Compiling Q for local use (non-MPI)

For your convenience, our base environment installation includes `gfortran=11.3.0`, which enables you to compile Q locally without MPI support. This is useful for testing purposes. To compile it, navigate to the `src/q6` folder in the Q repository and run:
```bash
cd src/q6
make all COMP=gcc
```

## Setting up HPC configurations

Currently, we require job configurations to be set in the `settings.py` file located in `src/QligFEP/settings/`. Check [here](/src/QligFEP/settings/settings.py) for a list of different HPCs we have successfully ran RBFE simulations on. To add your own HPC system, please follow the format used in the file. In the example, we show how to add a custom HPC configuration named `MY_HPC`:

```python
MY_HPC = {
    "NODES": "1",  # We recommend not to change this
    "NTASKS": "8",  # Number tasks (processes). Check the preprint for guidance on this value
    "TIME": "0-06:00:00",  # time for job execution; formatted as d-hh:mm:ss
        "MODULES": nljoin(
        [
            "module purge",  # Clear all loaded modules
            "module load OpenMPI/4.1.4-GCC-11.3.0",  # Load the MPI module used for compiling Q
        ]
    ),
    **Q_PATHS, # Keep this line as is; it passes the paths to Q executables
}

CLUSTER_DICT = {
    "CSB": CSB,
    # ...
    "MY_HPC": MY_HPC, # Make sure to add your HPC configuration here to use it on the CLI
}
```

When using the created configuration, make sure to pass the cluster name (e.g., `MY_HPC`) to the `qligfep` or to the `setupFEP` CLI using the `--cluster` argument.

## ⌨️ Command line interface (CLI)

Now you're set with the qligfep package. This includes the command-linde-interfaces (CLIs):

1. `qcog`: calculates the center of geometry (COG) of a ligand in a PDB/SDF file. If multiple ligands are found in sdf, the program will calculate the COG for all of them
2. `pdb2amber`: formats a PDB file to be used with Q's implementation of the AMBER forcefield;
3. `qprep_prot`: creates an input file for qprep (fortran program) and runs it to either: 1) solvate the protein structure; 2) create the water sphere.
4. `qparams`: used to generate ligand parameters;
5. `qlomap`: wraps `Lomap` to generate the `.json` perturbation mapping;
6. `qmapfep`: in-house developed method to generate the `.json` perturbation mapping, interactively visualize and add or remove edges.
7. `qligfep`: main CLI for running QligFEP simulations.
8. `setupFEP`: sets up all the the QligFEP files for a simulation, including protein and water systems.
9. `qligfep_analyze`: CLI to analyze the results of a QligFEP simulation.
10. `ligalign`: aligns a set of ligands to a reference ligand based on their maximum common substructure (MCS).

## Tutorials

We are working on the documentation and tutorials for QligFEP. In the meantime, please refer to the Tyk2 case study available in the [tutorials directory](/tutorials/Tyk2/README.md). In addition to that, you can check the [benchmarking section](#-benchmarking) below, which contains the link to our benchmarking repository with scripts to reproduce the results.

# 📊 Benchmarking

To check and reproduce QligFEP performance results, please refer to our [benchmarking repository](https://github.com/qusers/qligfepv2-BenchmarkExperiments).

For the preprint describing the benchmarking results, see:

> Alencar Araripe D, Díaz Holguín A, Poso A, van Westen GJP, Åqvist J, Gutiérrez-de-Terán H, et al. Doing More with Less: Accurate and Scalable Ligand Free Energy Calculations by Focusing on the Binding Site. ChemRxiv. 2025; [doi:10.26434/chemrxiv-2025-x3r3z](https://doi.org/10.26434/chemrxiv-2025-x3r3z-v3)

# 📚 Citations
Q6:       https://doi.org/10.1016/j.softx.2017.12.001

Q         https://doi.org/10.1016/S1093-3263(98)80006-5

QligFEP:  https://doi.org/10.1186/s13321-019-0348-5

QresFEP:  https://doi.org/10.1021/acs.jctc.9b00538

# ⏩ Q-GPU

**Q-GPU** is an adaptation of **Q** version 5.06 to run on GPUs.

## Note to the current version
The Qprep tool from **Q** is needed for the preparation of molecular topology files required by the MD engine Qdyn. Currently, Qprep is provided as fortran code, which is compiled upon installation. The workflow for a **Q-GPU** free energy simulation consists then in:

- An initial topology preparation stage that runs on a regular CPU  
- MD sampling using Qdyn, which runs on a CUDA-based GPU  
- The FEP analysis tool (qfep) provided in python (running both in GPU or CPU)  

> ⚠️ Integration with the QligFEP workflow is currently under active development and not yet available. For now, we only provide instructions to run our test cases. A publication describing Q-GPU and its performance is in preparation. Please refrain from using this version until the publication is out.

## Testing ##
**Q-GPU** includes various tests that compare the output of the original fortran code with the C/CUDA code. They are situated in the test folder and include:

1. interactions  
2. physical-properties  

The first folder includes test cases for the different type of interactions in **Q**, that is water-water (w-w), solute-solute (p-p) and Qatom-Qatom (q-q) interactions, and any mixture thereof. These tests run a single point energy calculation and are compared with the output from Q5.07.

To compile the code and run the tests, you must checkout the GPU feature branch and build the components individually in their respective folders.

**Prerequisites:**
*   Installation of QligFEP environment as described above.
*   A CUDA compiler (Tested on CUDA/10.1.243)

**Steps:**

1.  **Checkout the feature branch:**
    ```bash
    git checkout feature/qgpu
    ```

2.  **Build the test suite (q6):**
    Navigate to the `q6` folder and build:
    ```bash
    cd src/q6
    make test
    ```

3.  **Build the QGPU engine (core):**
    Navigate to the `core` folder and build:
    ```bash
    cd ../core
    make
    ```

4.  **Run the tests:**
    Once compiled, you can run the test scripts.
    ```bash
    python runTEST.py -a gpu -t 100
    ```