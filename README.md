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
  - [Compiling Q for HPC (MPI support)](#️compiling-q-for-hpc-mpi-support)
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

> 🖥️ **For HPC users only**: If you need MPI support for parallel simulations on HPC systems, you'll need to manually compile the MPI version of Q.

The basic `qprep` tool is automatically compiled during pip installation and is sufficient for local use. However, for running parallel simulations on HPC systems, you need to compile the MPI-enabled version:

On your HPC system, load the appropriate modules (system-dependent):

_Example for Snellius_:
```bash
module load 2021
module load gompi/2021a
```

Then compile the MPI version:
```bash
cd src/q6
make mpi COMP=gcc
```

Check [here](/src/QligFEP/settings/settings.py) for a list of different HPCs we have successfully ran RBFE simulations on. Module availability and loading is system-dependent, so please check with your system administrator if you have any issues.

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

# 📊 Benchmarking

To check and reproduce QligFEP performance results, please refer to our [benchmarking repository](https://github.com/qusers/qligfepv2-BenchmarkExperiments).

For the preprint describing the benchmarking results, see:

> Alencar Araripe D, Díaz Holguín A, Poso A, van Westen GJP, Åqvist J, Gutiérrez-de-Terán H, et al. QligFEP-2: an automated workflow for small molecule free energy calculations in Q. ChemRxiv. 2025; [doi:10.26434/chemrxiv-2025-x3r3z](https://doi.org/10.26434/chemrxiv-2025-x3r3z)

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

> ⚠️ Integration with the QligFEP workflow is currently under active development and not yet available. For now, we only provide instructions to run our test cases.publication describing Q-GPU and its performance is in preparation Please refrain from using this version until the publication is out.

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
    make
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
    python runTEST.py
    ```