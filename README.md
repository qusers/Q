# Q

Q is a monorepo for the Q molecular dynamics engine and accompanying free energy workflows. It currently contains Q, QligFEP, and QresFEP.

In this README, **Q repository** refers to the complete monorepo, while **Q engine** refers to the CPU molecular dynamics and free energy engine in `src/q6`.

## Components

| Component | Purpose | Status |
| --- | --- | --- |
| **Q engine** | CPU molecular dynamics and free energy engine for FEP, EVB, and LIE calculations | Available |
| **QligFEP** | Automated ligand relative binding free energy workflows | Available |
| **QresFEP** | Residue-level free energy workflows | Protein thermostability protocol available |
| **Q-GPU** | GPU implementation of the Q engine | Integration pending |

The Python workflows are currently distributed together in the `QligFEP`
package. Installing that package provides both the QligFEP and QresFEP commands.

## Choose a workflow

- Use the **Q engine** to run molecular dynamics, FEP, EVB, or LIE calculations directly.
- Use **QligFEP** to calculate ligand relative binding free energies.
- Use **QresFEP** to calculate mutation-induced changes in free energy. This repository currently includes only the protein thermostability protocol.

## Table of contents

- [⚙️ Installation](#️-installation)
  - [Linux](#linux)
  - [macOS](#macos)
  - [PyMOL for QresFEP mutant preparation](#pymol-for-qresfep-mutant-preparation)
  - [Compiling Q for HPC (MPI support)](#compiling-q-for-hpc-mpi-support)
  - [Compiling Q for local use (non-MPI)](#compiling-q-for-local-use-non-mpi)
  - [Setting up HPC configurations](#setting-up-hpc-configurations)
- [⌨️ Command line interface (CLI)](#️-command-line-interface-cli)
- [Tutorials](#tutorials)
  - [Protein thermostability with QresFEP](#protein-thermostability-with-qresfep)
  - [Non-equilibrium FEP (NEQ²)](#non-equilibrium-fep-neq2)
- [📊 Benchmarking](#-benchmarking)
- [📚 Citations](#-citations)
- [⏩ Q-GPU](#-q-gpu)
  - [Note to the current version](#note-to-the-current-version)
  - [Testing](#testing)

## ⚙️ Installation

We recommend `micromamba`, the lightweight version of `mamba`. Follow the
[micromamba installation instructions](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html),
then create the project environment:

### Linux
```bash
micromamba create -n qligfep_new python=3.11
micromamba activate qligfep_new
micromamba install gfortran=11.3.0 openff-toolkit=0.17.1 "openff-utilities>=0.1.12" openff-forcefields=2026.01.0 openmm=8.1.1 openff-nagl=0.5.4 openff-nagl-models=2025.9.0 lomap2 kartograf=1.0.1 michellab::fkcombu konnektor -c conda-forge --yes
```

After activating the environment, [clone the repository](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository), enter the `Q` directory with `cd Q/`, and install the `QligFEP` Python package. This package contains both the QligFEP and QresFEP workflows:
```bash
python -m pip install -e .
```

The `qprep` Fortran binary will be automatically compiled during installation.

<details>
<summary>To install everything in one line...</summary>

```bash
micromamba create -n qligfep_new python=3.11 gfortran=11.3.0 openff-toolkit=0.17.1 "openff-utilities>=0.1.12" openff-forcefields=2026.01.0 openmm=8.1.1 openff-nagl=0.5.4 openff-nagl-models=2025.9.0 lomap2 kartograf=1.0.1 michellab::fkcombu konnektor -c conda-forge --yes && micromamba activate qligfep_new && python -m pip install -e .
```
</details>

### macOS

As on Linux, [clone the repository](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository), enter the `Q` directory with `cd Q/`, create the environment, and install QligFEP:

``` bash
micromamba create -n qligfep_new python=3.11 gfortran=11.3.0 openff-toolkit=0.17.1 "openff-utilities>=0.1.12" openff-forcefields=2026.01.0 openmm=8.1.1 openff-nagl=0.5.4 openff-nagl-models=2025.9.0 lomap2 kartograf=1.0.1 davidararipe::kcombu_bss konnektor -c conda-forge --yes
micromamba activate qligfep_new
python -m pip install joblib scipy tqdm
python -m pip install -e .
```

The `qprep` Fortran binary will be automatically compiled during installation.

<details>
<summary>To install everything in one line...</summary>

```bash
micromamba create -n qligfep_new python=3.11 gfortran=11.3.0 openff-toolkit=0.17.1 "openff-utilities>=0.1.12" openff-forcefields=2026.01.0 openmm=8.1.1 openff-nagl=0.5.4 openff-nagl-models=2025.9.0 lomap2 kartograf=1.0.1 davidararipe::kcombu_bss konnektor -c conda-forge --yes && micromamba activate qligfep_new && python -m pip install joblib scipy tqdm && python -m pip install -e .
```
</details>

### PyMOL for QresFEP mutant preparation

`setup_resFEP` runs PyMOL's mutagenesis wizard to place each mutant side chain
on the wild-type backbone. Install PyMOL in the project environment:

```bash
micromamba install -n qligfep_new -c conda-forge pymol-open-source
```

`setup_resFEP` calls the `pymol` executable directly, so you do not need a
separate Python environment or the PyMOL Python package. If you already have
mutant residue PDBs, pass their directory with `setup_resFEP --mutant-pdbs` and
omit this dependency.

### Compiling Q for HPC (MPI support)

> [!IMPORTANT]
> The QligFEP and QresFEP cluster workflows currently generate Slurm submission
> scripts. Installation automatically compiles `qprep`, which is sufficient for
> preparing topologies and workflow inputs. To run the simulations, compile the
> serial Q tools and the MPI-enabled `qdynp` engine as described below.

On your HPC system, load the appropriate modules (system-dependent). We recommend using the GCC compiler suite and OpenMPI, as those are commonly available and compatible with `qdynp`. To check for module availability, use the command `module spider openmpi` or `module avail openmpi`.

In the output, look for a version compiled with GCC (e.g., `OpenMPI/4.1.4-GCC-11.3.0`) and load it using the `module load` command. After loading the module, navigate to the `src/q6` folder in the Q repository and compile both the serial and MPI versions of Q with the commands `make all` and `make mpi`. In the example, we show how to do this on the Snellius, the Dutch national supercomputer:

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

Cluster profiles are defined in
[`src/QligFEP/settings/settings.py`](/src/QligFEP/settings/settings.py). The file
contains profiles used for QligFEP and QresFEP calculations on several HPC
systems. To add another system, follow the existing format. This example defines
a profile named `MY_HPC`:

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

Pass the profile name, such as `MY_HPC`, to a QligFEP or QresFEP setup command
with `--cluster`.

## ⌨️ Command line interface (CLI)

### Q engine

Building Q provides these Fortran executables:

- `qprep` prepares and solvates molecular topologies.
- `qdyn` runs serial molecular dynamics calculations.
- `qdynp` runs MPI-parallel molecular dynamics calculations.
- `qfep` analyzes FEP energy files.
- `qcalc` analyzes Q trajectories and molecular properties.

### Structure and parameter preparation

Installing the Python package provides these shared preparation tools:

- `qcog` calculates the center of geometry of each structure in a PDB or SDF
  file.
- `pdb2amber` converts a PDB file for use with Q's AMBER force field.
- `qprep_prot` solvates a protein and prepares its spherical Q topology. It
  preserves crystallographic waters by default; pass `--strip-crystal-waters`
  to remove them. The command also records the sphere charge and the
  input-PDB-to-Q residue mapping in `prep.json` for downstream workflows such as
  QresFEP.
- `qparams` generates ligand parameters.

### QligFEP

- `qlomap` uses LoMap to generate a ligand perturbation mapping in JSON format.
- `qkonnektor` uses Konnektor to generate a ligand perturbation network.
- `qmapfep` generates and displays a perturbation mapping so that you can add or
  remove edges interactively.
- `qligfep` generates the input files for one ligand perturbation.
- `setupFEP` prepares a complete ligand FEP series, including its protein and
  water systems. Pass `--neq` to select the non-equilibrium NEQ² workflow.
- `qligfep_analyze` analyzes windowed ligand FEP results.
- `qligfep_neq_analyze` analyzes NEQ² switching-work results.

### QresFEP

- `qresfep` generates one protein or reference-peptide leg for an amino-acid
  mutation.
- `setup_resFEP` validates and prepares a mutation series, including one
  mutation-centered sphere and both thermodynamic-cycle legs per mutation.
- `qresfep_analyze` combines the protein and reference-peptide legs into
  mutation-induced folding free energy changes.

Run any Python command with `--help` to see its required inputs and available
protocol settings.

## Tutorials

- The [Tyk2 tutorial](/tutorials/Tyk2/README.md) covers the standard windowed
  ligand FEP workflow.
- The [Tyk2 NEQ² tutorial](/tutorials/Tyk2/neq2/README.md) covers the
  non-equilibrium workflow.
- The [T4 lysozyme QresFEP tutorial](/tutorials/T4L/README.md) covers
  mutation-induced changes in protein folding stability.

### Protein thermostability with QresFEP

The QresFEP implementation in this repository currently provides the protein
thermostability protocol. It calculates a mutation-induced change in folding
free energy from two simulations: the mutation in the folded protein and the
same mutation in a capped reference peptide. It subtracts the reference-peptide
free energy from the protein free energy:

```text
ddG_fold = dG_protein - dG_tripeptide
```

The integrated QresFEP-2 workflow provides the following features:

- Packaged OPLS-AA/M residue libraries and parameters.
- A hybrid topology with a shared backbone and separate wild-type and mutant
  side chains.
- Two consecutive FEP stages, with stage 2 starting from the final coordinates
  of stage 1.
- Mutation-centered spherical systems and a `prep.json` manifest that preserves
  the mapping between the input PDB and Q topology.
- AXA, GXG, single-residue, and native-neighbor reference peptide models.
- Charge-matched reference spheres for charge-changing mutations.
- Explicit controls for lambda windows, production steps, equilibration length,
  Q's `separate_scaling` option, trajectories, cleanup, and random seeds.
- A `--manuscript-settings` preset for the fixed protocol and seed vector used
  by the in-preparation charge-changing campaign.
- Validation of qdyn and qfep output before cleanup, with the standard EXPRESS
  LOG footer in each SLURM output file.

For a mutation list, `setup_resFEP` can generate mutant coordinates with PyMOL,
prepare one sphere per mutation, and write both thermodynamic-cycle legs:

```bash
setup_resFEP -i protein.pdb -M mutations.txt -mc A -FF OPLSAAM -r 25 \
    -c SNELLIUS -w 25 -s exponential -T 298 -R 10
```

After the cluster jobs finish, combine both legs with:

```bash
qresfep_analyze -p protein -t tripeptide -T 298
```

See the [T4 lysozyme tutorial](/tutorials/T4L/README.md) for mutant preparation,
single-mutation setup, batch setup, submission, and analysis.

<a id="non-equilibrium-fep-neq2"></a>

### Non-equilibrium FEP (NEQ²)

Alongside the standard windowed protocol, QligFEP supports the
**non-equilibrium NEQ²** workflow. NEQ² changes λ continuously from one end state
to the other over many short, independent switching trajectories. It recovers
ΔΔG with the Bennett Acceptance Ratio over the forward and reverse work
distributions. The independent switching trajectories can run concurrently on
a cluster.

To set up a non-equilibrium calculation, pass `--neq` to `setupFEP`. Analyze the
switching work with `qligfep_neq_analyze`. The serial `qdyn` engine selects
non-equilibrium mode when its input contains a `[lambda_scaling]` section. The
`make all` command builds this engine with the other Q binaries in `src/q6`.
See the [Tyk2 NEQ² tutorial](/tutorials/Tyk2/neq2/README.md) for the complete
workflow.

# 📊 Benchmarking

Use the [QligFEP benchmarking repository](https://github.com/qusers/qligfepv2-BenchmarkExperiments)
to reproduce the published performance results.

# 📚 Citations

To cite the latest version of QligFEP, cite:
```bibtex
@article{araripe2026qligfepv2,
  author  = {Alencar Araripe, David and Díaz-Holguín, Alejandro and Poso, Antti and van Westen, Gerard J. P. and Åqvist, Johan and Gutiérrez-de-Terán, Hugo and Jespers, Willem},
  title   = {Doing More with Less: Accurate and Scalable Ligand Free Energy Calculations by Focusing on the Binding Site},
  journal = {Journal of Chemical Information and Modeling},
  year    = {2026},
  volume  = {66},
  number  = {6},
  pages   = {3164--3172},
  doi     = {10.1021/acs.jcim.5c02932},
  url     = {https://doi.org/10.1021/acs.jcim.5c02932},
}
```
**Other relevant references:**

- [Q](https://doi.org/10.1016/S1093-3263(98)80006-5)
- [QligFEP](https://doi.org/10.1186/s13321-019-0348-5)
- [QresFEP](https://doi.org/10.1021/acs.jctc.9b00538)
- [QresFEP-2](https://doi.org/10.1038/s42004-025-01771-0)

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
