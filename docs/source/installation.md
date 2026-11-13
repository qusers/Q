# Installation

## Supported environment

The current package metadata requires Python 3.13. QligFEP is developed and
tested primarily on Linux/HPC systems and macOS. The Python workflow and the
Fortran Q engine have separate build requirements:

- **Python preparation tools:** Python 3.13 and the scientific/cheminformatics
  dependencies below;
- **serial Q programs:** a Fortran compiler, normally `gfortran`;
- **MPI dynamics:** an MPI installation with `mpif90` built against a compatible
  compiler;
- **cluster workflow:** Slurm and a configured cluster profile.

## Recommended environment

[Micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html)
is recommended because several molecular-modeling dependencies are distributed
through Conda channels.

### Linux

```bash
micromamba create -n qligfep python=3.13 \
  gfortran=11.3.0 \
  openff-toolkit=0.17.1 "openff-utilities>=0.1.12" \
  openff-forcefields=2026.01.0 openmm=8.2.0 \
  openff-nagl=0.5.4 openff-nagl-models=2025.9.0 \
  lomap2 kartograf=1.0.1 michellab::fkcombu konnektor \
  -c conda-forge --yes
micromamba activate qligfep
python -m pip install -e .
```

### macOS

```bash
micromamba create -n qligfep python=3.13 \
  gfortran=11.3.0 \
  openff-toolkit=0.17.1 "openff-utilities>=0.1.12" \
  openff-forcefields=2026.01.0 openmm=8.2.0 \
  openff-nagl=0.5.4 openff-nagl-models=2025.9.0 \
  lomap2 kartograf=1.0.1 davidararipe::kcombu_bss konnektor \
  -c conda-forge --yes
micromamba activate qligfep
python -m pip install -e .
```

The channel-qualified FKComBU package differs between the documented Linux and
macOS environments.

:::{note}
The dependencies declared directly in `pyproject.toml` do not include every
network-generation and ligand-parameterization package. Use the complete
environment above if you need `qparams`, `qlomap`, and `qkonnektor`.
:::

## What the Python installation builds

```bash
python -m pip install -e .
```

installs the QligFEP console scripts and attempts to compile `qprep`. If no
Fortran compiler is available, the installation can finish with a warning while
leaving preparation commands unusable. Check explicitly:

```bash
command -v qparams
command -v setupFEP
ls -l src/q6/bin/q6/qprep
```

## Compile the Q engine

### Serial build

The serial build includes `qprep`, `qdyn`, `qdum`, `qfep`, and `qcalc`:

```bash
cd src/q6
make all COMP=gcc
```

Binaries are moved to `src/q6/bin/q6/`.

### MPI build

On a cluster, first load a GCC-compatible MPI module. Module names vary by site;
ask the system administrator or inspect them with `module spider openmpi`.

```bash
module purge
module load <site-specific-GCC/OpenMPI-module>
cd src/q6
make all COMP=gcc
make mpi COMP=gcc
```

`make mpi` writes `src/q6/bin/q6/qdynp`. Use the same module environment at run
time as at link time.

### Clean rebuild

```bash
cd src/q6
make nuke
make all COMP=gcc
make mpi COMP=gcc
```

`make nuke` removes compiled objects and the binary directory. Do not run it in a
directory containing unarchived local build products.

## Verify the installation

Check the Python command surface:

```bash
for command in qcog pdb2amber qprep_prot qparams qlomap qkonnektor \
  qmapfep qligfep setupFEP qligfep_analyze qligfep_neq_analyze; do
  command -v "$command" || echo "missing: $command"
done
```

Check compiled programs:

```bash
ls -l src/q6/bin/q6/{qprep,qdyn,qdynp,qdum,qfep,qcalc}
```

Run source tests from the repository root:

```bash
pytest -m "not slow"
```

Build-level Q tests can be compiled with:

```bash
make -C src/q6 test COMP=gcc
```

A command being present is not a complete scientific validation. Before a
production network, run a small known preparation, serial MD job, MPI MD job,
and qfep analysis on the target machine.

## Common installation failures

### `src/q6/bin/q6` does not exist

Install `gfortran`, then compile at least qprep:

```bash
cd src/q6
make qprep-only COMP=gcc
```

For simulations, use `make all` and, where needed, `make mpi`.

### `mpif90` is not found

Load the cluster's MPI development module. A runtime-only MPI package is not
enough; the Fortran wrapper is required.

### Python import is successful but a CLI dependency is missing

Install the complete Micromamba environment. In particular, ligand
parameterization and network generation require packages beyond the base
scientific Python list.

### A binary fails with a shared-library error

Reload the compiler/MPI modules used during compilation or rebuild inside the
module environment that will be used by Slurm jobs.
