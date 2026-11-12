# Q manual

This is the manual of the Q programs in this repository. It is based on the
*Manual for Q6* (17 August 2017), split into pages and revised so that
commands, keywords, defaults and file formats are those of the source here.

```{toctree}
:maxdepth: 1
:caption: User guide

introduction
topology-preparation
running-qdyn
fep-files
analysis
scoring
tutorial-lie
```

```{toctree}
:maxdepth: 1
:caption: Reference

reference/force-fields
reference/qprep
reference/library-file
reference/parameter-file
reference/boundaries
reference/algorithms
reference/qdyn-input
reference/fep-file
reference/utilities
differences
bibliography
```

## How to read this manual

- [Concepts](../concepts/index.md) explains the spherical boundary, the
  electrostatics and the state-based FEP that the manual takes for granted.

- The **user guide** walks through the normal order of work: build a topology,
  run dynamics, define a perturbation, analyse.
- The **reference** pages hold the file formats, command lists and the
  equations behind the boundary conditions and algorithms.
- [Changes from the Q6 manual](differences.md) is for readers who know the
  2017 manual or have input files written for Q6.

Citations in square brackets, such as [1], refer to the
[bibliography](bibliography.md).

## Programs

| Program | Purpose | Invocation |
|---|---|---|
| `qprep` | Build and inspect topologies; convert coordinates | `qprep prep.inp > prep.log`, or interactive |
| `qdyn` | Molecular dynamics, FEP and EVB sampling | `qdyn md.inp > md.log` |
| `qdynp` | MPI version of `qdyn` | `mpirun -np 8 qdynp md.inp > md.log` |
| `qdum` | Read and check a `qdyn` input without running dynamics | `qdum md.inp > check.log` |
| `qfep` | Free-energy analysis of energy files | `qfep < qfep.inp > qfep.out` |
| `qcalc` | Structure, trajectory and scoring analysis | `qcalc < qcalc.inp > qcalc.out` |

Build them from `src/q6`:

```bash
make all COMP=gcc    # qprep, qdyn, qdum, qfep, qcalc
make mpi COMP=gcc    # qdynp
```

The executables are placed in `src/q6/bin/q6/`. See
[Installation](../installation.md).

## Files

| File | Format | Written by | Read by |
|---|---|---|---|
| `.pdb` | text, fixed columns | `qprep` | `qprep` |
| `.lib` | text, [fragment library](reference/library-file.md) | you, QligFEP | `qprep` |
| `.prm` | text, [parameter file](reference/parameter-file.md) | you, QligFEP | `qprep` |
| `.top` | topology | `qprep` | `qdyn`, `qprep`, `qcalc` |
| `.inp` | text, [qdyn input](reference/qdyn-input.md) | you, QligFEP | `qdyn` |
| `.fep` | text, [FEP file](reference/fep-file.md) | you, QligFEP | `qdyn` |
| `.re` | binary restart | `qdyn` | `qdyn`, `qprep`, `qcalc` |
| `.dcd` | binary trajectory | `qdyn` | `qprep`, `qcalc`, viewers |
| `.en` | binary state energies | `qdyn` | `qfep` |
