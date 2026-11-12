# Q manual

This is the Q manual, ported from the *Manual for Q6* (17 August 2017) and
split into pages. The text follows the original; program names, keywords and
defaults have been brought in line with the source in this repository, and
every place where the two differ is marked.

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

- The **user guide** walks through the normal order of work: build a topology,
  run dynamics, define a perturbation, analyse.
- The **reference** pages hold the file formats, command lists and the
  equations behind the boundary conditions and algorithms.
- [Differences from the Q6 manual](differences.md) lists what the 2017 manual
  describes that this source does not contain, and what this source adds.

Citations in square brackets, such as [1], refer to the
[bibliography](bibliography.md).

## Program names

The 2017 manual calls the programs `Qprep6`, `Qdyn6` and so on. This
repository builds them under shorter names, which this manual uses throughout.

| This repository | Q6 manual | Purpose |
|---|---|---|
| `qprep` | `Qprep6` | Build and inspect topologies; convert coordinates |
| `qdyn` | `Qdyn6` | Molecular dynamics, FEP and EVB sampling |
| `qdynp` | `Qdyn6p` | MPI version of `qdyn` |
| `qdum` | `Qdum6` | Read and check a `qdyn` input without running dynamics |
| `qfep` | `Qfep6` | Free-energy analysis of energy files |
| `qcalc` | `Qcalc6` | Structure, trajectory and scoring analysis |
