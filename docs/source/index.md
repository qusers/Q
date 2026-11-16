# Q

**Molecular dynamics and automated free-energy calculations in biomolecular
systems.**

Q is a molecular dynamics package designed for free-energy calculations:
free-energy perturbation (FEP), empirical valence bond (EVB) and linear
interaction energy (LIE). Two Python workflows automate FEP calculations with
Q, from system preparation to analysis:

- **QligFEP** calculates relative binding free energies of ligands.
- **QresFEP** calculates the free-energy change caused by amino-acid
  mutations.

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} Ligand FEP (QligFEP)
:link: qligfep/workflow
:link-type: doc
Prepare, run and analyse relative binding free-energy calculations for a
ligand series.
:::

:::{grid-item-card} Residue FEP (QresFEP)
:link: qresfep/workflow
:link-type: doc
Prepare, run and analyse free-energy calculations of amino-acid mutations.
:::

:::{grid-item-card} Q engine
:link: manual/index
:link-type: doc
Read the Q manual: programs, input files, file formats and algorithms.
:::

:::{grid-item-card} Command reference
:link: cli/index
:link-type: doc
Find every installed Python command and option.
:::
::::

:::{note}
This documentation describes the Fortran implementation of the Q engine,
which runs on CPUs. A GPU implementation, Q-GPU, is under development and is
not yet integrated with the workflows described here.
:::

```{toctree}
:maxdepth: 2
:caption: Getting started

about
installation
quickstart
cluster-configuration
```

```{toctree}
:maxdepth: 2
:caption: Concepts

concepts/index
```

```{toctree}
:maxdepth: 2
:caption: Q manual

manual/index
```

```{toctree}
:maxdepth: 3
:caption: QligFEP

qligfep/workflow
qligfep/preparation
qligfep/equilibrium
qligfep/neq
qligfep/mapping-format
```

```{toctree}
:maxdepth: 3
:caption: QresFEP

qresfep/workflow
```

```{toctree}
:maxdepth: 2
:caption: CLI reference

cli/index
cli/preparation
cli/networks
cli/setup
cli/analysis
cli/qresfep
```

```{toctree}
:maxdepth: 2
:caption: Tutorials and help

tutorials/tyk2
tutorials/t4l
troubleshooting
citations
contributing-docs
```
