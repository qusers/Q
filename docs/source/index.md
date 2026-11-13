# Q and QligFEP

**Molecular dynamics and automated free-energy calculations in biomolecular
systems.**

Q is a Fortran molecular-dynamics suite built for free-energy perturbation (FEP),
empirical valence bond (EVB), and linear interaction energy (LIE) calculations.
QligFEP is the Python command-line workflow that automates ligand
parameterization, perturbation-network construction, spherical system setup,
Q input generation, cluster scripts, and free-energy analysis.

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} New user
:link: quickstart
:link-type: doc
Install the software and prepare a small relative binding free-energy workflow.
:::

:::{grid-item-card} QligFEP user
:link: qligfep/workflow
:link-type: doc
Follow the complete preparation, simulation, and analysis workflow.
:::

:::{grid-item-card} Q engine user
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

:::{important}
The documentation follows the source in this repository. The Q manual is a
port of the 2017 2017 Q manual and marks where this source differs from it.
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
:maxdepth: 2
:caption: CLI reference

cli/index
cli/preparation
cli/networks
cli/setup
cli/analysis
```

```{toctree}
:maxdepth: 2
:caption: Tutorials and help

tutorials/tyk2
troubleshooting
citations
contributing-docs
```
