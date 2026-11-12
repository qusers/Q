# Utility programs

## qcalc

`qcalc` analyses a topology together with restart and trajectory files. It
reads its commands from standard input:

```bash
qcalc < qcalc.inp > qcalc.out
```

Run without redirected input, it is interactive: it asks for the topology,
shows a menu of calculations, asks for the parameters of each calculation
that is added, and then for the coordinate files to process. Record an
interactive session once and keep the answers as a script.

The calculations in this source are:

| Key | Calculation |
|---|---|
| `rmsd` | RMS coordinate deviation |
| `rmsf` | RMS coordinate fluctuation |
| `fit` | Least-squares fit of coordinates |
| `dist`, `angle`, `torsion`, `improper` | Geometry and its energy |
| `entropy` | Entropy from the covariance of coordinates |
| `nonbond` | Non-bonded interaction between two masks |
| `nb_prot_qatom` | Non-bonded interactions between residues and Q-atoms |
| `rdf` | Radial distribution function |
| `com`, `com_ke` | Centre of mass, and its kinetic energy |
| `xscore`, `chemscore`, `pmfscore` | [Scoring](../scoring.md) |

The 2017 manual lists only the three scoring functions for `qcalc`.

## qdum

| | |
|---|---|
| Function | Test `qdyn` input files quickly |
| Input | `qdyn` input file |
| Output | Log file, restart file |
| Usage | `qdum test.inp` |

`qdum` is `qdyn` without the dynamics loop. It reads all input, initialises
the simulation, writes a restart file and stops. It is built from the same
source as `qdyn`.

## External tools

The 2017 manual points to the graphical interface QGui [14]
(<https://github.com/qusers/qgui>) and the command-line tools QTools [15]
(<https://github.com/mpurg/qtools>) for set-up and analysis. In this
repository, QligFEP provides the set-up and analysis workflow for ligand
binding free energies; see [QligFEP workflow](../../qligfep/workflow.md).
