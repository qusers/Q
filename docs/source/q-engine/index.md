# Q molecular-dynamics engine

Q is a modular Fortran suite for molecular dynamics and free-energy calculations.
Its central design is the **Q-atom** model: a chosen subset of atoms can have a
different force-field description in each thermodynamic state, while the rest of
the environment retains one topology.

## Executables

| Executable | Role | Input style |
|---|---|---|
| `qprep` | Construct and inspect molecular topologies | Interactive or positional command file |
| `qdyn` | Serial MD, FEP/EVB, minimization, and lambda switching | Positional `.inp` file |
| `qdynp` | MPI fixed-lambda MD/FEP | Positional `.inp` file under `mpirun`/`srun` |
| `qdum` | Initialize a Qdyn input without the MD loop | Positional `.inp` file |
| `qfep` | Free-energy estimators and EVB profiles | Standard input |
| `qcalc` | Structural, geometric, interaction, and scoring analysis | Standard input |

## Build targets

From `src/q6`:

```bash
make all COMP=gcc  # qprep, qdyn, qdum, qfep, qcalc
make mpi COMP=gcc  # qdynp
make test COMP=gcc # test-flavored serial suite
```

Executables are moved to `src/q6/bin/q6/`.

## Model and force field

The potential combines bonded terms (bonds, angles, proper and improper
torsions), Coulomb interactions, and 6–12 Lennard-Jones interactions. Q reads
molecular fragments from `.lib` files and parameters from `.prm` files rather
than hard-coding one force field. The current repository supplies AMBER14sb,
CHARMM36, OPLS2005, and OPLS2015 protein files.

For free-energy calculations, Qdyn evaluates the Q-region energy under each state
and propagates coordinates using lambda-weighted state forces. State data are
written to binary `.en` files for qfep.

## Boundary models

### Spherical boundary

The simulation focuses on a sphere around a chosen center. Solute outside the
active region is restrained/excluded, solvent is kept inside by radial boundary
forces, and boundary polarization plus local reaction-field terms represent the
truncated environment. This model underpins the standard QligFEP workflow.

### Periodic boundary

A `[PBC]` section enables periodic simulation and optional Monte Carlo volume
moves. The nonbonded cutoff must remain below half the shortest box length.

## MD algorithms in the current source

The current engine uses the established leapfrog-style update, weak temperature
coupling, and SHAKE bond constraints. It can run an optional two-phase
steepest-descent minimization before MD. The free-energy engine offers standard
LJ softcore, softened LJ/Coulomb effective-distance treatment, and Gapsys-style
short-range force linearization.

For exact input and format tables, continue to [Running Q](running.md),
[Free-energy analysis](fep-analysis.md), and [File formats](file-formats.md), or
open the [complete maintained manual](manual.md).
