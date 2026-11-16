# About Q

The Q repository holds the Q molecular dynamics engine and two workflows built
on it, QligFEP and QresFEP.

## Q engine

Q is a molecular dynamics package specialised for free-energy calculations.
The engine consists of topology preparation (`qprep`), serial and MPI
dynamics (`qdyn`, `qdynp`), free-energy analysis (`qfep`), trajectory analysis
(`qcalc`), and an input checker (`qdum`).

Q supports spherical systems as well as periodic systems. The spherical
boundary model is especially useful when one region, such as a binding site or
a mutated residue, is the object of interest: atoms outside the sphere are
excluded from non-bonded evaluation and restrained, the solvent at the boundary
is controlled by radial and polarisation restraints, and a local reaction
field represents longer-range electrostatics. See [Concepts](concepts/index.md).

The engine described in this documentation is written in Fortran and runs on
CPUs. A GPU implementation, Q-GPU, is under development.

## QligFEP

QligFEP automates relative binding free-energy calculations of ligands. It
provides commands to:

- normalise protein atom names and prepare spherical systems;
- parameterise ligands with OpenFF and NAGL or AM1-BCC charges;
- build perturbation networks with LoMap or Konnektor;
- generate dual-topology FEP definitions and mapped-atom restraints;
- create equilibrium-window or non-equilibrium switching protocols;
- write Slurm submission scripts for configured clusters;
- combine water- and protein-leg results and report uncertainties.

## QresFEP

QresFEP automates free-energy calculations of amino-acid mutations. The
protocol available here is protein thermostability: the change in folding free
energy caused by a point mutation. It provides commands to:

- build the mutant side chain on the wild-type backbone;
- prepare a sphere centred on the mutated residue;
- generate the hybrid-topology residue and the two-stage FEP definition;
- set up the mutation in the folded protein and in a capped reference peptide;
- write Slurm submission scripts for configured clusters;
- combine both legs into the change in folding free energy.

The two workflows are distributed in one Python package, `QligFEP`. Installing
it provides the commands of both.

## Where to start

- [Installation](installation.md) covers the Python package and the engine.
- For ligands, read the [QligFEP workflow](qligfep/workflow.md), or follow the
  [quickstart](quickstart.md) and the [Tyk2 tutorial](tutorials/tyk2.md).
- For mutations, read the [QresFEP workflow](qresfep/workflow.md) and follow
  the [T4 lysozyme tutorial](tutorials/t4l.md).
- Use the [CLI reference](cli/index.md) for exact flags and output names.
- Use the [Q manual](manual/index.md) to run the engine programs directly.

## Support and development

Source code and issue tracking are hosted at
<https://github.com/qusers/Q>. When reporting a scientific or runtime problem,
include the commit or version, the generated `fep_config.json` or
`resfep_config.json`, compiler and MPI versions, input files, complete logs,
and the smallest case that reproduces the problem.
