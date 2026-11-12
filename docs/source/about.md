# About Q and QligFEP

## Q

Q is a molecular-dynamics package specialized for free-energy calculations. Its
current Fortran engine consists of topology preparation (`qprep`), serial and MPI
dynamics (`qdyn`, `qdynp`), free-energy analysis (`qfep`), trajectory analysis
(`qcalc`), and an input initialization checker (`qdum`).

Q supports local spherical systems as well as periodic systems. The spherical
boundary model is especially useful when a binding-site-centered region is the
scientific object of interest: atoms outside the sphere are excluded from
nonbonded evaluation and restrained, boundary solvent is controlled by radial
and polarization restraints, and a local reaction field can represent longer
range electrostatics.

## QligFEP

QligFEP automates relative ligand free-energy calculations with Q. It provides
commands to:

- normalize protein atom names and prepare spherical systems;
- parameterize ligands with OpenFF and NAGL or AM1-BCC charges;
- build perturbation networks with LoMap or Konnektor;
- generate dual-topology FEP definitions and mapped-atom restraints;
- create equilibrium-window or non-equilibrium switching protocols;
- write Slurm submission scripts for configured clusters;
- combine water- and protein-leg results and report uncertainties.

## Choosing the right documentation

- Start with the [quickstart](quickstart.md) for a first workflow.
- Use the [QligFEP workflow](qligfep/workflow.md) for protocol decisions.
- Use the [CLI reference](cli/index.md) for exact flags and output names.
- Use the [Q engine guide](q-engine/index.md) to run the Fortran programs
  directly.
- Use the [maintained manual](q-engine/manual.md) for the consolidated reference.

## Support and development

Source code and issue tracking are hosted at
<https://github.com/qusers/Q>. When reporting a scientific or runtime problem,
include the commit/version, generated `fep_config.json`, compiler/MPI versions,
input files, complete logs, and the smallest reproducible edge.
