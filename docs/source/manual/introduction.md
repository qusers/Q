# Introduction

## Background

Molecular dynamics (MD) simulations can be used to sample the thermally
accessible regions of conformational space using a microscopic model of the
molecular system. From the ensemble of sampled structures and their associated
potential energies (given by the force field) it is, in principle, possible to
calculate free energies. Binding free energies, solvation free energies and
activation free energies are particularly interesting to compute because they
are the direct result of thermodynamic or kinetic experiments. It is thus
possible both to verify calculated results against experimental data and to
make predictions that can be tested experimentally.

Q [1] is a set of tools tailored to the calculation of free energies by:

1. free energy perturbation (FEP) simulations [2, 3];
2. empirical valence bond (EVB) calculations [4, 5] of reaction free energies;
3. linear interaction energy (LIE) calculations [6–8] of receptor–ligand
   binding affinities.

The features that distinguish Q from other MD packages are:

- **The spherical boundary.** Q is intended for free-energy calculations in
  biomolecular systems solvated in a spherical droplet of explicit water. A
  spherical boundary [11–13] limits the size of the simulated system, so the
  simulation can focus on a region such as a binding site, and makes accurate
  treatment of long-range electrostatics inexpensive.
- **The choice of force field.** Force fields are defined in parameter files,
  separate from the program.
- **Ease of use.** The input and force-field files are organised in sections
  of keywords and are easy to read and modify. The programs give extensive
  diagnostics when they find a problem.
- **Dynamic memory allocation.** Q can simulate systems of any size without
  recompilation.

## History

Q was started in the 1990s by Johan Åqvist. Contributors to the code include
John Marelius, Karin Kolmodin, Isabella Feierberg, Martin Almlöf, Martin Ander,
Jens Carlson, Peter Hanspers, Anders Kaplan, Kajsa Ljungberg, Martin Nervall,
Johan Sund, Paul Bauer, Alexandre Barrozo, Masoud Kazemi, Irek Szeler and
Åke Sandgren.

## Overview of the programs

The structure of Q resembles that of many other MD programs.

- `qdyn` carries out the trajectory calculation. Besides its input file it
  needs a **topology**: a file that lists how atoms are bonded, together with
  all force-field parameters.
- `qprep` creates the topology from a PDB coordinate file and force-field
  files. It is interactive, and is also used to convert restart and trajectory
  files back to PDB.
- The **FEP file** holds the input for FEP, EVB and receptor–ligand
  simulations. It lists the atoms to be transformed, called **Q-atoms**, and
  their force-field parameters in the different states.
- `qfep` calculates free energies from the energy files written by `qdyn`.
- `qcalc` analyses trajectories and restart files.

```{mermaid}
flowchart TD
    lib["Fragment library .lib"] --> qprep
    prm["Parameters .prm"] --> qprep
    pdb["Coordinates .pdb"] --> qprep
    qprep["qprep"] --> top["Topology .top"]
    top --> qdyn
    inp["Input .inp"] --> qdyn
    fep["FEP file .fep"] --> qdyn
    qdyn["qdyn / qdynp"] --> log["Log"]
    qdyn --> re["Restart .re"]
    qdyn --> dcd["Trajectory .dcd"]
    qdyn --> en["Energies .en"]
    en --> qfep["qfep"]
    re --> qprep2["qprep / qcalc"]
    dcd --> qprep2
```

## Installation

See [Installation](../installation.md) for how to build the programs in this
repository. The system requirements and binary-release instructions of the
2017 manual refer to the separate `qusers/Q6` repository and are not repeated
here.

## How to cite

- Marelius, J.; Kolmodin, K.; Feierberg, I.; Åqvist, J. Q: a molecular
  dynamics program for free energy calculations and empirical valence bond
  simulations in biomolecular systems. *J. Mol. Graph. Model.* **1998**, *16*,
  213–225.

See [Citations](../citations.md) for QligFEP and the methods it uses.

## Licence

The 2017 manual reproduces the GNU General Public License, version 2. The
`LICENSE` file of this repository is the GNU General Public License,
version 3, and is the one that applies to this code.

## Useful tips

- To run FEP simulations of a ligand in water and bound to a protein with the
  same FEP file, use the `offset_name` keyword in the `[FEP]` section instead
  of renumbering all the atoms.
- Make a separate library file for your new molecules and leave the amino-acid
  library unchanged. Load both in `qprep`.
- You can add parameters to the parameter file without restarting `qprep`.
  Type `maketop` again and the updated file is used.
- For FEP simulations with dummy atoms, `qprep` warnings about missing
  parameters can be ignored if all of those interactions are redefined in the
  FEP file. It is possible, but in general not advisable, to write a topology
  with missing parameters and use it in `qdyn`.
- Use build rules in your fragment library entries to control the positioning
  of hydrogens.
- Improve scoring accuracy by averaging over, for example, every tenth frame
  of a short equilibration trajectory.
