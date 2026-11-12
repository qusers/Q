# Tutorial: binding affinity from LIE simulations

This tutorial describes a linear interaction energy (LIE) calculation of the
binding of stearate to a muscle fatty-acid binding protein (M-FABP). The
published example used the Q version of the GROMOS87 force field, which is
not among the force fields of this repository; the procedure is the same
with any force field, but the coefficients $\alpha$ and $\beta$ at the end
are specific to GROMOS87.

## Editing the PDB file

Download the structure of the M-FABP complex with stearate, PDB entry 1HMT.
The file needs editing before use.

- **Remove crystal waters.** In 1HMT, 17 waters were kept because they take
  part in binding the ligand or in other interactions. To decide which to
  keep, pick an atom of the ligand as the centre of the system and choose the
  size of the simulation sphere; here the radius is 18.0 Å. Keep the waters
  inside the sphere that appear to take part in interactions and lie inside
  the protein.
- **Keep only coordinate lines** of the protein, the ligand and the chosen
  waters. Remove blank lines and `TER` lines.
- **Insert a `GAP` line** between different molecules.
- **Rename cysteines** that form disulphide bridges to `CSS`.

Hydrogens need not be present; `qprep` adds them.

## Modelling ionic groups of the protein

In a spherical simulation, ionisable residues far from the centre are
modelled in a neutral form. In the Q-GROMOS87 library used for this example,
the default entries for ASP, GLU, ARG and LYS have the protonation state of
the ionic form but the net charge replaced by a dipole, and the charged side
chains are separate entries (`AS-`, `GL-`, `AR+`, `LY+`). Changing a residue
to its charged entry is referred to below as *turning on* its charge.

Decide which charged residues to turn on with the same rules as for the
waters. Residues near the ligand, forming a salt bridge or otherwise involved
in the function of the protein, and inside the 18.0 Å sphere, should be
charged. Residues closer than about 3–5 Å to the boundary should be neutral
unless they form an ion pair with a more central group. When the ligand has
an ionic group, as here, it is important to make the protein net neutral. In
this example residues 17, 72, 76, 77, 78, 106 and 126 were charged.

## Writing the library file

Write a library entry for the ligand, most easily by editing an existing one;
the amino-acid library is a good guide. The entry lists all atoms of the
ligand, the bonds and the charge groups. Each atom needs a name, a type and a
charge. The atom types are listed in the parameter file. Charges can
sometimes be taken by comparison with an amino acid. Name the atoms as in the
PDB file. Because the ligand is not connected to any other fragment, the head
and tail connections can be omitted.

You also need a PDB file of the ligand alone: copy the relevant lines from
the PDB file of the complex.

## qprep

Make two topologies: one for the ligand in the protein and one for the ligand
in water.

1. Read the libraries and the PDB file (`readlib`, `readpdb`). For the
   complex, read both the amino-acid library and the ligand library.
2. Create the disulphide bridges with `addbond atom atom`. Find the atom
   numbers of the sulphurs with `listres residue`. The numbers in the PDB file
   cannot be used, because atoms are renumbered when hydrogens are added.
3. Select the boundary and solvate.
4. Create the topology with `maketop`. `qprep` reports missing parameters
   here. Add them to the parameter file, with a comment after `!` so the
   change is visible, save the file and run `maketop` again.
5. Write the topology with `writetop`.

Check for high-energy bonds, angles, torsions and impropers with `checkbonds`,
`checkangs`, `checktors` and `checkimps`, each followed by an energy level. A
very high bond, angle or torsion energy may mean a wrongly connected
disulphide or a missing `GAP`. An improper with a very high energy may have
the wrong sign (180° instead of −180°); `changeimp` redefines these
automatically. Write the topology again after `changeimp`. Write a new PDB
file with `writepdb` to get the renumbered atoms with hydrogens.

## FEP files

In this example the FEP file only names the Q-atoms, that is, the ligand. In
the simulation without the protein these are all the atoms. In the
protein–ligand simulation, take the atom numbers of the ligand from the new
PDB file.

## Input files

The input file sets the number and length of steps, the temperature, the
topology and so on. Here the data collection was split into five identical,
consecutive runs to make it easier to restart after an interruption. That
gives five input files for each of the two simulations, and six more for the
equilibration of the complex.

Specific to each simulation are the centre of the simulation and water
spheres, the topology and FEP files, and the restraints. The simulation and
water spheres share the same centre: the ligand atom chosen earlier.

The equilibration warms the system from 0 K to 300 K in steps, with all heavy
atoms restrained, including those of the ligand. After equilibration,
5 × 50 000 steps of 1 fs are run for each simulation. All runs except the
first restart from the coordinates and velocities of the previous one. When a
new temperature is given, a random seed is needed as well; when the
temperature is unchanged, set the seed to zero.

## qdyn

Start a simulation with

```bash
qdyn filename.inp > filename.log
```

With many input files, put the commands in a script.

## Evaluating the simulation

Extract the ligand–surrounding van der Waals and electrostatic energies from
the log files and plot them against time, to see whether they drift or jump
during the simulation.

Look at the structures after each part of the simulation. Convert restart
files to PDB with `qprep`: read the topology, read the restart file with
`readx` and write a PDB file.

Average the energies, discarding the start of the run, and compare the
averages of the two halves of the remaining data as a check on convergence.

Then calculate the electrostatic free-energy contribution of the ligand's
interaction with ionic groups of the protein that were not turned on during
the simulation, $\Delta G_{\mathrm{onoff}}$.

The binding free energy follows from the LIE formula:

$$
\Delta G_{\mathrm{bind}} =
\alpha \left( \left\langle V^{\mathrm{LJ}}_{l-s} \right\rangle_{\mathrm{bound}} - \left\langle V^{\mathrm{LJ}}_{l-s} \right\rangle_{\mathrm{free}} \right)
+ \beta \left( \left\langle V^{\mathrm{el}}_{l-s} \right\rangle_{\mathrm{bound}} - \left\langle V^{\mathrm{el}}_{l-s} \right\rangle_{\mathrm{free}} \right)
+ \Delta G_{\mathrm{onoff}}
$$

where $V_{l-s}$ is the interaction energy between the ligand and its
surroundings, in the protein (*bound*) and in water (*free*).

Here $\beta$ = 0.5 (no deviation from electrostatic linear response for a
charged ligand) and $\alpha$ = 0.181 (from an earlier calibration with
GROMOS87). The result is $\Delta G_{\mathrm{bind}}$ = −8.0 kcal/mol, close to
the experimental value of −7.9 kcal/mol.

$\beta$ depends on the number of OH groups of the ligand. With GROMOS87 and
ligands without ionic groups, choose $\beta$ from the set of values in [8]
according to the composition of the ligand.
