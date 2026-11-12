# Preparing a molecular topology

The topology file made by `qprep` contains all the information about the
molecular system that `qdyn` needs. To make a topology you need:

- a **fragment library** file describing the atoms and connectivity of each
  protein residue, ligand and so on (the *library entries*);
- a **force-field parameter** file;
- **coordinates** for all non-solvent atoms, as a PDB file.

## Preparing coordinates

The PDB format [16] read by `qprep` is a fixed format: the number of spaces
between columns is significant and tab characters are not permitted. Only
`ATOM` and `HETATM` records are read. See
[Coordinate files](reference/qprep.md#coordinate-files) for the full list of
rules.

**Hydrogens** are not required. `qprep` generates their coordinates so that the
bond and angles to the hydrogen are at an energy minimum. To control the
placement further, a torsion can be defined in the library entry; see
[build rules](reference/library-file.md#build-rules). Hydrogen coordinates
that are present in the PDB file are not altered.

**Molecules.** `qprep` identifies the boundaries between molecules from the
type of fragment or from gap marker lines. Monomers of a polymer chain, such
as amino acids, have designated *head* and *tail* atoms in their library
entries, which define how they connect to their neighbours. Entries for
separate molecules such as solvent or ligands have no such linkage. To separate
two peptide chains, put a line containing only `GAP` after the last atom of
the first chain. A `TER` line has the same effect.

**Numbering.** Atom numbers in the PDB file are not significant. `qprep`
renumbers all atoms in one sequence starting at 1, which is necessary to
insert the hydrogens. Atom numbers in the PDB file read by `qprep` and in PDB
files written by `qprep` therefore differ. Residues are renumbered as well.
Residue numbers must be numeric; identifiers such as `60A` are not permitted.

**Crystallographic water.** Well-ordered waters that matter for structure or
function should be kept; most others can be removed. A large number of waters
on the protein surface disturbs the solvation algorithm and leads to an
inhomogeneous water density. `qprep` assumes that solvent molecules come after
solute molecules in the PDB file and keeps track of the last solute atom.
Solvent molecules are identified by residue name. Solvent added by the
solvation algorithm is placed at the end of the atom sequence. Solvent
molecules outside the simulation sphere (more than the water radius + 2 Å from
the solvent centre) are excluded, like solute atoms outside the sphere.

## Selecting a force field

A Q parameter file holds the interaction parameters and also overall
properties of the force field, such as the van der Waals combination rule.
Selecting a force field is therefore a matter of loading the matching fragment
library and parameter file into `qprep`. See
[Force fields](reference/force-fields.md) for the files in this repository.

## Adding library entries

Your simulations will probably include fragments that are not in the standard
library, and you need a library entry for each new residue, ligand, cofactor
or other molecule. Keep your entries in a separate file rather than adding
them to the standard library. Load the standard library first and then your
own. To modify a standard entry, add the modified entry under the same name to
your own file; it replaces the old one when loaded.

The best starting point for a new entry is a copy of an old one. See
[Fragment library file](reference/library-file.md).

## Running qprep

`qprep` is a command-oriented program designed to be run interactively. The
`help` command lists the available commands. The input parser prompts for the
parameters of each command, but also accepts them on the command line: you
may enter `readlib` and wait for the prompt, or type
`readlib my_ligands.lib` on one line.

`qprep` handles errors without stopping. When it finds a problem in a library,
parameter or PDB file it reports it and continues, so later errors in the file
are reported as well. After correcting the file, read it again; there is no
need to restart the program.

The steps below are in the order they are normally executed. See
[qprep reference](reference/qprep.md) for all commands.

### Load fragment libraries

Use `readlib`. Repeat it to load several libraries. If an entry name occurs
more than once, the latest definition takes precedence, and a warning is
issued. `clearlib` removes all library entries from memory.

### Load force-field parameters

Use `readprm`. Only one parameter file can be loaded. If you need to modify
it, edit and save it and load it again; it is not necessary to reload the
libraries or the structure.

### Load coordinates

Use `readpdb`. Before loading, the program verifies that all required library
entries are loaded and that the number of heavy atoms in each fragment is
correct. Only one file can be loaded; reading another clears the previous
structure.

### Choose the boundary condition

Use `boundary` to choose a sphere (`sphere` or `1`) or a box (`box` or `2`).
For a sphere you give the centre and radius; for a box, the centre and the
three side lengths. The centre is given as `x y z`, as
`residue_number:atom_name`, or as `mass` for the centre of mass. The boundary condition in the
topology must be consistent with the one used in `qdyn`.

### Add solvent

Use `solvate`. With a spherical boundary you give the centre of the solvent
sphere (as for `boundary`, or `boundary` to reuse the centre of the
simulation sphere), its radius, and one of these modes:

- **Grid** (`grid` or `1`). Randomly oriented molecules on a grid. This needs no solvent
  coordinate file; you only give the library entry to use. The density is
  taken from the library entry.
- **File** (`file` or `2`). Coordinates from a solvent file, which is similar to a PDB file;
  see [Solvent file](reference/force-fields.md#solvent-file-format). The file
  contains a sphere or a box of solvent. A box is replicated in all directions
  as needed, so it can solvate a system of any size.
- **Restart** (`restart` or `3`). Solvent coordinates from a restart file of a previous
  simulation of the same system.

With the first two methods `qprep` fills the sphere with solvent and then
deletes molecules whose heavy atoms are closer than a threshold to any heavy
atom of the loaded structure. The threshold is the preference `solvent_pack`.
Crystallographic solvent molecules more than 2 Å outside the solvent sphere
are excluded from the simulation, to avoid excessive radial restraining
forces.

With a periodic boundary, a solvent file containing a sphere cannot be used.
The threshold distance is then also applied to solvent molecules near the box
boundary, to avoid clashes between solvent molecules in neighbouring boxes.

The water model is defined by the library entry: TIP3P [18] in the libraries
of this repository (`HOH`).

:::{note}
Crystallographic waters and waters added by `solvate` must have the same
residue name. Otherwise `qprep` marks the topology as mixed-solvent, which
`qdyn` does not support.
:::

### Add cross-link bonds

Cross-links such as disulphide bridges are not in the fragment library and
are not generated automatically. To add them:

- `xlink` searches the solute for close, non-bonded atom pairs and asks you to
  confirm each bond;
- `addbond` adds a bond between two atoms, given as atom numbers or as
  `residue_number:atom_name`.

### Generate the topology

`maketop` creates the topology in memory. Its only input is a name. It:

- sets atom types and partial charges;
- generates the lists of bonds, angles, torsions and impropers;
- generates the neighbour-exclusion and 1–4 neighbour lists;
- generates coordinates for hydrogens that are not in the loaded structure;
- marks atoms outside the simulation sphere as excluded;
- calculates the *effective solvent radius* used to keep the solvent density
  correct. This radius is based on the number and distribution of solute and
  solvent atoms in the sphere and typically differs somewhat from the radius
  used in the solvation step.

### Verify the topology

A topology that was generated without errors is not necessarily correct. Check
it with:

- `checkbonds`: lists bonds with an energy above a threshold, which helps to
  find errors in the connectivity;
- `checkangs`: the same for angles;
- `checktors`: the same for torsions;
- `checkimps`: the same for impropers.

`checkimps` is important for force fields with harmonic impropers. With a non-periodic potential it makes a large difference whether
the angle is −179° or +179° if the minimum is at +180°, although the
structures are nearly the same. Impropers with the wrong sign give high
energies and strong forces that distort the molecule. The sign depends on the
order of the bonds in the library entry. Instead of editing library entries,
`changeimp` flips the sign of selected impropers, or of all impropers with an
energy above a threshold. This is not a problem with periodic impropers.

### Write topology and coordinate files

`writetop` saves the topology.

`writepdb` writes a PDB file of the atoms in the current mask (all atoms by
default). You need it to see the new atom numbers, which are the numbers to
use in restraints and in FEP files. `writemol2` writes a mol2 file.

### Preferences

Parameters that affect `qprep` but are not normally changed, such as those of
the solvation step, are held as preferences. `prefs` lists them and `set`
changes one. See [Preferences](reference/qprep.md#preferences).

## Scripted use

`qprep` also runs a command file given as its argument, which replaces the
interactive session:

```text
readlib ligand.lib
readlib protein.lib
readprm combined.prm
readpdb complex.pdb
boundary 1 -4.689 26.119 -30.570 25.0
solvate -4.689 26.119 -30.570 25.0 1 HOH
maketop complex
writetop complex.top
checkbonds 6
writepdb complex_prepared.pdb y
quit
```

```bash
qprep prep.inp > prep.log
```

Here `boundary 1` selects a sphere with the given centre and radius, and
`solvate ... 1 HOH` fills it from a grid with the library entry `HOH`. Read the
whole log; `qprep` continues after errors.
