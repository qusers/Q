# qprep reference

## Coordinate files

Coordinates are read as PDB files. Where nothing is said below, the PDB
standard applies.

- Only `ATOM` and `HETATM` records are read. Other records are ignored, with a
  warning.
- Residue numbers must be numeric; identifiers like `60B` are not allowed. The
  numbering need not start at 1, and `qprep` renumbers residues from 1.
- Molecules must be separated by a gap marker: a line containing only the
  word `GAP` in capitals, optionally preceded by white space. There should be
  no gap marker at the end of the file. `qprep` tries to detect gaps between
  molecules itself, but check that the molecules are properly separated in the
  final topology.
- Atom numbering is not significant. `qprep` renumbers atoms from 1, and the
  numbers change because hydrogens are inserted.
- Residue names must match library entry names.
- Atom names must match the atom names of the library entry.
- Temperature factors and occupancies are ignored and are optional.

## Atom masks

An atom mask defines a subset of atoms. It is used:

- in `qprep`, to select the atoms written to PDB and mol2 files;
- in `qdyn`, to select the atoms written to the trajectory;
- in `qcalc`, to select atoms for superposition and deviation calculations.

A mask line selects the atoms within a sequence that have a set of
properties:

| Property | Atoms |
|---|---|
| `solute` | all atoms except solvent |
| `heavy` | all atoms except hydrogens (atoms with a mass of at least 4 u) |
| `excluded` | atoms outside the simulation sphere, including excluded solvent |
| `restrained` | atoms in the restrained boundary shell and atoms outside the sphere |

- `not` before a property negates it.
- Several properties on one line are combined with a logical *and*.
- The line ends with the first and last number of an atom sequence. The word
  `residue` before the numbers selects a sequence of residues. The last number
  may be omitted to select one atom or residue.
- Several lines are combined with a logical *or*.
- Without a mask, all atoms are selected.

```text
heavy not excluded residue 1 104
residue 105 106
```

## Commands

| Command | Alias | Arguments | Description |
|---|---|---|---|
| `addbond` | `ab` | | Add extra bonds, such as S–S |
| `average` | `av` | | Average coordinates |
| `boundary` | `bc` | boundary condition | Set the boundary condition: sphere (1) or box (2) |
| `changeimp` | | | Redefine (specified) improper torsions |
| `checkangs` | `ca` | energy threshold | Check angle energies |
| `checkbonds` | `cb` | energy threshold | Check bond energies |
| `checkimps` | `ci` | energy threshold | Check improper torsion energies |
| `checktors` | `ct` | energy threshold | Check torsion energies |
| `clearbond` | `clearbonds` | | Clear extra bonds |
| `clearlib` | `cl` | | Clear all loaded libraries |
| `cleartop` | | | Clear the current topology |
| `help` | `h`, `?` | | Show the command list |
| `listres` | `lr` | residue number | List the atoms of a residue |
| `listseq` | `ls` | | List the residue sequence |
| `makeshell` | `ms` | | Fix the mask of the atoms in the restrained shell |
| `maketop` | `mt` | name | Generate the topology |
| `mask` | `ma` | mask definition, or `none` | Add to or clear the atom mask |
| `preferences` | `prefs` | | List the preferences |
| `quit` | `q` | | Quit |
| `readframe` | `rf`, `readtrajectory` | trajectory file, frame number | Load coordinates from a Q trajectory file |
| `readlib` | `rl` | library file | Read a library file. Repeat to load several. |
| `readnext` | `rn` | | Load the next frame of the current trajectory |
| `readpdb` | `rp` | PDB file | Read a PDB file |
| `readprm` | `rprm`, `readff`, `rff` | parameter file | Read a force-field parameter file |
| `readtop` | `rt` | topology file | Read a topology file |
| `readx` | `rx` | restart file | Load coordinates from a Q restart file |
| `set` | `s` | name, value | Set a preference |
| `solvate` | `so` | see below | Add solvent |
| `trajectory` | `tr`, `trj` | trajectory file | Open a trajectory file |
| `writemol2` | `wm` | molecule number, mol2 file, hydrogen flag, water flag | Write a SYBYL mol2 file |
| `writepdb` | `wp` | PDB file, gap flag | Write a PDB file of the atoms in the current mask |
| `writetop` | `wt` | topology file, title | Write the topology file |
| `xlink` | `xl`, `crosslink` | | Search for possible cross-links and make bonds |

`average` is in the source of this repository but not in the 2017 manual.

### solvate

For a **sphere**:

| Arguments | Effect |
|---|---|
| centre, radius, `grid`, solvent name | Add solvent molecules from a grid |
| centre, radius, `file`, file name | Add solvent from a file with a box or a sphere of solvent |
| centre, radius, `restart`, file name, solvent name | Add solvent from a restart file of the same solute with solvent |

For a **box**:

| Arguments | Effect |
|---|---|
| `grid`, solvent name, box centre, box size | Add solvent molecules from a grid |
| file name, box centre, box size | Add solvent from a file with a box of solvent |
| `restart`, solvent name, file name, box centre, box size | Add solvent from a restart file |

### writemol2

- *Molecule number* is the number of the molecule in the topology (molecules
  are separated by gaps). `0` writes all molecules.
- *mol2 file* is the file to create; an existing file is overwritten. `auto`
  generates the name from the template
  `coordfile.frame.molecule.mol2`.
- *Hydrogen flag* and *water flag* say whether hydrogens and water are
  written: `y` or `n`.

### writepdb

The *gap flag* says whether `GAP` lines are written between molecules: `y` or
`n`.

### xlink

For each pair of non-bonded heavy atoms closer than the cross-link limit
(2.1 Å), the program asks whether to add a bond.

## Preferences

Use with `set`, for example `set solvent_pack 2.6`.

| Name | Meaning | Default |
|---|---|---|
| `solvent_pack` | Minimum distance between solute and solvent heavy atoms when adding solvent | 2.4 Å |
| `solvent_names` | Residue names recognised as solvent | set from the loaded libraries |
| `solute_density` | Number density of the solute | 0.05794 Å⁻³ |
| `random_seed_solute` | Seed of the random sequence used to generate solute hydrogen coordinates | 179857 |
| `random_seed_solvent` | Seed of the random sequence used to generate solvent hydrogen coordinates | 758971 |

The 2017 manual also lists `max_xlink` (2.1 Å), the upper bond length used by
`xlink`. In this source it is a fixed value and cannot be set.
