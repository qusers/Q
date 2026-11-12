# Fragment library file

The fragment library is a text file that defines the molecular building
blocks: amino-acid residues, ligands and so on. `qprep` uses it to build a
topology from a PDB file. The layout follows the same conventions as the
parameter file.

Each entry starts with its name in curly braces, for example `{ALA}`. The
lists of atoms, bonds and so on are sections within the entry. Atoms are
referred to by name in all sections after `[atoms]`.

## Sections

### info

General information about the fragment. Optional.

| Keyword | Value | Comment |
|---|---|---|
| `SYBYLtype` | `RESIDUE` or `GROUP` | Optional. SYBYL substructure type, used only when writing mol2 files. |
| `PDBtype` | `ATOM` or `HETATM` | Optional, default `ATOM`. `CONECT` records are written for `HETATM` groups. |
| `Solvent` | `on` or `off` | If `on`, the entry is recognised as a solvent. |
| `Density` | number density (Å⁻³) | Solvents only. Used to solvate from a grid and to calculate the effective solvent radius. |

### atoms

Atom names, types and charges. The atom name must match the name in the PDB
files to be read. The order of the atoms is not important.

| Column | Description |
|---|---|
| 1 | sequence number, from 1 to the number of atoms |
| 2 | atom name (4 characters) |
| 3 | atom type |
| 4 | partial charge (*e*) |

### bonds

The bonds within the fragment.

| Column | Description |
|---|---|
| 1 | name of atom 1 |
| 2 | name of atom 2 |

### connections

The atoms that form bonds to neighbouring residues. Optional; do not define
it for entries that are whole molecules.

| Keyword | Value |
|---|---|
| `head` | name of the atom that is bonded to the tail of the preceding residue |
| `tail` | name of the atom that is bonded to the head of the next residue |

### build rules

Rules for generating hydrogen coordinates.

| Column | Description |
|---|---|
| 1 | type of rule. The only type is `torsion`. |
| 2 | name of the hydrogen to be generated |
| 3 | name of atom 2 of the torsion (the atom the hydrogen is bonded to) |
| 4 | name of atom 3 of the torsion |
| 5 | name of atom 4 of the torsion |
| 6 | target value of the torsion angle formed by the atoms in columns 2 to 5 |

### impropers

Improper torsions. Used only for force fields in which impropers are defined
explicitly rather than generated automatically (see the parameter file).

| Column | Description |
|---|---|
| 1 | name of atom 1. Put `+` or `-` before the name to refer to an atom in the next or previous residue. |
| 2 | name of atom 2, the central atom to which atoms 1, 3 and 4 are bonded |
| 3 | name of atom 3 (`+`/`-` as for atom 1) |
| 4 | name of atom 4 (`+`/`-` as for atom 1) |

### charge_groups

One charge group per line.

| Column | Description |
|---|---|
| 1 | name of the switching atom |
| 2… | names of the other atoms |

## Example

The CHARMM entry for alanine:

```text
{ALA}                  ! Alanine
[info]
SYBYLtype  RESIDUE     ! SYBYL substructure type
[atoms]                ! number, name, type, charge
 1  N    NH1   -0.470
 2  HN   H      0.310
 3  CA   CT1    0.070
 4  HA   HB     0.090
 5  CB   CT3   -0.270
 6  HB1  HA     0.090
 7  HB2  HA     0.090
 8  HB3  HA     0.090
 9  C    C      0.510
10  O    O     -0.510
[bonds]
CB  CA
N   HN
N   CA
O   C
C   CA
CA  HA
CB  HB1
CB  HB2
CB  HB3
[connections]          ! how to bond to the previous and next residue
head  N
tail  C
[impropers]
N  -C  CA  HN
C  CA  +N  O
[charge_groups]        ! switching atom first
N   HN   CA   HA
CB  HB1  HB2  HB3
C   O
```
