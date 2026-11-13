# Force-field parameter file

The parameter file is a text file with the same conventions as the FEP file.
It is divided into sections, in any order, each starting with a title in
square brackets. It holds the constants of
{eq}`ff-potential` for each combination of atom types.

## Sections

### options

| Keyword | Value |
|---|---|
| `vdw_rule` | Rule for combining Lennard-Jones parameters: `geometric` or `arithmetic` |
| `scale_14` | Scaling factor for electrostatic interactions between atoms in 1–4 positions |
| `switch_atoms` | Cut-off logic. `on`: use the switching atoms of the charge groups. `off`: include two charge groups if any pair of atoms is within the cut-off. |
| `name` | Name of the force field, printed in the log |
| `type` (or `force_field`) | Force-field family: `GROMOS` (default), `AMBER` or `CHARMM` |
| `improper_potential` | `harmonic` (default) or `periodic` |
| `improper_definition` | `automatic` (default) or `explicit`. With `explicit`, impropers are taken from the `[impropers]` sections of the library entries. |
| `coulomb_constant` | Coulomb constant, default 332 kcal·mol⁻¹·Å·*e*⁻² |

`vdw_rule` and `scale_14` are required.

### atom_types

A name, Lennard-Jones parameters and a mass for each atom type. There are
three sets of Lennard-Jones parameters: set 1 for normal non-bonded
interactions, set 2 for pairs of polar atom types (those listed in
`LJ_type2_pairs`) and set 3 for atoms in 1–4 position. The B parameter is
shared by sets 2 and 3.

| Column | Geometric rule | Arithmetic rule |
|---|---|---|
| 1 | atom type name, at most 8 characters | |
| 2 | $A_i$ for set 1 (kcal<sup>½</sup>·mol<sup>−½</sup>·Å⁶) | $R^*_i$ for set 1 (Å) |
| 3 | $A_i$ for set 2 | $R^*_i$ for set 2 |
| 4 | $B_i$ for set 1 (kcal<sup>½</sup>·mol<sup>−½</sup>·Å³) | $\varepsilon_i$ for set 1 (kcal·mol⁻¹) |
| 5 | $A_i$ for set 3 | $R^*_i$ for set 3 |
| 6 | $B_i$ for sets 2 and 3 | $\varepsilon_i$ for sets 2 and 3 |
| 7 | mass of the (extended) atom (u) | |
| 8 | SYBYL atom type (5 characters), optional | |

### atom_aliases

Alias names for atom types, for example to keep old numeric names usable in
library files.

| Column | Description |
|---|---|
| 1 | alias name, at most 8 characters |
| 2 | atom type name defined in `atom_types` |

### LJ_type2_pairs

Pairs of polar atom types that interact with Lennard-Jones set 2.

| Column | Description |
|---|---|
| 1 | atom type 1 |
| 2 | atom type 2 |

### bonds

One line for each pair of atom types. The pair 1–2 is the same as 2–1; list
only one.

| Column | Description |
|---|---|
| 1 | atom type 1 |
| 2 | atom type 2 |
| 3 | force constant (kcal·mol⁻¹·Å⁻²) |
| 4 | equilibrium length (Å) |
| 5 | SYBYL bond type (2 characters), optional |

### angles

One line for each triplet of atom types. 1–2–3 is the same as 3–2–1; list
only one.

| Column | Description |
|---|---|
| 1 | atom type 1 |
| 2 | atom type 2 |
| 3 | atom type 3 |
| 4 | force constant (kcal·mol⁻¹·rad⁻²) |
| 5 | equilibrium angle (°) |
| 6 | Urey–Bradley force constant, optional |
| 7 | Urey–Bradley equilibrium distance (Å), optional |

### torsions

A torsion is defined by four atom types, but often only the two middle ones
matter. That case is written with `0` or `?` in columns 1 and 4. A definition
with four types overrides a definition with two.

- The force constant of the cosine function is half the barrier height.
- The periodicity is the number of maxima in a full turn.
- The phase shift divided by the periodicity is the angle of the first
  maximum.
- The number of paths is the number of ways a two-atom torsion can be
  defined, that is, the product of the numbers of atoms bonded to the two
  middle atoms. It distributes the force over all the atoms involved.
- Several terms can be defined for the same atom types; they are added.

| Column | Description |
|---|---|
| 1 | atom type 1, or `0` or `?` to match any type |
| 2 | atom type 2 |
| 3 | atom type 3 |
| 4 | atom type 4, or `0` or `?` to match any type |
| 5 | force constant = ½ · barrier height (kcal·mol⁻¹) |
| 6 | periodicity. A minus sign means that more components follow on the next lines: all but the last component of a torsion are entered with negative periodicity. |
| 7 | phase shift (°) |
| 8 | number of paths |

### impropers

Force constant and equilibrium angle of improper torsions, for a quadruplet
of atom types. A field that is not a defined atom type acts as a wild card.

| Column | Description |
|---|---|
| 1–4 | atom types 1 to 4 |
| 5 | force constant (kcal·mol⁻¹·rad⁻² for harmonic impropers) |
| 6 | equilibrium angle (°) |

## Example

The start of each section of `AMBER14sb.prm`:

```text
[options]
name                           Q-Amber14SB
type                           AMBER
vdw_rule                       arithmetic !vdW combination rule (geometric or arithmetic)
scale_14                       0.8333 ! electrostatic 1-4 scaling factor
switch_atoms                   off
improper_potential             periodic

[atom_types]
* tac               R*1        R*2   epsilon1        R*3 epsilon2&3       mass
2C                1.908        0.0     0.1094      1.908     0.0547      12.01
3C                1.908        0.0     0.1094      1.908     0.0547      12.01

[bonds]
2C           2C                620.0      1.526
2C           3C                620.0      1.526

[angles]
2C           2C           C                 126.0      111.1
2C           2C           CO                126.0      111.1

[torsions]
?            2C           2C           ?                   1.4   3.0        0.0   9.0
?            2C           CA           ?                   0.0   2.0        0.0   1.0

[impropers]
?            C            O            ?                  10.5      180.0
?            C4           H4           ?                   1.1      180.0
```

With the arithmetic rule the Lennard-Jones columns hold $R^*$ and
$\varepsilon$; the 1–4 set has half the $\varepsilon$.
