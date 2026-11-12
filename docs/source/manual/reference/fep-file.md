# FEP file format

The FEP file (extension `.fep`) designates atoms of the topology as Q-atoms
and redefines the topology for these atoms in a number of states, between
which transformations can be made. See [FEP files](../fep-files.md) for
worked examples and [How Q does FEP](../../concepts/fep-in-q.md) for the
concepts.

## Layout

- The file is plain text, divided into sections that may appear in any
  order.
- A section title is enclosed in square brackets and must be the first item
  on its line. Titles are not case-sensitive.
- Within a section, data are either lines of values or keyword–value pairs.
  A keyword and its value must be on the same line. Keyword lines may appear
  in any order.
- White space is not significant.
- Comments start with `!`, `#` or `*` and may follow values or stand on
  their own line.

Only `[atoms]` is required.

## How the sections work together

- `[FEP]` gives the number of states. It can also give an offset that is
  added to all topology atom numbers, so that the same file serves, for
  example, the free and the bound ligand.
- New atom types are defined in `[atom_types]` and assigned to Q-atoms, per
  state, in `[change_atoms]`.
- Pairs of atoms between which a bond forms or breaks can use the
  exponential repulsion instead of Lennard-Jones by listing them in
  `[soft_pairs]`.
- Non-bonded interactions can be turned off per state in `[excluded_pairs]`.
- New bond, angle, torsion and improper types can be defined and used for
  **any** atoms of the topology, not only Q-atoms. These sections therefore
  refer to atoms by topology number, not Q-atom number. Definitions in the
  FEP file override those of the topology. A type number of 0 disables the
  interaction in that state.
- Angles, torsions and impropers that depend on a bond being formed or
  broken should be coupled to that bond: their energy is scaled by the ratio
  of the current Morse bond energy to the dissociation energy. See
  `[angle_couplings]`, `[torsion_couplings]` and `[improper_couplings]`.
- Extra SHAKE constraints between any pair of atoms are defined in
  `[shake_constraints]`.
- The off-diagonal Hamiltonian elements that mix states in EVB are defined
  in `[off_diagonals]`.
- With periodic boundary conditions a `[PBC]` section is needed. It names one
  switching atom for all Q-atoms, used to build the non-bonded pair lists
  between Q-atoms and their surroundings.

## Soft-core

If atoms vanish or appear, a soft-core potential may be wanted for the
Q-atoms concerned [43]. The soft-core potential of `qdyn` is

$$
V_{\mathrm{vdW}}(r_{ij}) = \frac{A_{ij}}{(r_{ij}^6 + \alpha)^2} - \frac{B_{ij}}{r_{ij}^6 + \alpha}
\qquad\text{or}\qquad
V_{\mathrm{vdW}}(r_{ij}) = \varepsilon_{ij} \left[ \frac{R_{ij}^{*12}}{(r_{ij}^6 + \alpha)^2} - 2\,\frac{R_{ij}^{*6}}{r_{ij}^6 + \alpha} \right]
$$ (softcore)

Depending on `softcore_use_max_potential` in `[FEP]`, the α values are read
directly from `[softcore]`, or calculated for each pair from the desired
potential at $r = 0$ given in `[softcore]`. Soft-core is available only for
Q-atoms: Q–Q, Q–solvent and Q–solute interactions. For a Q–Q interaction in
which both atoms have a soft-core value, the α that gives the lowest
potential at $r = 0$ is used.

[Soft-core and the endpoint problem](../../concepts/softcore.md) explains
what this potential does and does not do in practice.

## Sections

### atoms

Define the Q-atoms. Required.

| Column | Description |
|---|---|
| 1 | Q-atom number, counting from 1 |
| 2 | Topology atom number |

### PBC

| Keyword | Value | Comment |
|---|---|---|
| `switching_atom` | Topology atom number | Required with periodic boundary conditions. The atom used as switching atom for the non-bonded pair list of the Q-atoms. |

### FEP

General perturbation information.

| Keyword | Value | Comment |
|---|---|---|
| `states` | Number of FEP/EVB states | Optional, default 1 |
| `offset` | Topology atom number | Optional, default 0. Added to all topology atom numbers in the file. |
| `offset_residue` | Residue number | Optional. Sets the offset to the topology number of the first atom of that residue, minus one. |
| `offset_name` | Residue name | Optional. Sets the offset to the topology number of the first atom of the first residue with that name, minus one. |
| `qq_use_library_charges` | `on` or `off` | Optional, default `off`. With `on`, the charges of the topology are used for interactions among Q-atoms, so only the electrostatic interactions between Q-atoms and their surroundings change. For studies of, for example, electrostatic linear response. |
| `softcore_use_max_potential` | `on` or `off` | Optional, default `off`. With `on`, the values in `[softcore]` are the desired maximum potentials (kcal/mol) at $r = 0$, and `qdyn` calculates the pairwise $\alpha_{ij}$ of {eq}`softcore`. With `off`, the values are used directly as α. |

### change_charges

Redefine the charges of Q-atoms.

| Column | Description |
|---|---|
| 1 | Q-atom number |
| 2… | Charge (*e*) in state 1, state 2, … |

### atom_types

Define new atom types for Q-atoms: Lennard-Jones parameters and the
parameters of the exponential repulsion
$V_{\mathrm{soft}} = C_i C_j\, e^{-a_i a_j r_{ij}}$.

| Column | Description |
|---|---|
| 1 | Name, at most 8 characters |
| 2 | Lennard-Jones $A_i$ (kcal<sup>½</sup>·mol<sup>−½</sup>·Å⁶) for the geometric rule, or $R^*_i$ (Å) for the arithmetic rule |
| 3 | Lennard-Jones $B_i$ (kcal<sup>½</sup>·mol<sup>−½</sup>·Å³), or $\varepsilon_i$ (kcal·mol⁻¹) |
| 4 | Soft repulsion constant $C_i$ (kcal<sup>½</sup>·mol<sup>−½</sup>) |
| 5 | Soft repulsion distance parameter $a_i$ (Å<sup>−½</sup>) |
| 6 | Lennard-Jones $A_i$ or $R^*_i$ for 1–4 interactions |
| 7 | Lennard-Jones $B_i$ or $\varepsilon_i$ for 1–4 interactions |
| 8 | Atomic mass (u) |

### change_atoms

Assign Q-atom types to Q-atoms.

| Column | Description |
|---|---|
| 1 | Q-atom number |
| 2… | Q-atom type name in state 1, state 2, … |

### soft_pairs

Pairs that use the soft repulsion.

| Column | Description |
|---|---|
| 1 | Q-atom number of the first atom |
| 2 | Q-atom number of the second atom |

### excluded_pairs

Pairs to exclude from non-bonded interactions. Atoms that are not Q-atoms
can be excluded as well.

| Column | Description |
|---|---|
| 1 | Topology atom number of the first atom |
| 2 | Topology atom number of the second atom |
| 3… | Exclusion in effect (1) or not (0) in state 1, state 2, … |

### el_scale

Scale the electrostatic interaction of Q-atom pairs. Useful, for example,
when highly charged intermediates appear in FEP/EVB. The factor applies in
all states. Only Q-atom pairs can be scaled.

| Column | Description |
|---|---|
| 1 | Q-atom number of the first atom |
| 2 | Q-atom number of the second atom |
| 3 | Electrostatic scale factor (0 to 1) |

### softcore

Soft-core values of Q-atoms. Their meaning depends on
`softcore_use_max_potential`.

| Column | Description |
|---|---|
| 1 | Q-atom number |
| 2… | For state 1, state 2, …: the desired potential at $r = 0$ for all van der Waals interactions of this Q-atom, or the α of {eq}`softcore`. An α of 200 gives potentials at $r = 0$ of 10 to 50 kcal/mol for interactions between heavy atoms. 0 means no soft-core for this Q-atom in this state. |

### monitor_groups

Groups of atoms whose non-bonded interactions are to be monitored (printed
in the log).

| Column | Description |
|---|---|
| 1… | Topology atom numbers of the atoms of the group |

### monitor_group_pairs

Pairs of monitor groups whose total non-bonded interaction is calculated.

| Column | Description |
|---|---|
| 1 | Number of the first monitor group |
| 2 | Number of the second monitor group |

### bond_types

Q-bond types, with Morse or harmonic potentials:

$$
E_{\mathrm{Morse}} = D_e \left( 1 - e^{-\alpha (r - r_0)} \right)^2,
\qquad
E_{\mathrm{harmonic}} = \tfrac{1}{2} k_b (r - r_0)^2 .
$$

The two kinds can be mixed, but each bond type is one or the other. An entry
with four values is a Morse potential; an entry with three is harmonic.

| Column | Morse potential | Harmonic potential |
|---|---|---|
| 1 | Q-bond type number, from 1 | Q-bond type number, from 1 |
| 2 | Dissociation energy $D_e$ (kcal·mol⁻¹) | Force constant $k_b$ (kcal·mol⁻¹·Å⁻²) |
| 3 | Exponent $\alpha$ (Å⁻¹) | Equilibrium bond length $r_0$ (Å) |
| 4 | Equilibrium bond length $r_0$ (Å) | |

:::{note}
The 2017 manual gives the unit of the Morse exponent as Å⁻². The exponent
multiplies a length, so its unit is Å⁻¹.
:::

### change_bonds

Assign Q-bond types. SHAKE constraints of the redefined bonds are removed.
The order of the two atoms is not important.

| Column | Description |
|---|---|
| 1 | Topology atom number of the first atom |
| 2 | Topology atom number of the second atom |
| 3… | Q-bond type number in state 1, state 2, …, or 0 to disable the bond |

### angle_types

| Column | Description |
|---|---|
| 1 | Q-angle type number, from 1 |
| 2 | Harmonic force constant (kcal·mol⁻¹·rad⁻²) |
| 3 | Equilibrium angle (°) |

### change_angles

| Column | Description |
|---|---|
| 1 | Topology atom number of the first atom |
| 2 | Topology atom number of the middle atom |
| 3 | Topology atom number of the third atom |
| 4… | Q-angle type number in state 1, state 2, …, or 0 to disable the angle |

### torsion_types

| Column | Description |
|---|---|
| 1 | Q-torsion type number, from 1 |
| 2 | Force constant = ½ · barrier height (kcal·mol⁻¹) |
| 3 | Periodicity (number of maxima per turn) |
| 4 | Phase shift (°) |

### change_torsions

The order of the atoms (1, 2, 3, 4 or 4, 3, 2, 1) is not important.

| Column | Description |
|---|---|
| 1–4 | Topology atom numbers of the four atoms |
| 5… | Q-torsion type number in state 1, state 2, …, or 0 to disable the torsion |

### improper_types

| Column | Description |
|---|---|
| 1 | Q-improper type number, from 1 |
| 2 | Harmonic force constant (kcal·mol⁻¹·rad⁻²). Impropers defined here are always harmonic. |
| 3 | Equilibrium angle (°) |

### change_impropers

The order of the atoms (1, 2, 3, 4 or 4, 3, 2, 1) is not important.

| Column | Description |
|---|---|
| 1–4 | Topology atom numbers of the four atoms |
| 5… | Q-improper type number in state 1, state 2, …, or 0 to disable the improper |

### angle_couplings

Couple Q-angles to Q-bonds: the angle energy is scaled by the ratio of the
current Morse bond energy to the dissociation energy.

| Column | Description |
|---|---|
| 1 | Q-angle number (line number within `[change_angles]`) |
| 2 | Q-bond number (line number within `[change_bonds]`) |

### torsion_couplings

| Column | Description |
|---|---|
| 1 | Q-torsion number (line number within `[change_torsions]`) |
| 2 | Q-bond number (line number within `[change_bonds]`) |

### improper_couplings

| Column | Description |
|---|---|
| 1 | Q-improper number (line number within `[change_impropers]`) |
| 2 | Q-bond number (line number within `[change_bonds]`) |

### shake_constraints

Extra SHAKE constraints. The effective constraint distance is the sum of the
distances of the states, weighted by their λ. Constraints defined here do
not override those imposed by the `shake_*` keywords of the `qdyn` input; to
remove a constraint, redefine the bond as a Q-bond. The order of the atoms is
not important.

| Column | Description |
|---|---|
| 1 | Topology atom number of the first atom |
| 2 | Topology atom number of the second atom |
| 3… | Constraint distance (Å) in state 1, state 2, … |

### off_diagonals

Off-diagonal elements of the Hamiltonian,
$H_{ij} = A_{ij}\, e^{-\mu_{ij} r_{kl}}$, where $i$ and $j$ are states and
$k$ and $l$ are Q-atoms.

| Column | Description |
|---|---|
| 1 | State $i$ |
| 2 | State $j$ |
| 3 | Q-atom $k$ |
| 4 | Q-atom $l$ |
| 5 | $A_{ij}$ (kcal·mol⁻¹) |
| 6 | $\mu_{ij}$ (Å⁻¹) |
