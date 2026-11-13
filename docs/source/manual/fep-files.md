# FEP files

The FEP file defines a set of atoms as **Q-atoms** and redefines their
interaction parameters. Any force-field parameter of these atoms can be
controlled, and several **states** can be defined. The states may differ very
little, for example in the van der Waals parameters of one atom, or they may
represent different valence-bond structures: reactants and products of a
chemical reaction studied by EVB, or one state per elementary step of a
multi-step reaction. In such a model bonds, angles, torsions, charges and van
der Waals parameters may change for many atoms.

For each configuration, `qdyn` computes the energy of each state and writes
these energies to the energy file. The *mapping potential* that generates the
forces is a mixture of the states, set by the weight λ given to each state in
the `qdyn` input. `qfep` then calculates free-energy differences between the
states with the FEP formula, or with the umbrella-sampling approach to obtain
EVB ground-state free-energy profiles. See
[How Q does FEP](../concepts/fep-in-q.md) for the equations.

The FEP file has the same layout as the `qdyn` input: sections of data, most
of them optional. This page works through three examples. The complete format
is in [FEP file format](reference/fep-file.md).

## Example: charging a benzene molecule

This FEP file calculates the electrostatic contribution to the solvation free
energy of benzene. The molecule is defined in the topology; carbons have odd
atom numbers and hydrogens even.

```ini
[FEP]
states 2                 ! number of states

[atoms]                  ! Q-atom number, topology atom number
1   1
2   2
3   3
4   4
5   5
6   6
7   7
8   8
9   9
10  10
11  11
12  12

[change_charges]         ! Q-atom number, charge in state 1 and 2
1   -0.15  0.0
2   +0.15  0.0
3   -0.15  0.0
4   +0.15  0.0
5   -0.15  0.0
6   +0.15  0.0
7   -0.15  0.0
8   +0.15  0.0
9   -0.15  0.0
10  +0.15  0.0
11  -0.15  0.0
12  +0.15  0.0
```

- `states` in `[FEP]` is the number of states.
- `[atoms]` designates topology atoms as Q-atoms. The first column is the
  Q-atom number, used to refer to the atom later; the second is its number in
  the topology.
- `[change_charges]` gives the charge of Q-atoms in each state. All charges
  change here, but in general only those that change need to be listed.

Bonded and van der Waals parameters are untouched. The file defines two
charge states, one with C–H dipoles of ±0.15 *e* and one with no charges.

## Example: changing van der Waals parameters

This FEP file calculates the difference in hydration free energy between Na⁺
and K⁺. The ions have the same charge, so only two sets of Lennard-Jones
parameters are needed.

```ini
[FEP]
states 2

[atoms]                  ! Q-atom number, topology atom number
1  1

[atom_types]             ! new atom types
!Name   Ai      Bi    Ci   ai   Ai(1-4)  Bi(1-4)  Mass
Na      143.70  3.89  0.0  0.0  0.0      0.0      22.99
K       522.70  4.35  0.0  0.0  0.0      0.0      39.10

[change_atoms]           ! Q-atom number, type in state 1 and 2
1  Na  K
```

- `[atom_types]` gives the parameters of the atom types involved. The first
  column is the type name, followed by the Lennard-Jones parameters $A_i$ and
  $B_i$, or $R^*$ and $\varepsilon$, depending on the combination rule in the
  parameter file used for the topology. Columns four to seven are not used
  here: two parameters of the exponential repulsion and two Lennard-Jones
  parameters for 1–4 interactions. The last column is the mass.
- `[change_atoms]` says that Q-atom 1 is of type `Na` in state 1 and `K` in
  state 2.

No charges change, so `[change_charges]` is omitted.

:::{note}
The mass in `[atom_types]` is not used in the dynamics. The simulation uses
the masses of the topology.
:::

## Example: EVB states for a proton transfer

This example is from a protein tyrosine phosphatase. A proton is transferred
from a cysteine of the enzyme to the doubly charged phosphate group of the
substrate, phenyl phosphate [19]. In state 1 the proton (H) is on the sulphur
and the phosphate carries two negative charges. In state 2 the proton is on a
phosphate oxygen and one negative charge is on the sulphur. Charges, van der
Waals parameters, bonds, angles and torsions all change. EVB itself is
reviewed in [5].

The Q-atoms and their topology numbers are:

| Q-atom | Topology atom | Atom |
|---|---|---|
| 1 | 79 | Cβ of Cys |
| 2 | 80 | S |
| 3 | 81 | H (transferred) |
| 4 | 1542 | P |
| 5 | 1543 | O (ester) |
| 6–7 | 1544, 1545 | O |
| 8 | 1546 | O (acceptor) |

```ini
[FEP]
states 2

[atoms]                  ! Q-atom number, topology atom number
1  79
2  80
3  81
4  1542
5  1543
6  1544
7  1545
8  1546

[change_charges]         ! Q-atom number, charge in state 1 and 2
1   0.180   0.000
2  -0.450  -1.000
3   0.270   0.398
4   0.540   1.230
5  -0.360  -0.360
6  -0.860  -0.860
7  -0.860  -0.860
8  -0.860  -0.548

[atom_types]
!Type  Ai       Bi      Ci     ai     Ai(1-4)  Bi(1-4)  Mass
P      2303.00  59.35   0.0    1.581  2303.00  59.35    30.97
OE     600.00   23.25   70.0   1.581  600.00   23.25    16.00
OD     956.00   23.01   70.0   1.581  956.00   23.01    16.00
H      0.00     0.00    6.5    1.581  0.00     0.00     1.00
C2     2500.00  46.06   0.0    1.581  2500.00  46.06    14.00
SH     2001.57  44.74   165.0  1.581  2001.57  44.74    32.06
S-     2720.00  136.00  165.0  1.581  7200.00  136.00   32.06

[change_atoms]           ! Q-atom number, type in state 1 and 2
1  C2  C2
2  SH  S-
3  H   H
4  P   P
5  OE  OE
6  OD  OD
7  OD  OD
8  OD  OE

[soft_pairs]             ! Q-atom pairs with C*exp(-a*r) repulsion
2  3
3  8

[excluded_pairs]         ! atom i, j, exclusion flag in state 1 and 2
81  1544  0  1
81  1545  0  1

[bond_types]             ! number, De, alpha, b0 (Morse)
1  85.0   2.00  1.61
2  120.0  2.00  1.49
3  84.0   2.00  1.43
4  110.0  2.00  1.00
5  94.0   2.00  1.33
6  112.5  2.00  1.80
7  100.0  2.00  1.53

[change_bonds]           ! atom i, j, type in state 1 and 2 (0 = no bond)
1542  1546  2  1
80    81    5  0
1546  81    0  4

[angle_types]            ! number, force constant, equilibrium angle
1  95.0   109.6
2  140.0  120.0
3  115.0  120.0
4  110.0  109.6
5  0.0    0.0
6  110.0  113.0
7  95.0   96.0

[change_angles]          ! atom i, j, k, type in state 1 and 2 (0 = no angle)
1544  1542  1546  2  1
1545  1542  1546  2  1
1542  1546  81    0  4
79    80    81    7  0

[torsion_types]          ! number, force constant, multiplicity, phase
1  0.75  3.0  0.00
2  0.70  3.0  0.00

[change_torsions]        ! atom i, j, k, l, type in state 1 and 2 (0 = none)
1543  1542  1546  81  0  1
1544  1542  1546  81  0  1
1545  1542  1546  81  0  1
78    79    80    81  2  0

[angle_couplings]        ! Q-angle number, Q-bond number
3  3
4  2

[torsion_couplings]      ! Q-torsion number, Q-bond number
1  3
2  3
3  3
4  2

[off_diagonals]          ! state i, state j, Q-atom k, Q-atom l, Aij, muij
1  2  2  8  1.0  0.45
```

`[FEP]`, `[atoms]`, `[change_charges]`, `[atom_types]` and `[change_atoms]`
are used as before. Q-atom 2, the sulphur, changes its charge from −0.45 to
−1.00 and its type from `SH` to `S-`.

**Soft pairs.** For pairs of atoms that form or break a bond, an exponential
repulsion is more physical than $1/r^{12}$, which is too strong at very short
distances. The pairs listed in `[soft_pairs]` interact through

$$
V_{\mathrm{soft}} = C_i \cdot C_j \cdot e^{-a_i \cdot a_j \cdot r_{ij}}
$$

where $r_{ij}$ is the distance between the atoms. $C$ and $a$ are parameters
of the atom type, combined geometrically. There is no attractive $1/r^6$
term.

**Bonds.** Morse parameters for bonds that are broken or formed are given in
`[bond_types]`. `[change_bonds]` lists the bonds, as pairs of atoms, with the
bond type to use in each state. If the bond is already in the topology, its
harmonic potential is turned off. Type 0 means no bond.

**Angles and torsions** are redefined in the same way: types in
`[angle_types]` and `[torsion_types]`, assignments in `[change_angles]` and
`[change_torsions]`. Q-angles are harmonic, like normal angles. No impropers
change in this example.

**Couplings.** Angles, torsions and impropers depend on the existence of the
bonds between their atoms. They can be coupled to a Morse bond, in which case
their energy is scaled by the ratio of the current Morse bond energy to the
dissociation energy [20]. In `[angle_couplings]` the numbers are line numbers:
the third line of `[change_angles]` (P–O–H) is coupled to the third line of
`[change_bonds]` (O–H), and the fourth (C–S–H) to the second (S–H).
Torsions and impropers are coupled in the same way.

**Off-diagonal elements** of the Hamiltonian are defined in
`[off_diagonals]`:

$$
H_{ij} = A_{ij} \cdot e^{-\mu_{ij} \cdot r_{kl}}
$$

where $i$ and $j$ are states and $r_{kl}$ is the distance between Q-atoms $k$
and $l$. The record here mixes states 1 and 2 through Q-atoms 2 and 8, with
$A$ = 1.0 and $\mu$ = 0.45.

## Monitoring non-bonded interactions

To analyse, for example, receptor–ligand interactions, define groups of atoms
in the FEP file and let `qdyn` report the non-bonded energy between pairs of
groups. This example reports the interactions between the pteridine ring of a
dihydrofolate reductase inhibitor and some side chains and a cofactor amide.

```ini
[monitor_groups]         ! topology atom numbers of each group
266 267 268                            ! Glu 30 carboxylate
317 318 319 320 321 322                ! Phe 34 side chain
1897 1898 1899 1900 1901 1902 1903     ! pteridine ring 1 of MTX
1908 1909 1910 1911 1912 1913 1914     ! pteridine ring 2
1880 1881 1882 1883 1884 1885          ! amide of NADPH

[monitor_group_pairs]    ! pairs of group numbers
1 3
2 4
2 5
```

Five groups are defined and the interactions 1–3, 2–4 and 2–5 are
calculated. The energies are evaluated for each FEP state and printed with
the energy summaries in the log:

```text
= Monitoring selected groups of nonbonded interactions =
pair   Vwsum    Vwel   Vwvdw   1:Vel  1:Vvdw
   1  -58.48  -65.65    7.18  -65.65    7.18
   2   -1.72    0.00   -1.72    0.00   -1.72
   3   -0.08    0.00   -0.08    0.00   -0.08
```

The columns are the pair number, the λ-weighted total, the λ-weighted
electrostatic and Lennard-Jones energies, and the electrostatic and
Lennard-Jones energies of state 1. With a single state the weighted values
equal those of state 1.

`qcalc` has a similar function for saved trajectories.
