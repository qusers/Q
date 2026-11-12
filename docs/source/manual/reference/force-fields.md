# Force fields

Q is not tied to a force field. Force fields are defined in parameter files,
separate from the program, and any force field can be used as long as it has
the functional form

$$
\begin{aligned}
V_{\mathrm{pot}} ={}& \sum_{\mathrm{bonds}} \tfrac{1}{2} k_b (r - r_0)^2
+ \sum_{\mathrm{angles}} \tfrac{1}{2} k_\theta (\theta - \theta_0)^2
+ \sum_{\mathrm{dihedrals}} K_\varphi \bigl(1 + \cos(n\varphi - \delta)\bigr) \\
&+ \sum_{\substack{\mathrm{improper} \\ \mathrm{dihedrals}}} \tfrac{1}{2} k_\xi (\xi - \xi_0)^2
+ \sum_{\substack{\mathrm{atom} \\ \mathrm{pairs}\ i,j}} \left[ \frac{1}{4\pi\varepsilon_0} \frac{q_i q_j}{r_{ij}} + \frac{A_{ij}}{r_{ij}^{12}} - \frac{B_{ij}}{r_{ij}^{6}} \right]
\end{aligned}
$$ (ff-potential)

| Symbol | Meaning |
|---|---|
| $k_b$, $r_0$ | bond force constant and reference length |
| $k_\theta$, $\theta_0$ | angle force constant and reference angle |
| $K_\varphi$, $n$, $\delta$ | dihedral force constant, multiplicity (number of minima per turn) and phase shift ($\delta/n$ is the position of the first maximum) |
| $k_\xi$, $\xi_0$ | out-of-plane force constant and reference angle of an improper dihedral |
| $q_i$, $q_j$, $r_{ij}$ | partial charges of atoms $i$ and $j$ and their distance |
| $A_{ij}$, $B_{ij}$ | Lennard-Jones parameters of the pair |

## Lennard-Jones combination rules

Lennard-Jones parameters are defined per atom type and combined by one of two
rules.

**Geometric.** The atom-type parameters are $A_i = A_{ii}^{1/2}$ and
$B_i = B_{ii}^{1/2}$, and

$$
A_{ij} = A_i A_j, \qquad B_{ij} = B_i B_j .
$$

**Arithmetic.** Some force fields write the 6–12 potential as

$$
\varepsilon_{ij} \left[ \left( \frac{R^*_{ij}}{r_{ij}} \right)^{12} - 2 \left( \frac{R^*_{ij}}{r_{ij}} \right)^{6} \right]
$$

with the atom-type parameters combined as

$$
\varepsilon_{ij} = (\varepsilon_i \varepsilon_j)^{1/2}, \qquad R^*_{ij} = R^*_i + R^*_j .
$$

## Dihedrals and impropers

Several Fourier components, with different $K_\varphi$, $n$ and $\delta$, can
be added for the same dihedral to model rotation barriers more accurately. An
alternative, periodic form of the improper potential, with the same
trigonometric function as normal dihedrals, is also implemented.

## Charge groups and cut-offs

The fragments of a force field, such as amino-acid residues, are divided into
**charge groups**: groups of atoms whose partial charges add up to an integer.
Non-bonded cut-offs are applied to whole groups: either all pairwise
interactions between two groups are evaluated, or none. The size of the
charge groups differs between force fields, from a few atoms to whole
residues.

Some force fields designate a **switching atom** in each charge group and
apply the cut-off to the distance between switching atoms. Others include two
groups if any pair of their atoms is within the cut-off.

## Force fields in this repository

The files are in `src/QligFEP/FF/`:

| Force field | Library | Parameters |
|---|---|---|
| AMBER ff14SB | `AMBER14sb.lib` | `AMBER14sb.prm` |
| OPLS-AA/M | `OPLSAAM.lib` | `OPLSAAM.prm` |
| OPLS-AA (2005) | `OPLS2005.lib` | `OPLS2005.prm` |
| CHARMM36 | `CHARMM36.lib` | `CHARMM36.prm` |

`AMBER14sb` is the default of QligFEP. `OPLSAAM` is the OPLS-AA/M force field
of Robertson, Tirado-Rives and Jorgensen (2015), extracted from the original
parameter and topology files and curated in 2025; it is the force field of
QresFEP. The directory also holds an older translation named `OPLS2015`,
which is likewise labelled OPLS-AA/M; prefer `OPLSAAM`.

These are Q translations of the named force fields; identity with the
original cannot be guaranteed. Do not mix libraries and parameter files of
different force fields.

## What a force-field file encodes

Besides the parameters, a Q library and parameter file pair fixes these
properties of the force field:

- **Aliphatic hydrogens.** Explicit (all-atom) or part of an extended atom.
- **Lennard-Jones form**, and thus the combination rule: geometric for
  $A_{ij}, B_{ij}$, arithmetic for $\varepsilon_{ij}, R^*_{ij}$
  (`vdw_rule`). The treatment of 1–4 interactions is specific to each force
  field (`scale_14` and the 1–4 Lennard-Jones set).
- **Impropers.** Harmonic or periodic (`improper_potential`), generated
  automatically or listed in the library entries (`improper_definition`).
- **Charge groups** and the cut-off logic (`switch_atoms`).

| | `vdw_rule` | `scale_14` | `switch_atoms` | Impropers |
|---|---|---|---|---|
| `AMBER14sb.prm` | arithmetic | 0.8333 | off | periodic |
| `OPLSAAM.prm` | geometric | 0.5 | on | periodic, explicit |
| `CHARMM36.prm` | arithmetic | 1.0 | off | harmonic, explicit |

## Solvents

A solvent is a library entry with `solvent` set in its `[info]` section and
a number density. The libraries of this repository define TIP3P water as
`HOH`.

## Solvent file format

`qprep` can solvate a system by filling the empty space in the simulation
sphere or box with molecules from a solvent file. This is a PDB-like file
that contains a box or a sphere of solvent. The residue name in the file
selects the library entry used to generate bonds and so on.

A box is replicated in all directions, so a small box can solvate a large
sphere or box. A sphere cannot be replicated and must be larger than the
system. The box or sphere is translated to the solvent centre, so the origin
of the file is arbitrary. A sphere cannot be used for a system with a
periodic boundary.

| Line | Content |
|---|---|
| 1 | For a box: the side of the box (Å). For a sphere: the radius, followed by the word `sphere`. |
| 3*n* − 1 | Coordinates of the first atom of solvent molecule *n*, in PDB format |
| 3*n* | Coordinates of the second atom of molecule *n* |
| 3*n* + 1 | Coordinates of the third atom of molecule *n* |

For solvent molecules with more than three atoms, the remaining atoms follow
in the same format.
