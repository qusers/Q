# Electrostatics in the sphere

Electrostatic interactions decay slowly, and a simulation of a finite sphere
has to say what it does with them at three distances: within the cut-off,
beyond the cut-off but inside the sphere, and outside the sphere. This page
goes through the three in turn and then draws the consequences for charged
groups and for perturbations that change the net charge. It builds on
[The spherical boundary](spherical-boundary.md).

| Range | Treatment |
|---|---|
| Within the cut-off | explicit Coulomb interaction between atoms |
| Beyond the cut-off, inside the sphere | local reaction field (LRF): a multipole expansion |
| Q-atoms with everything in the sphere | explicit, no cut-off |
| Outside the sphere | not represented in the energy; assumed to be a uniform dielectric |

<!-- FIGURE TODO: electrostatics-ranges
What: schematic of the sphere seen from one charge group near the centre:
  an inner disc (radius = cut-off, 10 A) labelled "explicit", the rest of the
  sphere labelled "LRF expansion", and the outside shaded and labelled
  "continuum, epsilon = 80 (implicit)". A second, smaller inset for a Q-atom
  where the explicit region is the whole sphere.
Why: the table above is the central message of the page; a picture of the
  three ranges fixes it in memory.
How: a drawn schematic (matplotlib patches or a vector editor), not a
  rendering. Use the same colours for free/excluded regions as in the
  sphere-regions figure of the boundary page.
-->

## Charge groups and cut-offs

Atoms are grouped into **charge groups** in the fragment library. The partial
charges of a group add up to an integer, usually zero. Cut-offs are applied
to whole groups: either all atom pairs between two groups interact, or none
do.

The reason is the range of the interaction. Two neutral groups interact as
dipoles, with an energy that falls off as $1/r^3$. If the cut-off sphere of
an atom sliced through a group and included only some of its atoms, the atom
would see a net partial charge, with an interaction that falls off as $1/r$,
and the energy would jump as atoms cross the cut-off.

How the distance between two groups is measured depends on `switch_atoms` in
the parameter file:

- `on`: the distance between the designated **switching atoms** of the two
  groups (the first atom of each line in `[charge_groups]`);
- `off`: the groups interact if any pair of their atoms is within the
  cut-off.

There are separate cut-offs for solute–solute, solute–solvent and
solvent–solvent pairs, 10 Å each by default. The pair lists are rebuilt at
the `non_bond` interval.

## Q-atoms have no cut-off

The cut-off for interactions between Q-atoms and everything else is
`q_atom`, 99 Å by default, which is larger than any sphere. Q-atoms interact
among themselves without any cut-off. The atoms whose free energy is being
calculated therefore see every charge in the sphere explicitly, and their
state energies contain no truncation error.

## Beyond the cut-off: the local reaction field

For the other atoms, interactions beyond the cut-off are not dropped. With
`lrf on` they are included through the **local reaction field** method of
Lee and Warshel, a multipole expansion.

Consider a charge group $c$ with its centre at $\mathbf{R}_c$. All charges
$q_j$ of groups beyond the cut-off produce, near $\mathbf{R}_c$, an electrostatic
potential that varies slowly in space. It is expanded in a Taylor series
around the centre:

$$
\phi_c(\mathbf{R}_c + \mathbf{d}) \approx
\phi^{(0)} + \sum_{\alpha} \phi^{(1)}_{\alpha} d_\alpha
+ \tfrac{1}{2} \sum_{\alpha\beta} \phi^{(2)}_{\alpha\beta} d_\alpha d_\beta
+ \tfrac{1}{6} \sum_{\alpha\beta\gamma} \phi^{(3)}_{\alpha\beta\gamma} d_\alpha d_\beta d_\gamma
$$

where $\mathbf{d}$ is the displacement from the centre and the coefficients
are the potential and its derivatives at the centre, summed over the distant
charges:

$$
\phi^{(0)} = \sum_j \frac{q_j}{|\mathbf{R}_c - \mathbf{r}_j|}, \qquad
\phi^{(1)}_{\alpha} = \frac{\partial \phi^{(0)}}{\partial R_{c,\alpha}}, \qquad \dots
$$

The energy of atom $i$ of the group is $q_i\,\phi_c(\mathbf{r}_i)$ and the
force follows from the gradient. `qdyn` uses the expansion to second order
for the energy and to third order for the forces. Each pair is counted from
both sides, so the total LRF energy is
$E_{\mathrm{LRF}} = \tfrac{1}{2}\sum_i q_i\,\phi_{c(i)}(\mathbf{r}_i)$.

The practical points are:

- **The coefficients are computed when the pair lists are rebuilt**, every
  `non_bond` steps. They are then held fixed, but the expansion is evaluated
  at the current atom positions in every step. This is different from
  holding the long-range *forces* constant between updates: the field of the
  distant charges is frozen, but each atom still feels how that field changes
  as it moves within its group.
- **The cost is small**, because the expensive part is done only at list
  updates.
- **With the default `lrf` cut-off of 99 Å, every pair of groups in the
  sphere is covered**, either explicitly or through the expansion. No
  electrostatic interaction inside the sphere is neglected.
- **Q-atoms and excluded atoms are not part of the expansion.** Q-atoms do
  not need it, since they have no cut-off.

In `qdyn` the keyword `lrf` in `[MD]` is off by default. QligFEP turns it on
in every input it writes. Without it, interactions beyond the cut-off are
simply truncated, which gives large errors in electrostatic energies; do not
run free-energy calculations that way.

<!-- FIGURE TODO: lrf-expansion
What: one charge group (for example a peptide C=O) at the origin with its
  centre marked, a dashed circle for the cut-off, and a handful of distant
  charges outside it. Field lines or a smooth colour gradient of the potential
  of the distant charges across the group, with the group's atoms displaced
  from the centre by small vectors d.
Why: explains "Taylor expansion of the far field around the group centre"
  without the index notation.
How: matplotlib. Place 5-8 point charges at 12-20 A, compute their potential
  on a grid around the origin, show it as a contour plot within ~3 A of the
  centre, and overlay the exact potential vs the 2nd-order expansion along a
  line through the centre as a small inset to show how good the expansion is
  over the size of a group.
-->

## Outside the sphere

Nothing outside the sphere appears in the energy function: excluded solute
atoms have no non-bonded interactions, and there is no solvent. Physically,
the outside is taken to be a uniform dielectric continuum, bulk water. The
polarisation restraints make the surface water orient as it would next to
such a continuum, but the **energy** of polarising the continuum is not
computed.

For a sphere of radius $R$ with net charge $Q$ that energy is, to leading
order, given by Born's formula:

$$
\Delta G_{\mathrm{Born}} = -\frac{1}{4\pi\varepsilon_0}\,\frac{Q^2}{2R}\left(1 - \frac{1}{\varepsilon}\right) .
$$

It is not small. For $Q = \pm 1$, $R$ = 25 Å and $\varepsilon$ = 80 it is
−6.6 kcal/mol, and it grows with the square of the charge.

This has three consequences.

### The two legs must match

A binding free energy is a difference between a simulation in the protein
and a simulation in water. The missing contribution of the outside cancels
if it is the same in both, and it is the same if the two spheres have

- the **same radius**, and
- the **same net charge**.

This is the reason behind several rules of the Q workflow: the water leg
uses a sphere of the same size as the protein leg, and the net charge of the
protein sphere is controlled, as described next.

### Charged groups near the boundary are neutralised

An ionised side chain close to the edge of the sphere is a problem for two
reasons.

- In reality it is surrounded by solvent, and its field at the centre is
  screened by a dielectric constant of the order of 80. In the sphere it has
  water on one side only, or none if it is excluded, and acts on the ligand
  almost unscreened.
- It adds to the net charge of the sphere, which then differs from that of
  the water leg.

Ionisable groups far from the centre are therefore modelled in a neutral
form. Their interaction with the ligand is thereby neglected. For a neutral
ligand that interaction is small. For a charged ligand, it can be estimated
afterwards from Coulomb's law with a high dielectric constant (Marelius et
al.):

$$
\Delta G^{\mathrm{el}}_{\mathrm{corr}} = \frac{1}{4\pi\varepsilon_0}
\sum_{p}\sum_{l} \frac{q_p\,q_l}{\varepsilon\,r_{pl}}
$$

where $p$ runs over the neutralised ionic groups with formal charge $q_p$,
$l$ over the ligand atoms with partial charges $q_l$, and $\varepsilon$ is
typically 80. QligFEP does not apply this correction.

In QligFEP, `qprep_prot` does the neutralisation when it prepares the
protein:

- an ionisable group is neutralised if it is excluded from the sphere, or if
  it lies in the outer layer of the sphere, by default the outermost 3 Å
  (`--neutralize_boundary_offset`);
- ASP, GLU, ARG, LYS and doubly protonated HIS are replaced by their neutral
  forms (for example `ASH`, `GLH`, `ARN`, `LYN`, `HID` in AMBER naming);
  charged chain termini are neutralised as well;
- if a neutralised group forms a salt bridge with a group inside, the partner
  is neutralised too, so that an ion pair is not left as a lone charge;
- `--skip-neutralization` turns all of this off, which is strongly
  discouraged.

Check the log of `qprep_prot`: it lists the residues that were changed and
the total charge before and after.

<!-- FIGURE TODO: boundary-neutralisation
What: the protein inside a 25 A sphere with ionisable side chains drawn as
  sticks and coloured by treatment: charged and kept (inside 22 A), neutralised
  (in the outer 3 A layer or excluded), neutralised salt-bridge partner.
  Dashed circles at 22 A and 25 A. Ligand at the centre in a contrasting
  colour.
Why: users often do not realise that the prepared protein has different
  protonation states from their input; the picture shows which residues change
  and why.
How: PyMOL, from the protein PDB before and after `qprep_prot`. Select the
  residues whose name changed between the two files (for example ASP -> ASH)
  with a short Python script, colour the three classes, show the sphere
  outlines as in the sphere-regions figure.
-->

### Perturbations that change the net charge

If the two ligands of a perturbation have different net charges, the net
charge of the sphere changes along the path, and the Born term changes with
it. Two things are involved.

**The polarisation of the surface water.** The target orientation of the
surface water depends on the net charge of the Q-atoms at the current state
weights (see [A charged sphere](spherical-boundary.md#a-charged-sphere)), so
the surface water follows the changing charge correctly.

**The energy of the outside continuum.** It is not in the energy function
and therefore not in the free energy of either leg. It cancels between the
protein and the water leg only if the charge changes by the same amount, in
a sphere of the same size, in both; and even then the cancellation relies on
the two spheres responding to the charge in the same way, which a protein
sphere with frozen and neutralised parts and a pure water sphere do not do
exactly. Charge-changing perturbations are therefore harder to converge and
to validate than neutral ones.

`setupFEP` has an option for these edges, `-cm` (`--charge-method`):

| Value | What it does |
|---|---|
| `ion_match` (default) | adds Na⁺ or Cl⁻ counter-ions to the water leg so that its total charge equals that of the protein leg |
| `coalchemical_water` | turns a water molecule of the water leg into a Na⁺ or Cl⁻ along the perturbation, so that the net charge of the water leg does not change |
| `none` | no special handling |

Whichever is used, inspect both legs of a charge-changing edge separately,
check the total charges reported in the `qdyn` logs, and do not mix methods
within one network.

## Settings

| Setting | Where | `qdyn` default | QligFEP |
|---|---|---|---|
| `solute_solute`, `solute_solvent`, `solvent_solvent` | `[cut-offs]` | 10 Å | 10 Å |
| `q_atom` | `[cut-offs]` | 99 Å | 99 Å |
| `lrf` (on/off) | `[MD]` | off | on |
| `lrf` (cut-off) | `[cut-offs]` | 99 Å | 99 Å |
| `non_bond` (list and LRF update interval) | `[intervals]` | 10 steps | 25 steps |
| `charge_correction` | `[solvent]` | follows `polarization` | default |
| `switch_atoms`, `scale_14`, `coulomb_constant` | `[options]` of the parameter file | — | set by the force field |

## References

- Lee, F. S.; Warshel, A. A local reaction field method for fast evaluation
  of long-range electrostatic interactions in molecular simulations.
  *J. Chem. Phys.* **1992**, *97*, 3100–3107.
- Åqvist, J. Calculation of absolute binding free energies for charged
  ligands and effects of long-range electrostatic interactions.
  *J. Comput. Chem.* **1996**, *17*, 1587–1597.
- Åqvist, J.; Hansson, T. On the validity of electrostatic linear response
  in polar solvents. *J. Phys. Chem.* **1996**, *100*, 9512–9521.
- Marelius, J.; Kolmodin, K.; Feierberg, I.; Åqvist, J. *J. Mol. Graph.
  Model.* **1998**, *16*, 213–225.
