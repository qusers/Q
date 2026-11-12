# The spherical boundary

Most molecular-dynamics programs put the whole solute in a periodic box of
solvent. Q simulates a **sphere** cut out of the system: a droplet of explicit
water centred on the region of interest, with the part of the protein that
falls inside it. This page explains what is inside the sphere, what happens
at its edge, and what the model assumes about everything outside. The
treatment of electrostatics is on the [next page](electrostatics.md).

## The idea

A binding free energy or a reaction barrier is determined mostly by what
happens within some 10 to 20 Å of the ligand or the reacting atoms. The
spherical model spends the computer time there:

- atoms inside the sphere move freely and interact as in any simulation;
- atoms of the solute outside the sphere are kept in the system, fixed in
  place, to hold the structure together, but have no non-bonded interactions;
- solvent outside the sphere is not present. Its effect on the water at the
  surface is replaced by restraints that make the surface water behave like
  bulk water.

For a protein such as PTP1B, a periodic simulation needs a cubic box with
sides of about 85 Å; a sphere of 16 Å radius is enough to study its active
site (Marelius et al.). QligFEP uses a radius of 25 Å by default.

<!-- FIGURE TODO: sphere-vs-box
What: two panels of the same protein-ligand complex at the same scale.
  Left: a periodic box (cube outline, solvent as faint dots, whole protein).
  Right: Q's 25 A sphere centred on the ligand; protein inside the sphere in
  colour, protein outside in grey, water only inside the sphere.
Why: shows at a glance how much smaller the simulated system is, and that the
  protein is not truncated but partly "frozen".
How: PyMOL. Load a prepared QligFEP protein leg (the PDB written by qprep with
  water). Draw the sphere with a CGO sphere or a pseudoatom with
  `set sphere_scale`/`sphere_transparency`, centred on the centre used in
  qprep. Colour by distance from the centre: `select inside, all within 25 of
  centre`. For the left panel draw only a cube outline with CGO lines; no
  periodic simulation is needed. Label both panels with the atom counts.
-->

## Regions of the sphere

The sphere is defined in `qprep` by a centre and a radius, and is stored in
the topology. `qprep` then sorts the system into regions. The sorting is done
once, from the coordinates of the topology, and is never updated during the
simulation.

| Region | Which atoms | Treatment in `qdyn` |
|---|---|---|
| **Free** | solute and solvent inside the sphere | normal dynamics |
| **Excluded** | solute outside the sphere radius | no non-bonded interactions; harmonic restraint to the topology coordinates with 200 kcal·mol⁻¹·Å⁻²; bonded terms are still evaluated |

Two details matter in practice.

**Atoms are assigned by charge group, not one by one.** A
[charge group](electrostatics.md#charge-groups-and-cut-offs) is inside or
outside as a whole, depending on the position of its switching atom (or, for
force fields without switching atoms, of the group itself). A group is never
cut in two, so the boundary does not create fractional charges.

**Excluded atoms still hold the structure.** Bonds, angles and torsions that
cross the boundary are kept, so the free part of the protein remains attached
to a rigid scaffold with the shape of the starting structure.

:::{note}
`qdyn` can also restrain the solute in an outer shell of the sphere to its
starting coordinates (`shell_radius` and `shell_force` in `[sphere]`).
QligFEP does not use this: it writes `shell_radius` equal to the sphere
radius, so all solute inside the sphere is free and can adapt to the ligand.
See [Boundary conditions](../manual/reference/boundaries.md#solute-boundary-restraints).
:::

<!-- FIGURE TODO: sphere-regions
What: a cross-section through the simulation sphere of a real protein leg,
  with the two solute regions in two colours (free, excluded) and water shown
  as small spheres or a transparent surface. Annotate the sphere radius with a
  dashed circle and a label, and mark the centre.
Why: the table above is abstract; users need to see that "excluded" protein
  is still there and that water stops at the sphere.
How: PyMOL with a clipping slab (`clip slab, 6`) through the centre. Colour by
  distance of each residue from the centre (inside or outside R), or better by
  the excluded flag: write a PDB from qprep with `mask not excluded` to get the
  included atoms. The dashed circle can be added afterwards in a vector editor
  or drawn as CGO.
-->

## Solvent in the sphere

`qprep` fills the sphere with water and removes molecules that overlap the
solute. Because the solute occupies part of the volume, the radius that
gives the correct water density is not exactly the radius that was asked for.
`qprep` therefore calculates an **effective solvent radius** $r_w$ from the
number of water molecules and of solute heavy atoms in the sphere, and stores
it in the topology. $r_w$ is the radius that the boundary restraints use.

Crystallographic waters up to 2 Å outside the sphere are kept as part of the
solvent; those further out are excluded.

## Keeping water in: the radial restraint

Without a restraint, water evaporates from the surface of the droplet. A
simple wall is not enough: a half-harmonic wall alone gives a density that is
too low in the outermost few ångström, because nothing outside pulls on the
surface molecules as bulk water would.

Q therefore applies a potential with two parts to the oxygen of each water
molecule, as a function of its distance $r$ from the centre:

$$
V_{\mathrm{solvent}}(r) =
\begin{cases}
\tfrac{1}{2} K (r - r_0)^2 - D_e & r > r_0 \\[4pt]
D_e \left[ e^{2\alpha (r - r_0)} - 2\, e^{\alpha (r - r_0)} \right] & r \le r_0
\end{cases}
\qquad r_0 = r_w - \sqrt{\frac{k_B T}{K}}
$$

- Outside $r_0$, a harmonic wall with force constant $K$ (60 kcal·mol⁻¹·Å⁻²
  by default) pushes the molecule back in.
- Inside $r_0$, a Morse-type attraction of depth $D_e$ pulls molecules near
  the surface outwards. It decays within a few ångström and does nothing in
  the interior.
- $r_0$ is placed slightly inside $r_w$, by the thermal amplitude of a
  molecule in the harmonic wall.

$D_e$ and $\alpha$ were calibrated on pure water spheres of 12 to 30 Å so
that the density is uniform and equal to that of bulk water up to the
surface. They are computed from $r_w$; see
[Boundary conditions](../manual/reference/boundaries.md#radial-restraint)
for the functions. For a 15 Å sphere $D_e$ is about 1 kcal/mol and $\alpha$
about 0.5 Å⁻¹.

<!-- FIGURE TODO: radial-potential-and-density
What: two stacked panels sharing the x axis (distance from the centre, 10 to
  27 A for r_w = 25 A).
  Top: V_solvent(r) with the two branches in one line, r_0 and r_w marked by
  vertical dashed lines, D_e marked as the well depth. Add the half-harmonic
  wall alone as a thin grey line for comparison.
  Bottom: water number density in radial shells from a simulation, with and
  without the attractive term, and the bulk value 0.0335 A^-3 as a dotted
  line (the analogue of Fig. 1a of Marelius et al. 1999).
Why: connects the functional form to its purpose, a flat density profile.
How: top panel from the formula with matplotlib (K = 60, D_e and alpha from
  the functions on the boundaries page). Bottom panel needs two short water
  sphere simulations, one with `morse_depth 0`; compute the density from the
  trajectory with qcalc `rdf` or a script. If the simulations are not
  available, show only the top panel.
-->

## Orienting water at the surface: the polarisation restraint

The second problem of a surface is orientation. In bulk water with no field,
molecular dipoles point in all directions with equal probability. At the
surface of a droplet they do not: molecules turn their dipoles along the
surface to keep their hydrogen bonds, and the surface layer becomes ordered.
A water layer that is ordered in the wrong way produces an electrostatic
potential in the interior that is not there in reality.

Q restrains the orientation of the surface water to what it would be in
bulk. This is the surface-constraint idea of the SCAAS model of Warshel and
King, in the implementation described by Marelius et al.

For each water molecule, $\theta$ is the angle between its dipole vector and
the radial direction, pointing outwards from the centre.

1. **Shells.** The outermost 3 Å of the sphere is divided into three shells
   of 0.5, 1.0 and 1.5 Å thickness, counted from the outside.
2. **Target distribution.** For freely rotating dipoles, $\cos\theta$ is
   uniformly distributed, that is, $p(\theta) = \tfrac{1}{2}\sin\theta$. For
   a shell with $N$ molecules the target angles are therefore
   $\cos\theta^{0}_{k} = 1 - (2k - 1)/N$ for $k = 1, \dots, N$.
3. **Sorting.** At every step the molecules of a shell are sorted by their
   current $\theta$. The molecule with the $k$-th smallest angle is assigned
   the $k$-th target.
4. **Restraint.** Each molecule is turned towards its target by

   $$
   V_{\mathrm{pol}} = \tfrac{1}{2} K_{\mathrm{pol}} \left( \theta - \theta^{\mathrm{target}} \right)^2 ,
   $$

   with $K_{\mathrm{pol}}$ = 20 kcal·mol⁻¹·rad⁻² by default.

The restraint acts on the distribution of the shell, not on particular
molecules: a molecule is free to rotate, and when it does, it exchanges its
target with another molecule.

### A charged sphere

If the sphere contains a net charge, bulk water around it would not be
randomly oriented. The charge polarises it, and a continuum with dielectric
constant $\varepsilon$ has, at distance $r$ from a charge $Q$, the
polarisation

$$
P(r) = \frac{Q}{4\pi r^2}\left(1 - \frac{1}{\varepsilon}\right) .
$$

Dividing by the number density $\rho$ and the dipole moment $\mu$ of water
gives the average orientation that the molecules of a shell should have:

$$
\langle \cos\theta \rangle = \frac{P(r)}{\rho\,\mu} .
$$

To first order in this small number, the target distribution is shifted to

$$
\theta^{\mathrm{target}}_{k} = \theta^{0}_{k} - \tfrac{3}{2}\,\langle\cos\theta\rangle\,\sin\theta^{0}_{k} .
$$

`qdyn` applies this shift when `charge_correction` is on, which is the
default whenever the polarisation restraint is on. Three points about what it
uses:

- $Q$ is the net charge of the **Q-atoms**, weighted by the state weights of
  the window: $Q = \sum_s \lambda_s \sum_i q_i^{(s)}$. `qdyn` also adds up
  the charge of the other atoms in the sphere and prints it as the total
  charge of the system, but that total is not used for the target.
- The charge is treated as if it sat at the centre of the solvent sphere:
  only its magnitude and the radius of each shell enter.
- With `charge_correction off`, $Q$ is set to zero and the neutral target is
  used even for a charged ligand.
- $\varepsilon$ is 80, $\rho$ is 0.0335 Å⁻³, and $\mu$ is the dipole moment
  of the water model in the topology.

Finally, `qdyn` corrects for the difference between the restrained and the
achieved orientation. Every 100 steps it compares the average $\theta$ of
each shell with the target average and adds the difference to an offset that
is applied to all targets of that shell. The offsets are saved in the restart
file, so a continued run starts from the adapted values.

<!-- FIGURE TODO: polarisation-distribution
What: probability distribution of the angle theta (0 to 180 degrees) for
  water in the outermost 3 A of a pure water sphere: (1) with polarisation
  restraints, (2) without, (3) the ideal curve 0.5*sin(theta). This is the
  analogue of Fig. 1b of Marelius et al. 1999. Optional second panel: the
  same for a sphere with a +1 ion at the centre, with the shifted ideal
  curve, to show the Born shift.
Why: shows that an unrestrained surface is over-ordered (peaked near 90
  degrees) and that the restraint restores the bulk distribution.
How: two short simulations of a 15-25 A TIP3P sphere, `polarization on` and
  `polarization off`. Compute theta per water per frame from the DCD with a
  script (dipole = H1 + H2 - 2*O, radial = O - centre), keep waters with
  r > r_w - 3, histogram. Plot with matplotlib, three lines, direct labels.
-->

<!-- FIGURE TODO: surface-water-orientation
What: a rendering of the outer water layer of a sphere, each water drawn with
  an arrow for its dipole, coloured by cos(theta) on a diverging scale
  (blue = pointing in, red = pointing out, white = tangential).
Why: makes "polarisation of the surface water" concrete for readers who have
  never thought about it; pairs with the distribution plot above.
How: PyMOL with CGO arrows generated by a script from one frame of the
  simulations used above. A thin slab through the centre keeps it readable.
-->

## What the model assumes

- **The interesting events are local.** Motions that extend beyond the
  sphere are suppressed by the fixed scaffold. Loop rearrangements or domain
  motions that reach the boundary cannot be sampled.
- **The starting structure is right outside the sphere.** Excluded atoms
  never move.
- **The boundary is far enough away.** The surface water is restrained to
  behave like bulk, but the outer few ångström of any sphere are not a
  reliable environment. The ligand and the residues it interacts with should
  be well inside.
- **Volume is fixed.** There is no pressure coupling; the density is set by
  the effective radius and the radial restraint.
- **What is outside is a continuum.** This matters for charges and is the
  subject of the [next page](electrostatics.md).

For free-energy differences, a further practical rule follows: the systems
that are compared (for example the ligand in the protein and in water) should
use spheres of the same radius. See
[Outside the sphere](electrostatics.md#outside-the-sphere) for when this and
the net charge of the spheres matter.

## Settings

| Setting | Where | Default | QligFEP |
|---|---|---|---|
| Sphere centre and radius | `boundary` in `qprep` | none | centre given to `qprep_prot`, for example the centre of geometry of the ligand from `qcog`; radius 25 Å |
| `shell_radius` | `[sphere]` | 0.85 × radius | the sphere radius (no restrained shell) |
| `radial_force` ($K$) | `[solvent]` | 60 kcal·mol⁻¹·Å⁻² | 60 |
| `morse_depth`, `morse_width` ($D_e$, $\alpha$) | `[solvent]` | functions of $r_w$ | defaults |
| `polarization` | `[solvent]` | on | on |
| `polarization_force` ($K_{\mathrm{pol}}$) | `[solvent]` | 20 kcal·mol⁻¹·rad⁻² | 20 |
| `charge_correction` | `[solvent]` | follows `polarization` | default |

The keywords are described in
[qdyn input file](../manual/reference/qdyn-input.md) and the equations in
[Boundary conditions](../manual/reference/boundaries.md).

## References

- Warshel, A.; King, G. Polarization constraints in molecular dynamics
  simulation of aqueous solutions: the surface constraint all atom solvent
  (SCAAS) model. *Chem. Phys. Lett.* **1985**, *121*, 124–129.
- King, G.; Warshel, A. A surface constrained all-atom solvent model for
  effective simulations of polar solutions. *J. Chem. Phys.* **1989**, *91*,
  3647–3661.
- Essex, J. W.; Jorgensen, W. L. An empirical boundary potential for water
  droplet simulations. *J. Comput. Chem.* **1995**, *16*, 951–972.
- Marelius, J.; Kolmodin, K.; Feierberg, I.; Åqvist, J. Q: a molecular
  dynamics program for free energy calculations and empirical valence bond
  simulations in biomolecular systems. *J. Mol. Graph. Model.* **1998**,
  *16*, 213–225.
