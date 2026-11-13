# Boundary conditions

## Solute boundary restraints

Solute atoms outside the simulation sphere are excluded from non-bonded
interactions and tightly restrained to their initial coordinates by a
harmonic potential with a force constant of 200 kcal·mol⁻¹·Å⁻². Their bonded
interactions are evaluated as for atoms inside the sphere.

Solute atoms in the outermost shell of the simulation sphere are also
restrained to their initial coordinates with a harmonic potential, to avoid
distorting bonds across the boundary. The inner radius of this shell and the
force constant are set in the `qdyn` input (`[sphere]`, keywords
`shell_radius` and `shell_force`). The force constant defaults to
10 kcal·mol⁻¹·Å⁻². The shell radius does not change the size of the sphere; it only selects
the solute inside the sphere that is restrained. If the keyword is absent it
defaults to 85 % of the exclusion radius; a value of 1 or less is read as a
fraction of the exclusion radius. Setting it equal to the exclusion radius,
as QligFEP does, leaves all solute inside the sphere unrestrained. Atoms are assigned to the
excluded region and to the shell by charge group.

Atoms are assigned to the shell and to the excluded region before the
simulation starts, from the coordinates of the topology, and the assignment
is not updated. The coordinates used for it can instead be taken from a
separate restraint coordinate file (`[files]`, keyword `restraint`).

## Solvent boundary restraints

Solvent molecules near the boundary must be restrained to prevent
evaporation and to keep the density correct and uniform throughout the
sphere. [The spherical boundary](../../concepts/spherical-boundary.md)
explains the reasoning; this page gives the equations.

### Radial restraint

The central quantity is the **effective solvent radius** $r_w$, which is
close to but not equal to the radius of the simulation sphere. It is the
solution of

$$
V_{\mathrm{sphere}}(r_w) = \frac{4\pi}{3} r_w^3 = N_p(r_w)\,v_p + N_w\,v_w
$$

where $N_p(r_w)$ is the number of heavy solute atoms within $r_w$ of the
water centre, $v_p$ is the average volume per heavy atom in proteins
(17.3 Å³), $N_w$ is the fixed number of solvent molecules and $v_w$ is the
volume of a solvent molecule (29.9 Å³ for water).

The radial potential has a half-harmonic term that pushes solvent molecules
back into the sphere and a Morse-like term that pulls molecules from inside
towards the boundary [38]:

$$
V_{\mathrm{solvent}}(r) =
\begin{cases}
\tfrac{1}{2} K (r - r_0)^2 - D_e & \text{if } r > r_0 \\[4pt]
D_e \left[ \left( e^{\alpha (r - r_0)} \right)^2 - 2\, e^{\alpha (r - r_0)} \right] & \text{otherwise}
\end{cases}
$$ (solvent-radial)

where $r$ is the distance of the switching atom of the solvent molecule (the
oxygen of water) from the water centre, $K$ is the force constant
of the half-harmonic potential, $D_e$ is the depth of the Morse potential and
$\alpha$ is its exponent. $r_0$ is the effective solvent radius minus the
average deviation from the minimum of the half-harmonic potential at the
temperature $T$:

$$
r_0 = r_w - \sqrt{\frac{k_B T}{K}} .
$$

The values of $D_e$ and $\alpha$ depend on $r_w$. They are calculated from
empirical functions calibrated for water spheres of 12 to 30 Å:

$$
D_e(r_w) = 0.26\ \tfrac{\mathrm{kcal}}{\mathrm{mol}} \cdot e^{-0.19\,\text{Å}^{-1} (r_w - 15\,\text{Å})} + 0.74\ \tfrac{\mathrm{kcal}}{\mathrm{mol}}
$$ (solvent-de)

$$
\alpha(r_w) = \frac{0.20\ \text{Å}^{-1}}{1 + e^{0.4\,\text{Å}^{-1} (r_w - 25\,\text{Å})}} + 0.30\ \text{Å}^{-1}
$$ (solvent-alpha)

$K$, $D_e$ and $\alpha$ can be set in the input to override these values
(`[solvent]`, keywords `radial_force`, `morse_depth` and `morse_width`).

Water molecules of the topology that are outside the simulation sphere at the
start are excluded from non-bonded interactions and restraints.

### Polarisation restraint

Polar solvent molecules near the boundary are not randomly oriented as in
the bulk. A restraining force makes the surface molecules follow the
distribution of angles between the radial axis and the dipole vector that is
found in bulk solvent.

When the Q-atoms carry a net charge, it polarises the solvent and the
distribution of dipole angles changes. This shift of the average
polarisation, given by Born's formula, is taken into account unless
`charge_correction` is set to `off`. The charge used is the net charge of the
Q-atoms, weighted by the state weights; see
[A charged sphere](../../concepts/spherical-boundary.md#a-charged-sphere).

The restraints are applied in three thin shells, to minimise the radial
dependence of the polarisation that would occur in one thicker shell. The
outermost shell is 0.5 Å thick, the second 1.0 Å and the third 1.5 Å, so the
polarisation is restrained in the outermost 3 Å of the sphere.

In each shell the molecules are sorted by the angle between dipole vector
and radial axis. A force is applied to each molecule $i$ to move its angle
towards the angle of molecule $i$ in a sorted sequence that follows the
target distribution:

$$
V_{\mathrm{polarisation}}(\Theta_i) = \tfrac{1}{2} K_{\mathrm{pol}} \left( \Theta_i - \Theta_i^{\mathrm{target}} \right)^2
$$ (solvent-pol)

$K_{\mathrm{pol}}$ is set with `polarization_force` in `[solvent]`; the
default is 20 kcal·mol⁻¹·rad⁻². Every 100 steps the targets of each shell
are offset by the difference between the achieved and the target average
angle; the offsets are stored in the restart file.

## Periodic boundary conditions

With periodic boundary conditions (PBC) none of the restraints above are
applied and no atoms are excluded. This gives more flexibility in the system
but increases the run time.

- The box may be cubic or rectangular.
- The cut-off must never exceed half the shortest side of the box. This also
  holds during a constant-pressure run. If the cut-off is too large, the
  program stops.
- Particles are moved across box boundaries as whole molecules.
- With `rigid_box_center off`, the box is centred in each step on the
  geometric centre of the solute, or of the solvent if there is no solute.
  Use this for a large protein in water.
- With `rigid_box_center on`, the centre of the box stays as given in the
  topology. If solute is present, each solute molecule must then be a single
  charge group, which is done by changing the library file.

## Constant pressure algorithm

The constant-pressure algorithm combines molecular dynamics with Monte Carlo
sampling of the volume. A volume change is chosen at random,

$$
\Delta V = n_{\mathrm{rand}} \cdot \Delta V_{\max}
$$

where $n_{\mathrm{rand}}$ is a random number between −1 and 1 and
$\Delta V_{\max}$ is the maximum volume displacement of one move. The new
volume is $V' = V + \Delta V$. The side lengths $l_x, l_y, l_z$ are scaled
by

$$
l_i' = l_i \sqrt[3]{\frac{V'}{V}} ,
$$

so the proportions of the box are kept. The centre of mass of each molecule
is moved to

$$
r_i' = (r_i - c_i)\,\frac{l_i'}{l_i} + c_i
$$

where $c_i$ is the coordinate of the box centre, which need not be the
origin. Molecules are moved as a whole, so intramolecular distances stay
fixed.

The potential energy is then recalculated. Only the non-bonded interactions
are needed, because the interior of the molecules is unchanged. With

$$
\Delta W = (U'_{\mathrm{pot}} - U_{\mathrm{pot}}) + P_0 (V' - V) ,
$$

where $P_0$ is the target pressure, the new configuration is accepted with
probability

$$
P(\Delta V) =
\begin{cases}
e^{-\Delta W / k_B T} & \Delta W > 0 \\
1 & \Delta W \le 0 .
\end{cases}
$$ (npt-accept)

If $\Delta W$ is zero or negative the move is always accepted. Otherwise a
random number $n \in [0, 1]$ is drawn and the move is accepted if
$n \le e^{-\Delta W / k_B T}$. The acceptance ratio is controlled by
`max_volume_displ`, which is $\Delta V_{\max}$.
