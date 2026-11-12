# Units and algorithms

## Units

| Quantity | Unit |
|---|---|
| time | fs |
| temperature | K |
| length | Å |
| energy | kcal·mol⁻¹ |
| charge | *e* |

Other units are combinations of these, for example kcal·mol⁻¹·Å⁻² for a
harmonic distance force constant. Angles in input files are in degrees;
angle force constants are in kcal·mol⁻¹·rad⁻².

## Integration

`qdyn` integrates the equations of motion with the
leapfrog algorithm. It uses first-order expansions of positions $\mathbf{r}$
and velocities $\mathbf{v}$, evaluated at times half a step apart:

$$
\mathbf{r}(t + \Delta t) = \mathbf{r}(t) + \mathbf{v}(t + \Delta t/2)\,\Delta t
$$ (leapfrog-r)

$$
\mathbf{v}(t + \Delta t/2) = \mathbf{v}(t - \Delta t/2) + \frac{\mathbf{F}(t)}{m}\,\Delta t
$$ (leapfrog-v)

The offset between velocities and positions is set when the velocities are
generated, half a step back. Leapfrog is derived from the Verlet algorithm
and is symplectic and time-reversible. Its drawback is that positions and
velocities are known at different times, which introduces an error in the
energies; for small time steps the deviation is not significant.

## Temperature control

`qdyn` uses the Berendsen thermostat [39]. The system is
coupled to a heat bath of temperature $T_0$ with a relaxation time $\tau$
(`bath_coupling` in `[MD]`). The velocities are scaled so that the
temperature relaxes according to

$$
\frac{dT}{dt} = \frac{T_0 - T}{\tau} .
$$ (berendsen)

With `separate_scaling on`, solute and solvent are coupled to the bath
separately.

## Constraints

Bond lengths, and the geometry of water, are constrained with SHAKE,
controlled by the `shake_*` keywords of `[MD]`.

