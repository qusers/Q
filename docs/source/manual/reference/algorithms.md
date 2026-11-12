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

`qdyn` in this repository integrates the equations of motion with the
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

`qdyn` in this repository uses the Berendsen thermostat [39]. The system is
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

## Algorithms of Q6 that are not in this source

The 2017 manual also describes the velocity Verlet integrator and the
Nosé–Hoover chain and Langevin thermostats, selected with the `integrator`
and `thermostat` keywords. `qdyn` in this repository does not read these
keywords. Their equations are reproduced here for reference only.

**Velocity Verlet** evaluates positions and velocities at the same time:

$$
\mathbf{r}(t + \Delta t) = \mathbf{r}(t) + \mathbf{v}(t)\,\Delta t + \frac{\mathbf{F}(t)}{2m}\,\Delta t^2
$$

$$
\mathbf{v}(t + \Delta t) = \mathbf{v}(t) + \frac{\mathbf{F}(t + \Delta t) + \mathbf{F}(t)}{2m}\,\Delta t
$$

**Nosé–Hoover.** An extra degree of freedom $s$ with momentum $p_s$ and mass
$Q$ is coupled to the system. For $N$ particles with $g$ degrees of freedom
the Hamiltonian is

$$
\mathcal{H}_{\mathrm{Nos\acute{e}}}(\mathbf{p}, \mathbf{r}, p_s, s) =
\sum_{i=1}^{N} \frac{\mathbf{p}_i^2}{2 m_i s^2} + V(\mathbf{r}^N) + \frac{p_s^2}{2Q} + g k_B T \ln s
$$

with the equations of motion

$$
\frac{d\mathbf{r}_i}{dt} = \frac{\mathbf{p}_i}{m_i s^2}, \qquad
\frac{d\mathbf{p}_i}{dt} = \mathbf{F}_i, \qquad
\frac{ds}{dt} = \frac{p_s}{Q}, \qquad
s\,\frac{dp_s}{dt} = \sum_{i=1}^{N} \frac{\mathbf{p}_i^2}{m_i s^2} - g k_B T .
$$

With $\mathbf{p}' = \mathbf{p}/s$ and

$$
\mathcal{H}(\mathbf{p}', \mathbf{r}) = \sum_{i=1}^{N} \frac{\mathbf{p}_i'^2}{2 m_i} + V(\mathbf{r}^N) ,
$$

the microcanonical average of an observable in the extended system equals
its canonical average in the physical system (see [42]):

$$
\langle A(\mathbf{p}/s, \mathbf{r}) \rangle_{NVE} = \langle A(\mathbf{p}', \mathbf{r}) \rangle_{NVT} .
$$

A single thermostat variable does not sample the canonical ensemble properly
in some cases. The remedy is a chain of $C$ thermostats, each coupled to the
previous one [40]:

$$
\begin{aligned}
\frac{d\mathbf{r}_i}{dt} &= \frac{\mathbf{p}_i}{m_i} \\
\frac{d\mathbf{p}_i}{dt} &= \mathbf{F}_i - \frac{p_{s_1}}{Q_1}\,\mathbf{p}_i \\
\frac{ds_j}{dt} &= \frac{p_{s_j}}{Q_j} \\
\frac{dp_{s_1}}{dt} &= \left( \sum_{i=1}^{N} \frac{\mathbf{p}_i^2}{m_i} - g k_B T \right) - \frac{p_{s_2}}{Q_2}\,p_{s_1} \\
\frac{dp_{s_j}}{dt} &= \left( \frac{p_{s_{j-1}}^2}{Q_{j-1}} - k_B T \right) - \frac{p_{s_{j+1}}}{Q_{j+1}}\,p_{s_j} \\
\frac{dp_{s_C}}{dt} &= \frac{p_{s_{C-1}}^2}{Q_{C-1}} - k_B T
\end{aligned}
$$

:::{note}
The chain equations are given here in the standard form of Martyna, Klein
and Tuckerman [40]. The 2017 manual prints them with a factor $s_j$ on the
left-hand side and with $m_{j-1}$ in place of $Q_{j-1}$.
:::

**Langevin.** A friction force removes kinetic energy and a random force
adds it:

$$
\frac{d\mathbf{p}_i}{dt} = \mathbf{F}_i - \gamma\,\mathbf{p}_i + \sqrt{2 \gamma m_i k_B T}\;\boldsymbol{\Gamma}
$$

where $\gamma$ is the friction coefficient and $\boldsymbol{\Gamma}$ is
white noise with zero mean and unit variance. The stationary solution of the
corresponding Fokker–Planck equation is the Maxwell–Boltzmann distribution.

:::{note}
The 2017 manual prints the noise amplitude as
$\sqrt{2\gamma k_B T / m_i}$. For an equation in the momentum the amplitude
that satisfies the fluctuation–dissipation relation is
$\sqrt{2 \gamma m_i k_B T}$.
:::
