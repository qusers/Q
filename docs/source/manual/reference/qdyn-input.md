# qdyn input file

The name of the input file is the first argument on the `qdyn` command line.
The file is divided into sections, each starting with a heading in square
brackets. A section holds either lines of a keyword and a value, or lines of
values in a fixed order. Comments start with `!`, `#` or `*` and may follow
the data on a line or stand on their own line. The order of sections, and of
keywords within a section, is not significant. Many keywords have defaults;
`qdyn` reports in its log when it uses one. Sections without required
keywords are optional.

:::{warning}
`qdyn` ignores keywords it does not know, without an error. A misspelt
keyword silently leaves the default in place. Check the log for the values
that were actually used. Note the spelling of `polarization` and
`rigid_box_center`. See [Changes from the Q6 manual](../differences.md) if
your inputs were written for Q6.
:::

## Keyword sections

### MD

Basic simulation data. Required.

| Keyword | Value | Comment |
|---|---|---|
| `steps` | Number of MD steps | Required |
| `stepsize` | Step size $\Delta t$ (fs) | Required |
| `temperature` | Target temperature (K) | Required |
| `bath_coupling` | Relaxation time $\tau$ of the temperature bath (fs) [39] | Optional, default 10 fs. Must be at least `stepsize`. |
| `separate_scaling` | Separate temperature scaling of solute and solvent (`on`/`off`) | Optional, default `off` |
| `random_seed` | Integer seed of the random number generator | Optional, default 1. Used only for random initial velocities. |
| `initial_temperature` | Temperature (K) of the Maxwell distribution for random velocities | Required when no restart file is given. Give it only when velocities should be randomised. |
| `shake_solvent` | Constrain bonds and angles of the solvent (`on`/`off`) | Optional, default `on` |
| `shake_solute` | Constrain all solute bonds (`on`/`off`) | Optional, default `off` |
| `shake_hydrogens` | Constrain bonds to hydrogen atoms (`on`/`off`) | Optional, default `off` |
| `lrf` | Local reaction field: Taylor expansion of the electrostatic field beyond the cut-off (`on`/`off`) | Optional, default `off` |
| `force_rms` | Print the RMS force (`on`/`off`) | Optional, default `off` |
| `minimize` | Run an energy minimisation before the dynamics (`on`/`off`) | Optional, default `off` |
| `max_minimize_steps` | Maximum number of minimisation steps | Optional, default 1000 |
| `minimize_tolerance` | Convergence tolerance of the minimisation | Optional, default 0.1 |
| `minimize_step_size` | Initial step size of the minimisation | Optional, default 0.001 |

### PBC

Settings for periodic boundary conditions. The presence of the section
selects a periodic boundary.

| Keyword | Value | Comment |
|---|---|---|
| `rigid_box_center` | Keep the box centre fixed (`on`) or centre the box on the solute in each step (`off`) | Optional, default `off` |
| `constant_pressure` | Simulate in the isothermal–isobaric ensemble (`on`/`off`). The trial interval is set in `[intervals]`. | Optional, default `off` |
| `max_volume_displ` | Maximum volume change of one Monte Carlo step (Å³) | Required when `constant_pressure` is `on` |
| `pressure_seed` | Seed of the random number generator for the volume moves | Optional, default 3781. Change it when a constant-pressure simulation is split over several input files. |
| `pressure` | Target pressure (bar) | Optional, default 1.0. Used only with `constant_pressure on`. |
| `control_box` | Three box lengths (Å) to change the box to | Optional. Used only with `constant_pressure off`. |
| `put_solute_back_in_box` | Put solute molecules back in the box (`on`/`off`). Does not affect energies. | Optional. Turn off when measuring, for example, diffusion coefficients. |
| `put_solvent_back_in_box` | The same for solvent molecules | Optional |

### cut-offs

Cut-off radii for non-bonded interactions.

| Keyword | Value | Comment |
|---|---|---|
| `solute_solute` | Cut-off (Å) for solute–solute pairs | Optional, default 10 Å |
| `solvent_solvent` | Cut-off (Å) for solvent–solvent pairs | Optional, default 10 Å |
| `solute_solvent` | Cut-off (Å) for solute–solvent pairs | Optional, default 10 Å |
| `q_atom` | Cut-off (Å) for the interaction of Q-atoms with all other atoms | Optional, default 99 Å, that is, no cut-off. With PBC, compare with the box length. |
| `lrf` | Cut-off (Å) for the LRF expansion | Optional, default 99 Å. Read only with `lrf on`. Must not be smaller than the other cut-offs. |

### sphere

Not used with PBC.

| Keyword | Value | Comment |
|---|---|---|
| `shell_radius` | Inner radius (Å) of the restrained shell | Optional, default 85 % of the exclusion radius |
| `shell_force` | Force constant (kcal·mol⁻¹·Å⁻²) of the shell restraints | Optional, default 10 |
| `exclude_bonded` | Eliminate (`on`) or retain (`off`) bonded interactions between excluded atoms | Optional, default `off` |

`center` and `radius` are obsolete here; the sphere is defined in the
topology, and `qdyn` ignores them with a warning.

### solvent

Boundary conditions of the solvent sphere. Optional. Omit it for vacuum
simulations; it is not read with PBC. The old section name `[water]` is also
accepted.

| Keyword | Value | Comment |
|---|---|---|
| `radius` | Target solvent radius (Å) | Optional. Use only to override the target radius calculated from the topology. |
| `radial_force` | Force constant (kcal·mol⁻¹·Å⁻²) of the half-harmonic radial restraint | Optional, default 60 |
| `polarization` | Polarisation restraints in the outer solvent shell (`on`/`off`) | Optional, default `on` |
| `charge_correction` | Correct the polarisation restraints for the total charge of the Q-atoms by Born's formula (`on`/`off`) | Optional. Defaults to the value of `polarization`; `on` requires `polarization on`. |
| `polarization_force` | Force constant (kcal·mol⁻¹·rad⁻²) of the polarisation restraints | Optional, default 20 |
| `morse_depth` | Depth $D_e$ (kcal·mol⁻¹) of the Morse-type boundary attraction | Optional, default a function of the solvent radius, {eq}`solvent-de` |
| `morse_width` | Exponent $\alpha$ (Å⁻¹) of the Morse-type boundary attraction | Optional, default a function of the solvent radius, {eq}`solvent-alpha` |

`center`, `pack` and `model` are obsolete and ignored with a warning.

### intervals

Intervals, in steps, between saving data and updating the non-bonded lists.

| Keyword | Value | Comment |
|---|---|---|
| `non_bond` | Update interval of the non-bonded lists (and the LRF data) | Optional, default 10 |
| `output` | Interval for printing energy summaries | Optional, default 10 |
| `temperature` | Interval for printing the temperature | Optional, default: the `output` interval |
| `energy` | Interval for writing Q-atom energies to the energy file | Optional, default 0 (off) |
| `trajectory` | Interval for writing coordinates to the trajectory file | Optional, default 0 (off) |
| `volume_change` | Interval for Monte Carlo volume trials | Optional, default 10. Read only with `constant_pressure on`. |

### files

File names for input and output. Required.

| Keyword | Value | Comment |
|---|---|---|
| `topology` | Topology file | Required |
| `restart` | Restart file from which initial coordinates and velocities are read | Optional. Without it, coordinates are taken from the topology and random velocities are generated (see `[MD]`). |
| `final` | Restart file to which the final coordinates and velocities are written | Required |
| `trajectory` | Trajectory file | Required when the trajectory interval is above 0 |
| `energy` | Energy file | Required when the energy interval is above 0 |
| `fep` | FEP file | Required when λ values are given |
| `restraint` | Restart file with the coordinates to restrain to | Optional. Replaces the coordinates of the topology as the reference. The restraints themselves must still be specified, for example in `[sequence_restraints]`. |

### lambda_scaling

Selects non-equilibrium switching, in which the state weights change during
the run. Without this section the weights are constant.

| Keyword | Value | Comment |
|---|---|---|
| `scaling_parameter` | Shape of the switching schedule | Optional, default `sigmoidal` |
| `L_sigmoid` | Steepness of the sigmoidal schedule | Optional, default 8.0 |

See [Non-equilibrium FEP](../../qligfep/neq.md).

## Data sections

The following sections hold data in columns, without keywords.

### trajectory_atoms

An [atom mask](qprep.md#atom-masks) that selects the atoms written to the
trajectory. Several lines may be used.

### lambdas

| Column | Description |
|---|---|
| 1… | Weight of state 1, state 2, … Each $\lambda_i \in [0, 1]$ and $\sum_i \lambda_i = 1$. |

### sequence_restraints

Restrain sequences of atoms.

| Column | Description |
|---|---|
| 1 | Number of the first atom of the sequence |
| 2 | Number of the last atom of the sequence |
| 3 | Force constant (kcal·mol⁻¹·Å⁻²) of the harmonic potential |
| 4 | Restrain hydrogens too (0 = no, 1 = yes) |
| 5 | Reference: each atom to its initial coordinates (0 or absent), the sequence to its geometric centre (1), or the sequence to its centre of mass (2) |

### atom_restraints

Restrain the positions of single atoms.

| Column | Description |
|---|---|
| 1 | Atom number |
| 2, 3, 4 | Target x, y and z coordinates (Å) |
| 5, 6, 7 | Force constants (kcal·mol⁻¹·Å⁻²) in the x, y and z directions. Separate constants allow an atom to be restrained to a line or a plane. |
| 8 | FEP state in which the restraint is active (energies and forces are scaled by its λ), or 0 for all states |

### distance_restraints

Restrain atom–atom distances.

| Column | Description |
|---|---|
| 1 | Number of the first atom |
| 2 | Number of the second atom |
| 3 | Lower limit (Å) of the unrestrained region |
| 4 | Upper limit (Å) of the unrestrained region. Set lower = upper for a standard harmonic potential. |
| 5 | Force constant (kcal·mol⁻¹·Å⁻²) |
| 6 | FEP state in which the restraint is active (energies and forces are scaled by its λ), or 0 for all states |

Atoms can also be given as `residue:atom`.

### angle_restraints

Restrain the angle between three atoms.

| Column | Description |
|---|---|
| 1, 2, 3 | Numbers of the three atoms; the second is the central atom |
| 4 | Target angle (°) |
| 5 | Force constant (kcal·mol⁻¹·rad⁻²) |
| 6 | FEP state in which the restraint is active, or 0 for all states |

### wall_restraints

Half-harmonic (elastic wall) restraints of sequences of atoms. Not used with
PBC.

| Column | Description |
|---|---|
| 1 | Number of the first atom of the sequence |
| 2 | Number of the last atom of the sequence |
| 3 | Distance (Å) from the water centre beyond which the potential applies |
| 4 | Force constant (kcal·mol⁻¹·Å⁻²) of the harmonic potential |
| 5 | $D_e$ of the Morse potential, the depth of its minimum |
| 6 | Exponent $a$ of the Morse potential |
| 7 | Restrain hydrogens too (0 = no, 1 = yes) |
