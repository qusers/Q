# Running dynamics with qdyn

With a topology from `qprep` you can run MD, FEP and EVB simulations with
`qdyn`, with a spherical or a periodic boundary, serially or in parallel. Two
input files define a run:

1. the **`qdyn` input file**, which sets the time step, temperature, cut-offs,
   restraints and so on;
2. the **FEP file**, which redefines the topology for selected atoms. It is
   needed for FEP and EVB calculations and is described in
   [FEP files](fep-files.md).

## Simulation procedure

The MD simulation behind a free-energy calculation usually has several
stages. The first is run at a very low temperature with strong coupling to
the temperature bath, similar to an energy minimisation, to relax strain in
the initial structure. Then follow stepwise heating and equilibration at the
target temperature. After this comes the main simulation, during which
energies and structures are collected. For a perturbation, this phase is a
series of simulations on intermediate potentials defined by different sets of
weights for the FEP states.

Each sub-simulation has its own input file. The name of the input file is the
only command-line argument of `qdyn`; the input cannot be redirected with `<`.
A shell script that runs the stages in order looks like this:

```bash
qdyn relax.inp > relax.log
qdyn eq1.inp   > eq1.log
qdyn eq2.inp   > eq2.log
qdyn eq3.inp   > eq3.log
qdyn data1.inp > data1.log
qdyn data2.inp > data2.log
```

## Output

- **Log.** Progress, energy summaries and temperature go to standard output,
  normally redirected to a log file.
- **Restart file.** Final coordinates and velocities, used to start the next
  sub-simulation and, after conversion with `qprep`, to view the final
  structure. It is also written if the simulation stops because forces or
  velocities become too large, so it serves a diagnostic purpose.
- **Energy file.** Energies of the Q-atoms in each FEP state, written at a
  set interval.
- **Trajectory file.** Coordinates of all or a subset of atoms, written at a
  set interval.

The restart and energy files are Fortran binary files. The trajectory is in
DCD format, which many visualisation programs read. The trajectory-reading
routines of Q cannot read DCD files written by other programs.

## The input file

The input file is a text file divided into sections. Each section starts with
a heading in square brackets. There are two kinds of section:

- sections in which each line is a keyword and a value;
- sections in which all lines have the same layout and together form a data
  set, without keywords.

Only a few sections are required and they may appear in any order. Blank lines
are allowed. Comments start with `!`, `#` or `*` and may follow the data on a
line or stand on their own line. White space within a line is not significant,
but all data of a record must be on one line. Units are based on Å, K, fs and
kcal/mol; see [Units](reference/algorithms.md#units).

The sections are introduced below. The complete list of keywords and defaults
is in [qdyn input file](reference/qdyn-input.md).

### Dynamics control

`[MD]` defines the number of steps, the step size and the temperature.

```ini
[MD]
steps               10000
stepsize            2.0
temperature         300
shake_solvent       on
shake_solute        off
lrf                 on
```

### Periodic boundary

The presence of a `[PBC]` section selects a periodic boundary.

```ini
[PBC]
rigid_box_center    on
constant_pressure   on
max_volume_displ    65
pressure            1.5
```

`rigid_box_center on` keeps the box at fixed coordinates instead of centring
it on the solute. In this example the Monte Carlo constant-pressure algorithm
is used, with a target pressure of 1.5 bar and a maximum volume displacement
of 65 Å³.

### Non-bonded interactions

Cut-off radii for different categories of atoms are set in `[cut-offs]`.

```ini
[cut-offs]
solute_solute       10
solvent_solvent     10
solute_solvent      10
q_atom              10
lrf                 10
```

`q_atom` is the cut-off for interactions between Q-atoms and other atoms. No
cut-off is applied among Q-atoms. `lrf` is the cut-off of the local reaction
field treatment of long-range electrostatics. With a periodic boundary, all
cut-offs must be less than half the shortest box length.

### Sphere

`[sphere]` holds the parameters of the spherical boundary. The most used is
the shell radius: solute atoms outside this radius are restrained to their
coordinates in the topology.

```ini
[sphere]
shell_radius        18     ! inner radius of the restrained shell
shell_force         10     ! force constant
```

### Solvent

`[solvent]` controls the solvent boundary restraints of a spherical
simulation. It is omitted with a periodic boundary. The defaults are adequate
for most simulations.

```ini
[solvent]
polarization        on     ! solvent polarisation restraints
```

Polarisation restraints must be turned off for solvents other than
three-site water.

### Intervals

`[intervals]` sets how often the non-bonded pair lists are regenerated and
how often energies and coordinates are written.

```ini
[intervals]
non_bond            25
output              5
energy              0      ! no energy file
trajectory          100
```

Here the pair lists are regenerated every 25 steps, an energy summary is
written to the log every 5 steps, no energy file is written and coordinates
are written every 100 steps.

### Trajectory atoms

To store only a subset of atoms in the trajectory, give an
[atom mask](reference/qprep.md#atom-masks) in `[trajectory_atoms]`.

```ini
[trajectory_atoms]
heavy not excluded residue 1 104
residue 105 106
residue 109
```

This selects the heavy atoms of residues 1 to 104 that are inside the
simulation sphere, and all atoms of residues 105, 106 and 109.

### Files

`[files]` names the files to read and write. A topology and a name for the
final restart file are always required. A restart file starts the run from
the final coordinates and velocities of a previous run. A perturbation needs
an FEP file. Trajectory and energy files must be named if they are to be
written.

```ini
[files]
topology            molecule.top
final               data_01.re
trajectory          data_01.dcd
energy              data_01.en
fep                 molecule.fep
```

### State weights

For a perturbation, the weights of the FEP states are given on one line under
`[lambdas]`. For 70 % of state 1 and 30 % of state 2:

```ini
[lambdas]
0.70 0.30
```

See [How Q does FEP](../concepts/fep-in-q.md) for what the weights mean.

### Restraints

Several kinds of restraint are available.

**Sequence restraints** hold a sequence of atoms at their initial coordinates
(those of the topology) with a harmonic potential. They need the first and
last atom and a force constant, and can apply to heavy atoms only or to all
atoms. Instead of restraining each atom, the set can be restrained as a whole
to its initial geometric centre; the same force is then applied to all atoms.
With a low force constant, this keeps a small solute at the centre of the
sphere without hindering its rotation.

```ini
[sequence_restraints]
21  40  5.0  0
65  72  2.0  1  1
```

Atoms 21 to 40 are restrained to their initial positions with
5.0 kcal·mol⁻¹·Å⁻², hydrogens excepted (`0`). Atoms 65 to 72, including
hydrogens (`1`), are restrained as a group to their initial geometric
centre (`1`).

**Atom restraints** hold single atoms at a target position given in the
input. The restraint can apply in one FEP state, in which case the force is
scaled by the weight of that state, or in all states. Separate force constants
are given for x, y and z; setting one or two of them to zero restrains the
atom to a plane or a line.

```ini
[atom_restraints]
8  82.5  28.32  72.6  5.  5.  5.  0
```

Atom 8 is restrained to (82.5, 28.32, 72.6) with 5.0 kcal·mol⁻¹·Å⁻² along
all axes, in all FEP states (`0`).

**Distance restraints** restrain the distance between two atoms with a
harmonic or a flat-bottomed harmonic potential.

```ini
[distance_restraints]
13  20  4.5  5.0  10.0  1
```

Atoms 13 and 20 are held by a potential that is zero between 4.5 and 5.0 Å
and harmonic with 10.0 kcal·mol⁻¹·Å⁻² outside. It is active in FEP state 1
only.

**Wall restraints** apply a half-harmonic or Morse potential outside a given
distance from the solvent centre. They are used with a spherical boundary to
stop a molecule from moving towards the edge.

```ini
[wall_restraints]
80   99   14.0  5.0  0  0  0
102  102  14.0  5.0  0  0  0
```

Atoms 80 to 99 and atom 102 feel an inward harmonic force beyond 14 Å from
the sphere centre, with 5.0 kcal·mol⁻¹·Å⁻². The next two columns are the
depth $D_e$ and the exponent $a$ of a Morse potential, here unused. The last
column says whether hydrogens are included (`0` = no).

**Angle restraints** apply a harmonic potential to the angle between three
atoms. This is useful, for example, to stop an aspartate bound to a metal
centre from becoming bidentate.

```ini
[angle_restraints]
133  132  3400  180.0  3.0  2
```

The angle between atoms 133, 132 and 3400 is restrained to 180° with a force
constant of 3.0, in state 2.

## Examples

A minimal input, using defaults for everything optional:

```ini
[MD]
steps                 2000        ! number of steps
stepsize              1.0         ! step size (fs)
temperature           1           ! temperature (K)
initial_temperature   1           ! temperature for the initial velocities (K)

[files]
topology              molecule.top
final                 molecule.re
```

A more elaborate input:

```ini
[MD]
steps                 10000
stepsize              2.0
temperature           300
bath_coupling         10          ! temperature bath relaxation time (fs)
random_seed           57643       ! seed for the initial velocities
initial_temperature   300
shake_solvent         on          ! constrain bonds and angles of water
shake_hydrogens       on          ! constrain bonds to hydrogen
lrf                   on          ! local reaction field beyond the cut-off

[cut-offs]
solute_solute         10
solvent_solvent       10
solute_solvent        10
q_atom                10

[sphere]
shell_radius          18
shell_force           10.0

[solvent]
radial_force          60.0        ! radial restraint (kcal/mol/Å²)
polarization          on
polarization_force    20.0        ! polarisation restraint (kcal/mol/rad²)

[intervals]
non_bond              25
output                5
energy                10
trajectory            100

[trajectory_atoms]
heavy not excluded residue 1 104
residue 105 106
residue 109

[files]
topology              molecule.top
final                 data_01.re
trajectory            data_01.dcd
energy                data_01.en
fep                   molecule.fep

[lambdas]
0.70 0.30

[sequence_restraints]
21  40  5.0  0
65  72  2.0  1  1

[atom_restraints]
8  2.5  8.3  7.6  5.  5.  5.  0

[distance_restraints]
13  20  4.5  5.0  10.0  1

[wall_restraints]
80   99   14.0  5.0  0  0  0
102  102  14.0  5.0  0  0  0
```

## Running in parallel

`qdynp` is the MPI version of `qdyn`. It needs an MPI installation and is run
through the MPI launcher:

```bash
mpirun -np 8 qdynp md.inp > md.log
```

The log of the parallel version identifies itself as such in the version line
at the top. The parallel version can run on a single node but still needs
MPI. The 2017 manual recommends 2 to 12 processes, depending on the size of
the system.

## Checking an input

`qdum` is built from the same source as `qdyn` but has no dynamics loop. It
reads all input, initialises the simulation, writes a restart file and stops.
Use it to test an input quickly:

```bash
qdum md.inp > check.log
```
