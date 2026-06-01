# qdyn Process Walkthrough

This document follows the Fortran `qdyn` executable path in `src/q6/qdyn.f90`.
It does not describe the newer C++/CUDA `src/core` executable except where the
same concepts overlap.

## High-level Call Graph

```text
src/q6/qdyn.f90:18        program qdyn
  -> MPI setup or serial node setup
  -> signal handler setup
  -> src/q6/qdyn.f90:83   call startup
       -> src/q6/qdyn.f90:143  subroutine startup
       -> src/q6/md.f90:374    call md_startup
            -> src/q6/qatom.f90:180  call qatom_startup
                 -> src/q6/nrgy.f90:85  call nrgy_startup
            -> src/q6/trj.f90:56     call trj_startup
  -> master node setup:
       initialize, open files, topology, coordinates, FEP, dynamics prep
  -> optional MPI slave initialization
  -> distribute nonbonded work
  -> src/q6/md.f90:3917  call md_run
  -> close output files
  -> shutdown
  -> MPI finalization
```

## 1. Program Entry

The executable starts at `src/q6/qdyn.f90:18`:

```fortran
program qdyn
```

The main program body runs from `program qdyn` until `contains` at
`src/q6/qdyn.f90:138`. Routines after `contains`, such as `startup`,
`shutdown`, and `commandlineoptions`, are internal subroutines and only execute
when called.

At the top of the program, `qdyn` imports:

- `iso_fortran_env` for compiler/version environment support.
- `md` from `src/q6/md.f90`, which contains most MD state and algorithms.
- `mpiglob` from `src/q6/mpiglob.f90`, which contains MPI globals such as
  `nodeid`, `numnodes`, and `ierr`.

## 2. MPI or Serial Node Initialization

Code path: `src/q6/qdyn.f90:56-68`

If compiled with `USE_MPI`, `qdyn` calls:

```fortran
call MPI_Init(qdyn_ierr)
call MPI_Comm_rank(MPI_COMM_WORLD, nodeid, qdyn_ierr)
call MPI_Comm_size(MPI_COMM_WORLD, numnodes, qdyn_ierr)
```

This establishes:

- `nodeid == 0`: master node.
- `nodeid > 0`: slave/worker nodes.
- `numnodes`: total MPI process count.

If not compiled with `USE_MPI`, the run is serial:

```fortran
nodeid = 0
numnodes = 1
```

Most input parsing and output is master-only. Nonbonded force work can be
distributed to slaves in MPI builds.

## 3. Signal Handler Setup

Code path: `src/q6/qdyn.f90:70-77`, handlers at `src/q6/qdyn.f90:310-335`

`qdyn` installs handlers for:

- `SIGINT`
- `SIGKILL`
- `SIGABRT`

Each handler uses `use md` and calls `die(...)`, for example:

```fortran
call die('user request (control-C)')
```

`die` lives in `src/q6/md.f90:400`. It prints an error path, closes down as much
as possible, and exits/aborts.

## 4. qdyn Startup

Code path: `src/q6/qdyn.f90:83`

```fortran
call startup
```

The internal `startup` routine is defined at `src/q6/qdyn.f90:143`.

On the master node it prints the banner, compiler version, compiler options,
program version, and current date/time. The printed program name changes with
compile flags:

- `DUM`: qdum input checker
- `EVAL`: evaluation build
- default: qdyn

Then it initializes modules:

```fortran
call md_startup
```

## 5. md Module Startup

Code path: `src/q6/md.f90:374`

`md_startup` initializes modules used by MD:

```fortran
call qatom_startup
call trj_startup
pi = 4.0*atan(1.0)
deg2rad = pi/180.0
```

`qatom_startup` is at `src/q6/qatom.f90:180`. It calls:

```fortran
call nrgy_startup
```

`nrgy_startup` is at `src/q6/nrgy.f90:85` and is currently empty. It is only a
startup hook for the `nrgy` module.

`trj_startup` is at `src/q6/trj.f90:56` and initializes the trajectory module.

## 6. Master-only Input and System Setup

Code path: `src/q6/qdyn.f90:85-113`

Only `nodeid == 0` performs the full input and topology setup:

```fortran
if(.not. initialize()) call die('Invalid data in input file')
call open_files
call topology
call prep_coord
if ( nstates > 0 ) call get_fep
call prep_sim
call close_input_files
call init_shake
call make_nbqqlist
call shrink_topology
call nbmonitorlist
call init_trj
```

Slave nodes receive the necessary state later through `init_nodes`.

## 7. Reading the qdyn Input File

Code path: `src/q6/md.f90:2528`

```fortran
logical function initialize()
```

`initialize` reads the input file name from the command line. If no input file is
given, it dies:

```fortran
if (num_args .lt. 1) call die('no input file specified on the command line')
```

It then opens sections from the qdyn input file and fills global MD state.

Major input sections and variables:

- `PBC`: decides spherical boundary versus periodic box.
- `md`: steps, step size, target temperature, bath coupling, random seed,
  initial temperature, SHAKE flags, LRF flag, force RMS flag.
- `cut-offs`: `Rcpp`, `Rcww`, `Rcpw`, `Rcq`, and optionally `RcLRF`.
- `sphere`: restrained shell settings for spherical simulations.
- `solvent` or `water`: solvent radius, radial restraint, water polarization,
  charge correction, Morse parameters.
- `intervals`: nonbond update interval, output interval, temperature interval,
  trajectory interval, energy file interval, volume-change interval.
- `trajectory_atoms`: optional trajectory atom mask.
- `files`: topology, restart, final coordinate, trajectory, energy, FEP, and
  external restraint files.
- `lambdas`: FEP/EVB lambda values.
- restraint sections: sequence, atom, distance, angle, and wall restraints.

Important conversions:

- The input time step in femtoseconds is converted to internal units:

  ```fortran
  dt=0.020462*stepsize
  ```

- The temperature relaxation time is converted similarly:

  ```fortran
  tau_T=0.020462*tau_T
  ```

If the modern `md` section is missing, `initialize` falls back to
`old_initialize`, starting at `src/q6/md.f90:3370`, for compatibility with old
version-2 style input files.

## 8. Opening Input and Output Files

Code path: `src/q6/md.f90:1382`

`open_files` opens fixed Fortran units:

- Unit `2`: restart file, if `restart` is true.
- Unit `3`: final coordinates file, always opened for unformatted write.
- Unit `11`: energy output file, if `iene_cycle > 0`.
- Unit `12`: external implicit restraint file, if requested.

Errors call `die`.

Later, `close_input_files` at `src/q6/md.f90:1365` closes input-side units:

- `1`
- `2` if restart was used
- `12` if external restraints were used
- `13`

`close_output_files` at `src/q6/md.f90:1374` closes:

- `3`
- `10` if trajectory output was enabled
- `11` if energy output was enabled

## 9. Loading the Topology

Code path: `src/q6/md.f90:14705`

`topology` loads the topology file:

```fortran
if(.not. topo_load(top_file, require_version=4.15)) then
  call die('Failed to load topology.')
end if
```

It then:

- Sets `natom = nat_pro`.
- Computes the number of waters:

  ```fortran
  nwat = (natom - nat_solute) / 3
  ```

- Adds the final molecule-start sentinel:

  ```fortran
  istart_mol(nmol+1) = nat_pro + 1
  ```

- Converts angle/torsion/improper library data from degrees to radians.
- Inverts torsion path counts.
- Builds the `ljcod` matrix, including special Lennard-Jones pair codes.
- Applies arithmetic vdW preparation when `ivdw_rule == 2`.
- Checks that the boundary type in the topology agrees with the qdyn input.
- For PBC, computes inverse box lengths and verifies cutoffs are smaller than
  half the box length.

The topology code uses shared topology state from the `topo` module.

## 10. Preparing Coordinates

Code path: `src/q6/md.f90:13596`

`prep_coord` prepares positions and atom arrays before dynamics:

- If an external restraint coordinate file is used, it refreshes `xtop`.
- Distance restraints written as `residue:atom` are converted to atom numbers.
- For spherical boundary runs, it builds restrained shell atom lists using
  either `make_shell` or `make_shell2`.
- It calls `allocate_natom_arrays`, which allocates core arrays such as:
  `x`, `xx`, `v`, `d`, `winv`, and `iqatom`.
- If this is a restart, it reads coordinates and velocities from unit `2`.
- If PBC is active, it also reads `boxlength` and `boxcenter` from the restart.
- If not a restart, initial coordinates come from topology coordinates:

  ```fortran
  x(1:nat_pro*3) = xtop(1:nat_pro*3)
  ```

- It clears the Q-atom lookup array:

  ```fortran
  iqatom(:) = 0
  ```

## 11. Loading FEP/EVB Data

Code path: `src/q6/md.f90:1472`

`get_fep` runs only when `nstates > 0`.

It:

- Loads Q atoms from the FEP file using `qatom_load_atoms`.
- Maps topology atom numbers to Q-atom indices in `iqatom`.
- Allocates `qcrg(nqat,nstates)`.
- Copies topology charges into Q-atom charges as the initial state.
- Allocates and initializes the softcore lookup array `sc_lookup`.
- Loads the rest of the FEP file with `qatom_load_fep`.
- Adapts Q vdW parameters for arithmetic combination rules.
- Removes classical bonded terms redefined by Q bonded terms by setting their
  code to `0`.
- Applies special exclusions to the standard exclusion lists.

After this point, normal topology terms and Q-state-specific terms have been
separated enough for later energy evaluation.

## 12. Preparing Simulation Constants

Code path: `src/q6/md.f90:13787`

`prep_sim` prepares state that depends on topology, coordinates, and FEP data.

It:

- Prints the "Initializing dynamics" heading.
- Sets water oxygen/hydrogen charges for supported solvent models.
- For spherical boundary runs:
  - Calls `wat_sphere`.
  - Optionally calls `wat_shells` for water polarization restraints.
- For PBC:
  - Computes and prints total non-Q and total system charge.
- If LRF or PBC is active:
  - Allocates LRF arrays.
  - Builds `iwhich_cgp`, mapping each atom to its charge group.
- Builds inverse masses:

  ```fortran
  winv(:) = 1./iaclib(iac(:))%mass
  ```

- Optionally changes the box length from input and calls `put_back_in_box`.
- For PBC, computes molecule masses and atom masses.
- Scales all classical and Q charges by `sqrt(coulomb_constant)`:

  ```fortran
  crg(:) = crg(:) * sqrt(coulomb_constant)
  crg_ow = crg_ow * sqrt(coulomb_constant)
  crg_hw = crg_hw * sqrt(coulomb_constant)
  qcrg(:,:) = qcrg(:,:) * sqrt(coulomb_constant)
  ```

That charge scaling is why later Coulomb code can multiply scaled charges
directly instead of multiplying raw charges by the Coulomb constant every time.

## 13. Initializing SHAKE

Code path: `src/q6/md.f90:2326`

`init_shake` builds the SHAKE constraint data structure.

It:

- Allocates `shake_mol(nmol)`.
- Counts constrained bonds per molecule based on:
  - `shake_hydrogens`
  - `shake_solute`
  - `shake_solvent`
- Adds extra Q/FEP SHAKE constraints from `nqshake`.
- Allocates per-molecule constraint lists.
- Stores constrained atom pairs and target squared distances.
- Marks shaken classical bonds with code `-1`, so `shrink_topology` can remove
  them from the normal bonded list.
- Calculates degrees of freedom:
  - `Ndegf`
  - `Ndegfree`
  - solute and solvent variants
- Decides whether separate solute/solvent temperature scaling is possible.
- Removes angles where the end atoms are directly constrained by SHAKE.

If initial velocities are generated later, `initial_shaking` at
`src/q6/md.f90:2506` applies SHAKE to initial coordinates and velocities.

## 14. Building Q-Q Nonbond Lists

Code path: `src/q6/md.f90:1022`

`make_nbqqlist` prepares Q-atom to Q-atom nonbonded interactions that do not need
the normal charge-group update cycle:

```fortran
call make_qconn
nbqq_max = nbqq_count()
allocate(nbqq(nbqq_max, nstates), stat=alloc_status)
call nbqqlist
```

It prints the number of cutoff-independent nonbonded Q pairs for each state.

## 15. Shrinking the Topology

Code path: `src/q6/md.f90:14569`

`shrink_topology` removes interactions that should no longer be evaluated:

- Bonds with `cod <= 0`.
  - `0` means redefined or disabled.
  - `-1` means removed because SHAKE handles it.
- Angles with `cod == 0`.
- Torsions/impropers with `cod == 0`.
- If `exclude_bonded` is active, torsions and impropers whose atoms are all
  excluded.

This keeps later bonded loops from spending time on disabled interactions.

## 16. Nonbond Monitor Setup

Code path: `src/q6/md.f90:8769`

`nbmonitorlist` runs only when selected monitor group pairs exist.

It:

- Allocates `special_LJcod`.
- For each monitored atom pair and each FEP state, determines which LJ code
  should be used.
- Handles normal atom pairs, Q-involved pairs, 1-4 pairs, and Q-vdW special
  cases.

The actual monitor energy calculation happens inside `nonbond_monitor`.

## 17. Trajectory Setup

Code path: `src/q6/md.f90:13766`

`init_trj` runs only when `itrj_cycle > 0`.

It:

- Calls `trj_initialize` with frame count, interval, total steps, degrees of
  freedom, and topology file.
- Commits the trajectory atom mask with `trj_commit_mask`.
- Creates the trajectory file with `trj_create(trj_file)`.

Trajectory writes during the run go through `write_trj` at
`src/q6/md.f90:15442`.

## 18. Optional Initial Velocity Generation

Code path: `src/q6/qdyn.f90:106-112`

If `iseed > 0`, qdyn generates initial velocities:

```fortran
call maxwell
call initial_shaking
call stop_cm_translation
```

`maxwell` at `src/q6/md.f90:3787` generates Maxwellian velocities using atom
masses and `Tmaxw`.

`initial_shaking` applies SHAKE to the starting coordinates and velocities.

`stop_cm_translation` at `src/q6/md.f90:14670` removes center-of-mass
translation from velocities.

If `iseed == 0`, qdyn expects velocities to come from the restart file.

## 19. MPI Slave Initialization

Code path: `src/q6/qdyn.f90:116-119`

In MPI builds, if more than one node exists:

```fortran
if (numnodes .gt. 1) call init_nodes
```

`init_nodes` starts at `src/q6/md.f90:1901`.

It broadcasts from the master to slaves:

- Run parameters: atom counts, water count, steps, cutoffs, LRF flag,
  nonbond-list interval.
- PBC state: box center, box length, inverse box length, pressure flags.
- Topology state: atom types, charges, charge groups, exclusions, 1-4 lists,
  molecule starts, LJ parameters.
- Q-atom state: Q atom count, state count, Q atom mapping, Q charges, Q vdW
  parameters, lambda values, softcore lookup.
- Coordinates and velocities.

It allocates matching arrays on slave nodes and finally calls `allocate_mpi`.

Important limitation in the source comments: SHAKE is not parallelized here;
master handles the coordinate update and broadcasts coordinates afterward.

## 20. Distributing Nonbonded Work

Code path: `src/q6/qdyn.f90:124`

```fortran
call distribute_nonbonds
```

`distribute_nonbonds` starts at `src/q6/md.f90:1044`.

It:

- Counts PP, PW, QP, QW, and WW nonbonded pairs per charge group.
- Stores total pair counts in `totnbpp`, `totnbpw`, `totnbww`, `totnbqp`, and
  `totnbqw`.
- In serial mode, assigns all work to the master.
- In MPI mode, estimates bonded and nonbonded work, then assigns charge-group
  ranges to master and slaves for load balance.
- Allocates the nonbond arrays needed by each node.

This distribution defines what each node evaluates when `pot_energy_nonbonds`
runs.

## 21. Main Dynamics Entry

Code path: `src/q6/qdyn.f90:127`

```fortran
call md_run
```

The main MD routine starts at `src/q6/md.f90:3917`.

Its own comment is accurate:

```fortran
This subroutine has the main algorithms for the equations of motion.
```

## 22. md_run Setup Before the Loop

Code path: `src/q6/md.f90:3917-3988`

Before stepping:

- It sets profiling labels if profiling is enabled.
- It allocates MPI profiling arrays if needed.
- It sets `nat3 = natom*3`.
- On the master node:
  - Computes the hot-atom warning threshold energy `Ekinmax`.
  - Calls `temperature` to compute initial temperature and thermostat scaling.
  - Prints initial total/free temperatures, plus solute/solvent details when
    available.
  - Starts wall-clock timers.

`temperature` is at `src/q6/md.f90:3813`. It:

- Computes kinetic energy from velocities.
- Computes total, free, solute, solvent, and excluded temperatures.
- Updates `E%kinetic`.
- Calculates thermostat velocity scale factors.

## 23. Main MD Loop

Code path: `src/q6/md.f90:3995-4182`

Unless compiled with `DUM`, qdyn loops:

```fortran
do istep = 0, nsteps-1
```

Each time step performs the following sequence.

### 23.1 Nonbonded Pair-list Update

Code path: `src/q6/md.f90:3999-4047`

Every `NBcycle` steps:

- If PBC and requested, calls `put_back_in_box` for visualization and LRF
  consistency.
- Prints timing and estimated completion after the first update.
- Calls `make_pair_lists`.

`make_pair_lists` starts at `src/q6/md.f90:3716`. It selects the correct pair
list builders based on:

- PBC versus spherical boundary.
- LRF versus normal cutoff.
- Switch-atom mode versus charge-group-center mode.
- Geometric versus arithmetic handling is used later in energy evaluation.

It builds lists for:

- PP: protein/solute to protein/solute.
- PW: protein/solute to water.
- WW: water to water.
- QP: Q atoms to protein/solute.
- QW: Q atoms to water.

Q-Q lists are handled separately by `make_nbqqlist`.

### 23.2 Potential Energy and Force Derivatives

Code path: `src/q6/md.f90:4064`

```fortran
call pot_energy
```

`pot_energy` starts at `src/q6/md.f90:13278`.

It:

- Resets all energy accumulators.
- Resets the derivative/force array:

  ```fortran
  d(:) = 0.
  ```

- In MPI mode, prepares nonbond gather operations.
- Calls `pot_energy_nonbonds` on all nodes.
- On the master:
  - Calls `pot_energy_bonds`.
  - Applies shell restraints, sequence/position/distance/angle/wall restraints,
    solvent restraints, and water polarization restraints.
  - Computes Q-Q nonbonded interactions.
  - Computes Q bonded terms for each state.
- In MPI mode:
  - Slaves send their nonbonded energies and derivatives.
  - Master waits, sums `d`, `E`, and `EQ`.
- Builds Q-energy totals per state.
- Builds the total potential energy:

  ```fortran
  E%potential = ...
  ```

`pot_energy_nonbonds` at `src/q6/md.f90:13514` chooses the specific nonbonded
subroutines based on:

- PBC versus spherical boundary.
- `VDW_GEOMETRIC` versus `VDW_ARITHMETIC`.
- SPC optimized solvent routines versus general 3-atom solvent routines.
- Q-vdW mode.
- LRF mode.

`pot_energy_bonds` at `src/q6/md.f90:13477` chooses bonded terms based on force
field type:

- `FF_GROMOS`
- `FF_AMBER`
- `FF_CHARMM`

### 23.3 Optional Monitored Nonbond Energies

Code path: `src/q6/md.f90:4067-4070`

On output steps, if monitor groups exist:

```fortran
call nonbond_monitor
```

This computes special nonbonded summaries for selected atom groups.

### 23.4 Optional Off-diagonal EVB Terms

Code path: `src/q6/md.f90:4073`

If off-diagonal EVB terms exist:

```fortran
if ( noffd .gt. 0 ) call offdiag
```

`offdiag` starts at `src/q6/md.f90:12911`.

### 23.5 Leap-frog Coordinate and Velocity Update

Code path: `src/q6/md.f90:4076-4121`

On the master node, qdyn uses a Verlet leap-frog style update.

For solute atoms:

```fortran
v = (v - d*winv*dt) * Tscale_solute
x = x + v*dt
```

For solvent atoms:

```fortran
v = (v - d*winv*dt) * Tscale_solvent
x = x + v*dt
```

The old coordinate is saved in `xx` so SHAKE can correct the move and recompute
velocities.

### 23.6 SHAKE Constraint Correction

Code path: `src/q6/md.f90:4124-4127`

If constraints exist:

```fortran
niter=shake(xx, x)
v(:) = (x(:) - xx(:)) / dt
```

The corrected positions define corrected velocities.

### 23.7 Temperature and Thermostat Scaling

Code path: `src/q6/md.f90:4132`

After the coordinate update, master calls:

```fortran
call temperature(Temp,Tscale_solute,Tscale_solvent,Ekinmax)
```

This computes the temperatures for the updated velocities and prepares scale
factors for the next time step.

### 23.8 Constant-pressure Volume Move

Code path: `src/q6/md.f90:4140-4144`

If PBC and constant pressure are active, every `ivolume_cycle` steps after step
0:

```fortran
call MC_volume
```

`MC_volume` starts at `src/q6/md.f90:15615`.

### 23.9 Broadcast Updated Coordinates

Code path: `src/q6/md.f90:4146-4149`

In MPI builds, master broadcasts the updated `x` array to slaves:

```fortran
call MPI_Bcast(x, nat3, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
```

This is why slaves can evaluate nonbonded work for the next step even though the
coordinate integration itself is master-side.

### 23.10 Output During the Loop

Code path: `src/q6/md.f90:4152-4182`

Only the master writes output:

- Trajectory every `itrj_cycle` steps after step 0:

  ```fortran
  call write_trj
  ```

- Energy file every `iene_cycle` steps after step 0:

  ```fortran
  call put_ene(11, EQ, OFFD)
  ```

- Human-readable energy summary every `iout_cycle` steps:

  ```fortran
  call write_out
  ```

- Backup final-coordinate file every 1000 steps:

  ```fortran
  call write_xfin
  ```

- Temperature printout when the temperature changes by more than
  `TEMP_PRINT_THRESHOLD` or when `itemp_cycle` fires.

## 24. Final Energy and Coordinate Write

Code path: `src/q6/md.f90:4196-4208`

After the loop:

- If the final step is on a trajectory interval, qdyn writes a final trajectory
  frame.
- It rebuilds pair lists:

  ```fortran
  call make_pair_lists
  ```

- It recomputes final potential energy:

  ```fortran
  call pot_energy
  ```

- On the master:
  - Prints final `write_out`.
  - Writes final coordinates and velocities with `write_xfin`.

`write_out` starts at `src/q6/md.f90:15292`. It prints:

- Solute, solvent, solute-solvent, LRF, and Q-atom energy rows.
- Restraint energy rows.
- Total, potential, kinetic, and optional RMS force.
- Q-atom energy breakdowns by state.
- Optional monitored group nonbonded terms.
- Constant-pressure volume summary at the final step.

`write_xfin` starts at `src/q6/md.f90:15452`. It writes:

- `nat3` and coordinates.
- `nat3` and velocities.
- Water polarization shell data if present.
- PBC box length and center if PBC is active.

## 25. qdyn Shutdown

Code path: `src/q6/qdyn.f90:132`

```fortran
call shutdown
```

The internal `shutdown` routine starts at `src/q6/qdyn.f90:207`.

On the master node it prints normal termination text. Then it calls:

```fortran
call md_shutdown
```

`md_shutdown` starts at `src/q6/md.f90:388` and calls:

```fortran
call md_deallocate
call topo_deallocate
call qatom_shutdown
```

`md_deallocate` starts at `src/q6/md.f90:561`. It deallocates:

- Core atom arrays: `x`, `xx`, `v`, `d`, `winv`, `iqatom`.
- Nonbond lists: `nbpp`, `nbpw`, `nbww`, `nbqq`, `nbqp`, `nbqw`, `qconn`.
- Water polarization arrays.
- LRF arrays.
- Restraint arrays.
- MPI send/receive arrays in MPI builds.

Finally, if compiled with `USE_MPI`, `qdyn` calls:

```fortran
call MPI_Finalize(qdyn_ierr)
```

Code path: `src/q6/qdyn.f90:135-137`.

## 26. One-line End-to-end Summary

`qdyn` starts in `src/q6/qdyn.f90`, initializes MPI and MD modules, lets the
master parse input/topology/FEP/restart data, prepares constraints and nonbonded
work distribution, enters `md_run`, repeatedly rebuilds pair lists, computes
forces and energies, integrates coordinates with leap-frog plus SHAKE and
thermostat scaling, writes outputs, then deallocates MD/topology/Q-atom state
and finalizes MPI.

