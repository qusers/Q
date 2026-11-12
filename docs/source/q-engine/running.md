# Running Q directly

## Prepare a topology

Use a scripted qprep session for reproducibility:

```bash
src/q6/bin/q6/qprep prep.inp > prep.log
```

Example `prep.inp`:

```text
readlib ligand.lib
readlib protein.lib
readprm combined.prm
readpdb complex.pdb
boundary 1 -4.689 26.119 -30.570 25.0
solvate -4.689 26.119 -30.570 25.0 1 HOH
maketop complex
writetop complex.top
checkbonds 6
checkangs 6
checktors 6
checkimps 6
writepdb complex_prepared.pdb y
quit
```

The exact solvent entry and library/parameter files must match the chosen force
field. Check the complete qprep output and inspect the prepared PDB.

## Initialize-check an MD input

```bash
src/q6/bin/q6/qdum md.inp > qdum.log
```

`qdum` follows the Qdyn initialization path without running the dynamics loop. It
catches many syntax, topology, FEP, restraint, and file errors cheaply. It cannot
demonstrate runtime stability.

## Run serial dynamics

```bash
src/q6/bin/q6/qdyn md.inp > md.log
```

A minimal input needs `[MD]` values for `steps`, `stepsize`, and `temperature`,
plus `[files] topology` and `final`. Without `[files] restart`, provide
`initial_temperature`.

```text
[MD]
steps                50000
stepsize             2.0
temperature          298
initial_temperature  10
random_seed          12345
bath_coupling        100
shake_solvent        on
shake_hydrogens      on
lrf                   on

[cut-offs]
solute_solute        10
solute_solvent       10
solvent_solvent      10
q_atom               99
lrf                  99

[intervals]
non_bond             25
output               50
energy               10
trajectory           100

[files]
topology             complex.top
fep                  FEP1.fep
final                md.re
energy               md.en
trajectory           md.dcd

[lambdas]
1.000 0.000
```

## Run MPI dynamics

```bash
mpirun -np 8 src/q6/bin/q6/qdynp md.inp > md.log
```

On Slurm, follow site guidance; many systems prefer `srun`. QligFEP generates the
launch command from its cluster templates. Use MPI only after the same case
passes `qdum` and a short serial run.

## Continue from a restart

```text
[files]
topology  complex.top
restart   previous.re
final     continued.re
energy    continued.en
trajectory continued.dcd
fep       FEP1.fep
```

Keep the topology, FEP definition, lambda values, restraints, and atom ordering
consistent across a continuation. A restart contains coordinates and velocities;
it does not replace the control input.

## Run minimization before MD

Add to `[MD]`:

```text
minimize                on
max_minimize_steps      1000
minimize_tolerance      0.1
minimize_step_size      0.001
```

The minimizer is intended to remove severe local strain before dynamics, not to
repair incorrect chemistry. Inspect convergence and the resulting coordinates.

## Run non-equilibrium switching

The serial engine enters two-state switching mode when `[lambda_scaling]` is
present:

```text
[lambda_scaling]
scaling_parameter sigmoidal
L_sigmoid 8.0
```

Use equilibrated endpoint snapshots and run both directions. QligFEP automates
endpoint relaxation, snapshot spacing, switch naming, and BAR analysis; direct
use requires reproducing those controls carefully.

## Recognizing success

A scheduler status of `COMPLETED` is not sufficient. Check:

- Q's normal termination message;
- final restart and expected `.en`/`.dcd` files;
- temperature and energy stability;
- absence of SHAKE failures, overflows, or hot-atom cascades;
- correct number of steps and records;
- restart continuity between stages.
