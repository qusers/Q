# Shared probe restraint and completed-window accounting

Purpose: make the charge-only probe usable in the existing Q molecular dynamics
(MD) validation workflow, and check outputs before any future restart-chain
wrapper reuses them. This does not change the Surface Constraint All-Atom Solvent
(SCAAS) angular target, radial wall, forces, integrator or thermostat. No
high-performance computing (HPC) calculation was submitted for this checkpoint.

## Existing position restraint, bounded input support

The staged checker accepts the optional `[atom_restraints]` section only as
`atom x0 y0 z0 kx ky kz 0`. It requires unique mapped Q atoms, finite coordinates,
positive finite Cartesian force constants and state selector zero. Each complete
ordered definition enters the physical-system signature. Changing the reference,
constant or restrained atom in just one window/sign/direction fails consistency.
Other restraint sections and state-specific position restraints are not accepted.

This uses Q's existing potential

    V_position = (kx dx^2 + ky dy^2 + kz dz^2)/2,
    gradient = (kx dx, ky dy, kz dz),

with distances in angstrom and force constants in kilocalories per mole per
angstrom squared (kcal/mol/angstrom squared). Zero is the state selector for all
states, not a zero force constant or frozen atom. On the same configuration this
adds the same energy to each saved pure-state restraint/total and exactly one
copy to the lambda-weighted full potential. Its difference between charge states
is zero, although it affects the shared sampled ensemble and the probe's position.

The native fixed-coordinate regression tests this with a nonzero reference and
unequal Cartesian constants, both charge signs, seven lambda mixtures and four
angular/Born controls. It checks the whole-system energy, pure-state total and
restraint, serialized records, and the full Cartesian gradient difference. This
is an independent analytical check of the existing restraint's contribution, not
a fitted free-energy or dynamics result.

Q now emits `Q_BOUNDARY_AUDIT_V3`. Its new `RESTRAINT_COUNTS` record contains the
sequence, position, distance, angle and wall counts plus external restraint-file
flag. Zero or more `POSITION` records contain index, topology atom, reference
coordinates, force constants and state. The Python native checker requires exact
integer indices and the same shared restraint definitions as the staged input;
it rejects any unsupported extra restraint count. The engine change is read-only
reporting. Version-1/2 logs lack these checks and are deliberately rejected by the
current checker rather than silently treated as version 3.

## Completed-window check

After a run, use its saved preflight report (do not rerun the no-overwrite staged
check after outputs exist):

```sh
PYTHONPATH="$PWD/src" python -m QligFEP.charge_completion preflight.json \
  --series=-1-forward --window=0 --log /path/to/native.log
```

Series identifiers come from the retained manifest; window indices start at zero.
The command is read-only. Success emits JavaScript Object Notation (JSON) with
`gate: completed_window_consistency_passed` and `production_ready: false`.
Failure exits with status 2; it does not repair or discard a bad frame or restart.

The checker first verifies the native initialization against the actual staged
input/topology/charge/restart files. The selected series, window and engine
declarations must also match the retained manifest; editing a saved preflight
report cannot silently substitute another otherwise-valid window. It then checks:

- One normal native termination and final energy-summary marker.
- Exactly the expected number of saved frames for Q's current loop/output rule:
  `floor((steps-1)/energy_interval)`, with at least one frame.
- Every saved state index, lambda and energy value. The reader streams records
  rather than loading the whole energy file into memory.
- Every electrostatic and Lennard–Jones (LJ, repulsion/dispersion) subtotal and
  pure-state total, including the Born energy actually added by the native run.
- Finite final restart coordinates/velocities, matching atom dimensions and an
  unchanged frozen-offset record.

The supported binary dialect is little-endian, four-byte record markers and
state integers, fifteen double-precision numbers per state, and an empty
off-diagonal record after each two-state frame. No coupled-state records, other
endianness or silently truncated last frame are accepted. These restrictions
match the staged uncoupled two-state protocol; this is not a universal Q reader.

### What exactly is summed?

Q's serialized state fields contain redundant subtotals. If `qx` denotes all
nonbonded Q-atom interactions and `qq`, `qp`, `qw` denote Q–Q, Q–non-Q-solute and
Q–water interactions, respectively, then separately for electrostatics and LJ:

    qx = qq + qp + qw.
    E_state = bonded + qx_electrostatic + qx_LJ + restraint + B_applied.

Do not add `qx` and its three constituent pairs again. For this branch the
integrated exterior-Born constant is in the state **total**, not its restraint
field; omission or a second addition fails this identity. In declared post-hoc
or uncorrected control mode the native applied constant is zero, even though
the native checker also returns the independently derived prospective correction.
The completed-window check never applies a post-hoc correction to a free energy.

These saved totals contain the state-dependent part used for charge-perturbation
analysis, not all state-independent water and solute energies. Their difference
is the full state-energy difference for the restricted charge-only Hamiltonian.
The accounting tolerance is 1e-8 kcal/mol absolute with 1e-12 relative tolerance; residuals
are reported. This is floating-point consistency, not a physical error allowance.

The result records input/log/energy/final-restart fingerprints using SHA-256, a
cryptographic hash algorithm, along with the frame count and uncensored energy-gap
range. It does not clip large values, trim exact endpoints or infer that missing
endpoint caps cancel. A ghost-atom analysis is a separate, truncated protocol.

## What this does not establish

A finite last restart and normal exit do not prove that a trajectory equilibrated,
stayed stable, kept all interactions covered, or sampled the intended equilibrium
ensemble. Intermediate coordinates, temperature and solvent geometry still need
checks. Energy records do not contain a step identifier, so a frame-count check
alone cannot authenticate their temporal ordering or origin. The report does not
prove which source/build produced a supplied log;
an actual launch wrapper and immutable provenance record remain necessary.

Fresh-input integration tests run both signs and both directions at weights
0/0.5/1, with the shared position restraint, 298 kelvin and 100 MD steps per
window. They use one shared 20-step seed restart and are explicitly **not**
independent equilibrated replicas. Additional runs check nonintegrated Born modes;
adversarial tests alter saved totals, Born contributions, counts, offsets, record
contents, native restraint fields and retained-report declarations.

This completed-window check is now used by the separate
[isolated-build chain runner](BUILD_AND_CHAIN.md). It is not itself a scheduler,
statistical analysis or production readiness certificate. The remaining work is
specified in the [pilot protocol](PILOT_PROTOCOL.md).

Verification checkpoint: the focused suite passes **204 tests**, with two optional
integration skips and the same two documented archival partition expected failures,
in 26.51 seconds. The only native-engine change at this checkpoint is initialization
reporting; no production force, target or energy expression was changed.
