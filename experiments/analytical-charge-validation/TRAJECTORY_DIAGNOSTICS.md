# Native trajectory observations for the charge pilot

The diagnostic records observe Q's existing molecular dynamics (MD). They do not
change forces, coordinates, velocities, restraints, boundary targets or the
thermostat. They are enabled with the existing per-state spherical polarization
mode on the main process. No new sampler or trajectory propagation method is
introduced. High-performance computing (HPC) submission still requires approval.

## What is measured, and at what time?

`Q_CHARGE_TRACE_V1` is recorded at the energy-summary interval, including step zero
and the final step. It observes coordinates immediately after the force evaluation,
**before** the next position/velocity update. For the standard native MD loop this
is force geometry `istep`, not a post-update coordinate frame. The temperature is
the existing temperature value associated with that incoming state; the observer
does not call the temperature routine or perform another velocity rescaling.

Match observations to saved energy frames by native step number: energy records
exclude step zero and the final step, whereas this trace includes both. Applying
the same array-index discard to these different sequences would misalign them.

At **every** MD force geometry, including the final energy evaluation, the code
updates cumulative observations. A normal `N`-step run must therefore account for
`N+1` geometries, even though most are not printed individually:

- Largest atomic distance from the solvent center, including hydrogen and Q atoms.
- Largest Q-atom distance from that center.
- Minimum/maximum oxygen–hydrogen and hydrogen–hydrogen distances within waters.
- Minimum/maximum total and free temperatures; finite coordinates/velocities and
  finite, positive temperatures. “Free temperature” is Q's existing temperature
  for its free-atom degrees of freedom, not another thermostat.

The trace prints the current maximum Q-atom radius, current temperatures and the
cumulative extrema. For the one-atom pilot, Q-atom radius is probe displacement.
For a larger Q region it is its maximum radius, not a center-of-mass observable.

Periodic `QCT_DENSITY` records count oxygen positions in twenty equally spaced
radial bins between zero and the effective solvent radius. A twenty-first count
records waters at or outside that radius; no arbitrary exterior volume turns
that overflow count into a density. Interior number density divides each count
by its actual spherical-bin volume.

`QCT_SHELL` records each existing polarization shell's water population, sum of
radial dipole cosines and sum of squared cosines. The dipole direction points from
oxygen toward the hydrogen midpoint. Membership reproduces the engine's existing
hard inequalities, including waters beyond the outer surface. These are observation
bins and existing shells, **not smooth restraint windows or a new solvent model**.

## Coverage and failure checks

The reader requires exactly the expected snapshot sequence, evaluation counters,
density counts accounting for every water, ordered shell records and consistent
population/moment bounds. Cumulative maxima must not decrease and cumulative minima
must not increase. Missing, duplicate, malformed, nonfinite or failed records are
not accepted as a completed trajectory.

Q's existing temperature routine can zero a hot atom's velocity and remove its
kinetic energy from the reported temperature. The reader therefore rejects the
native `WARNING: hot atom` message; a subsequently benign temperature must not
hide that intervention. This checkpoint does not modify that existing behavior.

For the direct, spherical pilot, let `r_max` be the largest atomic radius over all
observed force geometries. The geometric inequality `distance(i,j) <= 2*r_max`
provides a conservative bound on pair separations. The gate requires this bound
to be strictly below **all four active pair cutoffs**, not just the Q-atom cutoff.
Unlike a check of sparse saved coordinates, the cumulative bound includes every
MD force geometry. This does not override the separate native checks for exclusions,
water-kernel compatibility or the restricted interaction settings.

The operational water-geometry alarm rejects a cumulative distance extremum that
drifts more than 0.005 angstrom beyond its initialized range. This is a declared
drift alarm, not a fit or a water-parameter calibration. It cannot certify that an
initially wrong geometry or force field is correct; actual prepared assets and
absolute distance ranges remain part of the evidence.

Temperature, density, shell orientation and probe-position summaries report extrema,
means and early/late means. The temperature range includes every force geometry;
the means use the printed snapshots. The summaries **do not establish stationarity,
canonical sampling, independent replicas or convergence**. Existing hard-shell,
radial-wall and thermostat limitations remain unresolved. Correlation-aware
multi-observable discard and block sensitivity analysis is still required.

## Consumption and provenance

The [completed-window checker](RESTRAINT_AND_COMPLETION.md) consumes a trace when
present. The [chain runner](BUILD_AND_CHAIN.md) additionally requires one when its
retained native source includes this diagnostic implementation. Removing an entire
trace cannot silently downgrade a new-build run to an old-build pass. Older build
records can still be inspected, but their `trajectory_diagnostics` is null; they
do not acquire all-evaluation coverage retroactively.

For a completed, verified window:

```sh
PYTHONPATH="$PWD/src" python -m QligFEP.charge_diagnostics \
  /path/to/plan.json --window 0
```

The command is read-only and rejects absent traces. It prints JavaScript Object
Notation (JSON) with `gate: native_trace_consistency_passed` and
`production_ready: false`. The raw native log retains the time series. The chain's
source fingerprints now include the diagnostic reader; freeze the complete
implementation before a run, rather than modifying a running chain's readers.

## Validation and remaining work

The native tests check both charge signs, all-evaluation counts, corrupted traces,
cutoff coverage failures and the geometry-drift alarm. Density, probe position
and shell moments are independently reconstructed from final restart coordinates.
Historical controls compared with the pre-observation engine at commit `4ad905d1` and
require byte-identical saved state energies and final restarts. These tests are
evidence for observational behavior on those controls, not physical correction
validation or universal performance neutrality.

After [constraint-solver integration](MODERNIZATION_INTEGRATION.md), the
observer-only test instead builds identical current source with trace calls
enabled and disabled. The old engine cannot isolate observation effects because
it also has different constraint dynamics. The current comparison retains the
byte-identical energy and restart requirement for both charge signs.

The extra per-step observation work changes runtime cost, so the earlier local
timings are historical estimates rather than timings of this build. Refresh the
bounded timing on the intended hardware before requesting the HPC allocation.
Campaign generation, endpoint-to-ladder transfer and multi-observable/replica
statistical qualification remain required. No long trajectory or HPC job is
authorized by this diagnostic checkpoint.

Verification checkpoint (2026-09-08): **302 passed, 2 optional integration skips,
2 documented archival partition expected failures**, in 93.64 seconds. Fresh
isolated builds used native source commit `53f2d015`, including the new observer.
The tests exercise its actual chain/endpoint/completion/analysis integration,
the diagnostics command, refusal of an entirely removed required trace, and
hot-atom warning rejection. The full focused boundary/accounting/endpoint suite
passed. No HPC calculation or long preparation was run.
