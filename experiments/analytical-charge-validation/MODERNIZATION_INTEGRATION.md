# Constraint-solver integration

This integration combines `feature/modernize-simulation` at `fec62d37`
with the clean analytical-correction branch at `6b1a96db`. Both original
branches remain unchanged. No Monte Carlo sampler work is imported.

LINCS means Linear Constraint Solver. SETTLE is the rigid-water constraint
solver supplied by the modernization branch. That branch also repairs SHAKE:
commit `90122307` rechecks all coupled bond constraints after corrections,
instead of accepting stale per-bond convergence flags.

The first comparison explicitly selects `constraint_algorithm shake shake`
(solute, solvent). This uses repaired SHAKE, retains the existing charge-only
Hamiltonian and timestep, and does not select new solvers implicitly. The
native log's solver selection is checked against the staged input and retained
in the parsed audit. Historical inputs without this key retain their legacy
SHAKE interpretation; absent historical solver evidence is not manufactured.

The broader modernization branch also contains bonded-force singularity handling,
an opt-in minimizer and changes to setup/output management. Therefore this is
not a SHAKE-only source comparison. Minimization and production-file cleanup
remain disabled in the restricted correction protocol; its raw energy and
restart evidence must be retained. No timestep or thermostat change is made
by this integration.

## Validation sequence

1. Compile serial Qdyn/Qprep and run the existing solver unit tests.
2. Check the unchanged 64-case single-water geometry perturbation against
   final squared-distance residuals, without relying on internal ready flags.
3. Verify pure-state energies, Born accounting and native initialization.
4. Compare tracing enabled/disabled using identical current source and solver.
   The test-only compile definition `Q_TEST_DISABLE_CHARGE_TRACE` suppresses
   only the two diagnostic calls. Production builds do not define it; campaign
   validation still requires traces from sources containing the trace routine.
   The historical `4ad905d1` executable is not a valid observer-only control
   because it used the defective solver.
5. Rerun the original 16-cell matrix without changing seeds, timestep,
   durations or geometry alarms. Remove expected-failure markers only when
   their underlying tests have been demonstrated to pass.
6. Check SETTLE separately before deciding to adopt it in a campaign.

Initial integration checks: native LINCS 29/29 and SETTLE 512/512 assertions;
100 template, minimization-option, isolated-constraint, state-energy and
diagnostic tests; 42 probe and native protocol tests passed. Full campaign
validation was pending at that checkpoint.

The subsequent restart-chain, endpoint and completion suite passed 79 tests.
All eleven campaign tests then passed with expected-failure handling disabled:
the original sixteen cells completed preparation, restart transfer, three charge
windows and analysis at the fixed aggregate 6.4 picosecond software-test budget.
No seeds or geometry alarms changed. The four obsolete campaign expected-failure
markers were removed only after this result.

The separate solver comparison uses a common neutral, repaired-SHAKE restart
at requested radius 14 angstrom (388 waters), orientation seed 758982 and
velocity seed 123. It runs 100 one-femtosecond steps for each charge sign and
each of SHAKE/SHAKE, SHAKE/SETTLE and LINCS/LINCS, retaining restart velocities.
All 25 tests passed, including independently reconstructed final water distances,
Born-inclusive saved-energy accounting, and rejection of missing, duplicated or
relabelled solver evidence. This does not assert identical trajectories between
different solvers or validate their equilibrium distributions. The campaign
remains explicitly SHAKE/SHAKE; SETTLE has not been adopted as its default.

These are software checks, not equilibrium or physical validation. No
high-performance computing (HPC) job has been submitted. The analytical
correction's statistical and physical validation requirements remain in force.

## Final regression checkpoint

On 2026-09-08, source commit `e90d1c429fe4417846257b289eebe83c640ea0e4`
passed **413 tests, with 2 skips and 2 expected failures**, in 147.67 seconds.
The remaining expected failures are the archival charge-partition checks,
not constraint or campaign checks. The four additional minimization/template
tests check that boundary overrides survive the merge.

The isolated 64-case SHAKE test's worst squared-distance relative residual is
`9.94479206325028748e-5`, below the unchanged nominal `1e-4` tolerance.
The original failing matrix cell (c11, negative reverse start, orientation seed
758982, velocity seed 123) now has water-distance extrema drift
`1.460484783222249e-6` angstrom, below the unchanged 0.005-angstrom alarm.
This is a replay on the integrated engine, not an attribution of every numerical
difference exclusively to SHAKE.

The successful full campaign and its isolated build are copied unchanged under
`runtime/integration-20260908/{campaign,build}`, excluded from version control.
The copied build passes archive/source/binary fingerprint validation. Original
temporary paths inside reports are retained: these are archival evidence copies,
not relocated runnable campaign manifests.

SHA-256 (a cryptographic hash algorithm) fingerprints:

- Build report: `99658ef6e0ef70de58fa8a3c84acd679170fbb0f43d069edddb8cd524d9fd6f8`.
- c11 seed log: `60f7ff2a151612528f448aacbab0201034b265e4314c0146ec7b83f98b683908`.

To repeat from this worktree with the test dependencies available in the active
Python environment:

```sh
make -C src/q6 qdyn qprep test-lincs test-settle FC=gfortran-11
mkdir -p src/q6/bin/q6
cp src/q6/qprep src/q6/qdyn src/q6/bin/q6/
PYTHONPATH="$PWD/src" PYTHONDONTWRITEBYTECODE=1 python -m pytest -q \
  test/q6 test/qligfep/test_charge_bar.py test/qligfep/test_charge_protocol.py \
  test/qligfep/test_endpoint_trim.py test/qligfep/test_historical_cmet_reproduction.py \
  test/qligfep/test_boundary_correction_cli.py test/qligfep/test_templates.py \
  test/qligfep/test_minimization_option.py
```

The test-only solver comparisons do not authorize a production solver switch.
Next work remains multi-observable discard/block sensitivity, between-replica
and direction/radius comparisons, bounded target-hardware timing and a reviewed
scheduler package. No production sampling result is implied by this checkpoint.
