# Charge-only validation: staged physical questions

Status: protocol under implementation, **not launch-ready or physically validated**.
This is a bounded plan for Q's existing molecular dynamics (MD), not a new sampler.
Surface Constraint All-Atom Solvent (SCAAS) forces, thermostat and radial wall remain
unchanged. High-performance computing (HPC) submission requires user approval.

## What would answer the original question?

Bookkeeping correctness and physical adequacy are different tests. The native
tests already examine the first; a small controlled free-energy campaign must
examine the second. Do not restart the large neutral-shell calibration to do this.

| Question | Comparison | Interpretation |
| --- | --- | --- |
| Are saved energies and correction constants consistent? | Same coordinates, both signs, integrated versus post-hoc Born | An exact software/accounting identity, not a physical validation |
| Does the exterior correction reduce finite-radius dependence? | Charge 0 to +1 and 0 to -1, separately, at multiple radii | A necessary consistency check; a larger finite sphere is not an infinite-system truth |
| Is a discrepancy caused by poor sampling? | Independently initialized forward/reverse ladders, time blocks, overlap and timestep checks | Unresolved sampling means inconclusive physics, not permission to tune the correction |
| Does retaining charged environment matter? | Matched zero, positive and negative non-Q backgrounds | Needed before transferring a central water-probe result to proteins |
| Is the monopole approximation sufficient for geometry? | Off-center probe and distributed background controls | Net charge alone cannot establish adequacy for an extended protein |

The non-Q background is the included solute charge outside the state-dependent
Q region. The proposed total-enclosed-charge angular target still requires
separate approval. **At zero background that candidate and the current target
coincide**, so preparation and initial testing need not wait for that decision.
Charged-background candidate comparisons remain conditional on it.

## Minimal system; no ghost endpoint

Use the repository's three-site transferable intermolecular potential water
model, TIP3P, from `AMBER14sb.lib`/`.prm`, unchanged. Prepare a one-atom synthetic
probe using the existing Amber `CT` steric parameters: radius parameter 1.908
angstrom, well depth 0.1094 kilocalories per mole (kcal/mol), mass 12.01 atomic
mass units. This is a declared diagnostic probe, **not a parameterized physical
ion or an experimental hydration target**. No charged result selects its parameters.

Only its charge changes: state 1 is 0, state 2 is +1 or -1 elementary charge.
Lennard–Jones (LJ, repulsion/dispersion), mass and all other terms are identical
in both states. Include both exact endpoints. Historical ghost-atom endpoint
trimming is a different analysis and is not applied here.

The existing `[atom_restraints]` option keeps the probe near the sphere center:
`1 0 0 0 10 10 10 0`. It means atom 1, reference position zero, three Cartesian
force constants 10 kcal/mol/angstrom squared, and state selector 0 (all states).
The potential is `0.5 k |r|^2`, not a fixed-position constraint. Its contribution
is identical in both pure states on any configuration. Keep it unchanged across
radii/signs/directions; record the probe's displacement distribution. A center
restraint does not prove the probe stays exactly at the center.

Use 298 kelvin, existing solvent/hydrogen bond constraints, no local reaction
field approximation, 99-angstrom pair cutoffs and pair-list rebuild every step.
The cutoff choice needs a whole-trajectory coverage check; it is not itself proof
of no truncation. Preserve the current angular/radial force constants 20/60,
dielectric 80 and common zero frozen angular offsets. Zero offsets define this
prospective Hamiltonian; they are not claimed to reproduce historically adapted
offsets or to be optimal. Do not adapt or fit them to charged results.

Qprep's requested grid radius is **not** the effective solvent radius. Retain
both, the water count and the actual native radius used by the dynamics and Born
term. For example, the fresh 10-angstrom preparation yields 146 waters and a
10.15-angstrom effective radius. These are not the old neutral-calibration droplets.

## Stage A: short local checks (implemented, limited)

The [probe tool](../../src/QligFEP/charge_probe.py) generates inputs with Qprep and
offers explicitly bounded Qdyn timing checks. Neither command submits jobs:

```sh
make -C src/q6 qprep qdyn FC=gfortran-11
PYTHONPATH="$PWD/src" python -m QligFEP.charge_probe prepare /tmp/my-probe/r10 \
  --qprep src/q6/qprep --radius 10 --seed 758971
PYTHONPATH="$PWD/src" python -m QligFEP.charge_probe timing /tmp/my-probe/r10 \
  /tmp/my-probe/r10-negative --qdyn src/q6/qdyn --sign=-1 --weight=1 --steps=1000
```

Use new output directories; existing directories are refused. Preparation records
inputs, topology, coordinates, both charge files and Qprep fingerprints using
SHA-256, a cryptographic hash algorithm. The timing check verifies normal exit,
finite saved state energies and lambda identities, the zero-background model and
unchanged zero offsets. It retains the input, log, energy file, restart and binary
fingerprint. JavaScript Object Notation (JSON) reports explicitly say
`production_ready: false`. The tool permits at most 2,000 MD steps per call and
60 seconds per native invocation; timeouts retain diagnostic files, not a pass.

Grid oxygens and seeded hydrogen directions are **not equilibrated water**.
Different orientation seeds alone do not establish independent equilibrium
replicas. A 0.1–1 picosecond (ps) smoke run is not an equilibration or timestep
assessment. Q's existing minimizer freezes solvent and Q atoms, so it must not
be represented as relaxing these generated waters.

The staged/native gates now verify the shared position restraint, and a
[completed-window check](RESTRAINT_AND_COMPLETION.md) verifies the saved state
accounting and final offset record. An [isolated-build chain runner](BUILD_AND_CHAIN.md)
now checks realized dependencies and records build/launch evidence. Before Stage B,
complete campaign generation, endpoint preparation, trajectory coverage, solvent geometry and
runtime/temperature diagnostics too.

## Stage B: capped feasibility pilot (not submitted)

Start with the nominal 10/14-angstrom preparations only. The prospective cap is:

- Two charge signs, two ladder directions, two independently prepared replicas
  per direction: eight ladders per radius, sixteen in total.
- Eleven state-2 weights: 0, 0.1, ..., 1. Forward traverses increasing weight;
  reverse decreases it without exchanging state labels.
- Up to 100 ps preparation per ladder to establish its starting endpoint; the
  exact heating/charging and equilibration schedule must be specified and audited
  before launch. No replica may be called independent merely by copying the same
  final restart. Failure to settle within the cap stops that ladder.
- Each window: 10 ps discarded settling plus 20 ps diagnostic sampling. These
  durations are a **cost cap**, not asserted sufficient equilibration or production.
  The sampling portion remains pilot data and is not used to select model parameters.
- Start with a candidate 1-femtosecond (fs) timestep. Compare against 0.5 fs over
  matched physical durations before qualifying production. The shorter timestep
  doubles step count, not requested physical sampling time.

This is `16 * (100 + 11 * 30) ps = 6.88 nanoseconds (ns)` of aggregate trajectory
at most. Time spent on failed preparation counts against the cap. At 1 fs this
is 6.88 million MD steps. A changed ladder, added replica or longer settling time
requires a revised budget, not an unrecorded extension. Production lengths are
chosen from pilot correlation/overlap/stationarity evidence, not from desired
free energies. More than two replicas per direction are required for the final
physical comparison; the initial two only assess feasibility.

Do not duplicate all trajectories for force-free exterior-Born on/off controls.
For a fixed state its contribution is a coordinate-independent constant, so
audited saved energies can yield both accounting views. Keep a short native
integrated/post-hoc equivalence test. An angular-target control is different:
it changes forces and generally needs its own sampled trajectories.

## Analysis and predeclared decisions

Let `w` be state-2 weight and `E1(x), E2(x)` the saved Q-state energy totals.
They are not the complete solvent-plus-solute potential: write the full pure-state
potential as `U_s(x) = U_common(x) + E_s(x)`, where the unsaved common part cancels
in `U2-U1 = E2-E1` for this restricted charge-only Hamiltonian. The sampled
potential is `U_common + (1-w) E1 + w E2`. For two adjacent weights `a,b`, the dimensionless forward
energy difference on a sample from `a` is
`beta (b-a) [U2(x)-U1(x)]`, with `beta = 1/(k_B T)`; the reverse sample uses its
negative evaluated in ensemble `b`. Use the topology/native unit convention and
record the Boltzmann constant `k_B` used in analysis. The adjacent-state estimator
is the Bennett acceptance ratio (BAR), not a one-sided average of energy gaps.
Its numerical solution succeeding does not establish phase-space overlap.
[The estimator documentation](https://pymbar.readthedocs.io/en/stable/other_estimators.html)
defines a normalized two-state overlap from 0 (none) to 1 (complete).

The analysis implementation must pass synthetic tests for known free-energy
differences, constant shifts, sign/direction reversal, unequal sample counts,
correlated data, poor overlap, truncated/corrupt files and Born double counting.
Pin its dependency versions before launch. Do not silently replace the saved
pure-state totals with one Coulomb component or a lambda-scaled output summary.

For each radius, `B_s = -k_e (1-1/80) (Q_env+q_s)^2/(2 R_eff)`. Report the raw
finite-sphere result and `Delta G_corrected = Delta G_raw + B2-B1`, **or** report
the integrated result and recover raw by subtracting that same difference.
Never add it to an already integrated result. These are fully correlated views
of the same data. In neutral background the Born shift is identical for the two
signs; it changes the charge-even average but not the charge-odd difference.
Positive and negative charging free energies need not be equal.

Proposed operational gates below are declared before a pilot free-energy result;
they are diagnostic choices, not universal guarantees or fitted physics:

1. Reject incomplete runs, nonfinite energies/coordinates, wrong state mappings,
   changed offsets/restraints, topology/asset mismatches, inconsistent water
   parameters, or uncovered pair separations. Preserve all failures.
2. Examine gap, temperature, probe displacement, radial density and shell
   population/orientation traces per window. Estimate discarded transient and
   statistical inefficiency (correlation-induced reduction of sample count).
   Use the most conservative observed timescale and test longer discards and
   block lengths. An automated equilibration detector is a diagnostic, not proof
   of stationarity or absence of slow modes. The
   [time-series documentation](https://pymbar.readthedocs.io/en/stable/timeseries.html)
   describes the relevant discarded-region and correlation estimates.
3. Require normalized adjacent-pair overlap at least 0.03, at least 100 effective
   observations per window and at least 20 blocks, each at least five times the
   estimated statistical inefficiency in saved-frame units, for any
   result presented as qualified. Otherwise mark insufficient sampling and propose
   a bounded longer pilot or finer ladder. Passing these minima is not sufficient.
4. Resample complete trajectory blocks jointly for all adjacent BAR terms using
   each window; do not independently resample the same frames for its two neighbors.
   Preserve covariance in the summed free energy. Analyze each independent ladder
   separately and retain between-replica uncertainty. Neighbor terms share samples;
   summing their individual variances as if independent is not generally justified.
   [The authors' analysis discussion](https://github.com/alchemistry/alchemical-best-practices/blob/main/paper/manuscript.tex)
   makes this dependence explicit. Compare block and replica-based intervals;
   disagreement or too few effective blocks means inconclusive, not a small error bar.
5. Report forward/reverse disagreement in a common 0-to-sign orientation,
   early/late production differences, and timestep differences with uncertainty.
   Work toward a 0.5 kcal/mol practical resolution for each primary contrast.
   A wide interval that includes zero does not demonstrate agreement.

Stage B primarily decides whether a sufficiently precise test is affordable; it
cannot qualify the general correction. For final radius comparisons require at
least three radii, at least four independent replicas per direction, converged
diagnostics and a frozen analysis plan. Use corrected differences from the largest
radius for each sign separately. For two smaller radii and two signs there are
four primary contrasts; use simultaneous 95% coverage (for example, individual
98.75% intervals via the Bonferroni adjustment for four predeclared tests).
Within the current 0.5 kcal/mol resolution target:

- An entire interval within [-0.5, +0.5] supports radius consistency at that
  resolution, **not exactness or a bulk reference**.
- An interval wholly outside that band is evidence of residual radius dependence
  after the sampling gates pass. Do not tune offsets, radius or dielectric to hide it.
- Any other interval is inconclusive at the declared resolution.

Compare raw and corrected radius trends, including uncertainty; a constant shift
alone cannot improve sampling. Even a neutral-background pass leaves charged
background, off-center/distributed geometry, hard shell crossings and the existing
thermostat/radial-wall equilibrium limitations unresolved. If those require new
forces or a new sampler, stop and ask for scope approval instead of redesigning Q.

## Budget and remaining handoff

See [the local timing evidence](LOCAL_TIMING.md). Translate measurements into a
candidate-HPC allocation only after a short same-build/same-input hardware check.
Keep serial aggregate compute and elapsed time with multiple concurrent jobs
separate. Avoid scaling the water count from the old neutral calibration: the new
Qprep systems and existing MD workload are different.

The next implementation work is campaign generation and endpoint preparation,
trajectory/runtime checks and tested BAR/uncertainty
analysis. A full launch manifest and scheduler script must not claim readiness
until those exist. After that, present the capped Stage B allocation for approval.
Neither this plan nor a passing native smoke test authorizes an HPC submission.
