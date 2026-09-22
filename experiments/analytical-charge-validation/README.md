# Charged-correction validation package — under development

Purpose: validate analytical charged-perturbation corrections in Q's existing
molecular dynamics (MD) with finite spherical Surface Constraint All-Atom Solvent
(SCAAS) boundaries. This directory is the production-package entry point, but
**it is not yet qualified for high-performance computing (HPC) production**.
No new sampler, integrator, thermostat, radial wall or solvent model is developed
here. The approved integration reuses the existing modernization branch's
constraint solvers and its SHAKE repair.

Implementation follow-up: [opt-in smooth angular boundary](SMOOTH_BOUNDARY_IMPLEMENTATION.md)
replaces hard membership/rank handling with a differentiable weighted angular
distribution and complete Cartesian derivatives. Legacy mode remains the
default. Scoped continuity/force, state-accounting and restart checks pass;
this is a mathematical repair, not a changed Born
convention or a qualified physical charged-perturbation correction.

The [paired smooth-boundary HPC pilot](SMOOTH_BOUNDARY_HPC_PILOT.md) tests both
targets and both legs against fresh legacy controls, gated by serial/16-rank
agreement and explicit model/restart checks. It uses one midpoint window per
case, not a full binding free-energy calculation.

The [completed pilot results and plots](SMOOTH_BOUNDARY_HPC_RESULTS.md) report
eight successful runs and 48 revalidated stages. Short-run numerical stability
passed on both targets; the physical charge-source/exterior-response question
remains open.

The [source/response follow-up](SOURCE_RESPONSE_CONVENTION.md) separates radial
field and potential, validates the distinction on 16 historical snapshots, and
adds a test-only reference-invariant exterior-response interface. No source
convention is silently changed in the native engine.

The [compact charging-derivative audit](CHARGING_DERIVATIVE_RESULTS.md) passes
81 native first/second derivative comparisons, including complete angular/Born
terms, explicit charge-roundoff accounting and charge-family versus endpoint
mixing checks. No dynamics or physical source change is implied by this result.

The [compact charging-response investigation](COMPACT_RESPONSE_GOAL.md) covers
autonomous protocol preparation, prerequisite local checks and a bounded HPC
pilot on matched finite droplets with both probe/background charge signs. Its
completion requires analyzed evidence, not merely job submission; it does not
authorize another binding campaign or a new sampler.

The [compact-response protocol](COMPACT_RESPONSE_PROTOCOL.md) defines a matched
two-particle/droplet comparison and the existing Born-completion hypothesis.
The [approved pilot](COMPACT_RESPONSE_HPC_PILOT.md) completed all 90 trajectories
and 540 stages. The [results, plots and next decision](COMPACT_RESPONSE_HPC_RESULTS.md)
show passing numerical checks and closely overlapping radius-response curves
for both background signs, but insufficient precision and failed qualification
gates. Retain the smooth repair; do not fit a correction or proceed to binding
production. A narrower charge-spacing/precision diagnostic is recommended.
The real production data were fully recovered and revalidated; 143 regression
tests pass. Total pilot plus prior gate allocation was 310.71 core-hours under
the approved 420 cap. No further jobs were submitted during analysis.

The approved [near-neutral autonomous follow-up](NEAR_NEUTRAL_RESPONSE_GOAL.md)
now has a [frozen, pilot-informed protocol](NEAR_NEUTRAL_PROTOCOL.md): zero
retained background, charges 0/±0.05/±0.10 e, both radii, three independent starts,
200 ps equilibration and 1 ns production. The smooth repair and correction are
unchanged. [Jobs 26533721/26533722 and recovery instructions](NEAR_NEUTRAL_HPC.md)
record the bounded 496 additional core-hour plan and passing 188-test suite.
The [completed results and inspected plots](NEAR_NEUTRAL_RESULTS.md) recover
all thirty trajectories: the 18 Å small-spacing response approaches the neutral
fluctuation estimate, but 12 Å density/centroid changes and replica disagreement
prevent physical qualification. Actual follow-up allocation was 387.79 core-hours.
No thresholds were relaxed. The [completion audit](NEAR_NEUTRAL_COMPLETION_AUDIT.md)
closes this focused investigation with an inconclusive physical outcome; next
diagnose the small-droplet redistribution using saved configurations, without
retuning the correction or launching another binding campaign.

The approved [autonomous redistribution diagnosis](REDISTRIBUTION_GOAL.md)
is [complete with an integrated diagnosis](REDISTRIBUTION_RESULTS.md): sustained
small-droplet repacking is dominated energetically by water–water interactions
and starts before the large centroid shift. All ninety scoped native force
checks pass, but preparation/history, finite-size bias and native dynamics are
not yet distinguished. The model and historical gates remain unchanged.
[Completion evidence](REDISTRIBUTION_COMPLETION_AUDIT.md) includes 77 passing
tests and inspected plots. A [neutral restart-history control](REDISTRIBUTION_NEXT_TEST.md)
was proposed but not submitted. It is now superseded by the approved
[matched solvent-baseline goal](SOLVENT_BASELINE_GOAL.md): audit preparation,
develop finite-droplet structural diagnostics and compare legacy versus smooth
boundary behavior with matched topology and dynamics. No new production
sampling is scheduled; a concrete bounded protocol precedes any allocation
request.

Design follow-up: [boundary design options and failure reproduction](BOUNDARY_DESIGN_OPTIONS.md)
documents the user-preferred continuous angular-boundary approach and the
complementary thermodynamic accounting approach. Unequal protein/water charges
and radii do not invalidate a corrected thermodynamic cycle; they prevent
assuming cancellation without justified corrections. The design note separates
those issues from demonstrated angular discontinuity and charge-label dependence,
and records how to reproduce each diagnostic. No replacement Hamiltonian is
implemented by this documentation update; the archived no-go evidence is unchanged.

Earlier: the six [matched-preparation compatibility jobs](MATCHED_COMPATIBILITY_RESULTS.md)
completed and passed independent verification of all 36 stages, using 36.68
allocated CPU-core-hours. A new fixed-coordinate native audit also exposes
finite angular-energy jumps at shell-membership changes. Together with the
reproduced charge-label dependence, this prevents physical qualification merely
from the passing bookkeeping/runtime checks. The full correction is not ready
for another production campaign; see [the decision](PHYSICAL_DECISION.md) and
[final requirement/evidence audit](INVESTIGATION_COMPLETION_AUDIT.md). The
investigation ends with an explicit rejection of current production use, not
with a qualified replacement model. The final diagnostic suite passes 205 tests.

Previously, all 24 standard-array runs in the [101-window campaign](STANDARD_101_RESULTS.md)
completed. The [follow-up audit](STANDARD_101_FOLLOWUP.md) fixes endpoint-trimmed
reading and passes native energy/force checks on actual protein configurations.
Replica sampling and the physical radius/charge convention remain unqualified;
completed simulations do not yet establish a validated analytical correction.

The [radius/background convention audit](RADIUS_BACKGROUND_CONVENTION.md) now
distinguishes the water radius, charge-group exclusion mask, and assumed
continuum interface. Sixteen read-only native snapshot exports confirm that
the included protein background is not geometrically enclosed by the water
radius. The existing convention is identified, but neither a simple radius
swap nor an angular-target replacement is physically qualified. The next step
is an explicit exterior-electrostatics derivation, not another protein campaign.

That [derivation and preparation-radius audit](EXTERIOR_RESPONSE_DERIVATION.md)
now confirms that Born already uses Qprep's **effective** radius in all eight
prepared setups. It separates exterior response to retained charges from
restoration of neutralized ionic groups, with tested spherical limiting cases.
The archived structures retain charged side chains in the effective outer
3 angstrom; the original neutralization mapping is the next required provenance
check before estimating a target-specific restoration. No correction was changed.

The [neutralization provenance and restoration audit](NEUTRALIZATION_RESTORATION_AUDIT.md)
has now recovered 25 c-Met and 53 Eg5 ionic-to-neutral changes from coordinate-
matched source structures. Their pattern matches a 25-angstrom, zero-offset
neutralization reused for the 20-angstrom campaign. Fixed-coordinate dielectric-80
restoration estimates, including separately excluded ionic groups, are about
-0.12 and -0.76 kcal/mol respectively; these are not revised binding free energies.
The next proposed control is explicit, matched-boundary preparation, not a Born
radius adjustment or an unqualified production campaign.

The user has authorized the [autonomous matched-boundary goal](MATCHED_BOUNDARY_GOAL.md),
including bounded compute. Its first preparation milestone reproduces all four
archived protein topology bodies and generates direction-consistent effective-
outer-shell controls with unchanged heavy-atom coordinates. Identity-based
FEP/restraint remapping and eight fresh-coordinate native accounting checks now
pass. The control backgrounds are both negative, and c-Met Arg54 changes its
fixed exclusion status when a proton is removed; these are recorded preparation
confounds, not concealed by adjusting charges or masks. Six fresh compatibility/
equilibration arrays completed on the cluster (26498482–26498487), using 36.68
of their capped 96 allocated CPU-core-hours. This is not a qualified production
campaign. The separate prepared-coordinate term ledger also confirms large
changes in the integrated Born gaps and much smaller screened ionic estimates;
neither is silently added to a free-energy result. See the goal document for
authoritative reports, paths, recorded observations, and the compute cap.

For visual inspection before further development, see the local
[seven-figure results gallery](runtime/standard-101-visual-review/README.md).
It includes replica estimates, separate protein/water legs, cumulative lambda
curves, adjacent-window overlap, split-time estimates, the Born constant shift,
and selected structural differences. The reproducible
[plotting script](plot_standard_results.py) reads the existing reports without
changing them. Graphics are available as PNG and scalable vector graphics (SVG);
exact endpoint windows remain excluded from free-energy analysis.

The [historical water-constraint blocker](CONSTRAINT_BLOCKER.md) is resolved in
the [modernization integration](MODERNIZATION_INTEGRATION.md): the unchanged
sixteen-cell [software campaign](CAMPAIGN.md) passes with repaired SHAKE.
Fresh campaign inputs explicitly select SHAKE/SHAKE and validate the native solver
report. Other solvers are checked separately, not silently adopted as defaults.
This resolves a software blocker, not the remaining statistical or physical questions.

The first [Snellius Eg5/c-Met smoke job](SNELLIUS_TARGET_SMOKE.md) also passed:
both protein-charge signs, water legs and setup directions retain correct
integrated/post-hoc Born accounting. This is a midpoint compatibility result,
not a converged binding free-energy comparison.

The user has now selected **no softcore**. The separate
[no-softcore validation batch](NO_SOFTCORE_VALIDATION.md) tests interior weights
0.0001, 0.5 and 0.9999 on both targets before longer sampling. Exact endpoints
are excluded; this short batch does not estimate free energies.

The [staged physical protocol](PILOT_PROTOCOL.md) now specifies the minimum
charge-only questions, a capped 10/14-angstrom feasibility pilot and proposed
uncertainty/failure criteria. A fresh Qprep probe generator and bounded existing-MD
timing tool are implemented. The [local timing evidence](LOCAL_TIMING.md) estimates
about 5.4 aggregate serial hours for that pilot before contingency or timestep
checks; it predates the native trajectory observer and is not a current-build HPC
allocation or convergence result. Shared position-restraint
support and a [completed-window check](RESTRAINT_AND_COMPLETION.md) now exist.
An [isolated build and one-window chain runner](BUILD_AND_CHAIN.md) now records
source/build/launch evidence and validates realized restart dependencies.
An [audited chain analysis](ANALYSIS.md) now reports raw/with-Born estimates,
overlap and conditional within-window block uncertainty. A bounded
[fixed-endpoint preparation schedule](ENDPOINT_PREPARATION.md) now retains restart
velocities and includes its grid-start seed in the compute cap. Campaign generation
and endpoint-to-ladder transfer are now tested at software-only durations.
Between-replica/direction/radius analysis, sampling qualification, current-build
hardware timing and the scheduler package still need work before production launch.
No angular target was changed.
Native [trajectory diagnostics](TRAJECTORY_DIAGNOSTICS.md) now cover every MD
force geometry for cutoff bounds and report temperature, water geometry, radial
density and shell orientation. Statistical qualification of these observables
remains incomplete.

## Staged-input preflight

The read-only [checker](../../src/QligFEP/charge_protocol.py) inspects actual
inputs and existing restart files before launch. From the clean worktree, use
the Python environment that imports this worktree's QligFEP:

```sh
PYTHONPATH="$PWD/src" python -m QligFEP.charge_protocol /path/to/manifest.json
```

Success prints a JSON (JavaScript Object Notation) report with
`gate: staged_input_consistency_passed` and **`production_ready: false`**.
Failure exits with status 2 and a diagnostic. It never starts Q, submits a job,
edits an input, copies a restart, or overwrites an output. The report records
actual file content fingerprints using SHA-256, a cryptographic hash algorithm.

The manifest has these exact top-level fields:

| Field | Required contents |
| --- | --- |
| `schema_version` | Integer `1` |
| `engine` | `binary` path, its pinned `sha256`, and full 40-character `source_commit` |
| `series` | Nonempty list of the series objects described below |

Each series declares `id` (unique), `system` (shared physical system identifier),
`sign` (-1 or +1), `direction` (`forward` or `reverse`), positive integer `replica`,
`born_mode`, Boolean `apply_born_posthoc`, and ordered `windows`.
Each window declares `input`, its pinned `sha256`, and `assets_sha256`, a mapping
with pinned `topology`, `fep`, and `restart` file hashes. Manifest paths resolve
relative to the manifest; Q's `[files]` paths resolve relative to each MD input's
directory, which must also be its launch working directory.

`system` means the same topology, boundary parameters, frozen offsets, Q-atom
mapping and state-1 charge vector. Do not use different names to bypass a
physical-system consistency check. Water/protein or different-radius systems
are distinct, but their relationship still needs the campaign-level audit.

### State and correction conventions

The free-energy perturbation (FEP) file accepts only `[FEP]` with `states 2`,
`[atoms]`, and `[change_charges]`. Every Q atom needs explicit charges in both
states. Q-region totals must be 0 and the declared sign. State identities remain
fixed across directions: forward increases state-2 weight from 0 to 1; reverse
decreases it from 1 to 0. Both exact endpoints are required in this **charge-only**
gate. Ghost-atom historical analysis remains a separate, explicitly truncated
Bennett acceptance ratio (BAR) workflow.

No atom-type, bonded, mass or softcore changes are accepted. This verifies the
absence of alchemical Lennard–Jones (LJ, repulsion/dispersion) changes, but native
inspection must still prove that the affected atoms have real, nonzero LJ
parameters in the topology. A zero-LJ dummy present in both states is not ruled
out by the restricted FEP syntax alone.

| `born_mode` | Q's `perstate_born_correction` | `apply_born_posthoc` |
| --- | --- | --- |
| `integrated` | `on` | `false` |
| `posthoc` | `off` | `true` |
| `control` | `off` | `false` |

These are explicit exterior-Born accounting modes, not three different angular
models. The checker does not calculate or apply a correction. The eventual
analysis must consume the audited mode, derive constants from native charge,
radius and Coulomb-unit definitions, and check saved energy records. Merely
declaring `posthoc` does not prove that a later analysis applies it correctly.

### Supported MD input contract

All simulation controls must be explicit. Unsupported sections or keys, repeated
sections/keys, coefficient overrides and ambiguous state mappings are rejected.
The required sections are `[MD]`, `[cut-offs]`, `[sphere]`, `[solvent]`,
`[intervals]`, `[files]`, and `[lambdas]`; the exact key allowlist is in the checker.
An optional `[atom_restraints]` section accepts only shared state-0 Cartesian
restraints on unique mapped Q atoms, with finite reference coordinates and strictly
positive force constants. Their complete ordered definitions enter the system
identity and are checked against the native record. Other extra restraints remain
unsupported; this does not make arbitrary restrained Q inputs acceptable.
This is intentionally not a general validator for arbitrary Q jobs.

- Use existing spherical MD, direct electrostatic interactions (`lrf off`, which
  disables the local reaction field approximation), solvent/hydrogen bond
  constraints on and solute constraints off. Declare all five cutoffs; native
  geometry/coverage checks are still required to prove no interactions are cut.
- Explicitly enable polarization, charge correction and per-state polarization;
  disable offset adaptation. Use dielectric 80 to match the current hardcoded
  angular dielectric factor, not a separately chosen exterior dielectric.
- Declare solvent radius and positive wall/force parameters, temperature,
  timestep, coupling, initial temperature, random seed, steps and output intervals.
  Seed zero explicitly retains restart velocities; a positive seed regenerates
  Maxwell velocities. Record one consistent strategy within each chain.
  Inputs in the same physical system must match except seed, paths, lambda and
  the declared exterior-Born on/off control.
- Every starting restart must already exist, contain finite coordinates and
  velocities with consistent dimensions, and have a valid finite offset record.
  Supported binary dialect: little-endian, four-byte record markers/integers,
  double-precision coordinates/velocities and single-precision offsets.
- Require identical offset-record hashes within and across series of one system.
  Reject reused output paths, existing outputs, or paths that overwrite inputs.

An as-yet-unwritten chained restart cannot pass this staged gate. The
[chain runner](BUILD_AND_CHAIN.md) instead validates planned dependencies and
rechecks realized files before each launch. Its planned-chain result is distinct
from this staged gate; no restart hash check is removed to make future files pass.

## Evidence and remaining work

### Native initialization audit

Q now emits a read-only `Q_BOUNDARY_AUDIT_V3` block during initialization when
per-state polarization is enabled (on the main process only). It records loaded
values before the internal Coulomb charge rescaling; it does not change any
force law, coordinate, velocity, or target. The
[native checker](../../src/QligFEP/boundary_native.py) compares that block against
one window from a **saved preflight report**, after the corresponding run:

```sh
PYTHONPATH="$PWD/src" python -m QligFEP.boundary_native preflight.json \
  --series positive-forward --window 0 --log /path/to/run.log
```

Use the actual series ID; window indices start at zero in the manifest's order.
Save the staged preflight report before launch: re-running the staged gate after
outputs exist intentionally fails its no-overwrite check. The native command
checks the retained manifest and engine hashes, rechecks the current input and
asset hashes, reconstructs the staged window, and rejects disagreement. It
prints JSON without rewriting any files. Success means
`native_initialization_consistency_passed`, still **not production readiness**.

The native block identifies the current angular convention explicitly as code
1, Q-region-only. No proposed total-charge option has been implemented. The
record payloads are:

| Record | Values in order |
| --- | --- |
| `META` | Atom, solute, water, Q-atom, state and shell counts; LJ combining-rule and solvent-type codes |
| `FLAGS` | Per-state polarization, integrated Born, adaptation, alchemical LJ changes, local reaction field, periodic boundary flags |
| `PARAMETERS` | Effective and requested water radii; topology Coulomb constant; dielectric; coefficient override; included/excluded non-Q solute charge; angular/radial force constants; Morse depth/width; maximum loaded atom distance from solvent center |
| `CENTER` | Solvent-center coordinates |
| `SOLUTE_BOUNDARY` | Effective inner restrained radius, exclusion radius, solute shell force constant |
| `RESTRAINT_COUNTS` | Sequence, position, distance, angle and wall restraint counts; external restraint-file flag |
| `POSITION` | Record index, topology atom, reference x/y/z, Cartesian force constants, state selector |
| `CUTOFFS` | Solute–solute, solute–water, water–water, Q-atom, local-reaction-field cutoffs |
| `WATER` | Number density and molecular dipole magnitude used in the target |
| `WATER_COMPATIBILITY` | Uniform water site types/charges; zero hydrogen LJ coefficients (Boolean flags) |
| `WATER_ATOM` | Site index, type, mass, charge, three A and three B topology LJ coefficients for the first water |
| `STATE` | State index, lambda, total/included/excluded Q-region charge, Born energy actually added |
| `QATOM` | Q index, topology index, exclusion flag, type; mass; three A and three B topology LJ coefficients; pure-state charges |
| `SHELL` | Shell index, outer radius, width, offset, pure-state angular strengths |

`CONVENTION` precedes these records; `BEGIN`/`END` delimit exactly one block.
Version 3 adds loaded position restraints and extra-restraint counts to version 2's
water-site records and compatibility checks. Versions 1 and 2 are rejected by the
current native gate; archived timing logs remain valid records of their older
checks, not passes of this stronger gate. No archived log is rewritten or upgraded.
Units are Q's existing angstrom, elementary-charge and kilocalorie-per-mole
conventions; the A/B coefficients retain the topology combining-rule convention,
not a universal sigma/epsilon interpretation. The checker compares atomic
charges and stored shell radii at their actual single-precision representation.
The global effective radius and LJ coefficients are double precision. It does
not enlarge a tolerance to conceal these representation differences.

The inactive local-reaction-field cutoff is reported as zero by convention;
Q does not read that cutoff key when the method is off, and its inactive variable
need not be initialized. Zero here is not an active zero-range interaction rule.
The initial cutoff bound uses twice the maximum loaded atom radius. A passing
geometric bound is not a proof of pair-list correctness or whole-trajectory
interaction coverage; the initial coordinates also precede initial bond-constraint
projection. The original/native logs and later trajectory coverage remain needed.

The checker rejects nonuniform water site types/charges, inconsistent neutral
water charges, and optimized SPC-like water flags with nonzero hydrogen LJ
coefficients. SPC means simple point charge; the optimized kernel label describes
its interaction assumptions, not proof of a particular named water model.
It also rejects excluded/non-solute Q atoms, nonpositive standard Q-atom LJ
coefficients, wrong loaded state charges/lambdas, mismatched frozen offsets and
inconsistent applied Born constants. It returns the independently recomputed
Born state constants even in post-hoc/control mode, distinguishing those values
from the zero Born energy actually added by the engine. It does not apply them
to a free-energy result. Source-build provenance, log-to-executable attribution,
full solvent-model identification, job completion, final output validation and
physical qualification are not proved by an initialization block.

The automated native smoke test constructs both signs and both ladder directions
on the existing Na/benzene/water fixture: 12 windows at state-2 weights 0, 0.5 and
1, four MD steps per window. It runs the preflight, launches the unchanged Q
workflow, checks finite saved state energies and exact saved lambda mappings,
and verifies frozen restart offsets. Shared starting coordinates are intentional
for this software test: it is **not** a set of independent equilibrated replicas,
a timestep assessment, or a physical endpoint/free-energy validation.

The current smoke test uses a generated copy selecting the general three-site
kernel, preserving every topology charge and LJ coefficient. The archive has an
incompatible optimized-water flag despite a small nonzero hydrogen LJ attraction;
the stronger native gate now rejects that original combination. The
[partition diagnosis](../../docs/analytical-charge-corrections/PARTITION_AUDIT.md)
documents the exact omitted term and the minimal no-softcore force-denominator
initialization fix needed by the spherical three-site charge-only kernel. This
copy is a software control, not an endorsed production water parameterization.
Earlier native passes did not check this water-kernel assumption and must not
be interpreted as passing the new gate. No archived topology or historical
result was changed, and no angular target was changed.

A separate fresh-probe integration test now exercises 12 actual windows with the
common position restraint: both signs, both directions and weights 0/0.5/1 at
298 K. Each is 100 existing-MD steps from a shared 20-step seed restart. The
staged/native and completed-window checks pass, including every saved pure-state
total and applied Born constant. These remain software tests with shared starts,
not independent equilibrated replicas. The chain wrapper revalidates completed
predecessors and their final-file fingerprints before reusing a restart.

Still required before production:

1. Resolve and document the angular background-charge convention exposed by the
   [partition audit](../../docs/analytical-charge-corrections/PARTITION_AUDIT.md).
   The [bounded opt-in proposal](CHARGE_CONVENTION_DECISION.md) is now documented
   and awaits approval; no target change has been implemented. Do not replace
   the physical target as an incidental preflight change.
2. Collect the now-implemented native initialization evidence for the actual
   campaign and tie it to an exact source build. Complete whole-trajectory
   interaction coverage, solvent-model identification and launch provenance;
   the local fixture's native pass cannot substitute for those checks.
3. A complete both-sign/control/radius/environment campaign, independent replicas
   and initializations, runtime/final-restart audits, and immutable provenance.
4. An analysis package with matched state/leg definitions, Born accounting,
   block/replica uncertainty, equilibration/stationarity and overlap checks, and
   predeclared physical acceptance/failure criteria. No fitting charged results.
5. A measured throughput benchmark and justified compute budget, followed by user
   approval before substantial HPC computation or submission. A short local
   baseline now exists; intended-hardware performance and the final allocation
   remain unverified.

The [approved goal](../../docs/analytical-charge-corrections/PROPOSED_GOAL.md) is
unchanged. The input gate and smoke test are components of that goal, not a
replacement completion criterion.

Verification checkpoint: **107 passed, 2 skipped, 2 expected failures in 13.47
seconds** across the focused native and analysis suite. This includes 34 new
input-contract/command-line tests and the native 12-window smoke test. The two
optional analysis-runtime skips and the two known non-angular partition failures
are unchanged. `git diff --check` passed. No HPC job was submitted.

Native-audit checkpoint: **123 passed, 2 skipped, 2 expected failures in 12.44
seconds**, after rebuilding serial Qdyn with GNU Fortran 11. This includes the
12-window native matrix, two additional nonintegrated-Born runs, 12 deliberately
altered native-record checks, and successful/failing native command-line cases.
All new native checks executed. The target-change proposal is still awaiting
approval. Distributed-execution and other compiler builds remain unverified.

Water-kernel checkpoint: **137 passed, 2 skipped, 2 expected failures in 29.08
seconds**, after rebuilding serial Qdyn. The archive's strict partition failures
remain explicit; their LJ and Coulomb causes are now reproduced quantitatively.
The general-water test copy passes the rounding-accounted partition, Cartesian
force/energy, state-bookkeeping and native MD checks. Version-2 compatibility
checks reject the original incompatible water flag and deliberately altered water
records. The angular target is unchanged and its proposed opt-in replacement
still awaits approval. No HPC jobs were submitted.
