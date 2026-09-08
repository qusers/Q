# Audited charge-only free-energy analysis

This is analysis of Q's existing molecular dynamics (MD), not a new sampler or
force model. The Bennett acceptance ratio (BAR) estimator consumes the saved
state-energy differences and estimates adjacent-window free energies. It is now
connected to the [audited restart-chain workflow](BUILD_AND_CHAIN.md). Neither a
finite estimate nor a confidence interval qualifies the physical correction.

## Command and data contract

```sh
PYTHONPATH="$PWD/src" python -m QligFEP.charge_analysis /path/to/plan.json \
  --discard-frames 1000 --bootstrap 1000 --seed 112
```

This command is read-only and submits no high-performance computing (HPC) job.
The example discard is 10 picoseconds only when energies are saved every 10
one-femtosecond steps; it is **not** an automatic equilibration prescription.
Use the actual run intervals and the preregistered discard choice. The argument
is mandatory; discarding all frames in any window fails. Exact charge-only
endpoints must remain in the complete ladder. No energy is clipped, no interior
window is silently omitted and no ghost-endpoint trimming is performed.

Analysis first revalidates the full chain's build/input/launch receipts, native
initialization, saved energy totals and final restarts. An incomplete chain or
changed output fails. It records fingerprints, discard counts, runtime versions
and its source-file fingerprints in JavaScript Object Notation (JSON). Fingerprints
use SHA-256, a cryptographic hash algorithm. The output always has
`production_ready: false`; `estimate_status` distinguishes insufficient sampling
from conditional gap-statistics evidence.

With a native observation trace, energy gaps are matched to temperature, probe
radius, oxygen density and shell-orientation records by **native step number**.
Step zero and the final trace snapshot are not energy frames and are excluded
from this matching. The same declared energy-frame discard then selects the
aligned observables. Missing coverage fails: the analysis does not interpolate
observables or silently drop energy samples. The report retains the matching
steps, intervals, observable names and native-log fingerprint.

The numerical dependencies are already part of this project's Python stack.
The [analysis requirements](analysis-requirements.txt) pin the versions observed
in this checkpoint's local environment: NumPy 2.5.2 and SciPy 1.18.0. Verify that
the intended HPC analysis environment can reproduce these versions and run the
tests before allocation. The command records actual versions; it does not install
packages or silently change the environment. PyMBAR is not a runtime dependency.

## State energies, temperature convention and BAR

Write the full state potential as `U_s(x) = U_common(x) + E_s(x)`, where `E_s` is
the serialized Q-state total. The unsaved common energy cancels from `U2-U1`.
For a state-2 interpolation weight `w`, the potential is
`U_common + (1-w) E1 + w E2`. Adjacent states `a,b` therefore use:

    W_forward = beta (b-a) [E2(x)-E1(x)]   on samples from a,
    W_reverse = -beta (b-a) [E2(x)-E1(x)]  on samples from b.

These are dimensionless energy differences, not recorded nonequilibrium switching
work. The program solves the BAR likelihood equation in a numerically stable
logarithmic form, retaining the actual forward/reverse sample-count ratio. This
equation and the distinction between estimates and overlap diagnostics are
described in the [estimator documentation](https://pymbar.readthedocs.io/en/stable/other_estimators.html)
and [primary estimation paper](https://doi.org/10.1063/1.2978177).

`beta = 1/(k_B T)`. The analysis reads `k_B`, the Boltzmann constant in Q's
kilocalorie-per-mole and kelvin units, from the retained compiled source snapshot.
The supported source declares it as default single precision, currently
`0.001986`; the code preserves that representation when promoting to double
precision. It does not silently substitute `0.001987` or a different precision.
An unrecognized source declaration fails rather than guessing the convention.
Temperature comes from the actual matched inputs. This matches the engine's
declared convention; it does not prove that its thermostat samples an exact
canonical equilibrium distribution.

The normalized two-state overlap is zero for no overlap and one for identical
ensembles. In the implementation, with `N_A,N_B` samples, pooled differences
`x = (W_forward, -W_reverse)` and BAR solution `f`, define

    p = 1/[1 + exp(-(f - log(N_A/N_B) - x))],
    overlap = (1/N_A + 1/N_B) sum p(1-p).

This is the two-state form of the multistate Bennett acceptance ratio (MBAR)
overlap scalar used by [PyMBAR's overlap routine](https://raw.githubusercontent.com/choderalab/pymbar/master/pymbar/other_estimators.py).
It is an overlap diagnostic, not a probability that the correction is right.
An exact discrete-distribution regression independently checks overlap 0.9 and
its known free energy with unequal sample counts. If overlap falls below the
implementation's 1e-12 resolution floor, the aggregate point estimate is withheld;
a numerical root alone must not turn disjoint samples into an apparent zero result.

All reported aggregate results use the canonical transformation **charge 0 to
the declared sign**, even for a decreasing-weight reverse chain. Pair-level
dimensionless estimates retain traversal direction; their normalized overlaps
do not depend on swapping forward and reverse roles.

## Born exactly once

The [completed-window checker](RESTRAINT_AND_COMPLETION.md) already verifies
the Born contribution in every saved pure-state total. Analysis independently
uses the native charge/radius/unit-derived state constants `B1,B2`:

- Integrated mode: subtract `B2-B1` from each saved state-energy gap to obtain
  the raw finite-sphere view; add that constant once to the final raw free energy.
- Post-hoc mode: saved gaps are raw; add `B2-B1` once to the final result.
- Control mode: the declared result remains raw. A separately labeled with-Born
  comparison view is available but is not applied to the control result.

The raw and with-Born views share exactly the same samples and uncertainty; their
intervals differ only by the known constant. They are not independent experiments.
Both signs in a neutral background receive the same exterior-Born shift; equal
positive/negative raw charging free energies are not required. The native test
runs all three accounting modes and confirms identical raw estimates and matching
integrated/post-hoc corrected estimates, without retuning or new force terms.

## Conditional uncertainty: what is and is not estimated

For each window the code estimates an autocorrelation factor `g`, the loss of
independent information due to temporal correlation. It uses the finite-series
autocorrelation, truncates nonpositive successive lag-pair sums and enforces
nonincreasing positive sums. `g` is clamped to at least one. This particular
diagnostic is not claimed to be PyMBAR's automated equilibration detector.
Correlation analysis and choosing a discarded transient are distinct operations;
the [time-series documentation](https://pymbar.readthedocs.io/en/stable/timeseries.html)
provides background on that distinction.

For native-trace builds, the block-selection factor is the **largest** measured
factor across the energy gap, total/free temperature, maximum Q-atom radius,
twenty radial oxygen-density bins, exterior water count, and each existing shell's
population and first/second radial-orientation sums. Here the Q region denotes
the atoms assigned state-dependent charge treatment; in this pilot it is the
single probe. The report includes every factor, the gap-only factor, the controlling
observable and constant traces. A constant empty density bin is not itself a
mixing failure or evidence of good mixing. A constant energy gap retains the
existing conservative failure rule.

Cumulative extrema from the geometry audit are not treated as stationary time
series for correlation estimation. Their separate safety checks remain in force.
Older builds without observations remain explicitly `gap_only` and do not gain
multi-observable qualification retrospectively.

Without `--block-length`, each window uses blocks at least `ceil(5*g)` frames
long. An explicit length is accepted but flagged if shorter than that diagnostic
minimum, using the largest observed factor when multiple observables are available.
Every window needs at least 100 estimated effective observations and
20 full blocks, and each adjacent overlap must be at least 0.03. Constant traces
cannot diagnose mixing. Failing any of these operational gates yields no
confidence interval, not an artificially small error bar.

When those gates pass, a circular moving-block bootstrap resamples consecutive
chunks of each window's existing frames. One resampled series per window is
reused for **both** neighboring BAR edges. The summed free energy consequently
retains shared-frame edge covariance; the report includes the covariance matrix
and conditional standard error. It does not sum independent edge variances.
If a bootstrap draw has inadequate overlap, the interval is withheld rather than
discarding only unfavorable draws. Draw counts and the random seed are reported.

The 95% interval is the 2.5/97.5 percentile interval of the jointly resampled
ladder totals. It is conditional on stationary sampling within the retained
regions and on the block-length choice. Circular blocks join the end and start
of a trace; they are inappropriate evidence of stationarity in a drifting trace.
An autoregressive test with known correlation timescale checks the estimator and
short-block rejection, not universal coverage of this interval in molecular systems.

Residual correlation **between chained windows** and **between replicas** is not
included. Independent-ladder replication, longer-discard/block sensitivity,
early/late and forward/reverse comparisons, other slow observables, and radius
contrasts still need analysis before a physical claim. The current
`statistical_gates_passed` field is deliberately not a production or general
equilibration gate. The older `gap_statistical_gates_passed` key remains an alias
for compatibility; `correlation_basis` identifies whether aligned observables
were included. A passing conditional analysis is labeled
`conditional_multi_observable_statistics` or `conditional_gap_statistics_only`,
not equilibrium sampling. Existing MD/thermostat and hard shell-boundary limitations
remain separate from energy accounting.

In particular, the [known SHAKE constraint defect](CONSTRAINT_BLOCKER.md) still
blocks production qualification, even if a short individual chain passes its
geometry alarm or a statistical minimum. This analysis update does not repair
that solver, change the alarm, authorize a new sampler, or validate the correction.

## Verification and next work

Tests cover known Gaussian and exact discrete free energies, unequal counts,
constant shifts, direction reversal, extreme/disjoint work values, invalid inputs,
correlation scale, joint edge covariance and absent intervals for inadequate
sampling. Native integration analyzes both charge signs/directions and all three
Born accounting modes. The short native chains are labeled insufficient sampling
and do not produce physical confidence intervals.

Verification checkpoint (2026-09-08): **252 passed, 2 optional integration skips,
2 documented archival partition expected failures**, in 51.46 seconds. This
includes a fresh isolated native build, the actual analysis command-line output,
rejection of incomplete chains and a missing discard choice, and all focused
boundary/accounting/endpoint regressions. The skips are the separately configured
native endpoint and historical reproduction integrations, not successes. Only
analysis, tests and documentation changed in this checkpoint; no native force,
energy expression or dynamics method changed.

Still required: between-replica/direction/radius comparisons, additional observable
diagnostics and sensitivity analyses, transfer from the
[endpoint preparation workflow](ENDPOINT_PREPARATION.md) and the complete pilot
campaign generator. The [physical protocol](PILOT_PROTOCOL.md) and its no-fitting
rules remain in force. No HPC experiment was submitted for this analysis work.

Multi-observable checkpoint (2026-09-08): **56 focused estimator and native-chain
tests passed in 28.07 seconds**. The new tests demonstrate that a slow nonenergy
observable can withhold an interval without changing the BAR point estimate,
reject blocks that are too short for that observable, and reject missing or
misaligned observations. Native chains cover both signs/directions and Born
accounting modes with actual step matching. This is a focused analysis check,
not a rerun or replacement of the previously recorded full suite. No constraint,
dynamics or native diagnostic code changed; the pending SHAKE repair remains
unapproved and the campaign's existing failures remain unresolved.
