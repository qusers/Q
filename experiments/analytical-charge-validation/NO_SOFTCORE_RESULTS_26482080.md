# No-softcore Eg5/c-Met results — job 26482080

Analysis date: 2026-09-09. Protocol: [no-softcore validation](NO_SOFTCORE_VALIDATION.md).

## Decision

**Pass the short-run numerical/bookkeeping gate on both targets. Do not yet
claim that the analytical correction gives accurate binding free energies.**
There is no observed charge-sign-dependent Born accounting failure in this
batch. Large near-endpoint energy differences remain, and their dominant
terms are Lennard–Jones (LJ, repulsion/dispersion) interactions, not Born.
This supports continuing the correction work while treating ghost-state
sampling and energy evaluation as a separate issue. It does not justify
fitting a new boundary model or changing boundary parameters to match targets.

## Completion and independent checks

Slurm reports `COMPLETED`, exit `0:0`, elapsed **46 minutes 33 seconds**, and
aggregate process CPU (central processing unit) time **2 hours 56 minutes
33 seconds**. Sixteen cores were allocated/billed; four simulations ran at
most concurrently. The build's LINCS (Linear Constraint Solver) tests passed
29/29 assertions; SETTLE (rigid-water constraint solver) passed 512/512.

All 24 paired checks passed: two targets, two legs, two setup directions,
and state-2 weights 0.0001, 0.5 and 0.9999. Each pair contains two 2-picosecond
molecular dynamics (MD) runs with identical initial conditions, differing only
in whether the Born energy is included during simulation. All 48 target runs
finished. The target budget was 96 picoseconds, plus 0.08 picoseconds of
synthetic probes. No exact endpoint was sampled.

The analysis downloaded the target inputs, reports and raw outputs, verified
**528 recorded output hashes**, and independently checked:

- All **95,952 saved energy frames** (1,999 per run) have the expected weights
  and finite values. These are correlated saved frames, not independent samples.
- Saved pure-state totals sum correctly, including the applied Born constant.
- Every integrated/post-hoc pair has byte-identical final restarts. Comparing
  all saved fields leaves a maximum Born-subtracted difference of
  **1.7621459846850485e-12 kilocalories per mole**, below the 1e-8 tolerance.
- Recomputing geometry/temperature diagnostics from all 48 raw logs reproduces
  the recorded reports. No hot-atom velocity resets occurred.
- The isolated build and all **883 pinned source files** still validate remotely.

The paired simulations are deterministic controls, not independent replicas.
Forward/reverse starts have different restart bytes in all four target/leg
combinations, but that alone does not establish independent equilibrium sampling.

## Short-run geometry and temperature

The mean free-degree-of-freedom temperature per run ranges from **297.865 to
298.415 kelvin**. All-force-evaluation extrema across runs are 278.044–309.435
kelvin. Water-distance extrema drift by at most **0.0000757 angstrom**, versus
the predeclared 0.005-angstrom alarm. These checks show no detected numerical
instability over these short trajectories; they do not establish equilibrium
or thermostat ensemble accuracy.

The conservative all-atom distance bound establishes pair coverage in every
run. However, Eg5 protein has only **1.165 angstrom** of minimum margin to the
99-angstrom cutoff, compared with about 18 angstrom for c-Met protein and
56 angstrom for water. This is a bound based on twice the maximum atom radius,
not a measured nearest omitted pair. Keep coverage checks in longer runs and
review the cutoff bound before extending their duration.

## Where the large energies come from

Use the saved post-hoc-mode pure-state difference `E2 - E1`, with Born absent.
The following are maxima of its absolute value across both setup directions,
in kilocalories per mole. They are **not free energies** or simulated total
potential-energy spikes.

| Target / leg | At weight 0.5 | At weights 0.0001 or 0.9999 |
| --- | ---: | ---: |
| c-Met / water | 52.2 | 34,858.5 |
| c-Met / protein | 100.4 | 50,062.1 |
| Eg5 / water | 40.3 | 35,444.1 |
| Eg5 / protein | 139.1 | 22,115.8 |

Saved component differences identify several distinct sources:

- c-Met reverse/protein, weight 0.0001, saved step 1942: total gap
  **50,062.14**, including **49,024.26** from perturbed-region–water LJ.
- c-Met forward/protein, weight 0.0001, saved step 156: total gap
  **7,871.88**, including **7,851.17** from LJ within the perturbed region.
- Eg5 forward/protein, weight 0.0001, saved step 989: total gap
  **22,115.75**, including **21,930.60** from perturbed-region–non-Q-solute LJ.
  In this protein leg that is the protein/environment channel, not water.

Here “Q” denotes atoms whose parameters are being perturbed. Component
attribution identifies interaction classes; identifying specific overlapping
atom pairs still requires coordinate-level inspection. The intraregion result
is not, by itself, proof of an incorrect exclusion table.

Large gaps also persist in the second half of the trajectories: for example,
the largest c-Met reverse/protein gap occurs at 1.942 picoseconds. They are not
confined to the first integration step. At small state weight, a large unscaled
energy in the weakly coupled state need not imply a comparably large sampled
Hamiltonian contribution. This explains why finite, numerically stable dynamics
can coexist with difficult cross-state energies; it does not establish overlap.

Thus **dropping exact endpoints avoids evaluating those endpoints but does not
remove near-endpoint tails**. The observations are consistent with the user's
ghost-overlap concern, extended beyond ligand–water interactions. There is no
evidence here that a Born bookkeeping spike causes those large LJ terms.

## Sign-sensitive Born accounting

For the actual forward state mapping, the computed constant change `B2 - B1`
is as follows, in kilocalories per mole:

| Target | Included non-Q protein charge, approximately (e) | Water constant change | Protein constant change |
| --- | ---: | ---: | ---: |
| c-Met | +4 | +8.12704 | +79.85806 |
| Eg5 | −6 | +8.11089 | −96.64527 |

These use actual unrounded charge tables, effective radii and the configured
dielectric 80. Included non-Q charges are not whole-protein charge labels.
Reversing the state mapping reverses the constant change. The opposite signs
in the protein legs are reproduced in the saved energy accounting.

These constants validate implementation of the selected formula, not its
physical adequacy. For a genuinely linear mixture sampled only from state-2
weight 0.0001 to 0.9999, the constant contribution is `0.9998 * (B2 - B1)`,
not the full endpoint change. Adding a constant cannot change trajectory
overlap; integrated and post-hoc modes must not both be counted as corrections.

## What to do next, staying within the original goal

1. **Audit fixed-coordinate energy and force evaluation for these actual
   no-softcore dual topologies.** Check that saved state energies reconstruct
   native energies at another lambda, and that force changes agree with the
   corresponding energy derivatives, including per-state polarization. Inspect
   the high-energy interaction classes/atom pairs and exclusions. This is a
   focused engine/analysis diagnostic, not a new sampler.
2. **Then run a replicated, equilibrated interior-lambda pilot on both targets
   and both legs**, with intermediate windows concentrated near the weakly
   coupled regions. Qualify adjacent-window overlap, time dependence and
   replica/direction agreement before a binding free-energy conclusion. Use
   appropriately requested higher concurrency; there is no need to duplicate
   every longer trajectory in both Born modes after this accounting gate.
3. **Keep the endpoint and boundary questions separate.** Do not tune the
   analytical correction to compensate for LJ tails. Label interior estimates
   truncated; do not assume the omitted endpoint contributions cancel between
   protein and water. Physical correction validation still needs the planned
   sign/radius/reference comparisons and uncertainty assessment.

No Bennett acceptance ratio (BAR) free-energy estimate was calculated from this
batch. Three widely separated weights, two picoseconds per window, preparation
under a different Hamiltonian, and unverified cross-lambda evaluation do not
support a qualified result. The existing charge-only BAR ladder requires exact
endpoints and cannot simply be fed this dual-topology interior dataset.

## Provenance and reproduction

Native/source release commit: `787e76ebf40c37b3a68dd63d1132df4ce78bba9f`.

Remote release:
`/projects/prjs2157/astra-charge-change-perturbation/releases/no-softcore-787e76eb-20260908`.

Local retrieval: `runtime/snellius-results-26482080/`. It includes complete
target inputs/outputs, but not the remote source/build directories or probe
artifacts; it is not a relocated runnable isolated build.

Remote and local `summary.json` SHA-256 (cryptographic fingerprint):
`90639ab817865bd7e7e390c71cd9f889e04077157baebc4c8f2d4c89cb7536c5`.

To reproduce the arithmetic, parse each `posthoc/states.en` with
`QligFEP.charge_completion.frames` and the corresponding `plan.json` weights.
Take column 1 of state 2 minus column 1 of state 1 for the gap. At the frame
with largest absolute gap, subtract the two state rows: columns 9, 11 and 13
are respectively intraregion, non-Q-solute and water LJ; column 6 is the total
electrostatic term, column 14 the restraint term, and columns 2–5 bonded terms.
Column indexing here is zero-based, matching the parser. Step numbers are
one-based saved-frame indices because energy is written every step, excluding
step zero and the final step. Recompute log diagnostics with
`QligFEP.charge_diagnostics.assess` at 2,000 steps and output interval 10.

This analysis changed no engine code, remote results or inputs, and submitted
no further calculations.
