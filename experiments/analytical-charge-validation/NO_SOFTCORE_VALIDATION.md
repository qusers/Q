# No-softcore Eg5/c-Met validation batch

The user explicitly selected no softcore on 2026-09-08. This batch follows the
[successful midpoint compatibility check](SNELLIUS_TARGET_SMOKE.md), whose
archived inputs used Gapsys softcore. No new sampler or molecular dynamics (MD)
engine change is introduced here.

## Fixed protocol

Use the same selected ligand transformations, both protein/water legs, both
forward/reverse setup directions, and archived replica 1: eight input cases.
For each case independently start at state-2 weights 0.0001, 0.5 and 0.9999.
The other state has weight `1 - w`. Exact weights 0 and 1 are forbidden.
Reverse setup files exchange chemical state identities; a common numerical
weight therefore does not imply a common chemical state across directions.

At each weight run 2,000 steps of 1 femtosecond at 298 kelvin, twice from the
same copied restart: once with integrated Born and once without Born for
post-hoc accounting. Thus there are 24 paired checks, 48 native target runs,
and 96 picoseconds of aggregate target dynamics. The pairs are deterministic
accounting controls, **not independent sampling replicas**. The three weights
are independent short checks, not a chained free-energy ladder.

Preserve topology, atom mapping, charge tables, atom-type definitions and
changes, and shared restraints. In the active FEP (free-energy perturbation)
copy only, remove the entire `[softcore]` section, set
`softcore_use_max_potential off` and `softcore_method standard`. Q does not
accept a method named `none`: its normal Lennard–Jones (LJ, repulsion/dispersion)
path is selected by absent softcore coefficients. Require the native message
`No softcore section found. Using normal LJ potentials.` Original input files
and their hashes are preserved separately and are never modified.

Retain the previous check's repaired SHAKE constraints, direct pair coverage,
radial/angular parameters, frozen zero offsets, per-state polarization,
dielectric 80 and native effective radii. Retain restart velocities. These
archived configurations were prepared under a different Hamiltonian
(potential-energy model); two picoseconds is **not** a claim of equilibration.

## Gates and limits

Every run must finish normally with finite saved state energies, expected
weights and flags, unchanged frozen offsets, complete energy/geometry records,
no hot-atom velocity resets, and water-distance drift no larger than the
existing 0.005-angstrom alarm. Actual distance bounds must establish full pair
coverage. Report state-total and state-gap ranges without clipping outliers.

Paired final restarts must match byte for byte. Each saved state total must
differ by precisely its force-free Born constant, within the existing absolute
accounting tolerance of 1e-8 kilocalories per mole. Other saved terms must match.
If very large unscaled ghost-state energies exhaust numerical precision, retain
the failure for diagnosis; do not loosen the tolerance or call it physical
evidence against the correction without investigating the cause.

`snellius_no_softcore.sbatch` requests four CPU (central processing unit) cores,
4 gigabytes of memory, and a 3.5-hour wall-time cap on `rome`. The partition
allocates/bills at least sixteen cores. Each native target run has a fifteen-
minute timeout; at most four cases run concurrently. No retry or requeue is
enabled. Build and existing native solver checks precede the targets, as do
four 20-step synthetic probes; the complete maximum is 96,080 MD steps,
or 96.08 aggregate picoseconds. Failed cases are recorded while other cases
finish; no longer sampling job is submitted automatically.

All source archives, binaries, copied references, adapted inputs, raw outputs
and provenance reports reside in a new immutable release under
`/projects/prjs2157/astra-charge-change-perturbation`. The existing home-directory
Python environment is only read. Original runs remain read-only.

## Interpretation and next decision

Passing allows a longer, replicated interior-lambda sampling pilot; it does
not establish that the correction improves binding free energies. Failure on
either target blocks that progression and must not be hidden by selecting the
other target. Inspect extreme-weight energies and geometry separately from
Born bookkeeping. No boundary parameters are fitted to target outcomes.

Before any BAR (Bennett acceptance ratio) free-energy analysis of these
type-changing dual topologies, verify fixed-coordinate energy/force linearity
and cross-state evaluation against the actual native path. The existing
charge-only full-endpoint analysis cannot simply be reused. These short runs
produce **no free-energy estimate**. A future interior result must be labeled
truncated: omitting exact endpoints does not prove that the missing endpoint
contributions vanish or cancel between protein and water. Apply a post-hoc
Born contribution only once and only over the actual sampled interval.
