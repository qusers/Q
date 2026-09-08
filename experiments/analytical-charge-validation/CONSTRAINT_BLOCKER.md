# Historical blocker: water-constraint convergence

Integration update (2026-09-08): the user approved reuse of the modernized
constraint solvers in a separate branch. Its existing SHAKE repair passes the
64-case final-residual test and the unchanged sixteen-cell software campaign,
including preparation, restart transfer, charge windows and analysis. See
[integration evidence and limitations](MODERNIZATION_INTEGRATION.md).
The geometry blocker is resolved for these tests; physical validation remains open.

The remainder of this document records the **pre-integration diagnosis**, not the
current implementation or approval status.

Historical status: diagnosed; solver repair required user approval. Neither the solver,
integration scheme, thermostat nor the geometry alarm was changed in this audit.
This is a problem in Q's existing dynamics, not evidence that the analytical
charged-perturbation correction has failed or that a new sampler is needed.

## What the complete seed matrix revealed

The prospective two-radius, two-sign, two-direction, two-replica matrix was prepared
at software-test durations. Every seed used twenty existing molecular dynamics
(MD) steps. The campaign check additionally consumed the seed's native diagnostics.
One declared cell failed the existing 0.005-angstrom water-geometry drift alarm:

| Setting | Value |
| --- | --- |
| Requested/effective solvent radius | 14 / 14.04 angstrom |
| Water count | 388 |
| Endpoint | Charge -1, reverse-ladder start, replica 2 |
| Solvent orientation / velocity seeds | 758982 / 123 |
| Timestep and duration | 1 femtosecond, 20 steps |
| Initial maximum oxygen–hydrogen distance | 0.9574682169 angstrom |
| Running maximum after 20 steps | 0.9642957648 angstrom |
| Increase beyond the initialized maximum | 0.0068275479 angstrom |

The observed free/total temperature also reached about 635 kelvin during this
unequilibrated grid start. That is a transient observation, not a converged ensemble
or evidence for charge-sign asymmetry. Changing seeds, dropping this cell or
relaxing the geometry alarm to produce a successful matrix is not justified.

## Independent, charge-free reproduction

The diagnostic [native helper](../../test/q6/shake_residual_audit.f90) calls Q's
unchanged SHAKE routine on one water, without forces, charge, solvent boundaries,
a thermostat or MD propagation. SHAKE is the conventional name of the iterative
algorithm that enforces prescribed interatomic distances; it is not introduced
here as a new sampling method.

The test uses oxygen–hydrogen distance 0.9572 angstrom and hydrogen–hydrogen distance
1.5136 angstrom, matching the prepared topology. It starts with both oxygen–hydrogen
lengths exact and perturbs the water angle by at most 0.03 radians in 64 fixed,
deterministic cases. Correcting the hydrogen–hydrogen distance moves shared atoms,
so the other two constraints must be rechecked on the resulting coordinates.

The unchanged native solver reports every constraint `ready`, yet the worst
relative **squared-distance** residual is `0.0147504212636415`, compared with its
actual stored tolerance `0.0000999999974738`: about **148 times the tolerance**.
This isolates a constraint-convergence defect independently of any charged result.

In `src/q6/md.f90`, `shake` clears readiness flags once before the iteration loop.
Once a constraint is marked ready it is skipped in subsequent iterations. Later
corrections to shared atoms can invalidate it, but termination tests the retained
flags rather than all final constraint residuals. The diagnostic regression states
the desired final residual bound and is a strict expected failure until repaired.
It does not implement or substitute another constraint solver.

## Retained evidence and reproducibility

The failing cell, source snapshot, build log and executables are copied unchanged
under `runtime/constraint-failure-20260908/` in the original
`analytical-charge-corrections` worktree, not this new integration worktree, excluded from version
control. The original build is source commit `bf1a8e06`; the runtime report retains
its original paths, not invented replacement provenance. Fingerprints use SHA-256,
a cryptographic hash algorithm:

| Retained artifact | SHA-256 |
| --- | --- |
| `c11/seed/native.log` | `7c31231b543e309a26a328f89e7fb7707a22588cd400de0478a5fa66e4e7685e` |
| `c11/seed/final.re` | `3bbd97f2122c01f9a865f0924912df04f9b7933747eb0ad84cb976e5ebfca7d0` |
| `c11/prepared/system.top` | `e912049c261687ad914872bb9f059d2edc68e0a965a669cef9e6c6d7a19a36aa` |
| `build/build.json` | `5a4ed8302aaab79f0465949bbe6b34f72fbee3e59c3694780fe4279d7bd7bfb7` |

The charge-free helper and its Python regression are tracked, so reproducing the
solver defect does not depend on those local runtime assets or HPC access.
No high-performance computing (HPC) job or long preparation was run.

## Decision needed

Recommended next action: approve a **minimal repair of the existing SHAKE
convergence check**, with native final-residual regressions and replay of the same
unchanged seed matrix. This affects shared constraint dynamics beyond the analytical
correction, so it is not being folded into the branch silently. It is not a proposal
for Monte Carlo, another integrator, another thermostat, softened water geometry
or fitted boundary parameters.

Until that decision, the campaign framework remains under development and the
positive end-to-end matrix tests are explicit expected failures at the geometry
gate. These are not successful preparation or production evidence. Bookkeeping
tests remain separate; the missing physical validation is not waived.

Verification checkpoint (2026-09-08): **309 passed, 2 optional integration skips,
7 expected failures**, in 114.59 seconds. The expected failures are two previously
documented archival partition checks, four campaign tests blocked by the native
geometry defect and one direct SHAKE final-residual regression. This is not seven
successful checks. The retained copied native build also passes its source/archive/
binary fingerprint validator. No production solver change was made.
