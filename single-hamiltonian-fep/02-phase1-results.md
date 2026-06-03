# Phase 1 — single-Hamiltonian dynamics: stability and soft-core firing

Date: 2026-06-03. Branch `feature/single-hamiltonian-fep`.
Test system: TYK2 `ejm_31 -> ejm_47`, protein leg, 15 Å sphere, AMBER14sb
(arithmetic vdW), the same edge used throughout `Q-softcore-verify`.
Harness: `Q-softcore-verify/experiments/single_h/` (`run_singleh.py`,
`sweep_lambda.py`), driving the worktree `qdyn`.

## Goal of Phase 1
Per David's decision: get a provably **numerically stable** single-Hamiltonian
trajectory with the **Gapsys soft-core actually firing**, before building the
free-energy estimator (Phase 2).

## What was built (commits)
1. `Add single-Hamiltonian (parameter-interpolation) FEP path to q6` — the
   `[FEP] single_hamiltonian on` flag, `single_h_init` (per-q-atom interpolated
   ε/σ/q + decoupling factor + ligand group), cross-ligand exclusion, and the
   three arithmetic-sphere routines `nonbon2_{qq,qp,qw}_singleh`.
2. `Fix single-Hamiltonian charge scaling and Coulomb soft-core radius` — the
   two bugs below.

## Two bugs found by running the first dynamics that exercises Gapsys
The Gapsys force path is **dormant in the multistate scheme** (the soft-core
multiplies parameters that are zero on the DUM side), so single-Hamiltonian mode
is the first time it actually drives an integrator. Two latent issues surfaced:

### Bug 1 — Coulomb constant missing from interpolated charge
`prep_sim` scales `qcrg *= sqrt(coulomb_constant)` **after** `get_fep`. The first
implementation precomputed `sh_qcrg` inside `get_fep`, so it missed the factor
and the Coulomb energy came out ~332× too small (el@0 = -2.28 instead of ~-138).
Fix: call `single_h_init` at the **end of `prep_sim`**, after the scaling.
Symptom→cause→fix, verified: el@0 went to -137.53.

### Bug 2 — mismatched LJ vs Coulomb soft-core radii → collapse
With the charge fixed, `sh_gap` crashed at step 54 (SHAKE) while the no-soft-core
single-H run was clean. The trace showed a q-atom collapsing onto a water:
`Q-wat el -> -354` (unbounded Coulomb attraction) while `Q-wat vdW -> +85`
(LJ bounded by the soft-core). Cause: `r_sc_q` used the Gapsys paper value
(~0.3 Å) while `r_sc_lj` is ~3 Å. Once the LJ wall is softened the atom is
penetrable, but the **hard, unbounded Coulomb** still pulls a partner straight
through → collapse. The no-soft-core run survives because its hard 1/r¹² wall
blocks penetration.

This matters specifically because David's scheme vanishes **vdW and Coulomb
simultaneously** (an appearing atom carries partial charge while its LJ is
already soft). The fix (GROMACS `sc-coul` style): **soften Coulomb over the same
range as the LJ wall**, `r_sc_q = r_sc_lj`. Wherever an atom is penetrable its
Coulomb is bounded too. `gapsys_alpha_q/sigma_q` are consequently unused by the
single-H Coulomb radius.

(The same `r_sc_q` formula is latent in the existing multistate Gapsys routines
but never bites there because q=0 on the DUM/soft-core side — left untouched to
avoid perturbing the validated multistate baseline.)

## Results

### Regression — flag off is non-breaking
`ms_gap` (multistate, flag off) with the worktree qdyn reproduces the documented
baseline Q-SUM(state 1, step 0) = **-278.55 kcal/mol** exactly.

### Stability at the midpoint (λ=[0.5,0.5], 100 steps, from the equilibrated restart)
| run | status | Q-SUM₁@0 | Q-SUM₁@end | el@0 | vdw@0 |
|-----|--------|---------:|-----------:|-----:|------:|
| ms_gap (multistate ref) | OK | -278.55 | -286.51 | -252.33 | -26.47 |
| sh_gap (single-H, gapsys) | OK | -161.12 | -234.02 | -137.53 | -23.84 |
| sh_std (single-H, no soft-core) | OK | -161.12 | -180.69 | -137.53 | -23.84 |

`sh_gap` and `sh_std` are identical at step 0 (soft-core dormant at an
unclashed config) and **diverge during dynamics** (Q-SUM₁ -234 vs -180) — the
Gapsys soft-core is doing real work, which is impossible in the multistate scheme.

### Stability across the whole λ range (chained sweep, 200 steps/window)
Chaining outward from the [0.5,0.5] restart, each window restarting from the
previous one's coordinates:

```
up:   0.6 0.7 0.8 0.9 1.0  -> all OK, no SHAKE failures
down: 0.4 0.3 0.2 0.1 0.0  -> all OK, no SHAKE failures
```

Endpoint consistency (the key physical check):
- λc = 0.0 → Q-SUM₁ = **-280.2** ≈ multistate state-1 (lig1 fully real, lig2 absent; baseline -278.6)
- λc = 1.0 → Q-SUM₁ = **-250.9** ≈ multistate state-2 (lig2 fully real; baseline -255.6)

The single Hamiltonian reduces to the correct pure-ligand state at each endpoint
and interpolates stably between them. No endpoint `.en` overflow (there is only
one Hamiltonian, no unsampled state to evaluate).

### Soft-core necessity on this edge
On this benign congeneric edge, the no-soft-core single-H also survives gradual
ramping and even a ghost→λc=0.3 "appearing into occupied space" jump — the
"kinetic protection" (small interpolated ε × steep real LJ) is enough here. This
matches the project's prior finding that soft-core is *cosmetic* for routine
perturbations. The new value is that soft-core **now fires and is stable**, so it
is available as a real safety net for hard cases (large size changes, tight
pockets, scaffold hops) — demonstrating that load-bearing benefit needs a harder
test system and is the recommended next validation.

## Status vs the two project goals
1. *Numerically stable scheme* — **stability achieved**; "comparable ΔΔG to
   multistate" needs the Phase-2 estimator (foreign-λ BAR).
2. *Full use of soft-core* — **achieved**: Gapsys fires and is stable in a real
   FEP trajectory for the first time in this codebase.

## Next (Phase 2)
- Foreign-λ BAR estimator: emit per-frame single-H energy at neighbouring λ
  windows + a BAR analyzer, to get a ΔΔG comparable to the multistate BAR.
- QligFEP `--single-hamiltonian` wiring so `setupFEP` emits the flag.
- Extend to the geometric/OPLS routines (`nonbond_*`).
- Harder test edge to demonstrate soft-core as load-bearing.
