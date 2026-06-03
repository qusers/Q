# Phase 2 — free-energy estimator and the ΔG-convergence finding

Date: 2026-06-03. TYK2 `ejm_31 -> ejm_47`, 15 Å sphere, AMBER14sb.
Harness: `Q-softcore-verify/experiments/single_h/` (`run_fep_singleh.py`,
`run_fep_chain.py`, `run_fep_water.py`, `bar_singleh.py`, `mbar_singleh.py`).

## What was built
- **Engine**: a `[foreign_lambdas]` MD-input section. At each energy-output frame
  the engine re-evaluates the single-Hamiltonian q-atom nonbonded energy at every
  listed coupling and writes `<ene_file>.fl`. Only the q nonbonded energy depends
  on λ (bonded = full-strength topology bonds at both ends; environment +
  restraints λ-independent), so `U(λ')−U(λ) = qx(λ')−qx(λ)` and these columns are
  exactly the reduced-potential differences BAR/MBAR need.
- **Analyzers**: pure-Python `bar_singleh.py` (adjacent-window Bennett BAR, with a
  forward/reverse overlap diagnostic) and `mbar_singleh.py` (all-states MBAR).
- Self-consistency verified: the self-coupling `.fl` column equals the force-path
  `Q-any` energy.

## The result: ΔG does not converge to the multistate value
Multistate baseline (gapsys, 25 sigmoidal windows): **ΔG_protein = 53.7 kcal/mol**.

Single-Hamiltonian estimates of the same ΔG_protein (same endpoints):

| schedule | sampling | BAR | MBAR |
|----------|----------|----:|-----:|
| 25 sigmoidal windows, multistate restarts | 1500 steps/window | 29.2 | 6.6 |
| 21 uniform windows, chained equilibration | 2500 steps/window | 30.3 | 30.5 |

The estimates are **schedule- and estimator-dependent (6.6–30) and all well below
54**. The uniform schedule is the most self-consistent (BAR ≈ MBAR ≈ 30), the
sigmoidal one is unreliable (its sparse middle — 0.156→0.25→0.5→0.75 — gives
almost no overlap for the nonlinear single-H path; MBAR collapses to 6.6).

## Root cause: a high-free-energy intermediate the estimators can't bridge
The single-H free-energy profile (uniform schedule, `f_k` from MBAR) **rises to
≈ +79 kcal/mol at λc = 0.5, then falls to +30 at λc = 1** — non-monotonic, unlike
the monotonic multistate profile (0 → +53.7). Per-window forward/reverse work
means are **+17…+24 kT vs −6…−16 kT**, i.e. overlap gaps of **5–15 kT per
window** (BAR wants ≈ 1–2 kT). With such poor overlap the exponential-averaging
estimators are biased low and do not recover the true ΔG.

Why the +79 kcal intermediate? At λc = 0.5 **both** ligands are simultaneously
present at half-ε (co-located, tethered, cross-ligand-excluded), so the
environment must accommodate two partial ligands at once — a strained,
high-free-energy state. ΔG is a state function, so the true single-H ΔG equals
the multistate 53.7 (the endpoints are bit-identical — at λc=0/1 the absent
ligand is pure DUM, so cross-ligand exclusion is a no-op there). The single-H
*path* simply detours over a large barrier that needs far more windows/sampling
to converge than the multistate linear-mixing path.

## Water leg + ΔΔG: not yet obtained
The plan was to test whether the offset cancels in ΔΔG = ΔG_prot − ΔG_water. The
chained single-H **water leg crashed at λc = 0.1** (SHAKE failure, step 413):
in dense solvent the appearing ligand clashes with many waters at once and the
per-pair-bounded soft-core cannot bound the *sum* over many pairs (the "Layer B"
limit from `CONTEXT.md §3.2`). The protein pocket has far fewer close waters, so
the protein leg chained cleanly. So ΔΔG is currently **inconclusive**.

## Honest status against the goals
- **Goal 2 (soft-core fires, stable):** met for the protein pocket (Phase 1).
- **Goal 1 (numerically stable + ΔG comparable to multistate):** *not yet*.
  - Stability: holds in the pocket across all λ; **fails in dense solvent**
    (water leg) at the appearing endpoint.
  - ΔG: the dual-topology single-H path has a high-barrier intermediate that
    makes ΔG convergence impractical with feasible sampling; estimates (6.6–30)
    are below the multistate 53.7 and not converged.

## Interpretation and options (for David)
The mechanics are correct (stable in pocket, soft-core fires, estimator
self-consistent, endpoints identical to multistate). The blocker is **the
dual-topology single-H path itself**: keeping both ligands co-present at
intermediate λ creates a strained, high-free-energy intermediate that is
expensive to sample and clash-prone in dense water.

1. **True single-topology** (one merged molecule; shared atoms morph A→B, unique
   atoms grow/vanish) removes the dual presence entirely — no two-ligands-at-once
   intermediate, no cross-ligand exclusion, a much lower barrier, and standard
   convergence. The engine already supports this (the interpolation is per-atom
   and `sh_group=0` shared atoms are handled); it needs a QligFEP single-topology
   builder. **This is the most promising path and what §9.2 ultimately points to.**
2. **More sampling / many more windows** for the dual-topology path (≥ the +79
   kcal barrier implies ~100+ windows for 1 kT/window). Likely impractical.
3. **Gentler water-leg protocol** (finer λ near the appearing endpoint, longer
   equilibration, stronger/decoupled soft-core in dense solvent) to fix the
   crash, then re-test ΔΔG cancellation — worth a try but does not address the
   convergence barrier.

Recommendation: pursue (1) single-topology, which is where the value is — it both
fixes convergence and is the regime where soft-core is genuinely load-bearing.
