# Single-H dual-topology — convergence root cause, soundness audit, and the RBFE case

Date: 2026-06-04. TYK2 `ejm_31 -> ejm_47`, protein leg, 15 Å sphere, AMBER14sb.
Builds on `03-phase2-results.md` (the ΔG-underconvergence finding). All numbers
below come from the existing 21-window `fep_chain/` foreign-λ data and the run
logs in `Q-softcore-verify/experiments/single_h/` — no new MD.

## 1. Convergence: confirmed underconverged, and 50 windows cannot fix it

Re-running the existing analyzers on the 21 uniform windows:
- BAR ΔG = 30.3, MBAR ΔG = 30.5 (multistate baseline 53.7).
- 19 of 20 adjacent windows flag POOR; overlap gaps run 5–13 kT across the whole
  path (largest on the flanks, smallest — 3.8 kT — at the λc≈0.5 apex).
- MBAR runs the full 2000-iteration cap **without converging** (`mbar_singleh.py`
  caps at line 76); the BAR≈MBAR≈30 "agreement" is therefore two unconverged
  numbers coinciding, not a convergence signal.
- Projecting gap ∝ Δλ from the measured 8.7 kT mean: ~90 windows for mean < 2 kT,
  ~140 for max < 2 kT. **50 windows lands at ~3–6 kT/window — still POOR.** More
  sampling per window does not help: overlap is set by Δλ, and debiasing a
  non-overlapping window costs ~e^gap more frames. The lever is more windows or a
  smoother path, not more time and not fewer windows.

## 2. Physical soundness audit — the implementation is correct

- **Cross-ligand exclusion fires exactly.** Runtime q–q pair count drops from 2291
  (multistate) to 1043 (single-H); the difference is 1248 = 32 × 39 = every
  ejm_31×ejm_47 cross pair, removed at neighbour-list build (`nbqqlist`). Intra-
  ligand and 1-2/1-3/1-4 exclusions are intact.
- **The boundary terms do not see the solute charge field**, so they cannot
  inflate the barrier and they are λ-independent: SCAAS `watpol` is a function of
  water dipole geometry only (never reads `qcrg`); LRF excludes q-atoms entirely
  (uses topology `crg`); radial/shell/geometric restraints are charge-free.
- **The estimator is unbiased here.** `single_h_qx_energy` (md.f90:11186) sums
  only q-atom nonbonded (q-q + q-prot + q-wat, LJ+Coulomb+Gapsys). It omits
  polarization/restraints/Born/LRF — all λ-independent at a fixed frame, so they
  cancel in `U(λ′)−U(λ)`. The Phase-2 assumption holds.
- **Config is byte-matched to multistate production** ([sphere], [solvent]
  polarisation force 20, lrf, cutoffs, shake, 20 distance_restraints fk=0.5,
  bottom 0–0.1 Å). Differences are sampling only: single-H ran 2500 steps/window
  vs production 5000, a continuous chain vs 5-stage per-window eq, fixed seed 42.
- **Charged-edge caveat (not this edge):** for ΔQ≠0 the Born `cstb`/`crgQtot`
  (md.f90:17099, 17220) is built on the combined two-half-ligand charge — it
  violates SCAAS's single-solute monopole assumption AND is frozen at the setup λ
  while the `.fl` omits it. For neutral TYK2 `crgQtot≈0`, so it is inert here, but
  it must be re-examined before trusting single-H on charge-changing edges.

## 3. Root cause of the +79 kcal barrier: convexity of charge interpolation

Per-window mean Q-energy components across the 21 windows (Qsum recovers the
multistate endpoints, −280 at λc=0 and −252 at λc=1):

```
lambda     QQ_el    QQ_vdw  Qprot_el Qprot_vdw   Qwat_el  Qwat_vdw      Qsum
  0.00   -216.76     14.45    -22.14    -43.85    -15.85     -2.53   -280.35
  0.50   -103.83     -0.53    -47.07    -12.87    -15.40      5.82   -170.33
  1.00   -188.72     15.66    -22.35    -47.72     -9.75     -6.27   -252.28
```

Going 0→0.5 the system loses ~96 kcal of stability, and one term carries it:
**Q–Q electrostatics, −217 → −104 = +113 kcal.** Q-prot and Q-wat barely move.

Q–Q is **intra-ligand** (cross-ligand pairs are excluded, so Q-Q = ejm_31-internal
+ ejm_47-internal). The mechanism is the convexity of Coulomb in charge:

- Single-H interpolates **charges** per atom, `q(λ)=(1−λ)q_A+λ q_B`. An intra-
  ligand pair has *both* charges scaled, so its Coulomb ∝ (1−λ)². At λc=0.5 each
  ligand sits at a **quarter** of its internal cohesion: ¼(−217)+¼(−189) ≈ −101 ≈
  the observed −104.
- Multistate mixes **energies**, `U=(1−λ)U_A+λU_B`, with each state's ligand at
  **full** charge. Its midpoint intra-Coulomb is ½(−217)+½(−189) = **−203**.

The gap, **−203 vs −104 = +99 kcal**, is essentially the entire ~+96 kcal barrier.
Multistate's energy-mixing is linear in λ and produces no hump; single-H's charge
interpolation is quadratic and produces one.

**Intra vs inter is the key distinction.** Only intra-ligand pairs scale *both*
charges (quadratic → convex → barrier). Ligand–environment pairs (q-prot, q-wat)
scale one charge against a fixed environment charge — linear in λ, no convexity,
no hump. The decomposition confirms it: Q-Q el humps by +113 while Q-prot el and
Q-wat el stay flat.

**Why multistate never has this:** its intermediate is `[full ligand] + [true
ghost]`. The off-ligand is pure DUM — zero charge (no intra-Coulomb to deflate)
and zero LJ (nothing to clash). Single-H's intermediate is `[two partially-real,
deflated ligands]`: partial charge gives the intra-Coulomb barrier (§3), and
partial LJ gives the dense-solvent clash that crashed the water leg
(`03-phase2-results.md`). Both symptoms come from the same choice — making a whole
vanishing ligand *partially* present instead of a clean ghost.

## 4. The RBFE case: the barrier should cancel in ΔΔG

The dominant barrier term is intra-ligand — a property of each ligand's own charge
interpolation, nearly independent of whether it sits in the pocket or in bulk
water. The same ~+113 kcal collapse therefore occurs in both legs and contributes
near-equally to ΔG_protein and ΔG_water, cancelling in ΔΔG = ΔG_prot − ΔG_water.
The poor-overlap *bias* is driven by the same intra-Coulomb swing, so to the
extent both legs are biased alike, the bias also cancels — the standard reason
RBFE tolerates imperfectly-converged legs. Caveats: cancellation is approximate
(pocket and bulk sample different ligand conformations and overlap quality), it
requires identical protocols, and the water leg must actually run.

## 5. Why single-topology is the clean fix, and a dual-topology alternative

- **Single-topology** keeps the shared core at near-full charge throughout (it
  morphs A→B rather than vanishing); only the few unique atoms scale. The
  quadratic intra-Coulomb collapse — which is maximal precisely because dual
  topology takes a *whole* ligand's charge to zero — mostly disappears, and the
  small changing region does not clash with dense solvent. This is the regime
  where single-H parameter interpolation + soft-core is genuinely the right tool.
- **Dual-topology, targeted:** scale each vanishing/appearing ligand's *intra*
  nonbonded by the energy (linear in λ, multistate-style) instead of by
  interpolated charge (quadratic), while keeping charge-interpolation + soft-core
  for the *inter*-molecular terms. Endpoints are unchanged (group-1 intra ×(1−λ),
  group-2 intra ×λ). This removes the §3 barrier and would prove the diagnosis,
  but it does not fix the inter-molecular water-leg clash (a separate symptom),
  and it edges dual-topology single-H back toward multistate for those terms.

## Status vs the two goals
- **Goal 2 (soft-core fires, stable in the pocket):** met (Phase 1).
- **Goal 1 (stable + ΔG comparable to multistate):** the dual-topology single-H
  path is correct and sound but pays a convexity penalty (intra-Coulomb barrier)
  and a clash penalty (water leg) that multistate avoids by using clean ghosts.
  These are intrinsic to making a whole ligand partially present, not bugs, and
  not schedule-fixable. Single-topology removes both at the source.
