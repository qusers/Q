# Single-Hamiltonian FEP in Q — design and verified code map

Worktree: `Q-single-hamiltonian-fep`, branch `feature/single-hamiltonian-fep`,
based on `feature/softcore-methods-v2` (HEAD `be1c3a97`).

Goal (from David): implement §9.2 of `future-prospect-softcore.md` — a
single-Hamiltonian perturbation where one coupling parameter λ vanishes vdW and
Coulomb **simultaneously** by interpolating per-atom parameters
ε(λ), σ(λ), q(λ) and wrapping them in the **Gapsys** soft-core, so soft-core
actually fires. Keep the dual-topology layout (both ligands present, DUM
endpoints), but design the engine so it equally supports single-topology
(one merged molecule whose atoms mutate type A→B or grow DUM→B).

The free-energy result must be **numerically stable** and **comparable or better
than the current multistate (dual-topology + DUM) scheme**, validated on the
TYK2 cycle in `/Users/davidararipe/projects/Q-softcore-verify`.

---

## 1. Why the existing soft-core is dormant (recap)

Q mixes *per-state Hamiltonians*: `U = Σ_s λ_s U_s(θ_s)` with **fixed** per-state
parameters θ_s. A vanishing atom is `DUM` (ε=σ=q=0) in its ghost state, so the
soft-core offset has nothing to act on there, and QligFEP's asymmetric
`(0,20)/(20,0)` puts α=0 on the real side. Result: standard ≡ beutler_coul
bit-for-bit, gapsys equivalent — soft-core never does work
(`Q-softcore-verify/experiments/sc_bitcompare`, `CONTEXT.md §7`).

§9.2 fixes this by replacing per-state energy mixing with **per-atom parameter
interpolation inside a single Hamiltonian**:

```
ε(λ) = (1-λ)ε_A + λ ε_B,   σ(λ),  q(λ) = (1-λ)q_A + λ q_B
V_LJ^sc(r,λ) = Gapsys/Beutler soft-core wrapping 4ε(λ)[(σ(λ)/r)^12 - (σ(λ)/r)^6]
```

A disappearing atom has θ_B = DUM (ε_B=0); it interpolates smoothly to nothing,
and the soft-core bounds the singularity while ε(λ) is still nonzero — exactly
the regime soft-core was designed for.

---

## 2. Verified code map (what we will touch)

All paths relative to the worktree `src/`. Line numbers verified against the
`feature/softcore-methods-v2` tree on 2026-06-03.

### 2.1 FEP-file parsing — `q6/qatom.f90`
- **`qatom_load_atoms`** (275–606) parses `[atoms]` **and the entire `[FEP]`
  header**. Header flag dispatch is an if-chain at **309–397** with a nested
  select-case for `softcore_method` (324–338). Reads `states`,
  `qq_use_library_charges`, `softcore_use_max_potential`, `softcore_method`,
  the four `gapsys_*` reals, `offset*`.
  **→ the new single-Hamiltonian flag is parsed here** via
  `prm_get_logical_by_key` (mirror line 320).
- SC constants: `SC_STANDARD=0`, `SC_BEUTLER_COUL=1`, `SC_GAPSYS=2` (38–41).
- **`qatom_load_fep`** (936–1725) parses the body sections. Relevant:
  `[change_charges]`→`qcrg(nqat,nstates)` (1003–1044), `[atom_types]`→
  `qavdw/qbvdw(nqlib,nljtyp=3)` (1046–1076), `[change_atoms]`→`qiac(nqat,nstates)`
  (1078–1108), `[softcore]`→`alpha_max`→`populate_sc_lookup` (1621–1675),
  **`[excluded_pairs]`→`exspec(:)` with per-state flags (1167–1196)**.
- Gapsys precompute **`softcore_init_gapsys`** (1883–1993): builds
  `gapsys_lj_rsc(nqat,natyps+nqat,nstates)` and
  `gapsys_q_lambda_factor(nqat,nstates)`. r_sc = `α_lj·((26/7)σ⁶(1-λ_s))^(1/6)`.
- Gapsys evaluators **`gapsys_eval_lj`** (1997–2026) and **`gapsys_eval_q`**
  (2030–2053) — pure, take (C12,C6,r_sc,r) / (qq,r_sc,r) → (V_lin, F_lin).
  **These are parameter-source-agnostic**: feed them interpolated C12(λ)/C6(λ)
  and r_sc from σ(λ) and they just work.

### 2.2 Orchestration & MPI — `q6/md.f90`
- **`get_fep`** (1715–1774): `qatom_load_atoms` → copy charges into `qcrg`
  (1746–1750) → alloc `sc_lookup` (1753) → `qatom_load_fep` (1757) → sqrt of
  arithmetic `qbvdw` (1763–1766) → `softcore_init_gapsys` **after λ set**
  (1773–1774). **→ a new single-Hamiltonian init routine is called here.**
- **`init_nodes` MPI_Bcast** (2539–2562): every new `[FEP]` scalar must be
  broadcast here **and** in the `md_test.f90` mirror (2288–2301).
- λ plumbing: `nstates` recomputed from MD `[lambdas]` field count (3395–3429);
  `EQ(1:nstates)%lambda` set there. `EQ(:)` is `type(q_energies)`
  (`nrgy.f90:51-58`), declared `md.f90:184`.

### 2.3 Nonbonded inner loops — `q6/md.f90`
Dispatched purely by `ivdw_rule` in **`pot_energy_nonbonds`** (15045–15119):
- `nonbond_*` = **geometric/OPLS** (ivdw_rule==1): avdw=√C12, bvdw=√C6, MULTIPLY.
- `nonbon2_*` = **arithmetic/AMBER+CHARMM** (ivdw_rule==2): avdw=R*=R_min/2,
  ADD then ⁶; bvdw=√ε, MULTIPLY.

Per loop: `do istate=1,nstates` → per-state LJ+Coulomb → force `dv = (… )·EQ(istate)%lambda`;
per-state energy accumulated **unweighted** into `EQ(istate)%qq/qp/qw`; energies
λ-weighted later (14974–14993). Gapsys branch gated on
`softcore_method==SC_GAPSYS .and. sc_lookup(...)>0`, fires per channel when
`dist < r_sc`.

**Routines for the AMBER14sb / 15 Å-sphere TYK2 test (ivdw_rule==2, no PBC):**
`nonbon2_qq` (9938–10093), `nonbon2_qp` (10254–10387), `nonbon2_qw` (10543–10780).
Full feature later adds geometric (`nonbond_qq`/`nonbond_qp_qvdw`/`nonbond_qw_*`),
`_box` PBC variants, `_lib_charges`, `_spc`/`_3atom` — 24 routines total.

### 2.4 Output / qfep contract — `q6/nrgy.f90`, `q6/qfep.f90`
- **`put_ene`** (nrgy.f90:94–119) writes one unformatted record **per state**
  (`q_energies`: λ + all components + total).
- **qfep** reconstructs window energies as `veff = Σ_s λ_s·V_s` from the stored
  per-state totals and does Zwanzig/BAR/OS between adjacent windows
  (qfep.f90:400–419). **This is linear in λ.**

### 2.5 QligFEP — `QligFEP/qligfep.py`
- `read_files` (176–246) merges lig1+lig2: lig1 atoms = (real, DUM), lig2 =
  (DUM, real). **The (state1,state2) columns already are the (A,B) endpoints.**
- `write_FEP_file` (351–446) writes `[FEP]`,`[atoms]`,`[change_charges]`,
  `[atom_types]`+DUM, `[softcore]` asymmetric, `[change_atoms]`.
- CLI flag precedent: `--softcore-method` → `qligfep_cli.py` →
  `write_FEP_file` → `[FEP]` line → engine. Follow this exact wiring.

### 2.6 Verification harness — `Q-softcore-verify`
- **Binaries hardcoded to the MAIN repo** `/Users/davidararipe/projects/Q/src/q6/bin/q6/`
  (`run_verify.py:32`, every `runLOCAL.sh:7-8`). Must build the worktree's qdyn
  and repoint, or false-pass on the old engine.
- Build: `cd src/q6 && touch qatom.f90 && make qdyn COMP=osx FC=gfortran-11`
  (and `qdynp`). gfortran-11 mandatory on this Mac.
- Python side runs via `micromamba run -n qligfep_new` importing QligFEP from
  the **main repo** `src` — repoint if QligFEP changes are exercised through setupFEP.
- Smallest probe: one window `qdyn md_0500_0500.inp` (~30–60 s) or the
  `experiments/alt_dualtop` 100-step λ=[0.5,0.5] pattern (seconds) — grep logs
  for `ABNORMAL TERMINATION` (serial qdyn returns rc=0 even on crash).
- Success: `analyze_closure.py` reads BAR sum(dG) at λ=0; ΔΔG_bind = prot−wat;
  cycle sum ≈ 0. Baselines: standard −1.07, beutler_coul −0.81, gapsys +2.85.

---

## 3. Design

### 3.1 The flag (decided)
A boolean `[FEP]` header key, default off → existing behaviour untouched.
Working name **`single_hamiltonian on`** (alt: `parameter_mixing single_hamiltonian`).
Orthogonal to `softcore_method`, so the user can pick `softcore_method gapsys`
*and* `single_hamiltonian on`. EVB / standard FEP keep the per-state path.

### 3.2 λ convention (decided)
Keep `nstates = 2` for storage/restart/MPI/output compatibility. Define a single
coupling **λ ≡ EQ(2)%lambda** (state-2 weight), λ: 0→1 maps endpoint A→B.
Interpolate per q-atom between the state-1 (A) and state-2 (B) columns that
already exist (`qcrg`, `qiac`→`qavdw/qbvdw`). No new FEP storage for endpoints.

### 3.3 Engine changes (decided)
1. Parse the flag (`qatom_load_atoms`), broadcast it (`init_nodes` + mirror).
2. New init routine (called from `get_fep` post-λ) that precomputes, per q-atom:
   interpolated ε(λ), σ(λ), q(λ) and the Gapsys r_sc from σ(λ) — analogous to
   `softcore_init_gapsys` but using interpolated σ(λ) instead of the endpoint σ.
3. A single-Hamiltonian branch in the nonbonded loops: when the flag is on,
   compute the interpolated pair C12(λ)/C6(λ)/qq(λ) **once** (no `do istate`
   double count for q-atoms) and call the existing `gapsys_eval_lj/_q`. Force is
   the full single-Hamiltonian force (no `·EQ(istate)%lambda` reweight).
4. **Cross-ligand exclusion** (dual-topology only): derive ligand groups from the
   endpoints — group 1 = atoms with θ_B=DUM (vanishing), group 2 = θ_A=DUM
   (appearing), group 0 = real-in-both (shared / single-topology). Skip
   nonbonded between group 1 and group 2 (they are co-located and tethered; at
   intermediate λ both carry partial ε and would clash — soft-core bounds a
   *pair* but not the *sum* over many pairs, see `ALTERNATIVE-DUALTOP.md` S3).
   This is auto-enabled by the flag and is a no-op for single-topology
   (no two disjoint vanishing groups). Reuses Option A's `q_ligand_group` idea,
   ported into this branch (Option A itself is not in this base tree).

### 3.4 QligFEP changes (decided)
- CLI flag + thread through `qligfep_cli.py` → `write_FEP_file` `[FEP]` header.
- Keep the existing dual-topology FEP content (DUM endpoints stay). The flag is
  what changes engine behaviour. No `[excluded_pairs]` writer needed if the
  engine auto-derives groups (3.3.4); kept as an option if we prefer explicit.

### 3.5 What is NOT changing
- Bonded terms (bonds/angles/torsions) stay per-state mixed — λ only vanishes
  **nonbonded** (vdW+Coulomb), matching the brief. Disappearing atoms keep their
  bonded topology (held by per-state bonds), same as today.
- `[distance_restraints]` ligand tether stays (co-locates the two ligands).

---

## 4. Open decision — the free-energy estimator

**Fact:** BAR through the *current* qfep is mathematically impossible for a true
single Hamiltonian. qfep reconstructs window energies as `Σ_s λ_s V_s` (linear in
λ) from two stored per-state totals. Single-Hamiltonian U(θ(λ)) is **nonlinear**
in λ (soft-core + interpolated σ), so two stored endpoint energies cannot
reconstruct it. We must choose an estimator path:

- **(A) TI via analytic dU/dλ** — accumulate ∂U/∂λ during MD (chain rule through
  ε,σ,q and r_sc; Gapsys SI Eqs 1.5–1.6), integrate over the 25 λ windows with a
  small analyzer. Correct, exploits the analytic Gapsys derivative, minimal new
  I/O. Different estimator from the BAR baseline but ΔΔG is comparable.
- **(B) Foreign-λ BAR/MBAR** — at each output frame, evaluate the q-atom energy at
  the neighbouring windows' couplings and store them; extend qfep to BAR on these.
  Apples-to-apples with the baseline estimator; more engine + qfep + file-format work.
- **(C) Stability-first** — ship the dynamics (interpolation + Gapsys + exclusion)
  and prove numerical stability + that soft-core fires (method-dependent
  trajectories on short runs), then pick (A) or (B) for the production number.

Recommendation: **(C) then (A)** — get a provably stable single-Hamiltonian
trajectory with soft-core firing first (fast to validate), then add TI for the
ΔΔG, with (B) as a later apples-to-apples cross-check.

### Decision (David, 2026-06-03)
- **Estimator = (B) foreign-λ BAR/MBAR**, because the production workflow uses
  **BAR summed per λ window** and we want something analogous. **Sequencing =
  stability-first is acceptable** ("stability first to later figure out the
  analysis part, that's fine").
- **Scope = arithmetic sphere first** (`nonbon2_qq/qp/qw`), then extend to
  geometric/OPLS. **PBC `_box` variants are out of scope** — the project does not
  use PBC (spherical boundary conditions are Q's edge); a PBC use-case would go
  to other software. So we drop all `_box` routines from the plan.

Concrete consequence for foreign-λ BAR: BAR needs, per sampled frame, the
single-Hamiltonian energy at the **previous, current, and next** λ window
(3 values) — the current 2-state `.en` record cannot hold this. Phase 2 will
emit per-frame foreign-λ energies to an auxiliary per-window file and run BAR
over adjacent windows (a small analyzer, or a focused `qfep` branch), preserving
the "sum dG across each λ window" shape. Designed after Phase-1 stability lands.

---

## 5. Risks / notes
- `qcrg` is `real(4)`; interpolate charges in `real(8)` locals (precision fine).
- Arithmetic vs geometric σ handling must mirror `softcore_init_gapsys`
  (R*→σ via /2^(1/6) for AMBER/CHARMM; (C12/C6)^(1/6) for OPLS).
- `md.f90` and `md_test.f90` are parallel copies — mirror every change.
- 24 nonbonded routines exist; first working version targets the 3 arithmetic
  sphere routines the TYK2 test exercises, then extend (no silent gaps — each
  unported routine should hard-error under the flag rather than silently mix).
- Untracked `md_baseline*.f90`, `qdyn_shake_tmp.f90` are not part of the build
  target; ignore unless the makefile references them.
