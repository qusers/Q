# Charge perturbations in Q — formula interpretations & slide flow

Companion to `formulas.tex`. Card IDs (A1, B4, C3, …) match the titled boxes in
the LaTeX sheet, so you can pull a formula and its explanation together.

Each entry gives: **what it says** (plain language), **the correct
interpretation** (what it actually means), and — where it's easy to get wrong —
**the trap** (the misreading to head off in the talk).

---

## Suggested slide flow (Foundation → bridge)

1. **What FEP computes — and that it's multistate.** `ΔG` is a *difference*
   between Hamiltonians, built as a chain of effective Hamiltonians (a mapping
   vector whose weights sum to 1) sampled over many small windows. → B1, B2
2. **Why a finite droplet at all.** Q simulates a spherical droplet, not a
   periodic box — cheaper, focused on the binding site. Price of admission: the
   boundary is artificial and must be *taught* to behave like bulk. → B2, B3
3. **SCAAS restraint #1 — density.** Keep the water density uniform to the edge.
   → A1, A3, B3
4. **SCAAS restraint #2 — polarization.** Keep the surface dipoles as disordered
   as bulk (`sin θ`), instead of over-orienting. → A2, A4, B4
5. **The catch for charged solutes.** The *correct* boundary polarization
   depends on the enclosed net charge (Langevin shift). → A5, B6
6. **The Born artifact.** A finite sphere misses the continuum polarization
   energy outside it; that missing energy scales as `Q²`. → B5
7. **How Q keeps neutral FEP honest.** Size- and charge-matching + neutralizing
   boundary residues so the `Q²` term cancels. → B7
8. **Charge-changing FEP breaks the cancellation.** `Q` differs between end
   states, so the Born term no longer cancels — an artifact of tens of kcal/mol.
   → C1, C2
9. **The per-state Born baseline.** Add the analytic Born self-energy of each
   *end state* to the energy BAR actually differences. This removes the leading
   exterior monopole term, not every charged-edge artifact. → C3, C4
10. **Use both legs and Q's actual radii.** Correct *both* legs, each with the
    effective solvent radius computed by `qprep` and used by `qdyn`. This is
    parameter-free; it strongly improves every target, but the residual range is
    about 1.5–17 kcal/mol rather than uniformly single-digit. → C5, C6
11. **Keep local and exterior terms distinct.** `phi_cog` is an in-sphere local
    potential diagnostic, not a measurement of missing exterior polarization.
    No transferable `Δq`/`φ` residual law survived the tests. → C7, C8
12. **Repair the sampled Hamiltonian.** Evaluate the charge-dependent water
    polarization restraint per alchemical state, mix its forces with `lambda`,
    and expose its endpoint energies to BAR; freeze its learned offset during
    production. Keep this independent of the Born switch. → C9

One-sentence spine of the talk: *A finite droplet needs both the missing exterior
Born self-energy and a thermodynamically consistent charge-dependent SCAAS
restraint; the first is a validated monopole baseline, while the remaining
artifacts must be separated by controlled boundary tests rather than fit away.*

---

## Part A — SCAAS polarization restraints (Warshel & King, 1985)

The problem SCAAS solves: with only a few hundred solvent molecules, the ones at
the surface don't feel the bulk they'd normally be embedded in. Left alone they
(i) thin out (low density) and (ii) over-orient their dipoles tangential to the
surface. SCAAS applies *restraints* that force the surface shell to reproduce the
average structure and polarization of an infinite system.

### A1 — Three-region partition of the potential *(Eq. 1)*
`V_tot = V^aa + V^as + V^ss + V^sb + V^ab + V^bb`

- **What it says.** Split the system into region (a) = solute + inner solvation
  shells, (s) = a thin surface shell, (b) = the bulk that extends to infinity.
  The total energy is all the intra- and inter-region interactions.
- **Interpretation.** Only (a+s) is simulated explicitly. The bulk (b) is
  replaced by *effective* terms: `V^ab` is a continuum model (this is where the
  Born-type energy lives), and **`V^sb` is the surface constraint** that makes
  region (s) behave as if (b) were still there. Everything in the talk hangs off
  this decomposition.
- **The trap.** `V^sb` is not a physical interaction with real bulk molecules —
  it's an *effective restraint* engineered to mimic their average effect.

### A2 — The surface-constraint force *(Eq. 4)*
`f_{j,n}^sb = −(K₀+ΔK_n)[Y_j − (⟨Y_j⁰⟩ + ΔY_j^n)]·h_j − μ(γ₀+Δγ_n)·Ẏ_j + A_n(t)`

- **What it says.** The force on a surface variable `Y_j` is a harmonic spring
  pulling it toward a reference value, plus Langevin friction (`−μγẎ`) and a
  random thermal kick (`A(t)`). `Y_j` is *either* a radial distance `R` *or* a
  polarization angle `θ`.
- **Interpretation.** This is the single mechanism behind both SCAAS restraints.
  The surface molecules are **sorted by magnitude** of `Y` each step, and index
  `j` labels the *j-th largest* — so the restraint acts on the *distribution*,
  not on tagged molecules (they diffuse in and out). `⟨Y_j⁰⟩` is the value that
  variable would have in an infinite system; the `Δ` corrections are tuned
  iteratively so the finite-system average matches it.
- **The trap.** It's a soft, ensemble-level restraint, not a hard geometric
  constraint on individual waters. Nothing is nailed down; the *average*
  structure is steered.

### A3 — Radial reference from the bulk RDF *(Eq. 11)*
`∫₀^⟨R_j⁰⟩ g(R)·4πρR² dR = (n₀−1) + j`

- **What it says.** Choose the target radius of the j-th surface molecule so that
  the number of molecules enclosed (obtained by integrating the *experimental*
  radial distribution function `g(R)` of bulk water) comes out right.
- **Interpretation.** This is the **density restraint**. It pins the surface
  shell at the radii that reproduce bulk water's number-density profile all the
  way to the edge, curing the surface thinning. In Q this role is played by the
  Morse + half-harmonic potential (B3).
- **The trap.** The reference comes from *bulk* structure (experimental `g(R)`),
  not from the finite droplet itself — that's what "teach it to behave like bulk"
  means concretely.

### A4 — Angular reference ⇒ the ideal `sin θ` law *(Eqs. 12–13)*
`(n_s/2)(1 − cos⟨θ_j⁰⟩) = j  ⟹  ⟨θ_j⁰⟩ = arccos(1 − 2j/n_s)`

- **What it says.** Integrating the ideal orientational distribution
  `P(θ) ∝ sin θ` gives the target angle for the j-th sorted surface molecule.
- **Interpretation.** This is the **polarization restraint**. For a *field-free*
  (neutral) system the correct distribution of the angle between a water's dipole
  and the radial direction is `sin θ` — i.e. orientations are isotropic, exactly
  as in bulk. The restraint pushes the surface toward this, killing the
  artificial tangential ordering.
- **The trap.** `P(θ) ∝ sin θ` looks "peaked at 90°," but that peak is just the
  solid-angle (`sin θ`) weighting of *random* orientations — it means *no net
  polarization*, not "waters prefer to lie flat."

### A5 — Charged-solute polarization target *(King Eq. 22)*
`θ_i⁰ = θ_i^{0,*} − (3/2)·L(X)·sin θ_i^{0,*}`

- **What it says.** When the solute carries a net charge, shift the target angle
  away from the field-free value by its Langevin response `L(X)`. Q implements
  this directly as `theta0 -= (3/2)·c·sin(theta0)`, with `c ∝ Q/r_shell²`.
- **Interpretation.** **This is the hook for the whole charge-perturbation
  story.** The physically correct boundary polarization is *not* the neutral
  `sin θ` when there's a net enclosed charge — it depends on `Q`. A finite model
  that ignores this gets the surface ensemble wrong for charged systems. The
  same monopole response also underlies the exterior Born term, but the angular
  restraint and the missing continuum self-energy are distinct Hamiltonian
  contributions.
- **The trap.** Don't present the `sin θ` restraint (A4) as "the" SCAAS
  polarization rule and stop there — for charged solutes the reference itself
  moves. A4 is the neutral limit of A5.

---

## Part B — The Born artifact & Q's boundary (Marelius et al., 1999)

### B1 — The multistate mapping potential *(Eqs. 2–3)*
`V_m = Σ_i λ_i^m · V_i,   Σ_i λ_i^m = 1`   →   `ΔG_initial→final = Σ_m ΔG_{m→m+1}`

- **What it says.** Each simulation window `m` runs on an *effective* Hamiltonian
  `V_m` that is a weighted sum of *all* `N` state potentials, with weights (the
  mapping vector `λ_m`) that sum to 1. The total `ΔG` is the sum of the small
  step-to-step FEP terms.
- **Interpretation.** Q doesn't jump `V_1 → V_2`; it walks a chain of blended
  Hamiltonians and adds up the pieces. The `Σλ = 1` constraint is what makes each
  window a genuine interpolation — endpoints recover pure states and the energy
  scale is preserved. A plain 2-state FEP is the special case `λ_m = (1−λ, λ)`,
  giving `V(λ) = (1−λ)V_1 + λV_2`. The same machinery mixes EVB resonance states,
  not just two ligands — which is why one program does both FEP and EVB.
- **The trap.** The superscript `m` labels the *window*; it is not an exponent.
  And the states are stored/evaluated **separately** — the mapping vector only
  weights them for the forces — which is exactly what makes a *per-state* term
  like the Born correction (C3) a well-defined object.

### B2 — The FEP master equation *(Eq. 1)*
`ΔG_{i→j} = −RT·ln⟨exp[−(V_j − V_i)/RT]⟩_i`

- **What it says.** The free-energy difference between two states is the ensemble
  average of `exp(−ΔV/RT)`, sampled on state `i`.
- **Interpretation.** Frame the whole talk here: everything about boundaries and
  Born terms is about making `V_i` and `V_j` comparable so the *only* thing that
  survives the subtraction is the physical change you meant to make. A boundary
  contribution that differs between `i` and `j` leaks straight into `ΔG`.
- **The trap.** Exactness of the formula is not the issue — comparability of the
  two potentials is.

### B3 — Q's radial water restraint *(Eq. 8)*
`V_water(r) = ½k_rad(r−r₀)² − D_e (r>r₀);  D_e(1−e^{a(r−r₀)})² − D_e (otherwise)`

- **What it says.** A piecewise radial potential on water oxygens: a half-harmonic
  wall outside `r₀`, and a Morse well inside that *attracts* waters toward the
  surface.
- **Interpretation.** Q's concrete implementation of the SCAAS **density**
  restraint (A3). Without the Morse attraction, surface density sags (their
  Fig 1a); the attraction pulls molecules back out and keeps density uniform to
  the boundary. `D_e` and `a` are calibrated per sphere radius.
- **The trap.** It's not just a wall to stop evaporation — the *attractive* term
  is what fixes density; a bare half-harmonic wall gives the wrong (low) surface
  density.

### B4 — The ideal boundary polarization Q restrains to
`p(θ) = ½ sin θ`

- **What it says.** The target distribution of the dipole–radial angle for
  boundary waters.
- **Interpretation.** Q's implementation of the SCAAS **polarization** restraint
  (A4): the 3 Å border is split into three subshells (0.5/1.0/1.5 Å), and within
  each, waters are sorted by angle and harmonically rotated toward their target
  in `½ sin θ`. Without it the surface over-polarizes (their Fig 1b). For a
  charged solute the target is corrected by the induced surface polarization
  (the A5 mechanism).
- **The trap.** Same as A4 — the `sin θ` shape is the *unpolarized* reference.

### B5 — The Born term — the finite-sphere artifact
`ΔG_Born = −(1/4πε₀)·(Q²/2r)·(1 − 1/ε)`

- **What it says.** The electrostatic free energy of polarizing a dielectric
  continuum (constant `ε`) *outside* a sphere of radius `r` that encloses net
  charge `Q`.
- **Interpretation.** This is the leading piece of `V^ab` (A1) — the energy the
  *explicit* droplet cannot contain, because there's no explicit solvent past the
  boundary. It's a **self-energy of the enclosed net charge**, and it scales as
  **`Q²`**. That `Q²` dependence is the crux: the term is large, but it **cancels
  exactly** between two systems that have the same radius `r` and the same net
  charge `Q`.
- **The trap.** This is not a per-atom interaction you can see in the force
  field; it's a continuum correction *outside* the simulated region. And it is
  only "harmless" *because* it cancels — that cancellation is a precondition, not
  a given (see C).

### B6 — Continuum polarization around a charge *(Fig. 1c)*
`P(r) = (Q/4πr²)·(1 − 1/ε)`

- **What it says.** The radial polarization density of a dielectric continuum at
  distance `r` from an enclosed charge `Q`.
- **Interpretation.** This is what the surface restraints are *for*: Q's Fig 1c
  shows the restrained droplet reproduces the continuum profile near its edge.
  It connects the **angular restraint** (A4/B4) to the same monopole response that
  produces the **Born energy** (B5). It does **not** mean that restraining the
  explicit shell automatically supplies the continuum energy outside the sphere.
- **The trap.** Do not collapse the restraint energy and the missing exterior
  self-energy into one term. Deviations at small `r` are also real discrete
  solvation structure, not necessarily restraint failure.

### B7 — Coulomb correction for neutralized ionic sites *(Eq. 9)*
`ΔG_corr^el = (1/4πε₀)·Σ_{p,l} q_p·q_l / (ε·r_{p−l})`

- **What it says.** Sum the Coulomb energy between ligand atoms `l` and the
  charged protein sites `p` that were *switched off* (neutralized), using a high
  dielectric `ε ≈ 80`.
- **Interpretation.** To make `Q` matchable and to respect the high-`ε` screening
  those distant sites would really feel, Q **neutralizes charged residues near
  the boundary**, then adds their interaction with the ligand back analytically.
  This is the practical "size- and charge-matching" bookkeeping that keeps the
  Born term (B5) cancelling in charge-conserving work.
- **The trap.** Neutralizing pocket residues that genuinely contribute to binding
  would corrupt the affinity — the neutralized sites must be far from the ligand,
  which is why the correction uses a high (screened) `ε`.

---

## Part C — The bridge: charge-perturbation corrections (this work)

The setup C inherits from B: for **charge-conserving** FEP the Born term (B5) is
identical in both end states and cancels — you never have to think about it. The
entire contribution of this work is what happens when it *doesn't* cancel.

### C1 — The finite-sphere Born self-energy, made concrete
`E_Born(Q) = −C·Q²,   C = k_e·(1 − 1/ε)/(2R)`

- **What it says.** B5 rewritten in Q's units: `k_e = 332.0637 kcal·Å/(mol·e²)`,
  `R` the sphere radius. Collapse the constants into a single coefficient `C`.
- **Interpretation.** Exactly the Born term of B4 — no new physics, just the
  operational form used in the engine. At `ε = 80, R = 25 Å` the continuum value
  is `C ≈ 6.56`.
- **The trap.** `C` here is the *continuum* coefficient. In practice each leg
  uses its own *effective radius* (C5), which is **not** the input sphere radius
  `rwat` — don't conflate the geometric input `R` with the effective `R` that
  enters `C`.

### C2 — The net charge enclosed in pure state *s*
`Q_s = Q_env + Σ_iq q_iq^(s),   Q_env = Σ (non-Q, non-excluded solute charges)`

- **What it says.** The total charge inside the sphere when the Q-region sits in
  a *pure* end state `s` (no λ-mixing): the fixed environment charge (protein +
  ions) plus the ligand's charge in that state.
- **Interpretation.** The correction is a function of this **per-state** enclosed
  charge. `Q_env` is computed once; the two end states differ only through the
  ligand term. For a charge-changing edge, `Q_1 ≠ Q_2` — that's the whole issue.
- **The trap.** It's the *net* enclosed charge that matters, environment
  included — not the ligand's charge alone. The `−2C·Q_env·Δq` cross term (C4) is
  precisely why the *environment* charge drives the artifact.

### C3 — The per-state Born addition — monopole baseline
`E_Q(s) += E_Born(Q_s) = −C·Q_s²`

- **What it says.** Add each end state's Born self-energy to *that state's*
  energy.
- **Interpretation.** The constant Born self-energy belongs in the **per-state**
  energy `E_Q(s)` that BAR differences. It is the leading, parameter-free
  missing-exterior monopole correction. The benchmark shows that it removes most
  of the 40–120 kcal/mol artifact but is not a complete residual model.
- **The trap.** Q's existing `charge_correction`/`wpol_born` path is not an
  explicit Born self-energy in the wrong accumulator. It changes the target and
  forces of the water-polarization restraint, so it changes sampling. The actual
  inconsistency is that this target depends on `lambda` while its corresponding
  pure-state restraint energies are absent from `E_Q(s)`.

### C4 — What the correction does to ΔΔG
`ΔΔG_Born = −C(Q_2² − Q_1²) = −2C·Q_env·Δq − C(q_2² − q_1²)`

- **What it says.** The correction's net contribution is the difference of the
  two end-state Born energies, which expands into a term linear in `Δq` plus a
  small quadratic term.
- **Interpretation.** Two properties fall straight out and are exactly what you
  want:
  - **Neutral edges** (`Q_1 = Q_2`) → **0**. The thousands of charge-conserving
    edges are untouched — safe to leave the correction on, and neutral golden
    tests stay byte-stable.
  - **Charge edges** → dominated by `−2C·Q_env·Δq`: linear in both the
    environment charge and the ligand charge change. Fully **parameter-free**
    (given `ε` and `R`).
- **The trap.** The headline term depends on `Q_env`, so the *sign and size* of
  the correction are pocket-specific — a net-negative pocket and a net-positive
  pocket get opposite corrections for the same `Δq`.

### C5 — The defensible analytic baseline: two legs, two radii
`ΔΔG_corr = ΔG_Born^protein(R_p) − ΔG_Born^water(R_w),   R_w ≈ 24.7 > R_p ≈ 23.0 Å`

- **What it says.** Apply the per-state Born (C3) to *both* FEP legs — the
  solvated complex and the free ligand in water — and subtract, each leg using its
  **own effective radius**.
- **Interpretation.** This is the form to use as the analytic baseline. Two
  independent facts support it:
  - **Two legs, subtracted.** Only the difference is physical; the bulk part
    common to both legs cancels. (Charge-conserving edges → both legs identical →
    zero, as in C4.)
  - **Use the radius Q actually uses.** `qprep` derives an effective solvent
    radius from the water volume and displaced protein heavy atoms, stores it in
    the topology, and `qdyn` restrains that radius. The protein leg is therefore
    usually smaller than the nominal sphere and the water leg. This is an engine
    geometry, not a fitted radius inferred from experimental affinities.
- **The trap.** A **protein-leg-only** variant (drop the water term) wins on some
  positive-`Q_env` targets but is **catastrophic on eg5** (`−5`) — it is not a
  law. The earlier "charge topology / φ-sign picks the winner" framing was
  **retracted**. Use the symmetric two-leg form for both signs.

### C6 — Validation: the leading artifact is removed, residuals remain
`raw artifact 40–120 kcal/mol  ──(two-leg per-leg Born)──▶  RMSE ~ 1.5–17`

- **What it says.** Post-hoc, two-leg per-leg Born strongly improves every tested
  target, but not uniformly: eg5 1.5, cmet_r35 2.7, tyk2 3.3, cmet 3.9, jnk1
  4.0, jak1 8.2, and ptp1b 17.0 kcal/mol RMSE.
- **Interpretation.** The evidence the Born *form* is right is that the **raw
  artifact itself is linear in `Q_env` at exactly the predicted slope `−2C`**
  (`R² = 0.9994`). The correction is parameter-free — no scale fit — yet lands
  the leading term close to experiment without fitting a scale.
- **The trap.** The clean `Q_env` scaling validates the leading monopole *form*;
  it does not establish that all remaining charged-edge error is another Born
  term or that the correction is ready to be called final.

### C7 — The `phi_cog` diagnostic: local, not exterior
`phi_cog = Σ_{j∈env} k_e·q_j / ‖r_cog − r_j‖²    ⟹    ΔG_monopole ≈ Δq·⟨phi_cog⟩_λ`

- **What it says.** The engine logs the screened electrostatic potential that the
  explicit in-sphere environment exerts at the ligand centroid. The `1/r²`
  kernel follows Q's distance-dependent dielectric `ε(r) = r`.
- **Interpretation.** It is a useful, force-free **local interaction diagnostic**.
  It is not the measure counterpart of C3: C3 concerns continuum polarization
  *outside* the sphere. The failed correlations are still informative because
  they rule out using this local potential as a transferable residual correction.
- **The trap.** The logged columns are pre-scaled by powers of `√k_e` internally;
  they are *proportional* diagnostics to be calibrated downstream, not
  directly-comparable absolute energies. Also note: this centroid `phi_cog`
  (screened, `1/r²`) is a different quantity from the *unscreened* `φ_pocket`
  (`1/r`) that was tested — and refuted — as the missing residual term.

### C8 — The residual is scatter, not a law (honest limits)
`residual ~ per-edge scatter + per-target offset   (no predictable Δq or φ term)`

- **What it says.** After two-leg per-leg Born, a target-dependent residual
  remains, and it has **no** reliable transferable functional form.
- **Interpretation.** Two-leg per-leg Born is the parameter-free **baseline**. The
  multi-session search for a "missing monopole / Galvani / `φ_pocket` term" ended
  **negatively**:
  - **No `Δq`-linear law.** tyk2 (7 well-sampled edges) shows ~8 kcal of residual
    scatter at a *single* `Δq = −1`; the earlier per-target "slopes" were
    single-`Δq=+1`-edge leverage (jak1's slope collapses 7.09 → 0 on dropping one
    edge).
  - **`φ_pocket` refuted** by two leverage-free arguments: eg5's environment is
    net `−5` (so `φ_pocket ≈ −135`) yet its residual slope is small and *positive*
    like every cation target; and for the identical cmet system, growing the
    sphere 25 → 35 Å makes `φ_pocket` *rise* while the residual *falls* — they move
    oppositely.
  - **Not established:** a universal residual sign or monotonic radius law. The
    current dataset is too heterogeneous to promote either observation into a
    correction formula.
  - **Mechanistic candidates to isolate:** higher boundary multipoles/off-centre
    response, a water-model surface potential, and the thermodynamic inconsistency
    in the charge-dependent polarization restraint. None justifies a fitted
    per-target `φ` or radius.
- **The trap — and the point of principle.** Report the residual as an
  **uncertainty band**; do **not** fold it into the mean as if corrected, and do
  **not** back-solve an "effective radius" (or effective `C`) and present it as a
  measured physical quantity. State the continuum `C`, the per-leg effective
  radii, and the residual gap honestly.

### C9 — Make the charge-dependent restraint thermodynamically visible
`H_m(x) = Σ_s λ_s^m [V_s(x) + W_pol^(s)(x) − C·Q_s²]`

- **What it says.** Treat the boundary polarization restraint just like every
  other alchemical Hamiltonian component: evaluate a pure-state restraint energy
  `W_pol^(s)` for every state, mix its forces with the window weights, and write
  each value into the per-state energy consumed by BAR.
- **Interpretation.** Current Q uses a lambda-mixed Q charge to move one global
  restraint target. The restraint therefore changes with lambda and affects the
  sampled ensemble, but BAR does not see the corresponding endpoint energy
  difference. Endpoint-resolved restraint targets remove that mismatch. The
  analytic Born self-energy stays independently switchable so its value can be
  tested rather than assumed.
- **The trap.** The learned `theta_corr` offset currently adapts every 100 steps.
  Production free-energy windows should load the equilibrated offset from the
  restart and freeze it; otherwise the Hamiltonian is explicitly time-dependent.
