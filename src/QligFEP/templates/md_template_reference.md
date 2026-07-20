# Q MD input reference (`md_template.py`)

This document explains every parameter that `md_template.py` writes into a Q
molecular-dynamics input (`.inp`) file, grounded in the engine that consumes it:
`src/q6/md.f90` (the `qdyn`/`qdynp` binary). For each parameter it gives the
meaning, units, default, the engine variable it populates, the `MDParameters`
field that feeds it, and — where it matters — the underlying physics.

The input file is parsed by the `initialize` function in `md.f90` (section parser
starting around `md.f90:3002`). Keys are read per `[section]`; a missing optional
key falls back to a hard-coded default (the `*_default` parameters near
`md.f90:167`).

`render_md_input()` lays the sections out in this fixed order:
`[MD] → [cut-offs] → [sphere] → [solvent] → [intervals] → [files] →
[correction]? → [trajectory_atoms] → [lambdas] → [sequence_restraints] →
[distance_restraints] → [wall_restraints]`.

---

## The simulation model in one paragraph

Q runs **spherical-boundary** (non-periodic) simulations using the SCAAS model
(Surface-Constrained All-Atom Solvent). Instead of periodic images, the solute
sits at the centre of a finite water droplet. Two outer radii matter: the
**exclusion radius** `rexcl_o` (set when the system is solvated in qprep, read
from the topology) is the hard edge of the droplet; the **inner restrained-shell
radius** `rexcl_i` (= `shell_radius`) is where the soft boundary restraints
begin. Atoms in the outer shell are gently restrained so the droplet edge mimics
bulk solvent: solute heavy atoms are pinned toward their crystallographic
positions, surface waters are held radially against evaporation, and (optionally)
surface-water dipoles are oriented to reproduce the dielectric response of bulk
water. Time is integrated with a Verlet leap-frog scheme and a Berendsen
weak-coupling thermostat. FEP is driven by `[lambdas]` over the alchemical
**Q-atoms** defined in the FEP file.

---

## `[MD]` — core dynamics

Parsed at `md.f90:3002–3129`. Integration happens in the main loop labelled
"Verlet leap-frog algorithm"; the thermostat is the `temperature` subroutine
(Berendsen weak coupling); initial velocities come from `maxwell`.

| `.inp` key | `MDParameters` field | Engine var | Units | Default | Meaning |
|---|---|---|---|---|---|
| `steps` | `steps` | `nsteps` | count | required | Number of MD steps to run. |
| `stepsize` | `stepsize` | `stepsize`→`dt` | fs | required | Integration time step. Converted once to Q internal time units: `dt = 0.020462 * stepsize` (1 internal unit ≈ 1/0.020462 ≈ 48.9 fs). |
| `temperature` | `temperature` | `Temp0` | K | required | Thermostat target temperature. `"T_VAR"` is a placeholder the job script substitutes. Must be > 0. |
| `bath_coupling` | `bath_coupling` | `tau_T` | fs | 10 | Berendsen bath relaxation time. Also scaled by 0.020462; must be ≥ `stepsize`. Larger = weaker coupling / slower thermalisation. |
| `random_seed` | (eq1 only) | `iseed` | int | 1 | Seed for Maxwell-Boltzmann velocity generation. Rendered (as `SEED_VAR`) only in `eq1.inp` (`is_eq1=True`). |
| `initial_temperature` | (eq1 only) | `Tmaxw` | K | — | Temperature for the initial Maxwell velocity draw. Rendered (as `1`) only in `eq1.inp`. If absent, `iseed` is forced to 0 and a restart file becomes mandatory. |
| `shake_solvent` | `shake_solvent` | `shake_solvent` | on/off | on | Constrain all solvent bonds (rigid water). Keep on — Q uses rigid water models. |
| `shake_hydrogens` | `shake_hydrogens` | `shake_hydrogens` | on/off | off | Constrain all bonds to hydrogen (lets you use the 2 fs step). |
| `shake_solute` | `shake_solute` | `shake_solute` | on/off | off | Constrain all solute bonds. |
| `lrf` | `lrf` | `use_LRF` | on/off | off | Enable the Local Reaction Field Taylor expansion beyond the cutoff (see `[cut-offs]`). |
| `separate_scaling` | `separate_scaling` | `separate_scaling` | on/off | off (template sends on) | Couple solute and solvent to the heat bath **separately** (independent scale factors) instead of as one combined bath. Avoids the "hot solute / cold solvent" artefact. |
| `minimize` | `minimize` | `do_minimize` | on/off | off | Run FIRE minimisation **before** MD. The `minimize_*` keys below are only emitted when this is on. |
| `max_minimize_steps` | `max_minimize_steps` | `max_min_steps` | count | 1000 | Max minimisation iterations per phase. |
| `minimize_tolerance` | `minimize_tolerance` | `min_tol` | kcal/mol/Å | 0.1 | Convergence threshold on the **RMS force**: stop when `sqrt(Σf²/N) < min_tol`. |
| `minimize_step_size` | `minimize_step_size` | `min_step` | Å | 0.001 | Steepest-descent step, used only by the `minimize_method sd` baseline; ignored by FIRE. |
| `minimize_freeze_qatoms` | `minimize_freeze_qatoms` | `min_freeze_qatoms` | on/off | off | Hold the Q-atoms (both dual-topology ligands) fixed during minimisation, relaxing only the environment. Off by default: the two ligands do not interact, so relaxing them is sound. FIRE only; the `sd` baseline always freezes Q-atoms regardless. |

**Thermostat detail.** The Berendsen scale factor each step is
`λ = sqrt(1 + (dt/tau_T)·(Temp0/T_inst − 1))`, applied to velocities in the
integrator update `v ← (v − f/m·dt)·λ`, `x ← x + v·dt`. With
`separate_scaling on`, `T_inst` and `λ` are computed independently for solute
and solvent degrees of freedom.

**Minimiser detail.** FIRE (Bitzek et al., Phys. Rev. Lett. 97, 170201 (2006)):
damped MD that steers the velocity towards the force and freezes the moment the
power `F·v` turns negative, adapting the time step to the local curvature with no
line search. Two phases: phase 1 relaxes the bonded terms only (so misplaced
protons cannot fall into Coulomb wells) and phase 2 minimises the full force
field with bonds to hydrogen constrained. Every atom moves, including the
ligands and solvent, unless `minimize_freeze_qatoms` is set. `minimize_method sd`
selects the older steepest-descent minimiser (which always froze the Q-atoms and
solvent) as a baseline. Both report RMS-force convergence against
`minimize_tolerance`.

---

## `[cut-offs]` — non-bonded ranges (Å)

Parsed at `md.f90:3132–3171`. All cutoffs are **hard spherical** cutoffs applied
when building pair lists (no switching function). The pair list is rebuilt every
`non_bond` steps (see `[intervals]`).

| `.inp` key | `MDParameters` field | Engine var | Default | Governs |
|---|---|---|---|---|
| `solute_solute` | `cutoff_solute_solute` | `rcpp` | 10 | Solute–solute atom pairs. |
| `solute_solvent` | `cutoff_solute_solvent` | `rcpw` | 10 | Solute–water pairs. |
| `solvent_solvent` | `cutoff_solvent_solvent` | `rcww` | 10 | Water–water pairs. |
| `q_atom` | `cutoff_q_atom` | `rcq` | 99 | Interactions involving alchemical **Q-atoms**. |
| `lrf` | `cutoff_lrf` | `rcLRF` | 99 | Outer radius of the LRF Taylor expansion (only used when `lrf on`). |

**Why `q_atom` is 99.** A near-infinite Q-atom cutoff makes every Q-atom see the
entire environment with no distance truncation. This removes cutoff-dependent
discontinuities from the alchemical energy and keeps the perturbation
thermodynamically consistent across all λ windows — essential for FEP/TI
correctness.

**LRF.** When `lrf on`, charges between the ordinary cutoff and `rcLRF` are not
dropped; their reaction field is captured by a Taylor expansion (monopole →
octupole, the `phi0..phi3` moments) about each charge group's centroid. The
engine **requires** `rcLRF ≥ max(rcpp, rcpw, rcww)` and aborts otherwise. Q-atoms
are always handled with explicit pairs, never LRF.

---

## `[sphere]` — solute boundary shell

Parsed at `md.f90:3177–3215` (only when not running PBC). Applied each step by
the `fix_shell` routine.

| `.inp` key | `MDParameters` field | Engine var | Units | Default | Meaning |
|---|---|---|---|---|---|
| `shell_radius` | `shell_radius` | `rexcl_i` | Å | 0.85·`rexcl_o` if omitted | Inner radius of the restrained boundary shell. |
| `shell_force` | `shell_force` | `fk_pshell` | kcal/mol/Å² | 10 | Harmonic force constant of that shell. |

**What it does.** Solute **heavy** atoms lying between `rexcl_i` and the exclusion
radius `rexcl_o` (the droplet edge, fixed at solvation time and read from the
topology) are flagged as "shell" atoms and harmonically restrained toward their
topology coordinates: `E = ½·fk_pshell·|x − x_top|²`. This keeps the protein
surface from drifting or unravelling at the boundary where solvation is
incomplete. QligFEP always sets `shell_radius` explicitly, so the 0.85·`rexcl_o`
fallback is not normally exercised.

---

## `[solvent]` — water boundary restraints & charge correction

Parsed at `md.f90:3218–3263`. The radial restraint is applied in
`restrain_solvent`, the polarization restraint in `watpol`; auto-defaults that
depend on the droplet radius are set in `wat_sphere`.

### Standard SCAAS restraints

| `.inp` key | `MDParameters` field | Engine var | Units | Default | Meaning |
|---|---|---|---|---|---|
| `radial_force` | `radial_force` | `fk_wsphere` | kcal/mol/Å² | 60 | Radial restraint holding surface waters against the droplet edge (a harmonic wall outside the target radius, a shallow Morse well just inside it). Prevents surface-water evaporation. |
| `polarization` | `polarization` | `wpol_restr` | on/off | on | Enable the SCAAS surface-water **dipole-orientation** restraint. |
| `polarization_force` | `polarization_force` | `fkwpol` | kcal/mol/rad² | 20 | Force constant of that orientational restraint: `E = ½·fkwpol·(θ − θ₀)²`, where θ is the angle between a surface water's dipole and the outward radial vector. Reproduces the dielectric ordering of bulk water at the droplet surface. |

(The engine also accepts optional `morse_depth`/`morse_width` for the radial
restraint; the template leaves them at their radius-dependent auto-values, so it
does not emit them.)

### Per-state finite-sphere Born correction

For **net-charge-changing** edges, a finite water droplet under-stabilises the
charged state relative to bulk. These opt-in keys add the analytical Born
self-energy back per alchemical state.

| `.inp` key | `MDParameters` field | Engine var | Units | Default | Meaning |
|---|---|---|---|---|---|
| `perstate_born_correction` | `perstate_born` | `perstate_born` | on/off | off | Add a per-state Born self-energy to each FEP state (see formula below). Requires a spherical boundary (aborts under PBC). |
| `born_coefficient` | `born_coefficient` | `born_C_override` | kcal/mol/e² | continuum value | Override the Born coefficient `C`. Omitted (`None`) ⇒ the engine derives `C` from continuum theory. |
| `charge_correction` | — (not exposed) | `wpol_born` | on/off | = `polarization` | A *separate*, older boundary-restraint term (Coulomb stabilisation of the water shells by the enclosed Q-charge). Not an analytical correction; not set by the template. Requires `polarization on`. |
| `born_dielectric` | — (engine default 80) | `born_eps` | — | 80 | Solvent dielectric ε used in the continuum `C`. Not emitted by the template. |

**Born formula** (functions in `qatom.f90:2057–2080`,
driven by `init_perstate_born` at `md.f90:16661`):

```
C            = k_e · (1 − 1/ε) / (2R)          ! continuum coefficient, or born_coefficient override
Q_s          = Q_env + Σ_iq qcrg(iq, s)        ! net droplet charge with Q-region in PURE state s
E_born(s)    = −C · Q_s²                        ! self-energy added to EQ(s)%total each step
```

with `k_e = 332.0637`, `ε = born_dielectric` (80 by default), `R = rwat` (the
droplet radius), and `Q_env` the non-excluded, non-Q solute charge. The term is
**constant per window**, so it exerts no force and perturbs no dynamics; it shifts
each state's total energy by a fixed amount. Its state-to-state difference
`−C·(Q₂² − Q₁²)` survives into Δ*G* for charge-changing edges and cancels exactly
for neutral edges (`Q₁ = Q₂`). Reference: King & Warshel 1989; Marelius et al.
1999.

---

## `[intervals]` — output cadence (steps)

Parsed at `md.f90:3267–3337`.

| `.inp` key | `MDParameters` field | Engine var | Default | Meaning |
|---|---|---|---|---|
| `output` | `interval_output` | `iout_cycle` | 10 | Steps between energy-summary print-outs to the log. |
| `energy` | `interval_energy` | `iene_cycle` | 0 (off) | Steps between writes to the energy file. Emitted only when set (production MD); requires the `energy` file key. |
| `trajectory` | `interval_trajectory` | `itrj_cycle` | 0 (off) | Steps between trajectory frames; requires the `trajectory` file key. |
| `non_bond` | `interval_non_bond` | `NBcycle` | 10 | Steps between non-bonded pair-list rebuilds. |

---

## `[files]` — filenames

Parsed at `md.f90:3355–3445`.

| `.inp` key | Source arg / field | Engine var | Meaning |
|---|---|---|---|
| `topology` | `params.topology` | `top_file` | Topology file (required). |
| `trajectory` | `trajectory_file` | `trj_file` | Trajectory output; required if `trajectory` interval > 0. |
| `restart` | `restart_file` | `restart_file` | Input restart (coordinates+velocities). Omitted for eq1; then velocities are generated from the topology + `random_seed`. |
| `energy` | `energy_file` | `ene_file` | Energy output; required if `energy` interval > 0. |
| `final` | `final_file` | `xfin_file` | Final coordinates/restart written at the end (required). |
| `fep` | `params.fep_file` | `fep_file` | FEP definition file (the Q-atoms and their per-state parameters). `"FEP_VAR"` is a job-script placeholder. |

---

## `[correction]` — charge-correction observable logger (optional)

Emitted only when `render_md_input` is given a `correction_file`. Parsed at
`md.f90:3448–3464`; written by `write_qcorr` (`md.f90:17115`). This is a **read-only
diagnostic** — it samples geometry/charges and logs electrostatic observables but
never touches forces, energies, dynamics, or BAR.

| `.inp` key | Source | Engine var | Default | Meaning |
|---|---|---|---|---|
| `interval` | `interval_energy` or `interval_output` | `iqcorr_cycle` | 100 | Steps between samples. |
| `file` | `correction_file` | `qcorr_file` | required | Output log path. |
| `kernel` | hard-coded `1` | `qcorr_kernel` | 1 | Screening kernel; only `1` (`ε(r)=r`, i.e. a `k_e/r²` distance-dependent dielectric) is supported. |
| `exclude_last_qatoms` | `correction_exclude_last` | `qcorr_exclude_last` | 0 | Drop the last N Q-atoms from the centroid. The co-alchemical counter-water is appended last in the FEP atom list and sits far from the ligand; excluding it keeps the centroid on the ligand. 0 keeps all. |

Per frame the logger writes `step`, the per-state λ values, and two observables,
both using the `ε(r)=r` kernel summed over charged non-Q solute atoms:

- `U_obs(s) = Σ_iq Σ_j qcrg(iq,s)·crg(j)·k_e/r_ij²` — the full charge-distribution
  integrand; `U_obs(2) − U_obs(1)` is the distributed correction.
- `phi_cog  = Σ_j crg(j)·k_e/r(cog,j)²` — the potential at the geometric centroid
  of the (ligand) Q-atoms; `Δq·⟨phi_cog⟩` is the net-charge monopole correction.

---

## `[trajectory_atoms]`, `[lambdas]`

- `[trajectory_atoms]` — the template emits the mask `not excluded`, i.e. write
  every non-excluded atom to the trajectory.
- `[lambdas]` — the two λ values for this window (`lambda1 lambda2`), one per FEP
  state. Parsed at `md.f90:3468`; stored as `EQ(:)%lambda`.

---

## Restraint sections

The template passes these as pre-formatted, whitespace-flush strings (one
restraint per line). Column layouts come from the parser
(`md.f90:3506–3645`):

- **`[sequence_restraints]`** — `atom_i atom_j fc H_flag to_center`. Harmonically
  restrains the atom range *i…j* with force constant `fc`; `H_flag` includes/excludes
  hydrogens; `to_center` restrains toward the geometric centre (1) rather than the
  topology coordinates (0). Used here to pin the counter-ion water oxygen.
- **`[distance_restraints]`** — `atom_i atom_j dist1 dist2 fc state`. Flat-bottom
  distance restraint between `dist1` and `dist2` with force `fc`, active in `state`
  (0 = all states). Atoms may be given as `res:atomnr`.
- **`[wall_restraints]`** — `atom_i atom_j dist fc aMorse dMorse H_flag`. Restrains
  the atom range *i…j* to within `dist`, with optional Morse parameters. Spherical
  boundary only (ignored under PBC).

---

## Key engine anchors

| Concept | Location |
|---|---|
| Input parser (`initialize`) | `md.f90:3002–3648` |
| Default constants (`*_default`) | `md.f90:167–184` |
| Verlet leap-frog main loop | `md.f90:~4599` |
| Berendsen thermostat | `temperature`, `md.f90:4202` |
| Maxwell velocity init | `maxwell`, `md.f90:4178` |
| Steepest-descent minimiser | `md.f90:4413` |
| Solute shell restraint | `fix_shell`, `md.f90:1705` |
| Solvent radial restraint | `restrain_solvent`, `md.f90:16527` |
| Water polarization restraint | `watpol`, `md.f90:16785` |
| Per-state Born setup | `init_perstate_born`, `md.f90:16661` |
| Born coefficient / self-energy | `qatom.f90:2057–2080` |
| Correction logger | `write_qcorr`, `md.f90:17115` |

Line numbers are anchors as of the `feature/charge-corrections` branch and may
drift; the routine names are stable.
