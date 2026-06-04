# Single-Hamiltonian FEP in Q — method, rationale, and results

This is the consolidated reference. The numbered docs `01`–`06` are the detailed
investigation journey; this file is the standalone summary of what the method is,
why it is built this way, and what it achieves.

- Worktree `Q-single-hamiltonian-fep`, engine changes in `src/q6/` only.
- Validation + builder prototype in `Q-softcore-verify/experiments/single_h/`.
- Test edge: TYK2 `ejm_31 -> ejm_47`, 15 Å SCAAS sphere, AMBER14sb (gapsys soft-core).

## 1. What it is

Conventional Q FEP mixes per-state Hamiltonians, `U = Σ_s λ_s U_s(θ_s)`, with fixed
per-state parameters. A vanishing atom is a `DUM` (ε=σ=q=0) in its ghost state, so
the soft-core has nothing to act on and never does real work — `standard`,
`beutler_coul`, and `gapsys` give bit-identical trajectories.

Single-Hamiltonian FEP (`[FEP] single_hamiltonian on`) replaces that with one
Hamiltonian whose per-atom parameters are λ-interpolated, ε(λ)/σ(λ)/q(λ), wrapped
in the Gapsys soft-core. vdW and Coulomb vanish together, and the soft-core bounds
the singularity while ε(λ) is still nonzero — the regime it was designed for. This
is the first scheme in which the Gapsys soft-core is load-bearing in a Q FEP
trajectory (`02-phase1-results.md`).

The engine change is per-atom interpolation, so it serves both topologies:
- **dual topology** — both ligands present, DUM endpoints, cross-ligand exclusion;
- **single topology** — one merged molecule; shared core atoms morph A→B, unique
  atoms grow/vanish (`sh_group=0` shared, `1` vanishing, `2` appearing).

## 2. Why single topology (the load-bearing finding)

The dual-topology single-H path does **not** converge: ΔG_protein 28–30 at 21
windows (true ≈ 53.7), with **every** window poorly overlapped and a free-energy
profile that climbs to **+79 kcal/mol at λ=0.5** (`03-phase2-results.md`).

Root cause (`04-convergence-rootcause-and-rbfe.md`): the barrier is the **intra-
ligand Coulomb collapse from charge interpolation**. A per-component decomposition
shows Q–Q electrostatics carrying essentially the whole hump (−217→−104→−189 kcal
across λ); Q-prot/Q-wat stay flat. Single-H interpolates *charges*, so an intra
pair (both atoms in the same vanishing/appearing group) scales as presence² →
quarter strength at λ=0.5. Multistate mixes *energies* of full-charge states
(linear) and has no hump (½·−217 + ½·−189 = −203 vs −104; the +99 kcal gap is the
barrier). The `single_h_linear_intra` flag confirms this — it flattens the
enthalpic hump exactly — but the free energy still does not converge (BAR/MBAR
≈ 72, still poor overlap): the cost is **ramping a whole ligand's interactions**
(~217 kcal, ~18 kT/window) on *any* path, not the barrier's shape. Multistate
converges because it compares two *intact* ligands — the per-window work is the
small congeneric difference, `Δλ·(U_B−U_A) ≈ Δλ·28`, ~8× less.

Single topology gets that same advantage: the shared core stays intact and morphs
(small difference, like multistate), only the few unique atoms ramp — and those
get the soft-core. It is the regime where single-H + soft-core is both correct and
load-bearing.

## 3. Result (`06-singletopology-ddg.md`)

`ejm_31 -> ejm_47`, single topology, 21 windows, 2500 steps, single replicate:

| leg | BAR | MBAR | POOR windows | MBAR iters |
|---|---:|---:|---:|---:|
| protein | 30.2 | 29.6 | 0 / 20 | converged |
| water   | 29.7 | 29.8 | 0 / 20 | 13 |

| Method | ΔΔG (kcal/mol) |
|---|---:|
| Experimental | −0.16 |
| Dual-topology multistate, gapsys (BAR) | −0.234 |
| **Single-topology single-H, MBAR** | **−0.19** |
| Single-topology single-H, BAR | +0.54 |

Both legs converge with excellent overlap; the MBAR ΔΔG sits within ~0.05 kcal of
both experiment and the established multistate result. This succeeds where the
dual-topology forms cannot: dual-topology single-H underconverges (above) and its
water leg crashes (a co-present whole ligand clashing in dense solvent). The
absolute leg ΔG (~30) is the single-topology morph value, not the dual
whole-ligand value (53.7); only ΔΔG is comparable across topologies.

## 4. The three single-topology build requirements

These are single-topology specifics that any builder must handle (resolved in the
prototype; being productionized as `SingleTopologyFEP`):

1. **Mapping** — rdFMCS connectivity map + an *element-safe positional
   post-process* pairing remaining same-element atoms within ~0.6 Å. rdFMCS alone
   leaves the substitution-point carbons (methyl C20 / R-group C21, 0.41 Å apart)
   split into an overlapping vanish+appear pair → degenerate bonded geometry →
   NaN force in dense water. kartograf's relaxed positional mode is *not* usable
   here — it maps a methyl H onto a ring C (element-wrong).
2. **Boundary parameters** — bonded terms at the core/R-group junction mix
   shared-core (ejm_31) and appearing (ejm_47) types, which the merged prm lacks.
   Copy shared↔appearing terms from their all-X equivalents; zero the
   vanishing↔appearing cross-terms (they couple mutually-exclusive states).
3. **Water-box equilibration** — solvate the pure water sphere with the
   pre-equilibrated `water.pdb` box, not the `1 HOH` grid: grid boundary waters
   are mis-oriented and the SCAAS water-polarization restraint explodes at the eq2
   heating step (overflowed restraint energy → water SHAKE failure). Equilibrate
   the water leg at λ=1 (big group real, small group buried).

Single topology also needs **no ligand tether** (one molecule; position restraint
only), unlike the dual-topology 20-pair distance restraint.

## 5. Engine changes (in `src/q6/`, committed)

- `[FEP] single_hamiltonian on` — per-atom ε/σ/q interpolation + Gapsys soft-core +
  auto cross-ligand exclusion (`qatom.f90:single_h_init/single_h_interp`,
  `md.f90:nonbon2_*_singleh`). Flag off reproduces the baseline bit-for-bit.
- `[foreign_lambdas]` → `.fl` per-frame foreign-λ q-energy for offline BAR/MBAR
  (`md.f90:write_foreign_lambda/single_h_qx_energy`).
- `single_h_linear_intra` — diagnostic flag that linearizes the intra-ligand
  Coulomb (proves the q² root cause; off in production).
- Arithmetic/SBC sphere path only (`nonbon2_*`); PBC `_box` and geometric/OPLS not
  yet ported. MPI broadcast of the flags is not added (validated serially).

## 6. Reproduce

Engine: `cd src/q6 && touch qatom.f90 && make qdyn COMP=osx FC=gfortran-11`.

Prototype builder + runs (`Q-softcore-verify/experiments/single_h/`):
```
python3 build_singletop.py            # map -> hybrid lib/prm/pdb/fep + complex.pdb + qprep.inp
# (qprep + fix_boundary_params.py iterated to 0 missing params; offset filled)
python3 run_singletop.py 21 2500      # protein leg: eq + single-H chain
python3 run_singletop_water.py 21 2500# water leg (eq at λ=1, pre-equilibrated box)
python3 bar_singleh.py singletop 20 ; python3 mbar_singleh.py singletop 20
python3 bar_singleh.py singletop_water 20 ; python3 mbar_singleh.py singletop_water 20
```

## 7. Status and next

- **Done**: engine path, estimator, root-cause analysis, single-topology builder
  prototype, validated ΔΔG on one edge.
- **In progress**: productionize as a `SingleTopologyFEP(QligFEP)` subclass behind
  `setupFEP --single-topology`; user-facing docs.
- **Later**: replicates + remaining TYK2 cycle edges + closure; geometric/OPLS and
  PBC engine routines; MPI broadcast of the flags for `qdynp`.
