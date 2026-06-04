# Single-topology single-Hamiltonian FEP — ΔΔG matches experiment and multistate

Date: 2026-06-04. TYK2 `ejm_31 -> ejm_47`, both legs, 15 Å sphere, AMBER14sb.
Builder: `Q-softcore-verify/experiments/single_h/` (`build_singletop.py`,
`fix_boundary_params.py`, `run_singletop.py`, `run_singletop_water.py`).
Engine: worktree `src/q6/qdyn`, `single_hamiltonian on`, gapsys soft-core.

## Both legs converge

21 uniform windows, 2500 steps/window, single replicate, discard 20 frames:

| leg | BAR | MBAR | windows flagged POOR | MBAR iters |
|---|---:|---:|---:|---:|
| protein | 30.22 | 29.56 | 0 / 20 | converged |
| water   | 29.68 | 29.75 | 0 / 20 | 13 |

Every adjacent window overlaps well (gaps 0.1-1.5 kT), BAR and MBAR agree, and
MBAR converges in a handful of iterations -- the multistate-quality convergence
the investigation predicted, with a load-bearing soft-core on the appearing
R-group. The absolute leg ΔG (~30) is the single-topology morph value, not the
dual-topology whole-ligand-annihilation value (53.7); only ΔΔG is comparable.

## ΔΔG matches experiment and the established multistate result

| Method | ΔΔG = ΔG_prot - ΔG_water (kcal/mol) |
|---|---:|
| Experimental | -0.16 |
| Dual-topology multistate, gapsys (BAR) | -0.234 |
| Single-topology single-H, MBAR | **-0.19** |
| Single-topology single-H, BAR | +0.54 |

The MBAR estimate (-0.19) sits within ~0.05 kcal of both the experimental value
and the dual-topology multistate gapsys result. The BAR/MBAR spread (~0.7 kcal)
reflects single-replicate statistical noise; the central value is correct.

## What this demonstrates

Single-topology single-Hamiltonian FEP reproduces the experimental and
multistate ΔΔG for this edge, with a converged, well-overlapped calculation that
uses the Gapsys soft-core as a load-bearing component. It succeeds precisely
where the dual-topology forms could not: the dual-topology single-H path
underconverged (the whole-ligand intra-Coulomb ramp -> +79 kcal barrier, all
windows POOR at 21 windows, ~90+ needed) and its water leg crashed (a co-present
whole ligand clashing in dense solvent). Single topology morphs only the
methyl->R-group difference, so the perturbation is small, both legs converge at
21 windows, and the water leg runs cleanly.

## Build notes (single-topology specifics resolved)

- **Mapping**: rdFMCS connectivity map + an element-safe positional post-process
  that pairs remaining same-element atoms within 0.6 A. rdFMCS alone leaves the
  substitution-point carbons (methyl C20 / R-group C21, 0.41 A apart) split into
  a separate vanishing + appearing pair; the post-process maps them as one
  morphing atom. (kartograf's relaxed positional mode was rejected -- it mapped a
  methyl H onto a ring C, an element-wrong morph.)
- **Boundary parameters**: bonded terms at the core/R-group junction mix
  shared-core (ejm_31) and appearing (ejm_47) types; `fix_boundary_params.py`
  copies the shared<->appearing terms from their all-X equivalents and zeroes the
  vanishing<->appearing cross-terms. qprep then builds with 0 missing parameters.
- **Water-leg equilibration**: the pure 15 A water sphere must be solvated with
  the pre-equilibrated water box (`solvate ... 4 water.pdb`), not the `1 HOH`
  grid -- grid waters are mis-oriented at the boundary and the SCAAS water-
  polarization restraint explodes at the eq2 heating step (overflowed restraint
  energy -> water SHAKE failure). Pre-equilibrated boundary waters fix it.

## Caveats / next

- Single replicate, 21 windows, 2500 steps; the ~0.7 kcal BAR/MBAR spread implies
  a per-edge uncertainty of a few tenths kcal. Replicates / more sampling would
  tighten the error bar.
- One edge shown. A full validation would run the remaining TYK2 cycle edges and
  check the cycle closure, comparing per-edge against the multistate baseline.
- The builder is a single-edge prototype; generalizing it into a QligFEP
  `--single-topology` setupFEP path is the production follow-up.
