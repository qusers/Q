# Initial correction and bookkeeping audit

Status: active investigation, not ready for high-performance computing (HPC)
production. This audit uses the existing Q molecular dynamics (MD) implementation.
No sampler, integrator, thermostat or solvent-boundary model was changed.

## What the correction actually contains

The Surface Constraint All-Atom Solvent (SCAAS) description distinguishes a
surface-solvent restraint from a continuum contribution representing the missing
exterior. They should not be conflated. The supplied
[King–Warshel paper](</Users/davidararipe/projects/Q/goat-Q-papers/warshel-king1985-polarization-constraints-in-molecular-dynamics-simulation-of-aqueous-solutions-the-surface-constraint-all-atom-solvent-scaas-model/hybrid_auto/warshel-king1985-polarization-constraints-in-molecular-dynamics-simulation-of-aqueous-solutions-the-surface-constraint-all-atom-solvent-scaas-model.md>)
distinguishes these contributions after its equation (1), and supplies an
angular target response to an applied field in equation (14).

The [Q methods paper](</Users/davidararipe/projects/Q/goat-Q-papers/Marelius-1999-Q-a-molecular-dynamics-program-for-free-energy-calculations/hybrid_auto/Marelius-1999-Q-a-molecular-dynamics-program-for-free-energy-calculations.md>)
describes the Born expression as the leading continuum term for the exterior,
not a complete correction for all properties of an arbitrary finite droplet.
Its spherical-boundary discussion also treats neglected, neutralized protein
charges separately in equation (9). The local paper transcriptions contain
optical character recognition errors; detailed equation claims should be checked
against original page images when ambiguity matters.

| Contribution | Current code | Consequence |
| --- | --- | --- |
| Charge-dependent angular restraint | `wat_shells` and `watpol` in [md.f90](../../src/q6/md.f90) | Changes water forces and state energies; controlled by `charge_correction` and `perstate_polarization` |
| Added exterior Born term | `init_perstate_born`, `pot_energy` and [boundary_corrections.f90](../../src/q6/boundary_corrections.f90) | Constant for a fixed state/radius/charge definition; changes state energy without changing coordinate forces |
| Historical endpoint overlap | [interior-only analysis](../endpoint-trim-analysis.md) | A separate counterfactual overlap problem, not repaired by either correction above |

Here a pure state is a complete endpoint parameter set evaluated at the current
coordinates. Lambda weights determine how those endpoint potentials combine.
The Q region is the set of atoms designated for state-dependent treatment; it
is not the entire solute or necessarily the entire enclosed charge.

## Born accounting: exact identities versus physical approximation

The current implementation uses

    B_s = -C (Q_env + q_s)^2
    C = k_e (1 - 1/epsilon) / (2 R).

`Q_env` is the sum of non-excluded, non-Q solute charges; `q_s` is the pure-state
Q-region charge. `R` is Q's effective/overridden solvent radius. The coefficient
now uses the topology's Coulomb conversion constant `k_e`. A user override of C
is available but is not a production calibration method.

For fixed definitions, B_s has no coordinate derivative. Its addition to a
state energy shifts that state's free energy by exactly B_s, independent of
sampling. Thus the state-1 to state-2 correction for **one leg** is

    delta G_B = B_2 - B_1
              = -C [2 Q_env (q_2-q_1) + q_2^2-q_1^2].

This is an algebraic consequence of the chosen term, not proof that it fully
represents the physical exterior. The protein-minus-water binding correction
requires subtracting the two leg corrections, using their own radii and charges.
An integrated Born term must not also be added post hoc.

For linear endpoint mixing over a truncated interval, the constant correction
is multiplied by the change in state-2 weight across that interval. This identity
does not restore omitted ghost-endpoint contributions. Bennett acceptance ratio
(BAR), the free-energy estimator used here, still requires a valid sampled
Hamiltonian and useful overlap; bookkeeping identities alone do not supply them.

The leading net-charge term cannot by itself establish accuracy for off-center
or spatially extended charge distributions, surface structure, or neglected
protein sites. Nor should positive and negative ligand perturbations in the
**same charged environment** be expected to have equal corrections: the
`2 Q_env delta q` cross term explicitly breaks that comparison's symmetry.

## New native state-energy checks

[state_energy_audit.f90](../../test/q6/state_energy_audit.f90) links the actual
Q energy implementation and evaluates fixed coordinates. It does not advance
dynamics. The test fixture is the existing Na/benzene/water topology, with one
real solute atom changed from charge zero to +1 or -1. Atom types, masses and
Lennard-Jones (LJ, repulsion/dispersion) interactions are identical in both states.
The remaining solute has nonzero charge, exercising the Born cross term.

Use a common nonzero angular offset of 0.03 radians in this diagnostic, frozen
across all evaluations. The existing radial wall is evaluated with its
temperature argument held at the snapshot value 298 kelvin; there is no change
to the wall or its temperature dependence in actual MD. Compare four explicit
diagnostic controls: both terms, angular term only, Born term only, and neither.
These subtraction controls are not proposed production Hamiltonians.

At state-1 weights 1, 0, 0.25, 0.5, 0.75, 0.0001 and 0.9999, the eight native
regression checks pass for both charge signs in 2.11 seconds total. They verify:

- Pure-state energies and restraints do not depend on the current mixing weight.
- The full potential and coordinate gradients obey the same endpoint mixture.
- The full endpoint potential difference matches the saved pure-state gap.
- Angular and Born contributions appear once in both the appropriate state
  records and the lambda-weighted total, including the angular output bucket.
- Born on/off leaves gradients bitwise unchanged; angular on/off does not.
- Actual binary state-energy records reproduce the audited totals/restraints.

The native Born serialization tests carried onto this branch additionally
compare actual short MD runs with Born on/off and unchanged final dynamics.
These checks establish selected bookkeeping identities on a fixture, not an
exhaustive force-gradient proof, production convergence or physical validation.

## Slide claims requiring qualification

The supplied [formulas.tex](</Users/davidararipe/projects/Q/charge-perturbation-slides/formulas.tex>)
is treated as a set of hypotheses; it has not been edited.

- **C1:** use the topology constant when describing the implemented coefficient;
  preserve the old constant only when reproducing historical calculations.
- **C3:** the old `wpol_born` variable controls the charge-dependent angular target.
  It is not the new, force-free `-C Q^2` state-energy addition. An omitted explicit
  term in BAR's stored gap does not imply a force-producing restraint has no
  effect on the sampled configurations or the resulting estimate.
- **C4/C5:** distinguish a one-leg free-energy correction from the two-leg binding
  correction; do not assume either radius or charge terms cancel.
- **C6/C8:** reported agreement and residual structure are empirical claims to
  audit from identified data with uncertainty. They do not prove the continuum
  approximation, and are not targets for parameter fitting.

## Next questions, within the approved scope

1. **Charge definition:** the Born term uses environment plus Q-region charge,
   but the angular target uses Q-region charge alone. Also, `born_dielectric` is
   configurable while the angular dielectric factor remains fixed at 0.98750
   (epsilon 80). These are verified code facts, not yet a diagnosis or permission
   to change the physical target. The subsequent [partition audit](PARTITION_AUDIT.md)
   preserves Born constants but reproduces angular energy/gradient dependence
   under fixed-charge relabeling, with a separate small non-angular control
   mismatch explicitly retained. Inspect the intended background-charge convention.
2. **Matched offsets:** the preserved historical reproduction reports different
   frozen angular offsets between forward and reverse directions. Its directional
   difference is therefore not same-Hamiltonian closure. Quantify this mismatch
   before attributing residuals to a correction formula; do not fit it away.
3. **Input/provenance controls:** distinguish effective versus requested radius,
   actual included versus excluded/neutralized charge, and integrated versus
   post-hoc corrections. Require matching definitions in the eventual protocol.
4. **Validation design:** use charge-only, matched-Hamiltonian tests with both
   signs and uncertainty/overlap checks. Do not interpret this fixture or the
   historical trimmed result as a full charged-protein validation.

No engine change is warranted from this first bookkeeping audit alone. The
historical reproduction's older suggested shell-crossing development direction
is superseded by the [approved current scope](PROPOSED_GOAL.md).

Verification checkpoint: **58 tests passed in 6.42 seconds; two optional
integration tests skipped** (the native QFEP endpoint check and the historical
raw-data reproduction have not been installed in this clean worktree). The new
eight native state-audit checks all ran. `git diff --check` passes.
