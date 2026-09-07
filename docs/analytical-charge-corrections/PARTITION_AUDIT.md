# Fixed-charge partition audit

This is a fixed-coordinate diagnostic of the existing Q molecular dynamics
(MD) implementation, not a free-energy calculation or a new sampling method.
The Q region is the set of atoms designated for state-dependent treatment.
Moving an unchanged atom into that set should not represent a physical charge
perturbation. Nevertheless, the current angular boundary target depends on that
designation, while the exterior Born term depends on the enclosed charge.

## Experiment and result

[The native audit](../../test/q6/test_q_region_partition.py) uses the same
Na/benzene/water coordinates and parameters as the
[state-energy audit](INITIAL_AUDIT.md). State 1 has charge zero on topology atom 1;
state 2 has charge -1 or +1. The additional comparison puts topology atom 13,
a sodium with unchanged charge +1 in both states, into the Q region. No atom
type, mass, coordinate or Lennard–Jones (LJ, repulsion/dispersion) parameter is
changed. Both runs have the same frozen angular offset, 0.03 radians in every
shell, and the same snapshot temperature argument for the existing radial wall.

The environment charge decreases by 1 and each Q-region charge increases by 1.
The enclosed charges and both Born constants remain identical. The angular
energies and coordinate gradients do not:

| Perturbation on atom 1 | Angular energy change, state 1 | Angular energy change, state 2 | Change in state-2 minus state-1 angular gap |
| --- | ---: | ---: | ---: |
| 0 to -1 | +1.756996 | +0.861267 | -0.895730 |
| 0 to +1 | +1.756996 | +2.652726 | +0.895730 |

Energies are in kilocalories per mole; changes mean relabeled minus original.
The largest isolated angular coordinate-gradient change is 0.445772
kilocalories per mole per angstrom. These are instantaneous values for one
configuration, not free-energy errors or estimates of the historical residual.

## Non-angular control did not pass invariance

The first version of the check required the entire non-angular potential to
agree to 1e-8 kilocalories per mole. It failed. At state 1, the relabeling changes
the aggregate LJ energy by +0.0039727933 and Coulomb energy by -0.0000005354
kilocalories per mole. The other audited energy components are unchanged at the
tested precision. The non-angular gradient changes too (maximum about 0.004176
kilocalories per mole per angstrom across these evaluations).

This is a separate numerical/control issue, not evidence that all the energy
difference is a boundary effect. The source has different arithmetic paths for
Q-water and non-Q-water interactions; for example, `nonbond_pw` multiplies stored
single-precision LJ coefficients before assigning the product to double
precision, whereas `nonbond_qw_spc` promotes coefficients before multiplication.
This is a candidate explanation, **not a demonstrated attribution of the full
control difference**. Neither kernel was changed for this audit.

The strict non-angular energy-invariance assertion remains an explicit expected
failure (`xfail`, with unexpected passes treated as failures) for each charge
sign. Component-accounting checks pass separately. No tolerance was enlarged to
declare the invariance gate passed.

To isolate the angular effect, subtract the Born-only control within each
labeling before comparing labelings:

    delta U_angular = (U_both - U_Born_only)_relabeled
                   - (U_both - U_Born_only)_original.

This agrees with the lambda-weighted change in directly accumulated pure-state
angular energies to 1e-8 kilocalories per mole at all seven tested weights. The
same control subtraction isolates a nonzero angular gradient difference that
obeys the endpoint mixture. Born-only and no-boundary controls are diagnostic
evaluations, not proposed production protocols.

## Implication for the correction

The Born bookkeeping passes this charge-partition check. The angular target's
dependence on Q-region labeling is independently reproduced after removing the
non-angular confound. This warrants resolving the intended background-field
convention before claiming a generally transferable analytical correction.
Simply replacing Q-region charge with total enclosed charge would change the
boundary forces and sampled ensemble, not merely repair a saved energy record.
It also would not prove the central-charge approximation appropriate for an
extended protein. This audit does not authorize or validate that replacement.

No engine implementation, sampler, thermostat, integrator or radial-wall change
was made. The targeted checks report **6 passed and 2 expected failures**; a
passing reproduction of the angular dependence is not a physical validation
pass. The next checks are matched frozen offsets in historical comparisons and
an explicit charge/background convention for the production protocol.

Combined checkpoint: **64 passed, 2 skipped, 2 expected failures in 8.65 seconds**
across the boundary-function, Born-serialization, native state-energy,
partition, command-line correction, endpoint-trimming and historical-reproduction
test modules. The skips remain the optional native QFEP free-energy analysis
and raw-data reproduction checks whose runtime assets are absent from this
clean worktree. All new native partition checks executed. `git diff --check`
also passed. No high-performance computing (HPC) experiment was submitted.
