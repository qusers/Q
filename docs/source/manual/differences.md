# Differences from the Q6 manual

This manual is a port of the *Manual for Q6* of 17 August 2017. The engine in
this repository identifies itself as Q 6.0.1 but does not contain everything
that manual describes, and it has features the manual does not mention. This
page lists the differences found while porting, by comparing the manual with
the keywords that `md.f90`, `qatom.f90`, `prep.f90`, `qfep.f90` and
`qcalc.f90` read.

## Described in the Q6 manual, absent from this source

| Feature | Q6 manual | In this source |
|---|---|---|
| Quantum corrections (BQCP) | Program `Qpi6`; `[QCP]` section of the `qdyn` input; `[qcp_atoms]` and `[qcp_mass]` sections of the FEP file | No `Qpi6` program and none of the sections |
| Group contributions | `[group_contribution]` section of the `qdyn` input | Not read |
| Thermostats | `thermostat` (`berendsen`, `langevin`, `nose-hoover`), `langevin_friction`, `nhchains` | Berendsen only; keywords not read |
| Integrators | `integrator` (`leap-frog`, `velocity-verlet`) | Leapfrog only; keyword not read |
| Constraint method | `constraint_method` (`SHAKE`, `LINCS`) | SHAKE only; keyword not read |
| Constraint selection | `shake_all`, `shake_heavy` | Not read |
| Temperature output | `verbose_temp` | Not read |
| Force fields | Amber95, Amber14, Amber/OPLS, CHARMM22, GROMOS87, GROMOS96, OPLS-AA, OPLS-AA/2015 in the `ff` folder | Not distributed; see [Force fields](reference/force-fields.md) |
| Organic solvent libraries | `CLF.lib`, `DCM.lib`, `MTH.lib`, `ETH.lib` | Not distributed |
| `qprep` preference `max_xlink` | Settable | Fixed at 2.1 Å |

The pages of this manual leave these out, except where noted. For the
original description, see the Q6 manual in the
[`qusers/Q6`](https://github.com/qusers/Q6) repository.

## Different in this source

| Item | Q6 manual | This source |
|---|---|---|
| Program names | `Qprep6`, `Qdyn6`, `Qdyn6p`, `Qfep6`, `Qcalc6`, `Qdum6` | `qprep`, `qdyn`, `qdynp`, `qfep`, `qcalc`, `qdum` |
| Spelling of keywords | `polarisation`, `polarisation_force`, `rigid_box_centre` | `polarization`, `polarization_force`, `rigid_box_center` |
| `bath_coupling` default | 100 fs | 10 fs |
| `separate_scaling` default | `on` | `off` |
| `shake_hydrogens` default | `on` | `off` |
| `lrf` default | `on` | `off` |
| `shell_radius` default | the outer radius (no restrained shell) | 85 % of the exclusion radius |
| `charge_correction` default | `on` | the value of `polarization` |
| Licence | GPL version 2 | GPL version 3 (`LICENSE`) |

Because `qdyn` ignores unknown keywords, an input written for Q6 with
`polarisation on` runs here without error, with the default applied instead.

## In this source, not in the Q6 manual

| Feature | Where |
|---|---|
| Energy minimisation before MD: `minimize`, `max_minimize_steps`, `minimize_tolerance`, `minimize_step_size` | `[MD]` |
| Non-equilibrium switching: `[lambda_scaling]` | `qdyn` input |
| `force_rms` | `[MD]` |
| `control_box` | `[PBC]` |
| `morse_width` | `[solvent]` |
| `improper_potential`, `improper_definition`, `coulomb_constant`, `name`, `type` | `[options]` of the parameter file |
| Thermodynamic integration, overlap sampling and BAR | `qfep` output |
| `rmsd`, `rmsf`, `fit`, geometry, `entropy`, `nonbond`, `nb_prot_qatom`, `rdf`, `com`, `com_ke` | `qcalc` |
| `average` | `qprep` |
| Additional soft-core options in `[FEP]` | Not documented here; see [Soft-core](../concepts/softcore.md) for the default potential |

## Corrections to the Q6 manual

| Place | Q6 manual | Corrected |
|---|---|---|
| Two-state ground-state energy (Eq. 7) | $(\varepsilon_1 + \varepsilon_2)^2$ under the square root | $(\varepsilon_1 - \varepsilon_2)^2$, as in `qfep.f90` |
| Free energy of the mapping potential (Eq. 5) | logarithm of a sum of averages | sum of logarithms of averages |
| Units of $R^*$ and $\varepsilon$ in the atom-type tables | kcal·mol⁻¹·Å¹² and kcal·mol⁻¹·Å⁶ | Å and kcal·mol⁻¹ |
| Unit of the Morse exponent in `[bond_types]` | Å⁻² | Å⁻¹ |
| Unit of the angle-restraint force constant | kcal·mol⁻¹·degrees⁻² | kcal·mol⁻¹·rad⁻² |
| Nosé–Hoover chain equations (Eqs. 29–31) | factor $s_j$ on the left; $m_{j-1}$ in the thermostat terms | standard form with $Q_{j-1}$ |
| Langevin noise amplitude (Eq. 32) | $\sqrt{2\gamma k_B T / m_i}$ | $\sqrt{2\gamma m_i k_B T}$ |
| Angle couplings in the proton-transfer example | text refers to angles 6 and 7 | the section holds line numbers 3 and 4 of `[change_angles]` |
