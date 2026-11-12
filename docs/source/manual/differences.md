# Changes from the 2017 Q manual

This page is for readers who know the Q manual of 17 August 2017 or
have input files written for it. The rest of this manual describes only what
the source in this repository does.

## Keywords that are read differently

`qdyn` ignores keywords it does not know, without an error. An input written
for that version therefore runs, but the keywords below do not have the effect they had
there.

| 2017 keyword | Here |
|---|---|
| `polarisation`, `polarisation_force` (`[solvent]`) | spelled `polarization`, `polarization_force` |
| `rigid_box_centre` (`[PBC]`) | spelled `rigid_box_center` |
| `thermostat`, `langevin_friction`, `nhchains`, `integrator`, `constraint_method`, `shake_all`, `shake_heavy`, `verbose_temp` (`[MD]`) | not read; the Berendsen thermostat, the leapfrog integrator and SHAKE are always used |
| `[QCP]`, `[group_contribution]` sections | not read |

## Defaults that differ

| Keyword | 2017 manual | Here |
|---|---|---|
| `bath_coupling` | 100 fs | 10 fs |
| `separate_scaling` | `on` | `off` |
| `shake_hydrogens` | `on` | `off` |
| `lrf` | `on` | `off` |
| `shell_radius` | the outer radius (no restrained shell) | 85 % of the exclusion radius |
| `charge_correction` | `on` | the value of `polarization` |

## Formats that differ

| Item | 2017 manual | Here |
|---|---|---|
| `[el_scale]` of the FEP file | one scale factor per pair | one scale factor per pair and state |
| `[impropers]` of the parameter file | two atom types | four atom types |
| `writepdb` in `qprep` | asks whether to write `GAP` lines | asks whether to write `TER` cards |
| `writemol2` in `qprep` | molecule number, file, hydrogen and water flags | file name only; writes the current mask |
| Program names | `Qprep6`, `Qdyn6`, `Qdyn6p`, `Qfep6`, `Qcalc6`, `Qdum6` | `qprep`, `qdyn`, `qdynp`, `qfep`, `qcalc`, `qdum` |

## Errors in the 2017 manual corrected here

| Place | 2017 manual | Here |
|---|---|---|
| Two-state ground-state energy (Eq. 7) | $(\varepsilon_1 + \varepsilon_2)^2$ under the square root | $(\varepsilon_1 - \varepsilon_2)^2$, as in `qfep.f90` |
| Free energy of the mapping potential (Eq. 5) | logarithm of a sum of averages | sum of logarithms of averages |
| Units of $R^*$ and $\varepsilon$ in the atom-type tables | kcal·mol⁻¹·Å¹² and kcal·mol⁻¹·Å⁶ | Å and kcal·mol⁻¹ |
| Unit of the Morse exponent in `[bond_types]` | Å⁻² | Å⁻¹ |
| Unit of the angle-restraint force constant | kcal·mol⁻¹·degrees⁻² | kcal·mol⁻¹·rad⁻² |
| Angle couplings in the proton-transfer example | text refers to angles 6 and 7 | the section holds line numbers 3 and 4 of `[change_angles]` |
