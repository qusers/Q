# Q Fortran Compatibility Notes

## Scope

This note summarizes the current differences between QGPU and Q Fortran for the native Q input path, based on the CDK2 benchmark:

- GPU case: `/home/mcpi-02/code/qligfepv2-BenchmarkExperiments/perturbations/cdk2/2.protein/FEP_1h1q_1oiu_gpu_test/FEP1/298/1`
- Fortran case: `/home/mcpi-02/code/qligfepv2-BenchmarkExperiments/perturbations/cdk2/2.protein/FEP_1h1q_1oiu_f90_test/FEP1/298/1`

LRF is intentionally excluded from this comparison because QGPU does not support it yet.

## Current Conclusion

QGPU is closer to Q Fortran after the restart and initialization fixes, but it is not yet Fortran-equivalent. The remaining differences are not only CUDA `atomicAdd` ordering effects. Some are semantic differences in the MD loop, FEP force terms, and input/output behavior.

Small numerical differences from `atomicAdd`, parallel reduction order, and CPU/GPU operation order are expected after the semantic gaps below are fixed. They should not explain large drift in restraint, bonded, or Q nonbonded energies by themselves.

## Already Aligned

The following items were fixed or verified during this debugging pass:

- Native `.re` restart input now accepts the Fortran water-polarisation third record.
- Native `.re` restart output now writes the water-polarisation third record when water shells exist.
- `wshell(:)%theta_corr` is preserved from Fortran restarts instead of being reset to zero.
- Old two-record restarts and fresh starts still initialize `theta_corr` to zero.
- Startup SHAKE / Maxwell velocity generation / center-of-mass removal now follows the Fortran trigger: it runs only when `initial_temperature` is explicitly present with a positive random seed.
- Restart runs without `initial_temperature` should preserve restart coordinates and velocities instead of applying initial SHAKE again.
- The `fix_shell` force constant matches the Fortran fixed-atom restraint constant for the current case.

## Remaining Semantic Differences

| Priority | Area | QGPU behavior | Q Fortran behavior | Impact |
| --- | --- | --- | --- | --- |
| High | Final step handling | `main.cu` currently runs `handler.run(ctx.md.steps + 1)` | Fortran runs `nsteps` dynamics steps, then performs a final energy/output pass without another integration | QGPU can write final energies/restarts from a different state. This should be replaced with a Fortran-style final energy-only pass. |
| High | Softcore FEP terms | Softcore records are parsed and stored, but no force/runtime use was found | Fortran applies softcore lookup in Q nonbonded interactions | This affects the current CDK2 FEP because its `FEP1.fep` contains `[softcore]` entries. Q electrostatic/vdW energies may differ for semantic reasons. |
| High | CUDA Q bonded terms | CUDA Q angle, bond, and torsion calls are commented out | Fortran evaluates Q bond, angle, torsion, and improper terms for each state | Current CDK2 logs show zero Q bonded energy, but general FEP runs will not match. |
| High | Q improper terms | Q impropers are parsed, but no CPU/CUDA force implementation is active | Fortran evaluates `qimproper` | Any FEP using Q impropers will diverge. |
| High | Offdiagonal/FEP couplings | Offdiags, soft pairs, excl pairs, torsion couplings, and improper couplings appear parsed but not applied in runtime force paths | Fortran has active logic for these FEP features | Runs that use these features are not Fortran-equivalent. |
| Medium | Shell force input | `shell_force` is parsed, but pshell force uses the hard-coded constant `10.0` | Fortran uses the input shell force | Harmless for the current input because it uses `10.0`; wrong if the input changes. |
| Medium | Polarisation input | Native parsing effectively enables water polarisation | Fortran follows the input setting | Harmless for the current polarisation-on case; wrong for polarisation-off inputs. |
| Medium | Lambda state count | QGPU rejects more than two lambda states | Fortran supports more states | Current case has two states, so this does not affect the CDK2 benchmark. |
| Medium | CHARMM Urey-Bradley | No active Urey-Bradley support was found | Fortran evaluates Urey-Bradley terms for CHARMM-style inputs | Not expected to affect the current AMBER-style case. |
| Medium | Box/PBC/pressure features | QGPU native path is sphere-focused and does not match the full Fortran box/PBC/pressure feature set | Fortran supports more ensemble and boundary-condition modes | Separate from LRF; only relevant for inputs using those modes. |
| Low | Hot-atom diagnostic | QGPU prints `Ekin / Boltz / 3` | Fortran prints `2 * Ekin / Boltz / 3` | Diagnostic value differs by a factor of two; dynamics are unaffected. |
| Low | Stdout format | QGPU prints a bracketed debug-style summary | Fortran writes its native `write_out` format | Human comparison requires parsing/mapping fields; not byte-compatible output. |

## Needs Further Audit

These areas may be compatible enough for the current large-cutoff CDK2 case, but they still need a focused audit before claiming full Fortran parity:

- Nonbonded cutoff and pair-list behavior. Fortran uses pair lists updated on `NBcycle`; QGPU uses a different backend-specific scheduling model.
- Native `.en` energy-file byte compatibility with Fortran/qfep expectations.
- DCD trajectory header/frame compatibility beyond the frame-by-frame coordinate checks already done.
- Thermostat mode handling. The current case uses the expected simple temperature coupling path, but other Fortran thermostat modes are not confirmed.

## Current CDK2 Benchmark Impact

For this specific CDK2 case:

- Restart polarisation state should now match Fortran.
- Restart initialization should no longer disturb `eq2`-`eq5` / `md_*` restarts by applying initial SHAKE again.
- `shell_force = 10.0`, two lambda states, SPC water, and polarisation-on settings match QGPU's current assumptions.
- The current FEP does use softcore entries, so Q nonbonded forces/energies can still differ for a real semantic reason.
- Existing GPU logs in the benchmark directory may predate the restart and initialization fixes, so they should not be used as final evidence without rerunning.

## Recommended Fix Order

1. Replace the `steps + 1` workaround with a real Fortran-style final energy/output/restart pass.
2. Implement softcore in QGPU Q nonbonded force and energy kernels.
3. Enable and validate CUDA Q bond, angle, and torsion terms.
4. Implement Q improper terms on CPU and CUDA.
5. Implement or explicitly reject unsupported FEP features such as offdiags, soft pairs, excl pairs, and FEP couplings.
6. Make input-controlled values such as `shell_force` and polarisation mode drive the runtime behavior.
7. Validate `.en`, `.re`, and DCD output compatibility with Fortran tools.

