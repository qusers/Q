# Local native timing: cost evidence, not production validation

Measured 2026-09-08 on macOS 26.5.2, arm64, using serial Qdyn compiled with GNU
Fortran 11.5.0. Molecular dynamics (MD) used 298 kelvin, 1-femtosecond (fs)
steps, zero frozen offsets, state-2 weight zero, and the settings in the
[probe tool](../../src/QligFEP/charge_probe.py). Each run starts from unequilibrated
grid water; changing its velocity seed does not make this a physical replica test.

The machine-readable [evidence record](LOCAL_TIMING.json) contains exact times,
engine/topology/input/output fingerprints, water counts and loaded radii.
Fingerprints use SHA-256, a cryptographic hash algorithm. Original inputs, native
logs, binary energies, final restarts and preparation/timing reports are retained
locally under `runtime/local-timing-20260908/`, excluded from version control.
The committed record is not a replacement for those raw files; archive both for
any later reproducibility handoff. No historical or neutral-calibration data was
modified. The runs do not establish clean source-to-executable build provenance.

## Results

Wall times include process startup, pair-list rebuild every step, energy output
every 10 steps and text diagnostics every 10 steps. Throughput is simulated
nanoseconds per day (ns/day), not effective independent sampling per day.

| Grid radius (angstrom) | Loaded radius | Waters | Steps per run | Wall seconds, three runs | Median ns/day |
| --- | --- | --- | --- | --- | --- |
| 10 | 10.15 | 146 | 1,000 | 0.928, 1.098, 1.041 | 82.96 |
| 14 | 14.04 | 388 | 1,000 | 4.439, 4.731, 4.647 | 18.59 |
| 22 | 22.04 | 1,502 | 200 | 13.497, 13.771, 14.447 | 1.25 |

The 1,000-step 22-angstrom attempt hit the tool's 60-second safety timeout. It
has no passing timing report and is excluded from the completed-run medians;
its failed-attempt log is retained and fingerprinted. The three shorter checks
completed normally. The timeout is not a measured physical instability or proof
of an engine failure. An earlier negative-endpoint attempt completed natively
but exposed an error in the new Python checker: Q's empty off-diagonal energy
record has zero payload bytes. The checker was corrected, its rejection tests
were added, and the endpoint rerun passed. That first attempt is not counted as
a passing report either. No Q energy format or force was changed to satisfy it.

These timings are too short to remove startup/noise effects or measure correlated
sampling. No error bar on long-run throughput is claimed. A loaded computer,
different compiler, constraints, output schedule or parallel configuration may
change them. Charged endpoints and equilibrated configurations also need timing
on the intended hardware. The fresh topology's native Coulomb constant is 332.0,
not the archived fixture's 332.0716; derive the correction from the actual loaded
value rather than copying a historical constant.

## Implication for the capped pilot

The [Stage B plan](PILOT_PROTOCOL.md) caps aggregate trajectory at 3.44 ns per
radius: eight ladders, each at most 100 picoseconds (ps) preparation plus eleven
windows of 10 ps discarded settling and 20 ps diagnostic sampling. At the median
local rates, `serial hours = trajectory_ns / (ns/day) * 24` gives:

| Included radius | Aggregate trajectory | Estimated local serial hours |
| --- | --- | --- |
| 10-angstrom grid | 3.44 ns | 1.00 |
| 14-angstrom grid | 3.44 ns | 4.44 |
| Initial two-radius total | 6.88 ns | 5.44 |
| Optional 22-angstrom addition | 3.44 ns | 65.80 additional |

A factor-of-two **planning contingency**, not a statistical upper confidence
bound, gives roughly 11 serial hours for the two-radius 1-fs pilot. Repeating
the entire physical duration at 0.5 fs would approximately double that repeated
run's compute to 22 hours with contingency; running both complete timestep
campaigns would sum to about 33 hours. A smaller targeted timestep assessment
should be explicitly budgeted separately instead of silently adding it.

These are **local estimates, not requested or approved allocations**. Before a
high-performance computing (HPC) request, perform a short same-input measurement
on the intended build/hardware and convert aggregate steps to the scheduler's
actual allocation units. A central processing unit (CPU) core-hour allocation,
whole-node charge and elapsed time for concurrent jobs are not interchangeable.
Do not assume linear parallel speedup or price the run from unrelated sampling
code. No HPC job was submitted.

The cost separation argues for the smaller feasibility pilot first. It does not
make two radii sufficient to validate an exterior correction; the larger-domain
test remains necessary before the corresponding physical claim. If even the
small pilot cannot provide stable, overlapping data at the required precision,
report that limitation and revise the experiment with the user. Do not recover
affordability by hiding endpoint intervals or tuning charged results.

## Verification checkpoint

The expanded focused suite passes **156 tests**, with two optional integration
skips and the same two explicitly documented expected failures on the archival
water-kernel partition fixture. Nineteen added tests exercise new-input preparation,
seed reproducibility, both charge signs at weights 0/0.5/1, finite native energies,
unchanged real Lennard–Jones (LJ, repulsion/dispersion) probe parameters, saved state
mapping, frozen offsets, input tampering, overwrite refusal and runtime limits.
The suite took 30.55 seconds at this checkpoint. Passing these software tests is
not evidence of equilibrium convergence or physical charge-correction accuracy.

The existing campaign preflight still rejects the new shared position restraint;
its bounded support/native verification must be implemented before this probe
can enter a production launch manifest. All generated reports continue to say
`production_ready: false`.
