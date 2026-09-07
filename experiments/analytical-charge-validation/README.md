# Charged-correction validation package — under development

Purpose: validate analytical charged-perturbation corrections in Q's existing
molecular dynamics (MD) with finite spherical Surface Constraint All-Atom Solvent
(SCAAS) boundaries. This directory is the production-package entry point, but
**it is not yet qualified for high-performance computing (HPC) production**.
No new sampler, integrator, thermostat, radial wall or solvent model is included.

## Staged-input preflight

The read-only [checker](../../src/QligFEP/charge_protocol.py) inspects actual
inputs and existing restart files before launch. From the clean worktree, use
the Python environment that imports this worktree's QligFEP:

```sh
PYTHONPATH="$PWD/src" python -m QligFEP.charge_protocol /path/to/manifest.json
```

Success prints a JSON (JavaScript Object Notation) report with
`gate: staged_input_consistency_passed` and **`production_ready: false`**.
Failure exits with status 2 and a diagnostic. It never starts Q, submits a job,
edits an input, copies a restart, or overwrites an output. The report records
actual file content fingerprints using SHA-256, a cryptographic hash algorithm.

The manifest has these exact top-level fields:

| Field | Required contents |
| --- | --- |
| `schema_version` | Integer `1` |
| `engine` | `binary` path, its pinned `sha256`, and full 40-character `source_commit` |
| `series` | Nonempty list of the series objects described below |

Each series declares `id` (unique), `system` (shared physical system identifier),
`sign` (-1 or +1), `direction` (`forward` or `reverse`), positive integer `replica`,
`born_mode`, Boolean `apply_born_posthoc`, and ordered `windows`.
Each window declares `input`, its pinned `sha256`, and `assets_sha256`, a mapping
with pinned `topology`, `fep`, and `restart` file hashes. Manifest paths resolve
relative to the manifest; Q's `[files]` paths resolve relative to each MD input's
directory, which must also be its launch working directory.

`system` means the same topology, boundary parameters, frozen offsets, Q-atom
mapping and state-1 charge vector. Do not use different names to bypass a
physical-system consistency check. Water/protein or different-radius systems
are distinct, but their relationship still needs the campaign-level audit.

### State and correction conventions

The free-energy perturbation (FEP) file accepts only `[FEP]` with `states 2`,
`[atoms]`, and `[change_charges]`. Every Q atom needs explicit charges in both
states. Q-region totals must be 0 and the declared sign. State identities remain
fixed across directions: forward increases state-2 weight from 0 to 1; reverse
decreases it from 1 to 0. Both exact endpoints are required in this **charge-only**
gate. Ghost-atom historical analysis remains a separate, explicitly truncated
Bennett acceptance ratio (BAR) workflow.

No atom-type, bonded, mass or softcore changes are accepted. This verifies the
absence of alchemical Lennard–Jones (LJ, repulsion/dispersion) changes, but native
inspection must still prove that the affected atoms have real, nonzero LJ
parameters in the topology. A zero-LJ dummy present in both states is not ruled
out by the restricted FEP syntax alone.

| `born_mode` | Q's `perstate_born_correction` | `apply_born_posthoc` |
| --- | --- | --- |
| `integrated` | `on` | `false` |
| `posthoc` | `off` | `true` |
| `control` | `off` | `false` |

These are explicit exterior-Born accounting modes, not three different angular
models. The checker does not calculate or apply a correction. The eventual
analysis must consume the audited mode, derive constants from native charge,
radius and Coulomb-unit definitions, and check saved energy records. Merely
declaring `posthoc` does not prove that a later analysis applies it correctly.

### Supported MD input contract

All simulation controls must be explicit. Unsupported sections or keys, repeated
sections/keys, coefficient overrides and ambiguous state mappings are rejected.
The supported sections are `[MD]`, `[cut-offs]`, `[sphere]`, `[solvent]`,
`[intervals]`, `[files]`, and `[lambdas]`; the exact key allowlist is in the checker.
This is intentionally not a general validator for arbitrary Q jobs.

- Use existing spherical MD, direct electrostatic interactions (`lrf off`, which
  disables the local reaction field approximation), solvent/hydrogen bond
  constraints on and solute constraints off. Declare all five cutoffs; native
  geometry/coverage checks are still required to prove no interactions are cut.
- Explicitly enable polarization, charge correction and per-state polarization;
  disable offset adaptation. Use dielectric 80 to match the current hardcoded
  angular dielectric factor, not a separately chosen exterior dielectric.
- Declare solvent radius and positive wall/force parameters, temperature,
  timestep, coupling, initial temperature, random seed, steps and output intervals.
  Inputs in the same physical system must match except seed, paths, lambda and
  the declared exterior-Born on/off control.
- Every starting restart must already exist, contain finite coordinates and
  velocities with consistent dimensions, and have a valid finite offset record.
  Supported binary dialect: little-endian, four-byte record markers/integers,
  double-precision coordinates/velocities and single-precision offsets.
- Require identical offset-record hashes within and across series of one system.
  Reject reused output paths, existing outputs, or paths that overwrite inputs.

An as-yet-unwritten chained restart cannot pass this staged gate. A future launch
wrapper must run it when the required assets exist, or explicitly validate the
planned dependency graph and recheck the realized files. Do not remove restart
hash checks to make a prospective job pass.

## Evidence and remaining work

### Native initialization audit

Q now emits a read-only `Q_BOUNDARY_AUDIT_V2` block during initialization when
per-state polarization is enabled (on the main process only). It records loaded
values before the internal Coulomb charge rescaling; it does not change any
force law, coordinate, velocity, or target. The
[native checker](../../src/QligFEP/boundary_native.py) compares that block against
one window from a **saved preflight report**, after the corresponding run:

```sh
PYTHONPATH="$PWD/src" python -m QligFEP.boundary_native preflight.json \
  --series positive-forward --window 0 --log /path/to/run.log
```

Use the actual series ID; window indices start at zero in the manifest's order.
Save the staged preflight report before launch: re-running the staged gate after
outputs exist intentionally fails its no-overwrite check. The native command
checks the retained manifest and engine hashes, rechecks the current input and
asset hashes, reconstructs the staged window, and rejects disagreement. It
prints JSON without rewriting any files. Success means
`native_initialization_consistency_passed`, still **not production readiness**.

The native block identifies the current angular convention explicitly as code
1, Q-region-only. No proposed total-charge option has been implemented. The
record payloads are:

| Record | Values in order |
| --- | --- |
| `META` | Atom, solute, water, Q-atom, state and shell counts; LJ combining-rule and solvent-type codes |
| `FLAGS` | Per-state polarization, integrated Born, adaptation, alchemical LJ changes, local reaction field, periodic boundary flags |
| `PARAMETERS` | Effective and requested water radii; topology Coulomb constant; dielectric; coefficient override; included/excluded non-Q solute charge; angular/radial force constants; Morse depth/width; maximum loaded atom distance from solvent center |
| `CENTER` | Solvent-center coordinates |
| `SOLUTE_BOUNDARY` | Effective inner restrained radius, exclusion radius, solute shell force constant |
| `CUTOFFS` | Solute–solute, solute–water, water–water, Q-atom, local-reaction-field cutoffs |
| `WATER` | Number density and molecular dipole magnitude used in the target |
| `WATER_COMPATIBILITY` | Uniform water site types/charges; zero hydrogen LJ coefficients (Boolean flags) |
| `WATER_ATOM` | Site index, type, mass, charge, three A and three B topology LJ coefficients for the first water |
| `STATE` | State index, lambda, total/included/excluded Q-region charge, Born energy actually added |
| `QATOM` | Q index, topology index, exclusion flag, type; mass; three A and three B topology LJ coefficients; pure-state charges |
| `SHELL` | Shell index, outer radius, width, offset, pure-state angular strengths |

`CONVENTION` precedes these records; `BEGIN`/`END` delimit exactly one block.
Version 2 requires the water-site records and compatibility checks absent from
version 1; old native logs must not be passed off as this stronger gate.
Units are Q's existing angstrom, elementary-charge and kilocalorie-per-mole
conventions; the A/B coefficients retain the topology combining-rule convention,
not a universal sigma/epsilon interpretation. The checker compares atomic
charges and stored shell radii at their actual single-precision representation.
The global effective radius and LJ coefficients are double precision. It does
not enlarge a tolerance to conceal these representation differences.

The inactive local-reaction-field cutoff is reported as zero by convention;
Q does not read that cutoff key when the method is off, and its inactive variable
need not be initialized. Zero here is not an active zero-range interaction rule.
The initial cutoff bound uses twice the maximum loaded atom radius. A passing
geometric bound is not a proof of pair-list correctness or whole-trajectory
interaction coverage; the initial coordinates also precede initial bond-constraint
projection. The original/native logs and later trajectory coverage remain needed.

The checker rejects nonuniform water site types/charges, inconsistent neutral
water charges, and optimized SPC-like water flags with nonzero hydrogen LJ
coefficients. SPC means simple point charge; the optimized kernel label describes
its interaction assumptions, not proof of a particular named water model.
It also rejects excluded/non-solute Q atoms, nonpositive standard Q-atom LJ
coefficients, wrong loaded state charges/lambdas, mismatched frozen offsets and
inconsistent applied Born constants. It returns the independently recomputed
Born state constants even in post-hoc/control mode, distinguishing those values
from the zero Born energy actually added by the engine. It does not apply them
to a free-energy result. Source-build provenance, log-to-executable attribution,
full solvent-model identification, job completion, final output validation and
physical qualification are not proved by an initialization block.

The automated native smoke test constructs both signs and both ladder directions
on the existing Na/benzene/water fixture: 12 windows at state-2 weights 0, 0.5 and
1, four MD steps per window. It runs the preflight, launches the unchanged Q
workflow, checks finite saved state energies and exact saved lambda mappings,
and verifies frozen restart offsets. Shared starting coordinates are intentional
for this software test: it is **not** a set of independent equilibrated replicas,
a timestep assessment, or a physical endpoint/free-energy validation.

The current smoke test uses a generated copy selecting the general three-site
kernel, preserving every topology charge and LJ coefficient. The archive has an
incompatible optimized-water flag despite a small nonzero hydrogen LJ attraction;
the stronger native gate now rejects that original combination. The
[partition diagnosis](../../docs/analytical-charge-corrections/PARTITION_AUDIT.md)
documents the exact omitted term and the minimal no-softcore force-denominator
initialization fix needed by the spherical three-site charge-only kernel. This
copy is a software control, not an endorsed production water parameterization.
Earlier native passes did not check this water-kernel assumption and must not
be interpreted as passing the new gate. No archived topology or historical
result was changed, and no angular target was changed.

Still required before production:

1. Resolve and document the angular background-charge convention exposed by the
   [partition audit](../../docs/analytical-charge-corrections/PARTITION_AUDIT.md).
   The [bounded opt-in proposal](CHARGE_CONVENTION_DECISION.md) is now documented
   and awaits approval; no target change has been implemented. Do not replace
   the physical target as an incidental preflight change.
2. Collect the now-implemented native initialization evidence for the actual
   campaign and tie it to an exact source build. Complete whole-trajectory
   interaction coverage, solvent-model identification and launch provenance;
   the local fixture's native pass cannot substitute for those checks.
3. A complete both-sign/control/radius/environment campaign, independent replicas
   and initializations, runtime/final-restart audits, and immutable provenance.
4. An analysis package with matched state/leg definitions, Born accounting,
   block/replica uncertainty, equilibration/stationarity and overlap checks, and
   predeclared physical acceptance/failure criteria. No fitting charged results.
5. A measured throughput benchmark and justified compute budget, followed by user
   approval before substantial HPC computation or submission.

The [approved goal](../../docs/analytical-charge-corrections/PROPOSED_GOAL.md) is
unchanged. The input gate and smoke test are components of that goal, not a
replacement completion criterion.

Verification checkpoint: **107 passed, 2 skipped, 2 expected failures in 13.47
seconds** across the focused native and analysis suite. This includes 34 new
input-contract/command-line tests and the native 12-window smoke test. The two
optional analysis-runtime skips and the two known non-angular partition failures
are unchanged. `git diff --check` passed. No HPC job was submitted.

Native-audit checkpoint: **123 passed, 2 skipped, 2 expected failures in 12.44
seconds**, after rebuilding serial Qdyn with GNU Fortran 11. This includes the
12-window native matrix, two additional nonintegrated-Born runs, 12 deliberately
altered native-record checks, and successful/failing native command-line cases.
All new native checks executed. The target-change proposal is still awaiting
approval. Distributed-execution and other compiler builds remain unverified.

Water-kernel checkpoint: **137 passed, 2 skipped, 2 expected failures in 29.08
seconds**, after rebuilding serial Qdyn. The archive's strict partition failures
remain explicit; their LJ and Coulomb causes are now reproduced quantitatively.
The general-water test copy passes the rounding-accounted partition, Cartesian
force/energy, state-bookkeeping and native MD checks. Version-2 compatibility
checks reject the original incompatible water flag and deliberately altered water
records. The angular target is unchanged and its proposed opt-in replacement
still awaits approval. No HPC jobs were submitted.
