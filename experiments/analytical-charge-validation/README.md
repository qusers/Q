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

The automated native smoke test constructs both signs and both ladder directions
on the existing Na/benzene/water fixture: 12 windows at state-2 weights 0, 0.5 and
1, four MD steps per window. It runs the preflight, launches the unchanged Q
workflow, checks finite saved state energies and exact saved lambda mappings,
and verifies frozen restart offsets. Shared starting coordinates are intentional
for this software test: it is **not** a set of independent equilibrated replicas,
a timestep assessment, or a physical endpoint/free-energy validation.

Still required before production:

1. Resolve and document the angular background-charge convention exposed by the
   [partition audit](../../docs/analytical-charge-corrections/PARTITION_AUDIT.md).
   The [bounded opt-in proposal](CHARGE_CONVENTION_DECISION.md) is now documented
   and awaits approval; no target change has been implemented. Do not replace
   the physical target as an incidental preflight change.
2. Native initialization evidence for effective radius, included/excluded charge,
   topology Coulomb constant, shell geometry, real LJ parameters and interaction
   coverage, tied to the exact source build and input hashes.
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
