# Isolated native build and one-window restart-chain execution

These tools connect the existing Q molecular dynamics (MD) executable to recorded
source/build evidence and the actual input/restart/output files used by a charge
validation chain. They do not change any force, sampler, thermostat or angular
target. They do not submit jobs to high-performance computing (HPC) systems.
Production physical and statistical qualification remains incomplete.

## Build an immutable native source snapshot

From the clean implementation worktree:

```sh
PYTHONPATH="$PWD/src" python -m QligFEP.charge_build "$PWD" /tmp/q-charge-build \
  --compiler gfortran-11 --commit HEAD
```

Use a new destination. The tool resolves the requested Git commit, archives its
tracked `src/q6` tree and compiles serial Qdyn and Qprep in the fresh extracted
directory. Uncommitted native changes and existing object files are **not** used.
It rejects prebuilt native artifacts in the archive and unsafe extraction paths.
No checkout, reset, cleanup or edit of the current source worktree is performed.

The exact source archive and extracted files, build log, compiler identity,
explicit make command, controlled build environment and executable fingerprints
are retained. Fingerprints use SHA-256, a cryptographic hash algorithm. The
JavaScript Object Notation (JSON) build report is `build.json`; it declares
`gate: isolated_native_build_recorded` and `production_ready: false`.

The build deliberately does not inherit `MAKEFLAGS`, `FC`, `FFLAGS` or library
injection environment variables. It uses the repository's default portable
Fortran options and the explicitly selected compiler, not machine-specific
performance flags or a parallel message-passing build. The compiler and make
executables and source files are fingerprinted and checked for changes during
the build. This is reproducibility evidence, not a signed attestation, a proof
about the compiler's internals, or a complete archive of system libraries.

Each build has a 60-second cap. On timeout the tool terminates the process group
it created, including compiler children, and preserves the incomplete directory
and log. There is no successful build report. An existing directory is never
reused or cleaned automatically. Compiler/make must already be installed; this
command does not install software or access the network.

Build again on the intended HPC operating system and architecture before timing
or running there. A macOS binary is not a Linux production executable. Archive
the complete build directory with the run package, not just its small report.
Freeze and retain the Python implementation as well; native source provenance
does not by itself identify a later analysis or launch-script version.

## Two different preflights

The original staged-input gate requires every starting restart to exist and pins
its contents before launch. That remains useful for independently staged windows.
A chain cannot satisfy that condition for a future predecessor output. The new
chain tool therefore separates:

1. **Planned chain consistency:** check all existing input/topology/charge files,
   state identities, fixed settings and exact restart dependency paths. Future
   restart contents are explicitly unverified; no placeholder hash is invented.
2. **Realized window consistency:** immediately before launching one window,
   inspect its real restart and require its hash to match the preceding window's
   revalidated final output, or the pinned initial restart for the first window.

No nonexistent asset is passed off as a staged-input pass. These gates are not
interchangeable, and neither proves equilibrium sampling.

## Chain plan contract

The plan is a JSON object with these exact top-level fields:

| Field | Contents |
| --- | --- |
| `schema_version` | Integer 1 |
| `engine` | Binary path, `sha256`, full `source_commit`; all must match the isolated build |
| `build_report`, `build_report_sha256` | Path and fingerprint of the retained `build.json` |
| `initial_restart`, `initial_restart_sha256` | Existing initial restart and fingerprint |
| `series` | One series with `id`, `system`, `sign`, `direction`, `replica`, `born_mode`, `apply_born_posthoc`, `windows` |

Series metadata follows the [staged protocol](README.md): charge 0 to declared
sign, fixed pure-state identities, explicit Born accounting and complete monotonic
weight ladder including both endpoints. Each window declares `input`, `sha256`,
and `assets_sha256` containing **only** existing `topology` and `fep` hashes.
Here FEP means free-energy perturbation; its file contains only the allowed
charge-state definitions. Future restart hashes are recorded when they exist.

Plan-level paths resolve relative to the plan file. Q `[files]` paths resolve
relative to the corresponding input's directory, which is the native launch
directory. Each window needs a distinct directory; its energy and final restart
must be distinct paths inside that directory. The first restart must be the
declared initial asset; every later restart must be the immediately preceding
window's final file. Skipping a dependency, changing a frozen offset, overwriting
an input/build artifact, or changing within-chain settings is rejected.

All windows retain the same topology/charge file fingerprints and Hamiltonian
(potential-energy model) settings apart from lambda (the charge-state interpolation
weight) and the random seed. The chain validates a single series; it does
not replace the future campaign-level checks across signs, radii and replicas.
The source and binary hashes in an isolated build report are rechecked, not merely
copied as an unverified commit label beside an unrelated executable.

An explicit schema-version-2 `endpoint_preparation` purpose supports the
[fixed-endpoint preparation schedule](ENDPOINT_PREPARATION.md). It additionally
pins the grid-start origin and total preparation budget, keeps a fixed endpoint
weight and retains restart velocities. Only its segment lengths may differ.
It does not relax the complete monotonic ladder requirement for version 1 and
cannot be passed to free-energy analysis. Both purposes remain unqualified for
production. Driver fingerprints now also include the endpoint and probe modules;
freeze the full implementation before starting a chain and do not mix old and
new driver versions within an existing run.

## Inspect or execute exactly one next window

```sh
PYTHONPATH="$PWD/src" python -m QligFEP.charge_chain inspect /path/to/plan.json
PYTHONPATH="$PWD/src" python -m QligFEP.charge_chain run-next /path/to/plan.json \
  --max-steps 2000 --timeout 60
```

`inspect` is read-only and reports the planned total step count.
`run-next` runs **at most one** Qdyn process, directly and serially. The default
per-invocation limits are 2,000 MD steps and 60 seconds. They are safety limits for
local checks, not a sampling protocol. Any larger physical calculation or HPC
submission still requires user approval and a justified aggregate budget. There
is no scheduler script or automatic full-chain loop in this command.

Before launching, completed predecessors are rechecked against their actual
inputs, native logs, all saved energy records and final restarts. Receipt and
file fingerprints must still agree. The next input's restart must match the
verified predecessor exactly. The current validation/launch module fingerprints
must match those used by earlier completed windows; do not edit code mid-chain.

Within each window the wrapper reserves these new files, using exclusive creation:

- `charge-started.json`: exact binary command, launch directory, engine identity,
  validation/launch source fingerprints, time, host, wrapper process identifier
  (PID), limits and available Slurm job/array identifiers. Slurm is the scheduler;
  recording its environment identifiers is not querying or authenticating a job.
- `charge-preflight.json`: the actual realized input/restart descriptor.
- `charge-native.log`: the invoked program's captured standard output and errors.
- `charge-completed.json`: only after zero native exit status and successful
  [completed-window checks](RESTRAINT_AND_COMPLETION.md), including saved state
  totals/Born accounting and unchanged frozen offsets.
- `charge-failed.json`: an observed launch, timeout or validation failure; the
  other attempt files remain intact for diagnosis.

Existing attempt or output files are never overwritten. An interrupted attempt
without a completion receipt is **not assumed dead or safe to rerun**. Inspect the
actual process or scheduler state before deciding how to recover it. This tool
does not remove a claim file, retry a failed window, salvage partial energies or
continue an orphaned calculation. A failure blocks successor launch. If all
windows already have verified completion receipts, the command only rechecks
them and reports `chain_outputs_consistent`; it starts no new process.

An operating-system crash may leave a partial receipt or output. Such files do
not become success evidence. Hashes and launch records protect against accidental
mixing and make the workflow traceable; they are not tamper-proof attestations.
Normal completion remains distinct from trajectory stability and equilibration.

## Evidence and remaining work

Native integration tests build fresh Qdyn/Qprep from the committed source, prepare
fresh probe water with that Qprep and run three 100-step windows per chain for
both signs and both directions. The initial restart is a 20-step software seed,
not an independently equilibrated endpoint. Tests verify realized restart hashes,
no-op completion rechecks, failure/timeout preservation, budgets, dependency
errors, protected outputs and changed build artifacts/commit labels.

No charged result is fitted and no long trajectory is used for these tests. The
remaining production work includes endpoint preparation/equilibration schedules,
campaign generation and cross-series identity/independence checks, trajectory and
temperature/geometry diagnostics, statistical analysis, target-hardware timing
and an approved HPC allocation. The optional total-charge angular target still
requires separate approval. All reports continue to say `production_ready: false`.

Verification checkpoint: **227 passed, 2 optional integration skips, 2 documented
archival partition expected failures**, in 49.41 seconds. This includes a fresh
isolated native build and the chain/failure tests. Only Python workflow tools,
tests and documentation changed here; no native force or energy expression changed.
