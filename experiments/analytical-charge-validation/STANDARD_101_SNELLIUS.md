# Standard Snellius 101-window target subset

User-authorized follow-up to [job 26482080](NO_SOFTCORE_RESULTS_26482080.md).
This uses the [TYK2 tutorial submission layout](../../tutorials/Tyk2/README.md),
not one encompassing Python worker-pool job.

## Selection and submission

Keep the same outcome-independent edge selection as the previous checks:

- c-Met: `CHEMBL3402742_23 -> CHEMBL3402744_300` and its reverse setup.
- Eg5: `CHEMBL1085666 -> CHEMBL1089056` and its reverse setup.

Each target/direction has `1.water/FEP_<lig1>_<lig2>` and
`2.protein/FEP_<lig1>_<lig2>`. Each edge directory contains the standard
`FEP_submit.sh`, which submits `inputfiles/runSNELLIUS.sh` as a three-replica
Slurm scheduler array. Thus eight arrays produce 24 replica jobs. Each replica
uses 16 Message Passing Interface (MPI) ranks on `rome`, account `ugsei19097`,
with a 24-hour wall-time cap and no automatic requeue. Runtime predictions must
be based on actual longer-window timing, not the startup-dominated MPI check.

The staging command explicitly selects the cluster:

```sh
PYTHONPATH="$PWD/src" python experiments/analytical-charge-validation/standard_101.py stage \
  experiments/analytical-charge-validation/runtime/standard-101-references \
  experiments/analytical-charge-validation/runtime/standard-101-stage \
  --remote-root /projects/prjs2157/astra-charge-change-perturbation/releases/standard-101-20260909-v2 \
  -c SNELLIUS
```

This calls the same `QligFEP.write_runfile` and `write_submitfile` renderers used
by `qligfep --cluster SNELLIUS`. It reuses already prepared reference topologies,
parameter tables and restraint mappings instead of repeating ligand preparation
or changing the atom mapping. It is not a claim that `setupFEP` was rerun.

Submission is from **within each edge directory**, as in the tutorial:

```sh
cd /projects/prjs2157/astra-charge-change-perturbation/releases/standard-101-20260909-v2/runs/cmet-fwd/2.protein/FEP_CHEMBL3402742_23_CHEMBL3402744_300
bash FEP_submit.sh
```

Do not repeat submission of an existing edge. Each replica refuses to overwrite
an existing run directory. Protein arrays are submitted before water arrays.

## Sampling and analysis distinction

Preserve the original 101-point sigmoidal ladder, midpoint initialization and
two restart branches. Preserve 10 picoseconds per molecular dynamics (MD)
window and 131 picoseconds of preparation, including the original minimization
and heating schedule. Convert the original 2-femtosecond phases to 1-femtosecond
steps while doubling their step counts; keep the initial 0.2-femtosecond phase.
Each replica therefore runs 1,141 picoseconds: **27.384 nanoseconds aggregate**
over 24 replicas. Seeds are 2924, 25360 and 21448, the reference protocol's set.

Sample all 101 windows, retaining endpoint data. For qfep, list only the
**99 interior energy files**: exclude `md_1000_0000.en` and `md_0000_1000.en`.
For this archived schedule the retained state-2 interval is **0.001–0.999**,
not the 0.0001–0.9999 weights used in the separate short checks. Preserve the
original qfep histogram/settings header except its energy-file count. The
100-frame discard corresponds to 1 picosecond at the new energy-write interval.

Results are truncated interior estimates. Do not assume missing endpoint
contributions cancel. Cross-lambda reconstruction and sampling/overlap checks
remain necessary before physical interpretation of Bennett acceptance ratio
(BAR) free energies; automatic qfep output is not itself qualification.

## Hamiltonian and safety settings

- No softcore: remove its coefficient section and disable max-potential mode;
  retain all other free-energy perturbation (FEP) tables and topology bytes.
- Explicit repaired SHAKE/SHAKE constraints with original hydrogen selections.
- Per-state polarization enabled, adaptation disabled throughout, fresh zero
  angular offsets; retain the original radial and angular force constants.
- Integrated Born enabled with dielectric 80 and actual effective solvent
  radii. Do not apply Born a second time in postprocessing. Duplicated
  integrated/post-hoc trajectories are no longer needed for every window.
- Direct electrostatics, local reaction field disabled, pair-list updates every
  step. Raise all pair cutoffs to 120 angstrom to provide margin beyond the
  earlier Eg5 conservative bound of approximately 98 angstrom. The earlier
  99-angstrom checks already covered all pairs, so this does not intentionally
  introduce a new physical interaction model.
- Keep standard output directories and sequential restart dependencies.
  Disable trajectory-file output, but **do not clean energy files, inputs,
  logs, topologies or restarts**.

After every native invocation, validate normal completion, no hot-atom resets,
solver/boundary settings, geometry and pair-coverage diagnostics, and unchanged
zero offsets. Interior energy records must be finite and sum correctly. Exact
endpoint energy records are retained but intentionally excluded from this
energy-analysis gate. Every run is checked before its restart is used by the
next stage. Failed checks stop that replica; no automatic retries occur.

## MPI build qualification

The engine is rebuilt from native commit
`787e76ebf40c37b3a68dd63d1132df4ce78bba9f`, using GNU Fortran 12.3.0 and
OpenMPI 4.1.5. Job **26487594** passed the 16-rank checks for both protein
targets in 69 seconds, including compilation. The largest serial/MPI saved
energy difference over the short checks was **3.780087354243733e-12**
kilocalories per mole. The largest integrated/post-hoc Born residual was
**5.684341886080802e-14**; paired final restart bytes matched.

An earlier attempt, job 26487538 in `standard-101-20260909`, stopped because the
checking script copied a report file and then refused to overwrite it. This
was a validation-wrapper file collision, not a Qdyn failure. Its files are
preserved; the corrected attempt used the new `-v2` release directory.

Each production replica verifies pinned executable, driver and input hashes
before starting. Native build and MPI qualification receipts remain in the
new project-space release; nothing is installed into the home environment.

## Locations

Remote release and driver/build provenance:

`/projects/prjs2157/astra-charge-change-perturbation/releases/standard-101-20260909-v2/`

Remote standard-layout runs:

`/projects/prjs2157/astra-charge-change-perturbation/releases/standard-101-20260909-v2/runs/`

Local staged inputs (not live remote results):

`/Users/davidararipe/projects/Q-rebase/.worktrees/analytical-corrections-modernize/experiments/analytical-charge-validation/runtime/standard-101-stage/`

`protocol.json` in each edge records duration, execution order and input hashes;
`analysis-scope.json` records the exact endpoint exclusions and Born convention.
