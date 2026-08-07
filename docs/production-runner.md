# Platform-independent production FEP runner

Windowed FEP setup supports a lean production mode designed for campaigns that
are limited by filesystem metadata or inode quotas.

```bash
setupFEP -FF OPLS2015 -c SNELLIUS -S sigmoidal -r 25 -w 100 --production
```

`qligfep ... --production` provides the equivalent mode for one edge/leg.

## Platform contract

The production runner works directly in each replicate directory. It does not
inspect or use `$TMPDIR`, assume node-local storage, or perform a stage-out.
Scheduler and MPI details are kept in a thin generated launch script; the
scientific workflow is implemented by the platform-independent
`qligfep-run` Python command.

The execution host must provide:

- the installed QligFEP package (including `qligfep-run`);
- the configured `qdyn`/`qdynp` and `qfep` executables;
- the optional launcher named by the cluster adapter, such as `mpirun`.

## File lifecycle

For each replicate, the runner:

1. Discovers the restart dependency graph from the generated `eq*.inp` and
   `md_*.inp` definitions.
2. Renders a single `current.inp`; topology and FEP files are referenced at
   their shared locations rather than copied.
3. Disables trajectories for both equilibration and production stages.
4. Appends all qdyn output to `run.log`.
5. Alternates two checkpoint slots per active branch. A newly written
   checkpoint is recorded atomically before its consumed parent is retired.
6. Runs `qfep` against the native per-window `.en` files.
7. Converts every `.en` record into one human-readable `energies.csv`, including
   the original filename, frame, state, lambda, all Q energy components, and any
   off-diagonal data.
8. Removes native `.en` files only after qfep succeeds, CSV validation succeeds,
   and the conversion metadata has been persisted.

A completed sequential run normally retains:

```text
energies.csv
md_0000_1000.re
qfep.out
run.log
run-state.json
```

Centre-start protocols retain both terminal endpoint restart files.

## Resume and failure behaviour

`run-state.json` is atomically replaced after each successful transition. Running
the same generated launch script again resumes from the first incomplete stage.
The runner rejects a state file created with different temperature, seed,
launcher, qdyn executable, trajectory policy, input directory, or stage list.

A qdyn process is successful only if it both returns zero and prints
`terminated normally`. A failed qdyn stage does not promote its checkpoint.
Failures in qfep or CSV conversion preserve all native `.en` files.

## Direct runner usage

Generated scripts normally invoke the runner, but it can also be called directly:

```bash
qligfep-run \
    --input-dir /path/to/FEP_edge/inputfiles \
    --work-dir /path/to/FEP_edge/FEP1/298/1 \
    --temperature 298 \
    --seed 12345 \
    --fep-file FEP1.fep \
    --qdyn /path/to/qdynp \
    --qfep /path/to/qfep \
    --launcher "mpirun -n 16"
```

Useful diagnostic options:

- `--keep-en`: retain the native binary energies after conversion.
- `--trajectories`: retain DCD output instead of the production default.
- `--energy-csv NAME`: change the consolidated CSV filename.

The standalone `qen2csv` command can convert existing energy files without
running dynamics or deleting its inputs:

```bash
qen2csv md_*.en -o energies.csv --summary-json energies-summary.json
```
