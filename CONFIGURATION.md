# QligFEP configuration: profiles and protocols

> Temporary home for this guide until the documentation branch is merged.

QligFEP reads its settings from YAML. Two concepts are kept apart:

- A **profile** says *where and how* jobs run: Slurm array or local, which engine,
  which modules, how many cores. It must never change the simulated physics.
- A **protocol** says *what* is simulated: method, lambda windows, stage lengths,
  any qdyn input keyword.

The same protocol can therefore run on a laptop and on a cluster, and one machine
can run many protocols.

## Quick start

```bash
qligfep_config init      # write a commented profile for this machine
qligfep_config check     # test it: binaries, mpirun, libraries, sbatch
qligfep_config show      # list profiles and protocols and where they come from

setupFEP -FF AMBER14sb -j mapping.json          # uses default_profile
cd 1.water/FEP_lig1_lig2 && ./FEP_submit.sh     # sbatch array, or local run
```

Nothing has to be configured to keep working as before: the historical clusters
(`SNELLIUS`, `TETRA`, `HABROK`, `LOCAL`, ...) ship as packaged profiles, so
`-c SNELLIUS` behaves as it did. `-c` is now optional once `default_profile` is set.

## Where settings come from

Lowest to highest priority:

| Layer | Location | May contain |
|---|---|---|
| Packaged defaults | `src/QligFEP/config/defaults.yaml` | profiles, the `default` protocol |
| User file | `~/.config/qligfep/config.yaml` (or `$QLIGFEP_CONFIG`) | profiles |
| Project file | `./qligfep.yaml` (or `--config PATH`) | profiles and protocols |
| Command line | flags of `qligfep` / `setupFEP` | everything they expose |

A name defined in several layers is merged key by key, so a higher layer states
only what it changes. To set your account on a packaged profile:

```yaml
# ~/.config/qligfep/config.yaml
default_profile: TETRA
profiles:
  TETRA:
    slurm: {account: my-project}
```

Protocols are deliberately not read from the user file. They define the science
of a project, so they live in the project's `qligfep.yaml`, next to the results.

Unknown keys are errors, with the file, the key and the closest valid name, so a
typo cannot silently fall back to a default.

## Profiles

```yaml
profiles:
  habrok:
    execution: slurm-array        # slurm-array | local
    engine: qdynp                 # qdynp (MPI) | qdyn (serial)
    modules: [OpenMPI/4.1.4-GCC-11.3.0]   # or: from-build
    setup: module purge           # shell lines run before the modules are loaded
    resources:
      ntasks: 8                   # cores per replicate (auto = all cores, local only)
      time: "0-08:00:00"
      mem_per_cpu: auto           # MB; auto = 512 (equilibrium) / 1700 (NEQ); null omits it
      max_parallel_replicates: 1  # local execution only
      neq_cores: auto             # see "Non-equilibrium cores"
    slurm:
      account: null
      partition: null
      export: ALL                 # job inherits the submit environment; null leaves it to the site
      extra_directives: ["--qos=short"]
    systems:                      # per-leg overrides
      water: {resources: {ntasks: 4}}

  habrok-softcore:                # same site, another engine build
    extends: habrok
    q_bin: /home/me/clones/Q-softcore/src/q6/bin/q6
    modules: from-build
```

### Execution types

The generated run script is the same for both types. It takes its task index and
core count from `QLIGFEP_TASK_ID` / `QLIGFEP_NTASKS` and falls back to the Slurm
array variables. Only `FEP_submit.sh` differs:

- `slurm-array`: `./FEP_submit.sh` calls `sbatch`; `./FEP_submit.sh 1 3` submits
  only tasks 1 and 3.
- `local`: `./FEP_submit.sh` runs the tasks on this machine,
  `max_parallel_replicates` at a time, each writing `slurm.run<task>.local.0.out`.
  Those logs carry the same footer as Slurm logs, so `qligfep_analyze` and
  `qligfep_neq_analyze` read local and cluster runs alike.

Local execution supports non-equilibrium FEP when the profile uses `engine: qdynp`
and OpenMPI's `mpirun` is available.

### Job environment

Jobs are submitted with `#SBATCH --export=ALL`, so they inherit the environment
you submit from. Activating your (micro)mamba environment before running
`./FEP_submit.sh` is normally all that is needed. The forwarded environment
contains variables such as `PATH`, not shell functions, so `micromamba activate`
is not available *inside* the job; nothing in the generated scripts needs it.
For a site that ignores the export, put the activation lines in `setup`.

Every job first checks that the Q binaries exist and are executable on the node it
landed on, and stops with a clear message otherwise.

### Which engine is used

In order: `--q-bin`, `$QLIGFEP_Q_BIN`, the profile's `q_bin`, and finally the
binaries compiled inside the QligFEP installation that runs the command. The last
one is the default and keeps several clones independent: each checkout runs its
own engine without any configuration.

`modules: from-build` loads the modules recorded next to the binaries. Record
them once, in the environment the binaries were compiled in:

```bash
qligfep_config build-info          # writes <q_bin>/build-info.json
```

With Lmod only the modules you loaded yourself are recorded, not the dependencies
they pulled in.

Clones built with different toolchains then get the right modules from one profile.

### Non-equilibrium cores

A non-equilibrium task fires `2 x switch_pairs` switches concurrently.

- `neq_cores: auto` under Slurm requests enough cores for all of them in one wave,
  and never fewer than `ntasks`. This is the previous behaviour. For local execution
  `auto` is the profile's `ntasks`.
- `neq_cores: 16` is used as is. With fewer cores than switches the script runs
  them in successive waves; choose a divisor of the switch count to avoid idle cores.

The choice is logged when the inputs are generated.

## Protocols

```yaml
# ./qligfep.yaml
default_protocol: fast-ligand-simulation

protocols:
  fast-ligand-simulation:
    description: Shorter sampling for a first look.
    # extends: default            # implied for every protocol
    method: equilibrium           # equilibrium | nonequilibrium
    timestep: 2fs
    temperature: 298              # or a list: [298, 310]
    replicates: 5
    sphere_radius: 25
    equilibrium:
      windows: 50
      sampling: sigmoidal
      start: 0.5
    stages:
      eq4: {time: 10ps}
      eq5:
        time: 100ps
        inp:
          MD: {bath_coupling: 5}
          intervals: {trajectory: 500}
      production:
        time: 10ps                # per lambda window
        inp:
          cut-offs: {solute_solvent: 12}

  neq-screen:
    method: nonequilibrium
    nonequilibrium:
      switch_pairs: 5
      schedule: sigmoidal
      steepness: 8.0
    stages:
      relax: {time: 10ps}         # one-time endpoint relaxation
      spacing: {time: 2ps}        # equilibration between successive switches
      switch: {time: 100ps}       # the lambda switch
```

Select one with `--protocol NAME`, or set `default_protocol`.

### Methods

`method` chooses between windowed equilibrium FEP (analysed with qfep) and
non-equilibrium fast switching (analysed with BAR over the work distributions).
Each method has its own options block and its own stages:

| Method | Options block | Stages |
|---|---|---|
| `equilibrium` | `equilibrium:` | `eq1`-`eq5`, `production` |
| `nonequilibrium` | `nonequilibrium:` | `eq1`-`eq5`, `relax`, `spacing`, `switch` |

Setting options or stages that belong to the other method is an error rather than
being ignored.

### Stage overrides

Each stage accepts:

- `time`: a duration with a unit (`fs`, `ps`, `ns`), converted with that stage's
  stepsize. It stays correct if the timestep changes.
- `steps`: the length in MD steps instead.
- `inp`: any qdyn input keyword, as `section: {key: value}`.

Keywords QligFEP already writes are replaced in place. Other keywords are added
to their section. A keyword this version does not recognise is passed to qdyn with
a warning at setup time. qdyn itself ignores keywords it does not know, so that
warning is the only notice of a typo.

The sections QligFEP generates from the system (`files`, `lambdas`, the restraint
sections, `lambda_scaling`) and the per-replicate keys (`random_seed`,
`initial_temperature`) cannot be overridden.

### Command-line flags

Flags given on the command line win over the protocol. `-w 40` with a protocol
that sets 50 windows runs 40. Flags that are not given take the protocol's value.

## What was used for a run

Each `FEP_*/inputfiles/` directory contains:

- `resolved_config.yaml`: the fully resolved profile and protocol, including the
  engine directory and where it came from;
- `fep_config.json`: the equivalent command line, the protocol name and the git
  commit of the QligFEP checkout.

## Reference: `qligfep_config`

| Command | Purpose |
|---|---|
| `init` | Detect scheduler, cores, modules and Slurm accounts; write a commented user file. An existing file is never modified; the proposal is printed instead. |
| `show` | List profiles and protocols. `--profile NAME` / `--protocol NAME` print one fully resolved. |
| `check` | Validate all files, then test a profile here: binaries, `sbatch`, `mpirun` and the engine's shared libraries after the profile's setup lines and modules. |
| `build-info` | Record loaded modules and compiler/MPI versions next to the binaries. |

## Current limits

- `qresfep` and `setup_resFEP` accept profile names (including your own
  Slurm profiles) but do not use protocols, per-system overrides, `partition`,
  `export` or `extra_directives` yet. QresFEP still has no local execution.
- The switching phase of non-equilibrium runs uses OpenMPI's `mpirun` directly;
  `mpi_command` only changes the launch of the equilibration stages.
- `qligfep_config check` runs on the node you call it from. Compute nodes can
  differ, so a first short job is the real test.
