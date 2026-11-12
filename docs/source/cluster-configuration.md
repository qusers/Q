# Cluster configuration

QligFEP writes Slurm scripts from Python dictionaries in
`src/QligFEP/settings/settings.py`. Existing profile names are exposed as the
allowed values of `qligfep -c` and `setupFEP -c`.

## Existing profiles

The current source contains `CSB`, `SNELLIUS`, `ALICE`, `HEBBE`, `KEBNE`,
`STALLO`, `UPPMAX`, `TETRA`, `DARDEL`, `LOCAL`, and `LOCALP`. These profiles are
examples of environments used during development—not a guarantee that historical
module names and accounts remain valid.

## Add a profile

Define a dictionary and add it to `CLUSTER_DICT`:

```python
MY_CLUSTER = {
    "NODES": "1",
    "NTASKS": "8",
    "TIME": "0-06:00:00",
    "MODULES": nljoin(
        [
            "module purge",
            "module load OpenMPI/<site-version>-GCC-<site-version>",
        ]
    ),
    "ACCOUNT": "my-project",  # omit if the site does not require it
    **Q_PATHS,
}

CLUSTER_DICT = {
    # existing profiles...
    "MY_CLUSTER": MY_CLUSTER,
}
```

The Q paths resolve to this checkout's `src/q6/bin/q6/`. Compile the binaries on
a filesystem visible from the compute nodes.

## Serial and MPI local profiles

- `LOCAL` uses serial `qdyn` and one task.
- `LOCALP` uses `qdynp` and the detected CPU count.

The non-equilibrium setup rejects both local profiles. Use a Slurm cluster
profile for that workflow.

## Validation procedure

Before a full network:

1. compile Q while the intended compiler/MPI modules are loaded;
2. inspect `mpif90 --showme` or the site equivalent;
3. generate one edge;
4. read the generated `FEP_submit.sh` completely;
5. verify executable paths, modules, account/partition directives, task count,
   memory, and wall time;
6. submit one replicate;
7. confirm that both serial and MPI executables can load their runtime libraries;
8. check all stage logs for normal termination.

:::{warning}
Generated scripts may include a recursive wildcard removal when `-clean` is
used. Do not enable cleanup until output paths and retention requirements have
been verified.
:::

## Why `qdyn` and `qdynp` are both needed

Equilibrium-window production normally uses `qdynp` on clusters.
Non-equilibrium workflows use MPI dynamics for fixed-lambda equilibration and
serial `qdyn` for independent switching trajectories. Therefore a cluster used
for the latter workflow needs both binaries.
