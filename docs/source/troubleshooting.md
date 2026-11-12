# Troubleshooting

## Start with provenance

For every problem record:

```bash
python --version
gfortran --version
mpif90 --version  # on MPI systems
git rev-parse HEAD
```

Also retain `fep_config.json`, the complete command, input directory listing,
qprep/Qdyn/qfep logs, and scheduler output.

## Python and command errors

### Command not found

Activate the installation environment and reinstall editable entry points:

```bash
micromamba activate qligfep
python -m pip install -e .
command -v setupFEP
```

### Import failure in ligand/network commands

The base `pyproject.toml` dependency list does not cover every OpenFF and network
package. Install the full environment from [Installation](installation.md).

### `src/q6/bin/q6` is missing

Install `gfortran`, then:

```bash
make -C src/q6 qprep-only COMP=gcc
```

Build `all` and `mpi` before simulation.

## Preparation errors

### Unknown residue or atom

- match PDB names to the selected force-field library;
- remove/resolve alternate locations;
- verify protonation and terminal residue naming;
- separate disconnected molecules correctly;
- load ligand/cofactor libraries before reading coordinates.

### Sphere center is wrong

Confirm all SDF conformers share a coordinate frame. Remove misaligned outliers
or calculate the center from the bound reference ligand. Visualize the final
sphere.

### Qprep reports high bonded energies

Inspect missing/duplicated atoms, cross-links, atom order, ligand alignment,
coincident coordinates, and parameter assignments. A topology file being written
does not make the structure valid.

### Unexpected total charge

Sum ligand library charges, review protonation and boundary neutralization, and
compare the charge reported by qprep with the intended model. For charge-changing
edges, inspect both water and protein setup logs.

## Simulation errors

### SHAKE failure

Find the first failed stage and first constrained bond in the log. Common causes
are steric clashes, bad geometry, inappropriate timestep, incorrect bond
parameters, nearly decoupled overlapping atoms, or incompatible restarts.
Shorten the timestep only after fixing incorrect input.

### Temperature or energy explosion

Check initial bonded/nonbonded components, units, restart/topology consistency,
constraint settings, solvent overlaps, and atom mappings. Run `qdum`, then a
very short serial `qdyn` case with frequent output.

### MPI job fails immediately

Ensure the run-time compiler/MPI modules match the build environment, compute
nodes can see all binary/input paths, and the site's preferred launcher and task
layout are used.

### Slurm job says completed but output is incomplete

Check Q's own normal termination and expected restart/energy records. Shell
scripts can terminate without propagating an internal program failure if not
configured to stop on errors.

## Equilibrium analysis errors

### Missing edge

Check that the edge name and direction match the mapping and both roots contain
`FEP_<from>_<to>`. Use `-amiss` only for intentionally partial analysis.

### qfep output is unstable

Check ordered windows, missing/truncated `.en` files, equilibration discard,
neighboring-window overlap, endpoint energies, forward/reverse disagreement,
and replicate outliers. Add sampling based on diagnostics, not uniformly by
habit.

## Non-equilibrium analysis errors

### No shared edges found

The analyzer considers only identically named `FEP_*` directories under both leg
roots. Verify `-p` and `-w` paths and names.

### Single-edge flags appear ignored

Supply both `-pe` and `-we`; one by itself falls back to batch mode.

### Poor BAR overlap

Increase switching duration and/or sample count, improve endpoint equilibrium,
inspect failed/outlier switches, and verify work units and direction. Very low
overlap means uncertainty from a naive estimate is unreliable.

## Reporting a bug

Create the smallest reproducible case, remove private structures if necessary,
and open an issue at <https://github.com/qusers/Q/issues>. Include expected and
observed behavior, exact commands, versions, traceback/log, and whether the
problem reproduces in a clean environment.
