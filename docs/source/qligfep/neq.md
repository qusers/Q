# Non-equilibrium protocol

The non-equilibrium workflow drives lambda continuously between two end states.
It runs many forward and reverse switches from endpoint equilibrium snapshots and
uses Bennett acceptance ratio (BAR) to estimate the free-energy difference.

## Setup

```bash
setupFEP -FF AMBER14sb -c MY_CLUSTER -j mapping.json \
  -r 25 -T 298 -R 10 -ts 2fs -rs 42 \
  --neq --neq-reps 5 --neq-steps 50000 \
  --neq-eq-steps 1000 --neq-relax-steps 5000 \
  --neq-steepness 8 --neq-schedule sigmoidal
```

`--windows` and `--sampling` are not used in this mode.

## Controls

| Option | Meaning | Default |
|---|---|---|
| `--neq` | Enable switching workflow | off |
| `--neq-reps` | Forward/reverse switch pairs per replicate | 5 |
| `--neq-steps` | MD steps in each switch | 50,000 |
| `--neq-eq-steps` | Endpoint sampling steps between switches | 1,000 |
| `--neq-relax-steps` | One-time endpoint relaxation | 5,000 |
| `--neq-steepness` | Sigmoid steepness `L` | 8.0 |
| `--neq-schedule` | `sigmoidal` or `linear` | `sigmoidal` |

These are implementation defaults. Switch duration, endpoint decorrelation, and
sample count must be validated for the studied transformations.

## Generated files

Instead of the fixed-window `md_xxxx_xxxx.inp` series, each edge contains:

- ordinary `eq1`–`eq5` preparation/equilibration inputs;
- `relax_0` and `relax_1` endpoint-relaxation templates;
- `eq6_0` and `eq6_1` endpoint-spacing templates;
- `neq_0` and `neq_1` switching templates;
- a run script that repeats both directions.

The exact `[lambda_scaling]` section is appended to switch inputs. Fixed-lambda
equilibration can run with `qdynp`; switches run with serial `qdyn` and are
distributed as independent tasks.

## Output and failures

Forward and reverse switch logs are named `neq_1_*.log` and `neq_0_*.log` in the
current workflow. Only completed switches with parseable work contribute to BAR.
Preserve failed logs; asymmetric failure can bias the apparent distributions.

## Analyze a network

```bash
qligfep_neq_analyze -p 2.protein -w 1.water \
  -T 298 -u kcal -o neq_results.csv \
  -j mapping.json -exp ddg_value -t TARGET
```

Work units are `kcal` by default. Alternatives are `kJ` and `kT`; the selection
changes the BAR beta conversion. It must match the engine output.

The result CSV contains edge and per-leg free energies, protein-minus-water
ddG, uncertainty, overlap, and switch counts. With `-j`, an annotated
`<mapping-stem>_ddG.json` is written. With `-exp`, experimental comparisons and
`<target>_neq_ddG_plot.png` are also produced. Unless `-norun` is passed, Slurm
run diagnostics are collected when available.

## Analyze one edge

Supply both edge paths:

```bash
qligfep_neq_analyze \
  -pe 2.protein/FEP_ligandA_ligandB \
  -we 1.water/FEP_ligandA_ligandB \
  -o one_edge.csv
```

Both `-pe` and `-we` are required to activate single-edge mode. If only one is
supplied, the current implementation falls back to searching the root leg
directories.

## Quality control

Check:

- endpoint trajectories are equilibrated;
- snapshots used for switches are sufficiently decorrelated;
- forward and reverse work distributions both have adequate samples;
- overlap is not dominated by a few outliers;
- results stabilize with more switches and/or longer switching time;
- failed-switch rates are similar and understood in both directions;
- protein and water legs use the same thermodynamic edge and charge method;
- replicate estimates agree before pooling.
