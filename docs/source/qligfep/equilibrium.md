# Equilibrium-window protocol

The equilibrium workflow samples a set of fixed lambda windows for every edge,
leg, temperature, and replicate. Qfep then combines state-energy records from the
windows.

## Setup

```bash
setupFEP -FF AMBER14sb -c MY_CLUSTER -j mapping.json \
  -r 25 -T 298 -R 10 -w 100 -S sigmoidal \
  -l 0.5 -ts 2fs -rs 42
```

The shown values are current CLI defaults except `-rs`, which is unset by
default, and the site-specific cluster name. Supplying them explicitly improves
readability of archived commands.

## Lambda controls

| Option | Meaning | Default |
|---|---|---|
| `-l`, `--start` | Start from midpoint (`0.5`) or endpoint (`1`) | `0.5` |
| `-w`, `--windows` | Total windows | `100` |
| `-S`, `--sampling` | Lambda spacing | `sigmoidal` |
| `-R`, `--replicates` | Independent repeats | `10` |
| `-T`, `--temperature` | Temperature or comma-separated temperatures | `298` |

Available spacing values are `linear`, `sigmoidal`, `exponential`, and
`reverse_exponential`. The correct number/spacing depends on transformation and
overlap; defaults are a protocol starting point, not universal convergence proof.

## Dual topology and restraints

Both ligand topologies remain present. State-dependent nonbonded interactions
turn each ligand on/off while mapped atoms are restrained together. The default
restraint method is `heavyatom_p`.

A method string combines:

1. comparison: `heavyatom`, `aromaticity`, `hybridization`, or `element`;
2. surrounding-atom mode: `p`, `ls`, or `strict`;
3. optional Kartograf maximum distance, for example `_1.2`.

Example:

```bash
setupFEP ... -rest element_strict_1.2 -drf 0.5
```

`-drf` is the production distance-restraint force constant in kcal mol⁻¹ Å⁻².
Inspect generated `[distance_restraints]`; automatic mapping is not infallible.

## Softcore

```bash
setupFEP ... -sc standard
```

`standard` is the default. With it, the `[softcore]` section that QligFEP
writes has no effect on the calculation; see
[Soft-core and the endpoint problem](../concepts/softcore.md).

## Charge-changing edges

`-cm` selects:

- `ion_match` (default): add counter-ions to the water leg so its total charge
  matches the protein leg;
- `coalchemical_water`: transform a water into a co-alchemical sodium/chloride
  representation;
- `none`: disable this handling and retain the raw finite-sphere charge artifact.

Charge-changing transformations require explicit validation of both generated
legs, system charges, and convergence. Do not combine methods within a network
without a defined thermodynamic rationale.

## Submission

Each edge directory contains `FEP_submit.sh`, normally launching a Slurm array of
replicates. `setupFEP` does not submit it. Inspect the script, then submit a small
pilot according to site policy.

The `-clean` option inserts wildcard deletion at the end of generated jobs. For
example, `-clean dcd` removes matching trajectory files. Delay cleanup until
normal termination and data-retention checks are complete.

## Analysis

```bash
qligfep_analyze -p 2.protein -w 1.water \
  -j mapping.json -t TARGET -m ddGbar
```

Methods are `ddG`, `ddGf`, `ddGr`, `ddGos`, and `ddGbar`. `ddGbar` is the default.
Use `-v` for verbose q-energy, BAR, and run-data CSVs; `-norun` skips Slurm log
diagnostics; `-amiss` allows absent mapped edges.

Primary outputs:

- `<target>_FEP_results.json` — detailed per-edge/leg data;
- `<mapping-stem>_ddG.json` — mapping copy with `Q_ddG_avg`, `Q_ddG_sem`, and
  `Q_ddG_std` on edges;
- `<target>_ddG_plot.png` when an experimental key is supplied;
- verbose CSVs when `-v` is supplied.

## Convergence review

Review window completeness, adjacent-window overlap, forward/reverse estimates,
replicate spread, leg cancellation, outlier edges, and cycle closure where the
network contains cycles. `-amiss` is for intentionally partial analysis; report
which edges were omitted.
