# FEP setup commands

`qligfep` and `setupFEP` share protocol options. `qligfep` prepares one edge and
one system; `setupFEP` iterates a mapping and invokes it for both legs.

## `qligfep`

```text
qligfep -l1 LIGAND1 -l2 LIGAND2 -FF FORCEFIELD
         -s {water,protein,vacuum} -c CLUSTER [OPTIONS]
```

Required options:

| Option | Description |
|---|---|
| `-l1`, `--lig1` | Starting (`from`) ligand name |
| `-l2`, `--lig2` | Ending (`to`) ligand name |
| `-FF`, `--forcefield` | Bundled force-field name or path base |
| `-s`, `--system` | `water`, `protein`, or `vacuum` |
| `-c`, `--cluster` | Profile name from `CLUSTER_DICT` |

The command creates `FEP_<lig1>_<lig2>/inputfiles`, combines ligand/system
parameters and coordinates, runs qprep, creates dual-topology atom mappings and
restraints, writes FEP/equilibration/production inputs, writes run/submission
scripts, and stores provenance in `inputfiles/fep_config.json`.

`-pq`, `--protein_charge` is an advanced integer value used by `setupFEP` to pass
the protein-leg charge into water-leg charge matching.

## `setupFEP`

```text
setupFEP -FF FORCEFIELD -c CLUSTER [-j MAPPING_JSON]
         [-po | -wo] [OPTIONS]
```

Additional options:

| Option | Default | Description |
|---|---|---|
| `-j`, `--json_map` | auto-discover | Network JSON |
| `-po`, `--protein_only` | off | Generate only `2.protein` |
| `-wo`, `--water_only` | off | Generate only `1.water` |

Without `-j`, discovery checks `mapping.json`, then `lomap.json`, then exactly one
other JSON. Protein setup runs before water setup. For an edge explicitly marked
`same_charge: false`, water setup can read the prepared protein charge; therefore
water-only charge-changing setup expects the corresponding protein edge to
already exist.

`-po` and `-wo` are not parser-level mutually exclusive flags. Do not pass both;
the current implementation prioritizes protein-only behavior.

## Shared general options

| Option | Default | Description |
|---|---|---|
| `-r`, `--sphereradius` | `25` | Sphere radius, Å |
| `-b`, `--cysbond` | `auto` | Automatic or explicit `atom_atom,...` disulfides |
| `-l`, `--start` | `0.5` | Start at midpoint or endpoint (`0.5`, `1`) |
| `-T`, `--temperature` | `298` | One temperature or comma-separated values |
| `-R`, `--replicates` | `10` | Independent replicates |
| `-S`, `--sampling` | `sigmoidal` | `linear`, `sigmoidal`, `exponential`, `reverse_exponential` |
| `-w`, `--windows` | `100` | Equilibrium-window count |
| `-ts`, `--timestep` | `2fs` | `1fs` or `2fs` |
| `-clean`, `--files-to-clean` | none | One or more suffix/pattern fragments for post-run removal |
| `-rs`, `--random_state` | random | Seed for reproducible FEP seed generation |
| `-log`, `--log-level` | `info` | Logging threshold |

Although help text mentions an AMBER default for `-FF`, the parser marks it
required; always provide it.

## Restraint and water options

| Option | Default | Description |
|---|---|---|
| `-rest`, `--restraint_method` | `heavyatom_p` | Dual-topology atom-equivalence method |
| `-drf`, `--distance_restraint_force` | `0.5` | Production force constant, kcal mol⁻¹ Å⁻² |
| `-wath`, `--water-thresh` | `1.4` | Water-clash threshold, Å |
| `-wath-ligo`, `--wath-ligand-only` | off | Apply threshold only around ligand atoms |

The restraint method combines a comparison (`heavyatom`, `aromaticity`,
`hybridization`, `element`), a surrounding mode (`p`, `ls`, `strict`), and an
optional Kartograf maximum distance, e.g. `element_strict_1.2`.

## Hamiltonian controls

| Option | Default | Choices |
|---|---|---|
| `-sc`, `--softcore-method` | `standard` | `standard`, `beutler_coul`, `gapsys` |
| `-cm`, `--charge-method` | `ion_match` | `none`, `ion_match`, `coalchemical_water` |

These options alter the scientific protocol, not only performance. Record and
validate them.

## Non-equilibrium controls

| Option | Default | Description |
|---|---|---|
| `-neq`, `--neq` | off | Generate switching rather than windowed inputs |
| `-neqr`, `--neq-reps` | `5` | Forward/reverse pairs per replicate |
| `-neqs`, `--neq-steps` | `50000` | Switching steps |
| `-neqes`, `--neq-eq-steps` | `1000` | Endpoint spacing steps |
| `-neqrs`, `--neq-relax-steps` | `5000` | Initial endpoint relaxation steps |
| `-L`, `--neq-steepness` | `8.0` | Sigmoidal schedule steepness |
| `-neqsched`, `--neq-schedule` | `sigmoidal` | `sigmoidal` or `linear` |

Non-equilibrium mode is rejected for the `LOCAL` and `LOCALP` profiles. The
number of switch pairs must be at least one.

## Examples

```bash
# One protein edge
qligfep -l1 ligandA -l2 ligandB -FF AMBER14sb \
  -s protein -c SNELLIUS -r 25 -rs 42

# Whole equilibrium network, both legs
setupFEP -FF AMBER14sb -c SNELLIUS -j mapping.json \
  -r 25 -R 10 -w 100 -S sigmoidal -rs 42

# Protein leg only
setupFEP -FF AMBER14sb -c SNELLIUS -j mapping.json -po

# Whole non-equilibrium network
setupFEP -FF AMBER14sb -c SNELLIUS -j mapping.json \
  --neq --neq-reps 10 --neq-steps 50000
```
