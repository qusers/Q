# Network commands

## `qlomap`

Wrap LoMap to split ligands, score pairs, and write a mapping.

```text
qlomap -i SDF_OR_DIRECTORY [-o OUTPUT_DIRECTORY]
       [-t SECONDS] [-v VERBOSITY] [-exp PROPERTY]
```

| Option | Required | Default | Description |
|---|---:|---|---|
| `-i`, `--input` | yes in practice | — | Multi-molecule SDF or SDF directory |
| `-o`, `--output` | no | input stem/directory | Output directory |
| `-t`, `--time` | no | `30` | Pair/MCS time limit in seconds |
| `-v`, `--verbose` | no | `info` | LoMap verbosity value |
| `-exp`, `--exp-key` | no | none | SDF property containing experimental dG |

For a multi-molecule file, the output directory is created beside the input and
contains one SDF per ligand. The command writes `<output>/lomap.json`. Nodes
contain names, SMILES, formal charge, and SDF properties. Edges contain mapping
score/weight, direction, strictness, and `same_charge`. With `-exp`, nodes receive
`dg_value` and edges receive to-minus-from `ddg_value`.

```bash
qlomap -i ligands.sdf --exp-key r_exp_dg
```

## `qkonnektor`

Generate a configurable perturbation network using Konnektor and Kartograf.

```text
qkonnektor -i SDF_OR_DIRECTORY [-o DIR]
  [-n {mst,rmst,star,nedges,cyclic}]
  [-s {restraint,kartograf,combined}] [-rest METHOD]
  [-p PROCESSES] [-cl LIGAND]
  [--n_redundancy N] [--connectivity N]
  [--separate_charges] [--charge_changes_score SCORE]
  [-exp PROPERTY] [--self-solve] [-log LEVEL]
```

| Option | Required | Default | Description |
|---|---:|---|---|
| `-i`, `--input` | yes | — | Multi-molecule SDF or SDF directory |
| `-o`, `--output` | no | input stem/directory | Output directory |
| `-n`, `--network` | no | `mst` | `mst`, `rmst`, `star`, `nedges`, or `cyclic` |
| `-s`, `--scorer` | no | `combined` | `restraint`, `kartograf`, or product-like combined score |
| `-rest`, `--restraint_method` | no | `heavyatom_p` | Mapping method for restraint scoring |
| `-p`, `--processes` | no | `1` | Parallel scoring processes |
| `-cl`, `--central_ligand` | no | none | Central node for a star network |
| `--n_redundancy` | no | `2` | Redundant MST iterations for `rmst` |
| `--connectivity` | no | `3` | Minimum degree target for `nedges` |
| `--separate_charges` | no | off | Build separate formal-charge networks |
| `--charge_changes_score` | no | `0.0` | Combined-score multiplier for charge changes |
| `-exp`, `--exp-key` | no | none | Experimental dG SDF property |
| `--self-solve` | no | off | Re-embed misaligned outliers and rewrite SDF coordinates |
| `-log`, `--log-level` | no | `info` | Logging threshold |

At least two molecules are required. Output is `<output>/mapping.json` with nodes
and setup-compatible edges.

`--charge_changes_score 0.0` blocks charge-changing edges in the combined scorer;
`0.5` halves their score and `1.0` removes the penalty. The parser does not
currently constrain the numeric range, so use values from 0 through 1.

`--self-solve` can rewrite ligand coordinates in place. Back up the input and
review the re-embedded structures.

```bash
qkonnektor -i ligands.sdf -n cyclic -s combined -p 8
qkonnektor -i ligands.sdf -n star -cl reference_ligand
```

## `qmapfep`

Generate an in-house mapping or turn an existing mapping into browser assets.

```text
qmapfep -i INPUT_SDF [-l INPUT_JSON]
        [-m {MFP,Tanimoto,MCS,SMILES}]
        [-o OUTPUT_TEXT] [-wd DIRECTORY]
```

| Option | Required | Default | Description |
|---|---:|---|---|
| `-i`, `--insdf` | yes in practice | — | Multi-molecule SDF |
| `-l`, `--load_from_json` | no | none | Existing mapping to adapt/render |
| `-m`, `--metric` | no | `MFP` | `MFP`, `Tanimoto`, `MCS`, or `SMILES` |
| `-o`, `--otxt` | no | input string before first `.`, lowercased | Output JSON name/stem |
| `-wd`, `--workdir` | no | input string before first `.` | Working directory |

The command writes `<workdir>/<output>.json`, ligand drawings under
`<workdir>/img/`, and `<workdir>/qmapfep.html`. Charge groups are represented as
separate graph components. `-l` still requires the SDF for molecule rendering.
The current default is derived with a simple split of the full input string at
the first period, not with path-aware stem handling. For path-qualified or
dotted filenames, pass both `-wd` and `-o` explicitly.

Open `qmapfep.html` locally, import the generated JSON, inspect/edit the graph,
and export the revised JSON for `setupFEP`.

```bash
qmapfep -i ligands.sdf -l ligands/lomap.json -wd .
```
