# Analysis commands

(qligfep-analyze-cli)=
## `qligfep_analyze`

Analyze equilibrium-window water and protein legs.

```text
qligfep_analyze -j MAPPING_JSON -t TARGET
  [-p PROTEIN_DIR] [-w WATER_DIR] [-norun]
  [-exp KEY] [-experr KEY]
  [-m {ddG,ddGf,ddGr,ddGos,ddGbar}]
  [-v] [-lamb N] [-amiss] [-log LEVEL]
```

| Option | Required | Default | Description |
|---|---:|---|---|
| `-j`, `--json-file` | yes | — | Mapping used for setup |
| `-t`, `--target` | yes | — | Target/output name |
| `-p`, `--protein-dir` | no | `2.protein` | Protein-leg root |
| `-w`, `--water-dir` | no | `1.water` | Water-leg root |
| `-norun`, `--no-run-data` | no | off | Skip Slurm output parsing |
| `-exp`, `--experimental-key` | no | none | Experimental edge-value key |
| `-experr`, `--experimental-error-key` | no | none | Reserved parsed option; currently not applied downstream |
| `-m`, `--method` | no | `ddGbar` | Result estimator |
| `-v`, `--save-verbose` | no | off | Write detailed CSV tables |
| `-lamb`, `--n-lambdas` | no | infer | Number of lambda windows to retain |
| `-amiss`, `--allow-missing` | no | off | Continue when mapped edge directories are missing |
| `-log`, `--log-level` | no | `info` | Logging threshold |

The analyzer reads water then protein edge trees, obtains qfep data, subtracts
legs, and writes:

- `<target>_FEP_results.json`;
- `<mapping-stem>_ddG.json` with `Q_ddG_avg`, `Q_ddG_sem`, and `Q_ddG_std`;
- `<target>_ddG_plot.png` when `-exp` is provided;
- `<target>_qEnergies_verbose.csv`, `<target>_dgBar_verbose.csv`, and
  `<target>_run_data.csv` with `-v`.

```bash
qligfep_analyze -j mapping.json -t Tyk2 -m ddGbar \
  -exp ddg_value -v
```

(qligfep-neq-analyze-cli)=
## `qligfep_neq_analyze`

Apply BAR to completed non-equilibrium switching logs.

```text
qligfep_neq_analyze
  [-p PROTEIN_ROOT] [-w WATER_ROOT]
  [-pe PROTEIN_EDGE -we WATER_EDGE]
  [-T KELVIN] [-u {kT,kcal,kJ}]
  [-o CSV] [-j MAPPING_JSON] [-exp KEY]
  [-t TARGET] [-norun] [-log LEVEL]
```

| Option | Default | Description |
|---|---|---|
| `-p`, `--protein-dir` | `2.protein` | Protein-leg root for batch mode |
| `-w`, `--water-dir` | `1.water` | Water-leg root for batch mode |
| `-pe`, `--protein-edge` | none | One protein edge, used with `-we` |
| `-we`, `--water-edge` | none | One water edge, used with `-pe` |
| `-T`, `--temperature` | `298.0` | Temperature for beta/conversion |
| `-u`, `--work-units` | `kcal` | `kT`, `kcal`, or `kJ` |
| `-o`, `--output` | `neq_results.csv` | Main output CSV |
| `-j`, `--json-file` | none | Mapping to annotate |
| `-exp`, `--experimental-key` | none | Experimental edge-value key |
| `-t`, `--target` | `neq` | Plot/output label |
| `-norun`, `--no-run-data` | off | Skip Slurm diagnostics |
| `-log`, `--log-level` | `info` | Logging threshold |

Batch mode analyzes only `FEP_*` directories present in both leg roots. To
activate single-edge mode, supply **both** `-pe` and `-we`.

The main CSV contains calculated edge/leg values, uncertainties, overlap, and
forward/reverse sample counts. With `-j`, the command writes
`<mapping-stem>_ddG.json`; adding `-exp` computes correlation metrics and writes
`<target>_neq_ddG_plot.png`. Unless `-norun` is selected, available Slurm logs
are summarized in `<output-stem>_run_data.csv`.

```bash
# Batch network
qligfep_neq_analyze -p 2.protein -w 1.water \
  -T 298 -u kcal -o neq_results.csv \
  -j mapping.json -exp ddg_value -t Tyk2

# One edge
qligfep_neq_analyze \
  -pe 2.protein/FEP_A_B -we 1.water/FEP_A_B \
  -o A_B_neq.csv
```
