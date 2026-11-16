# Residue FEP commands

`qresfep` prepares one leg of one mutation; `setup_resFEP` prepares both legs
for a list of mutations; `qresfep_analyze` combines the results. See the
[QresFEP workflow](../qresfep/workflow.md).

(qresfep-cli)=
## `qresfep`

```text
qresfep -m MUTATION -mc CHAIN -S {protein,tripeptide} -c CLUSTER [OPTIONS]
```

Run it in a directory that holds the output of `qprep_prot` for a sphere
centred on the mutated residue. It creates `FEP_<MUTATION>/`.

Required options:

| Option | Description |
|---|---|
| `-m`, `--mutation` | `<wild-type><position><mutant>`, such as `LEU39ALA` or `L39A`; the position is the residue number in the PDB given to `qprep_prot` |
| `-mc`, `--mutchain` | Chain of the mutated residue |
| `-S`, `--system` | Leg to set up: `protein` or `tripeptide` |
| `-c`, `--cluster` | Profile name from `CLUSTER_DICT` |

(resfep-protocol-options)=
Protocol options, shared with `setup_resFEP`:

| Option | Default | Description |
|---|---|---|
| `-FF`, `--forcefield` | `OPLSAAM` | The only supported force field |
| `-t`, `--tripeptide` | `A` | Residues next to the mutated one in the reference peptide: `A` (Ala), `G` (Gly), `X` (none), `Z` (native neighbours) |
| `-w`, `--windows` | `25` | λ windows per stage; a mutation has two stages |
| `-s`, `--sampling` | `exponential` | `linear`, `sigmoidal`, `exponential`, or `reverse_exponential` |
| `-l`, `--start` | `1` | λ endpoint the simulations start from: `1`, `0.5`, or `0` |
| `-ts`, `--timestep` | `2fs` | `1fs` or `2fs` |
| `-T`, `--temperature` | `298` | Temperature in K; several as `T1,T2,...` |
| `-r`, `--replicates` | `10` | Independent replicates (`-R` in `setup_resFEP`) |
| `-eqs`, `--eq5-steps` | 2.5 ns at 2 fs | Length of the final equilibration, in steps |
| `-ps`, `--production-steps` | `10000` | Production steps per λ window |
| `-cof`, `--cofactors` | none | Cofactor base names; each needs a `.pdb`, `.lib` and `.prm` |
| `-clean`, `--files-to-clean` | keep all | Suffixes of files to delete when a run finishes, such as `dcd` |
| `--no-trajectories` | off | Write no DCD trajectories |
| `--separate-scaling` | `on` | Value of `separate_scaling` in the `qdyn` input |
| `--manuscript-settings` | off | Apply the fixed settings of the charge-changing protocol; see the [workflow](../qresfep/workflow.md#charge-changing-mutations) |
| `-rs`, `--random_state` | none | Seed from which the replicate seeds are drawn |
| `--seeds` | none | One explicit seed per replicate; excludes `-rs` |
| `-log`, `--log-level` | `info` | Logging threshold |

```bash
qresfep -m LEU39ALA -mc A -S protein -FF OPLSAAM -c SNELLIUS \
  -w 25 -s exponential -l 1 -ts 2fs -T 298 -r 10
```

(setup-resfep-cli)=
## `setup_resFEP`

```text
setup_resFEP -i PROTEIN_PDB (-m MUTATION [...] | -M FILE) -mc CHAIN -c CLUSTER [OPTIONS]
```

For every mutation it builds the mutant residue, prepares a sphere centred on
the mutated residue, and sets up the requested legs. It accepts the protocol
options of `qresfep` (with `-R` for the number of replicates) and the options
below.

| Option | Required | Default | Description |
|---|---:|---|---|
| `-i`, `--input_pdb_file` | yes | — | Prepared protein PDB; mutations are numbered as in this file |
| `-m`, `--mutations` | one of | — | Mutations on the command line, such as `LEU39ALA L39G` |
| `-M`, `--mutations-file` | one of | — | File with one mutation per line; blank lines and `#` comments are ignored |
| `-mc`, `--mutchain` | yes | — | Chain of the mutated residues |
| `-c`, `--cluster` | yes | — | Profile name from `CLUSTER_DICT` |
| `-r`, `--sphereradius` | no | `25` | Sphere radius, Å |
| `-cog`, `--center_of_geometry` | no | per mutation | One sphere centre `x y z` for the whole series |
| `-nbo`, `--neutralize_boundary_offset` | no | `3.0` | Width, in Å, of the outer layer in which ionisable residues are neutralised |
| `-mp`, `--mutant-pdbs` | no | build with PyMOL | Directory of mutant residue PDBs named `<MUTANT><POSITION>.pdb` |
| `--legs` | no | both | `protein`, `tripeptide`, or both |
| `-o`, `--output-dir` | no | current directory | Where the series is written |
| `-n`, `--validate-only` | no | off | Check the structure and every mutation, then stop |

Output:

```text
mutants/                    mutant residue PDBs built by PyMOL
work/<MUTATION>/            the prepared sphere of each mutation
protein/FEP_<MUTATION>/
tripeptide/FEP_<MUTATION>/
```

```bash
setup_resFEP -i 2LZM_prep.pdb -M mutations_neutral.txt -mc A \
  -FF OPLSAAM -r 25 -c SNELLIUS -w 25 -s exponential -T 298 -R 10
```

(qresfep-analyze-cli)=
## `qresfep_analyze`

```text
qresfep_analyze [-p PROTEIN_DIR] [-t TRIPEPTIDE_DIR] [-T TEMPERATURE]
  [-e {dG,dG_forward,dG_reverse,dG_bar}] [-exp CSV] [-o CSV] [-log LEVEL]
```

| Option | Default | Description |
|---|---|---|
| `-p`, `--protein` | `protein` | Directory with the `FEP_*` directories of the protein leg |
| `-t`, `--tripeptide` | `tripeptide` | Directory with the `FEP_*` directories of the reference leg |
| `-T`, `--temperature` | `298` | Temperature subdirectory to read |
| `-e`, `--estimator` | `dG_bar` | `dG_bar` (Bennett acceptance ratio), `dG` (Zwanzig), `dG_forward`, or `dG_reverse` |
| `-exp`, `--experimental` | none | CSV with `mutation` and `ddG_exp` columns; RMSE, MUE and R are logged |
| `-o`, `--output` | `resfep_results.csv` | Results table |
| `-log`, `--log-level` | `info` | Logging threshold |

The results table has one row per mutation, with the columns `mutation`,
`wild_type`, `position`, `mutant`, `dG_protein`, `dG_protein_sem`,
`n_protein`, `dG_tripeptide`, `dG_tripeptide_sem`, `n_tripeptide`,
`ddG_fold`, `ddG_fold_sem` and `crashes` (the number of failed replicates in
both legs). Free energies are in
kcal/mol. A mutation with a missing leg has `NaN` for `ddG_fold`.

```bash
qresfep_analyze -p protein -t tripeptide -T 298 -exp experimental.csv
```
