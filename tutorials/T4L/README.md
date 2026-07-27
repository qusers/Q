# QresFEP — amino-acid mutation free energies

This tutorial sets up, runs and analyses free energy calculations for point
mutations in T4 lysozyme, following the dual-topology protocol of
[Koenekoop et al. (2025)](https://www.nature.com/articles/s42004-025-01771-0).

What a mutation FEP measures is the shift in folding free energy:

```
ddG_fold = dG_protein − dG_tripeptide
```

Each mutation is therefore run twice — once in the folded protein, once in a
small capped reference peptide — and the two are subtracted. Both legs are
needed before a mutation has a result.

The workflow is:

1. Build the mutant side chains
2. Prepare the protein
3. Generate the FEP input files
4. Submit
5. Analyse

`setup_resFEP` does steps 1–3 for a whole list of mutations:

```bash
setup_resFEP -i 2LZM_prep.pdb -M mutations_neutral.txt -mc A -FF OPLSAAM -r 25 -c SNELLIUS \
             -w 25 -s exponential -l 1 -ts 2fs -T 298 -R 10 -clean dcd
```

Each mutation gets a directory of its own, holding a water sphere centred on the
residue being mutated. That is not a convenience: a sphere describes one place,
and its centre decides both whether the mutated residue is simulated at all and
which *other* charges are neutralised as out-of-sphere. Sections 1–3 below
describe what the command does, one mutation at a time.

## Prerequisites

- Q6 compiled (`cd src/q6 && make all && make mpi`)
- QligFEP installed (`pip install -e .`)
- **PyMOL** on your `PATH`. It builds the mutant side chains, and it caps the
  reference peptide, so the tripeptide leg cannot run without it.

## The system

[2LZM](https://www.rcsb.org/structure/2LZM) is a 1.7 Å structure of T4
lysozyme, one of the systems the QresFEP protocol was benchmarked on.

`2LZM.pdb` is the raw RCSB download. Raw PDB files are not directly usable: they
carry parts that should not be simulated, are often missing side chains or whole
loops, and have no hydrogens. `2LZM_-_hbond-opt.pdb` is the output of
Schrödinger's [Protein Preparation Wizard](https://www.schrodinger.com/life-science/learn/white-papers/protein-preparation-workflow/),
which adds hydrogens and picks protonation states from a predicted hydrogen-bond
network. `2LZM_prep.pdb` is that file after renaming histidines to the
protonation-specific names OPLS-AA/M uses (`HID`/`HIE`/`HIP`) and trimming
terminal atoms. It is the file this tutorial starts from.

## 1. Build the mutant side chains

`mutations_neutral.txt` and `mutations_charged.txt` each list ten mutations from
the published benchmark — charge-maintaining and charge-changing respectively.

QresFEP needs each mutant residue as its own PDB, positioned on the wild-type
backbone. PyMOL's mutagenesis wizard does that, and `setup_resFEP` drives it —
writing the residue PDBs to `mutants/`. To build them yourself instead:

```bash
python make_mutants.py 2LZM_prep.pdb mutations_neutral.txt
pymol -c mutagenesis.pml
```

Pass the directory holding the result to `setup_resFEP -mp`, which is also how to
supply residues PyMOL cannot build — its wizard knows the standard twenty, not
Q's protonation variants (`ASH`, `GLH`, `HIP`, `LYN`).

This writes one file per mutation — `ALA39.pdb`, `GLY25.pdb`, and so on:

```
ATOM      1  N   ALA A  39      38.865  31.620  20.999  1.00  0.00           N
ATOM      2  CA  ALA A  39      39.713  30.684  21.730  1.00  0.00           C
...
ATOM      8 1HB  ALA A  39      41.547  31.477  20.801  1.00  0.00           H
```

PyMOL numbers branched hydrogens as `1HB`, `2HB`, `3HB` while the force field
libraries use `HB1`, `HB2`, `HB3`. `qresfep` rewrites these names as it reads the
file, so no manual `sed` step is needed.

## 2. Prepare the protein

```bash
qprep_prot -i 2LZM_prep.pdb -FF OPLSAAM -r 25 -cog 41.088 31.308 21.828
```

This solvates the system under spherical boundary conditions and writes:

| File | Contents |
| --- | --- |
| `protein_processed.pdb` | the protein as handed to `qprep` |
| `water.pdb` | the water sphere |
| `prep.json` | sphere centre and radius, enclosed charge, disulfides, and the PDB→Q residue numbering |
| `top_p.pdb`, `dualtop.top`, `complexnotexcluded.pdb` | `qprep` output |

`prep.json` is what `qresfep` reads. It exists because a mutation is requested in
the numbering of *your* PDB (`LEU39ALA` on chain A) while every atom index in a
FEP file is in Q's own numbering, and nothing in the PDB files records how the
two relate.

> **Centre the sphere on the residue you are mutating.** `qprep_prot` neutralises
> charged residues that lie outside the boundary — an out-of-sphere `GLU` is
> prepared as `GLH`. If the residue you name has been neutralised, `qresfep` stops
> with an error rather than silently building a hybrid residue from the wrong
> library entry. The coordinates above are the CB of residue 39; `setup_resFEP`
> recomputes them for each mutation.
>
> To keep one centre for a whole series — mutating residues around a bound ligand
> that has to stay in the same sphere — pass `setup_resFEP -cog x y z`. Every
> mutation is then checked against that sphere before anything is written, and the
> series is refused if a residue reaches beyond the radius or into the restrained
> shell where charges are neutralised.

## 3. Generate the FEP input files

```bash
qresfep -m LEU39ALA -mc A -S protein    -FF OPLSAAM -c SNELLIUS -w 25 -s exponential -l 1 -ts 2fs -T 298 -r 10
```

The main options:

- `-m` the mutation, as `<wild-type><position><mutant>`, in your input PDB's
  numbering. `L39A` works too.
- `-mc` the chain of that residue.
- `-S` which leg: `protein` or `tripeptide`.
- `-t` what flanks the mutated residue in the reference peptide: `A` (Ala, the
  default), `G` (Gly), `X` (none), or `Z` (the native neighbours). `Z` is less
  stable for charge-changing mutations.
- `-w` lambda windows **per stage**. A dual-topology mutation runs two stages, so
  `-w 25` is 50 windows in total.
- `-s` lambda spacing. `exponential` gives stage 1 exponential and stage 2
  reverse-exponential spacing, concentrating windows at the end of each stage,
  where that stage's own topology is switching and the free energy moves fastest.
- `-l` the lambda endpoint to start from.
- `-ts`, `-T`, `-r` timestep, temperature and number of replicates.
- `-c` the cluster profile, from `settings.py`.

This creates `FEP_LEU39ALA/`:

| File | Contents |
| --- | --- |
| `inputfiles/L2A.lib` | the hybrid residue: wild-type and mutant side chains in one library entry, the mutant's atoms lower-cased |
| `inputfiles/FEP1.fep` | stage 1 — discharge the wild-type side chain, mutant present as a soft-core ghost |
| `inputfiles/FEP2.fep` | stage 2 — charge and grow in the mutant, turn the wild type into a ghost |
| `inputfiles/OPLSAAM_merged.prm` | force field plus zero-force terms for the bonded terms that span the two topologies |
| `inputfiles/eq*.inp` | equilibration |
| `inputfiles/md_{1,2}_*.inp` | production, per stage and lambda window |
| `inputfiles/qfep.inp` | free energy analysis input |
| `inputfiles/runSNELLIUS.sh` | the SLURM script driving both stages |
| `inputfiles/resfep_config.json` | how this directory was generated |
| `FEP_submit.sh` | submits the replicate array |

Now the reference leg. **Move the first leg out of the way** — both write to
`FEP_LEU39ALA`:

```bash
mkdir -p protein tripeptide
mv FEP_LEU39ALA protein/
qresfep -m LEU39ALA -mc A -S tripeptide -t A -FF OPLSAAM -c SNELLIUS -w 25 -s exponential -l 1 -ts 2fs -T 298 -r 10
mv FEP_LEU39ALA tripeptide/
```

### All mutations at once

`setup_resFEP` does all of the above for a whole list — builds the mutant
residues, rebuilds the sphere around each mutated residue, runs both legs, and
sorts them into `protein/` and `tripeptide/`:

```bash
setup_resFEP -i 2LZM_prep.pdb -M mutations_neutral.txt -mc A -FF OPLSAAM -r 25 -c SNELLIUS \
             -w 25 -s exponential -l 1 -ts 2fs -T 298 -R 10 -clean dcd
```

leaving:

```
mutants/                    the mutant residue PDBs PyMOL built
work/<MUTATION>/            the prepared sphere for each mutation
protein/FEP_<MUTATION>/
tripeptide/FEP_<MUTATION>/
```

The whole series is checked before any of it is set up — residue names against
the force field library, then every mutation against the structure — and all
problems are reported together. `-n` runs that check and stops.

`run_mutations.sh` is the same thing as a shell loop, kept as a worked example of
what the command does.

## 4. Submit

```bash
for d in protein/FEP_* tripeptide/FEP_*; do (cd "$d" && bash FEP_submit.sh); done
```

That is 10 mutations × 2 legs × 10 replicates = 200 jobs. Each job runs both FEP
stages in sequence: stage 2 restarts from stage 1's final coordinates, so they
cannot be split apart or reordered.

## 5. Analyse

```bash
qresfep_analyze -p protein -t tripeptide -T 298
```

Mutations are discovered from the `FEP_<WT><POS><MUT>` directory names — the name
fully describes the perturbation, so unlike the ligand workflow there is no
mapping file to keep in step.

The output, `resfep_results.csv`, carries one row per mutation: each leg's dG with
its standard error over replicates, the number of replicates that produced a
result, `ddG_fold` with propagated error, and a count of failed replicates.

Options:

- `-e` which estimator to report: `dG_bar` (Bennett acceptance ratio, the
  default), `dG` (Zwanzig), or the one-directional `dG_forward` / `dG_reverse`.
  A large gap between the forward and reverse sums is the clearest sign that a
  stage is underconverged.
- `-exp results.csv` scores against experiment; the file needs a `mutation`
  column and a `ddG_exp` column. RMSE, MUE and R are logged.
- `-o` where to write the table.

Mutations missing a leg are reported with a `NaN` ddG and listed in a warning, so
an unfinished set is visible rather than quietly dropped.
