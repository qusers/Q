# QresFEP: amino-acid mutation free energies

This tutorial shows you how to set up, run, and analyze point-mutation free
energy calculations for T4 lysozyme. The workflow follows the QresFEP-2
hybrid-topology method described by
[Koenekoop et al. (2025)](https://www.nature.com/articles/s42004-025-01771-0).

## What QresFEP calculates

For protein stability, QresFEP models the mutation in two environments: the
folded protein and a capped reference peptide that represents the unfolded
state. Calculate the change in folding free energy with this thermodynamic
cycle:

$$ \Delta \Delta G_{\text{fold}} = \Delta G_{\text{calc(protein)}} - \Delta G_{\text{calc(tripeptide)}} $$

Use the same wild-type-to-mutant transformation, hybrid topology, and
restraints in both legs so that you can subtract their free energies
consistently. You need both legs to calculate a result.

QresFEP-2 is a **hybrid-topology** method. The wild-type and mutant residues
share one set of backbone coordinates, but their side chains have separate
coordinates and parameters. QresFEP restrains equivalent side-chain heavy atoms
that overlap initially. These restraints prevent the side chains from moving
into neighboring atoms ("flapping"). QresFEP also disables non-bonded
interactions and bonded terms that would directly couple the two side chains.

The workflow is:

1. Build the mutant side chains
2. Prepare the protein
3. Generate the FEP input files
4. Submit
5. Analyze

`setup_resFEP` performs steps 1 through 3 for a list of mutations:

```bash
setup_resFEP -i 2LZM_prep.pdb -M mutations_neutral.txt -mc A -FF OPLSAAM -r 25 -c SNELLIUS \
    -w 25 -s exponential -l 1 -ts 2fs -T 298 -R 10 -clean dcd
```

Each mutation gets its own directory and a water sphere centered on the mutated
residue. The sphere center determines whether the target residue is inside the
simulated region and which ionizable residues in the outer boundary layer are
neutralized. Sections 1 through 3 describe this process for one mutation.

## Prerequisites

- Q6 compiled (`cd src/q6 && make all && make mpi`)
- QligFEP installed (`pip install -e .`)
- **PyMOL** on your `PATH` to build the mutant side chains. Install it in the
  QligFEP environment with
  `micromamba install -n qligfep_new -c conda-forge pymol-open-source`.
  QresFEP calls the executable directly. It builds the reference peptide from a
  packaged fragment, so the tripeptide leg does not require PyMOL. If you provide
  ready-made residue PDBs with `setup_resFEP -mp`, you can omit PyMOL. This option
  is useful on clusters.

## The system

[2LZM](https://www.rcsb.org/structure/2LZM) is a 1.7 Å structure of T4
lysozyme, one of the systems the QresFEP protocol was benchmarked on.

`2LZM.pdb` is the raw RCSB download. Raw PDB files are not directly usable: they
may contain components that should not be simulated, lack atoms, and do not
contain the hydrogens required for MD. The published benchmark structures were
prepared with Schrödinger's [Protein Preparation Wizard](https://www.schrodinger.com/life-science/learn/white-papers/protein-preparation-workflow/):
asparagine/glutamine flips and histidine protonation at pH 7 were considered,
titratable states were assigned using PropKa, non-protein heterogroups were
removed, and only crystal waters directly hydrogen-bonded to the protein were
retained.

`2LZM_-_hbond-opt.pdb` is the prepared, hydrogen-bond-optimized structure.
`2LZM_prep.pdb` additionally uses the protonation-specific histidine names
expected by OPLS-AA/M (`HID`/`HIE`/`HIP`) and removes incompatible terminal
atoms. It is the input structure used in this tutorial.

## 1. Build the mutant side chains

`mutations_neutral.txt` lists ten charge-maintaining examples selected from the
published T4L benchmark. `mutations_charged.txt` contains ten charge-changing
examples for the extended protocol developed for a separate manuscript that is
currently in preparation.

QresFEP needs a separate PDB for each mutant residue, positioned on the
wild-type backbone. `setup_resFEP` runs PyMOL's mutagenesis wizard to place each
mutant side chain and writes the resulting residue PDBs to `mutants/`. To
generate these files manually, run:

```bash
python make_mutants.py 2LZM_prep.pdb mutations_neutral.txt
pymol -c mutagenesis.pml
```

If you already have mutant residue PDBs, pass their directory to
`setup_resFEP -mp`. Use this option for residues that PyMOL cannot build. Its
wizard recognizes the 20 standard amino acids, but not Q's protonation variants
(`ASH`, `GLH`, `HIP`, and `LYN`).

PyMOL writes one file per mutation, such as `ALA39.pdb` or `GLY25.pdb`:

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

`qprep_prot` prepares the 25 Å-radius water sphere used by Q. The published
method uses SCAAS spherical boundary conditions with the local reaction field,
centers the sphere on Cβ (or a side-chain hydrogen for glycine), and neutralizes
ionizable residues in the outer 3 Å boundary layer. The command writes:

| File | Contents |
| --- | --- |
| `protein_processed.pdb` | the protein as handed to `qprep` |
| `water.pdb` | the water sphere |
| `prep.json` | sphere center and radius, enclosed charge, disulfides, and the PDB→Q residue numbering |
| `top_p.pdb`, `dualtop.top`, `complexnotexcluded.pdb` | `qprep` output |

`qresfep` reads `prep.json` to translate residue numbers from your input PDB to
Q's topology. For example, you specify `LEU39ALA` using the input PDB's
numbering, but the FEP files use Q's atom and residue numbers. The prepared PDB
files do not preserve this mapping.

> **Center the sphere on the residue you are mutating.** `qprep_prot` neutralizes
> charged residues outside the boundary. For example, it prepares an
> out-of-sphere `GLU` as `GLH`. If the requested residue has been neutralized,
> `qresfep` reports an error instead of building a hybrid residue from the wrong
> library entry. The coordinates above are the Cβ coordinates of residue 39.
> `setup_resFEP` calculates the appropriate center for each mutation.
>
> If you need one center for a whole series, such as when mutating residues around
> a bound ligand that must remain in the same sphere, pass
> `setup_resFEP -cog x y z`. Before writing any files, `setup_resFEP` checks every
> mutation against that sphere. It rejects the series if a residue extends beyond
> the radius or enters the outer boundary layer where charges are neutralized.

## 3. Generate the FEP input files

```bash
qresfep -m LEU39ALA -mc A -S protein -FF OPLSAAM -c SNELLIUS -w 25 -s exponential -l 1 -ts 2fs -T 298 -r 10
```

The main options:

- `-m` the mutation, as `<wild-type><position><mutant>`, in your input PDB's
  numbering. `L39A` works too.
- `-mc` the chain of that residue.
- `-S` which leg: `protein` or `tripeptide`.
- `-t` controls the reference peptide: `A` uses Ala flanks (AXA, the current
  default), `G` uses Gly flanks (GXG), `X` uses only the capped mutable residue,
  and `Z` preserves its native neighbors (ZXZ). The 2025 paper used ZXZ as its
  pragmatic reference and found no statistically significant effect when these
  four models were compared. This tutorial and the in-preparation
  charge-changing protocol use AXA.
- `-w` gives the number of lambda windows **per stage**. The hybrid-topology
  transformation has two stages, so `-w 25` produces 50 windows in total.
- `-s` controls lambda spacing. With `exponential`, stage 1 uses an exponential
  ladder and stage 2 its reverse. Together they concentrate sampling toward the
  relevant ends of the two-stage pathway. The 2025 paper describes the overall
  endpoint-enriched strategy as sigmoidal; the names refer to different views of
  the complete versus per-stage schedule.
- `-l` the lambda endpoint to start from.
- `-ts`, `-T`, `-r` timestep, temperature and number of replicates.
- `-c` the cluster profile, from `settings.py`.

The transformation is split into two consecutive stages. Stage 1 removes the
wild-type side-chain charges while introducing soft-core treatment of the
side-chain van der Waals interactions. Stage 2 introduces the mutant charges
and removes the soft-core treatment, leaving the mutant side chain interacting
normally and the wild-type side chain as a ghost. Stage 2 starts from the final
coordinates of stage 1.

### Charge-changing mutations

Mutating to or from a charged residue changes the net charge represented by the
alchemical side chain. For these transformations, QresFEP adds SOD or CLA ions
to the reference-peptide sphere so that its total charge matches the prepared
protein sphere. QresFEP places the ions inside the solvent boundary and
restrains them away from its outer layer. The in-preparation manuscript uses
this charge-matched reference leg.

`qresfep` creates `FEP_LEU39ALA/`:

| File | Contents |
| --- | --- |
| `inputfiles/L2A.lib` | the hybrid residue: wild-type and mutant side chains in one library entry, the mutant's atoms lower-cased |
| `inputfiles/FEP1.fep` | stage 1: remove the wild-type charges and introduce soft-core interactions |
| `inputfiles/FEP2.fep` | stage 2: introduce the mutant charges and remove its soft-core interactions |
| `inputfiles/OPLSAAM_merged.prm` | force field plus zero-force terms for the bonded terms that span the two topologies |
| `inputfiles/eq*.inp` | equilibration |
| `inputfiles/md_{1,2}_*.inp` | production, per stage and lambda window |
| `inputfiles/qfep.inp` | free energy analysis input |
| `inputfiles/runSNELLIUS.sh` | the SLURM script that runs both stages |
| `inputfiles/resfep_config.json` | how this directory was generated |
| `FEP_submit.sh` | submits the replicate array |

Next, set up the reference leg. Both commands write to `FEP_LEU39ALA`, so move
the protein leg before you run the second command:

```bash
mkdir -p protein tripeptide
mv FEP_LEU39ALA protein/
qresfep -m LEU39ALA -mc A -S tripeptide -t A -FF OPLSAAM -c SNELLIUS -w 25 -s exponential -l 1 -ts 2fs -T 298 -r 10
mv FEP_LEU39ALA tripeptide/
```

### All mutations at once

`setup_resFEP` performs the same steps for a mutation list. It builds the mutant
residues, prepares a sphere around each mutated residue, generates both legs,
and sorts them into `protein/` and `tripeptide/`:

```bash
setup_resFEP -i 2LZM_prep.pdb -M mutations_neutral.txt -mc A -FF OPLSAAM -r 25 -c SNELLIUS \
    -w 25 -s exponential -l 1 -ts 2fs -T 298 -R 10 -clean dcd
```

The command above retains the 25-window setup from the legacy tutorial. The
published protocol H used 50 unevenly spaced windows in each stage, 20 ps per
window and 10 independent replicas. At a 2 fs timestep, 20 ps corresponds to
10,000 production steps. Use `--production-steps` to change that length, and
`--separate-scaling on|off` to map directly to Q's `separate_scaling` setting.

`--manuscript-settings` refers to the separate manuscript currently in
preparation, not to protocol H in the 2025 paper. It selects the settings for
that campaign in one option: a 25 Å sphere, AXA reference peptide, 50 windows
per stage, exponential/reverse-exponential spacing, 2 fs timestep,
2,500,000-step eq5, 10,000 steps per production window, `separate_scaling off`,
and its ten fixed random seeds. Charge matching is automatic for every
charge-changing perturbation and does not depend on selecting the preset.
Operational options such as `-clean` and `--no-trajectories` remain independent.
For example:

```bash
setup_resFEP -i 2LZM_prep.pdb -M mutations_charged.txt -mc A -c SNELLIUS \
    --manuscript-settings -clean dcd
```

`setup_resFEP` creates this directory layout:

```
mutants/                    the mutant residue PDBs PyMOL built
work/<MUTATION>/            the prepared sphere for each mutation
protein/FEP_<MUTATION>/
tripeptide/FEP_<MUTATION>/
```

Before setup begins, `setup_resFEP` validates the residue names against the
force field library and each mutation against the structure. It reports all
problems together. Use `-n` to run these checks without generating any files.

`run_mutations.sh` provides a shell-loop example of the same batch workflow.

## 4. Submit

```bash
for d in protein/FEP_* tripeptide/FEP_*; do (cd "$d" && bash FEP_submit.sh); done
```

This example submits 10 mutations × 2 legs × 10 replicates, for a total of 200
jobs. Each job runs both FEP stages in sequence. Stage 2 restarts from stage 1's
final coordinates, so you cannot split or reorder the stages.

## 5. Analyze

```bash
qresfep_analyze -p protein -t tripeptide -T 298
```

`qresfep_analyze` discovers mutations from the `FEP_<WT><POS><MUT>` directory
names. Each directory name fully describes its perturbation, so this workflow
does not require the mapping file used by the ligand workflow.

The output file, `resfep_results.csv`, contains one row per mutation. Each row
reports the dG and standard error for both legs, the number of completed
replicates, `ddG_fold` with propagated error, and the number of failed
replicates.

Options:

- `-e` which estimator to report: `dG_bar` (Bennett acceptance ratio, the
  default), `dG` (Zwanzig), or the one-directional `dG_forward` / `dG_reverse`.
  A large gap between the forward and reverse sums is the clearest sign that a
  stage is underconverged.
- `-exp results.csv` scores against experiment; the file needs a `mutation`
  column and a `ddG_exp` column. The command logs RMSE, MUE, and R.
- `-o` where to write the table.

If a mutation is missing a leg, `qresfep_analyze` reports a `NaN` ddG and lists
the mutation in a warning. This behavior keeps unfinished calculations visible.
