# T4 lysozyme tutorial

This tutorial calculates the change in folding free energy for ten point
mutations of T4 lysozyme with QresFEP. The files are in `tutorials/T4L/`. Run
the tutorial in a copy of that directory, because the commands write many
files next to the inputs. Read the [QresFEP workflow](../qresfep/workflow.md)
for the method.

## Prerequisites

From the repository root, install the Python package and build Q:

```bash
python -m pip install -e .
make -C src/q6 all COMP=gcc
make -C src/q6 mpi COMP=gcc
```

You also need PyMOL on your `PATH`
([installation](../installation.md#pymol-for-qresfep)) and a valid
[cluster profile](../cluster-configuration.md).

## The files

| File | Contents |
|---|---|
| `2LZM.pdb` | [2LZM](https://www.rcsb.org/structure/2LZM), a 1.7 Å structure of T4 lysozyme, as downloaded |
| `2LZM_-_hbond-opt.pdb` | the structure after protein preparation: hydrogens added, protonation states assigned, hydrogen-bond network optimised |
| `2LZM_prep.pdb` | the same, with the histidine names of OPLS-AA/M (`HID`, `HIE`, `HIP`); **the input of this tutorial** |
| `mutations_neutral.txt` | ten mutations that keep the net charge |
| `mutations_charged.txt` | ten mutations that change the net charge |
| `make_mutants.py`, `run_mutations.sh` | scripts for doing the steps by hand |

A structure from the PDB cannot be simulated as it is: it lacks hydrogens and
may lack atoms or contain molecules that should not be simulated. Prepare your
own structures in the same way before using them.

The mutations of this tutorial are those of `mutations_neutral.txt`:

```text
LEU39ALA
SER44ALA
MET120ALA
ASN55GLY
GLN105GLY
ILE3THR
PHE153LEU
THR26SER
VAL75THR
TYR25GLY
```

## 1. Check the inputs

```bash
cd tutorials/T4L
setup_resFEP -i 2LZM_prep.pdb -M mutations_neutral.txt -mc A -c SNELLIUS -n
```

With `-n`, `setup_resFEP` checks the structure and every mutation and writes
nothing. Use your cluster profile instead of `SNELLIUS`.

## 2. Set up the series

```bash
setup_resFEP -i 2LZM_prep.pdb -M mutations_neutral.txt -mc A -FF OPLSAAM -r 25 \
    -c SNELLIUS -w 25 -s exponential -l 1 -ts 2fs -T 298 -R 10 -clean dcd
```

For each mutation this builds the mutant side chain with PyMOL, prepares a
25 Å sphere centred on the mutated residue, and sets up both legs with 25 λ
windows per stage and 10 replicates. `-clean dcd` deletes the trajectories
when a run finishes; omit it for a first run, so that the trajectories remain
available for inspection.

The result:

```text
mutants/                 ALA39.pdb, ALA44.pdb, ...
work/LEU39ALA/           the prepared sphere of each mutation
protein/FEP_LEU39ALA/    protein leg
tripeptide/FEP_LEU39ALA/ reference-peptide leg
```

## 3. Inspect one mutation

```bash
less work/LEU39ALA/prep.json
less protein/FEP_LEU39ALA/inputfiles/resfep_config.json
less protein/FEP_LEU39ALA/inputfiles/FEP1.fep
less protein/FEP_LEU39ALA/inputfiles/FEP2.fep
```

Open the prepared sphere in a molecular viewer and confirm that residue 39 is
at its centre.

## 4. Submit

Submit one mutation in both legs first:

```bash
(cd protein/FEP_LEU39ALA && bash FEP_submit.sh)
(cd tripeptide/FEP_LEU39ALA && bash FEP_submit.sh)
```

When both have finished normally, submit the rest:

```bash
for d in protein/FEP_* tripeptide/FEP_*; do (cd "$d" && bash FEP_submit.sh); done
```

The whole series is 10 mutations × 2 legs × 10 replicates = 200 jobs. Each
job runs both stages of the transformation.

## 5. Analyse

```bash
qresfep_analyze -p protein -t tripeptide -T 298
```

`resfep_results.csv` has one row per mutation, with the free energy of each
leg, `ddG_fold` and its standard error, and the number of replicates. To
check convergence, compare the two directions:

```bash
qresfep_analyze -p protein -t tripeptide -T 298 -e dG_forward -o forward.csv
qresfep_analyze -p protein -t tripeptide -T 298 -e dG_reverse -o reverse.csv
```

## Setting up one mutation by hand

`setup_resFEP` runs three commands for each mutation. To run them yourself
for `LEU39ALA`, first build the mutant residues:

```bash
python make_mutants.py 2LZM_prep.pdb mutations_neutral.txt
pymol -c mutagenesis.pml
```

Then prepare the sphere. The centre is the Cβ atom of residue 39:

```bash
qprep_prot -i 2LZM_prep.pdb -FF OPLSAAM -r 25 -cog 41.088 31.308 21.828
```

Then set up the two legs. Both commands write to `FEP_LEU39ALA`, so move the
first leg away before the second command:

```bash
qresfep -m LEU39ALA -mc A -S protein -FF OPLSAAM -c SNELLIUS \
    -w 25 -s exponential -l 1 -ts 2fs -T 298 -r 10
mkdir -p protein tripeptide
mv FEP_LEU39ALA protein/

qresfep -m LEU39ALA -mc A -S tripeptide -t A -FF OPLSAAM -c SNELLIUS \
    -w 25 -s exponential -l 1 -ts 2fs -T 298 -r 10
mv FEP_LEU39ALA tripeptide/
```

## Charge-changing mutations

`mutations_charged.txt` holds ten mutations to or from a charged residue. Set
them up in a separate copy of the tutorial directory:

```bash
setup_resFEP -i 2LZM_prep.pdb -M mutations_charged.txt -mc A -c SNELLIUS \
    --manuscript-settings
```

See [Charge-changing mutations](../qresfep/workflow.md#charge-changing-mutations)
for what `--manuscript-settings` selects.

## Expected outcomes

Before treating the tutorial as complete, you should have:

- inspected the prepared sphere of at least one mutation;
- run at least one mutation through both legs and the analysis;
- accounted for every failed replicate; and
- compared the forward and reverse estimates of each mutation.
