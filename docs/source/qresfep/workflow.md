# QresFEP workflow

QresFEP prepares free-energy calculations of amino-acid mutations for Q. The
protocol available in this repository is **protein thermostability**: the
change in folding free energy caused by a point mutation. It follows the
QresFEP-2 hybrid-topology method of
[Koenekoop et al. (2025)](https://doi.org/10.1038/s42004-025-01771-0).

For a worked example, see the T4 lysozyme [tutorial](../tutorials/index.md).

## The thermodynamic cycle

The mutation is simulated in two environments: in the folded protein
(**protein leg**) and in a capped reference peptide that represents the
unfolded state (**tripeptide leg**). The change in folding free energy is the
difference between the two:

$$
\Delta\Delta G_{\mathrm{fold}} = \Delta G_{\mathrm{protein}} - \Delta G_{\mathrm{tripeptide}}
$$

Both legs use the same wild-type-to-mutant transformation, hybrid topology and
restraints, so that their free energies can be subtracted. A result needs both
legs.

## Hybrid topology and the two stages

The wild-type and the mutant residue share one set of backbone atoms. Their
side chains have separate coordinates and parameters. QresFEP

- restrains equivalent side-chain heavy atoms that overlap in the starting
  structure, which keeps the side chains from swinging into neighbouring
  atoms, and
- removes the non-bonded interactions and bonded terms that would couple the
  two side chains to each other.

The transformation runs in two consecutive stages:

| Stage | FEP file | What changes |
|---|---|---|
| 1 | `FEP1.fep` | the wild-type side-chain charges are removed and the side-chain van der Waals interactions become soft-core |
| 2 | `FEP2.fep` | the mutant charges are introduced and the soft-core treatment is removed |

After stage 2 the mutant side chain interacts normally and the wild-type side
chain is a non-interacting ghost. Stage 2 starts from the final coordinates of
stage 1, so the stages cannot be run separately or in another order.
[Soft-core potentials](../concepts/softcore.md) explains the soft-core
potential of Q.

## Inputs

| Input | Purpose |
|---|---|
| prepared protein PDB | protonated structure with the residue names of the force field (`HID`, `HIE`, `HIP`, …) |
| list of mutations | `<wild-type><position><mutant>`, such as `LEU39ALA` or `L39A`, numbered as in the input PDB |
| chain identifier | chain of the mutated residues |
| cluster profile | executable paths and Slurm settings; see [Cluster configuration](../cluster-configuration.md) |
| PyMOL, or mutant residue PDBs | coordinates of each mutant side chain on the wild-type backbone |

QresFEP uses the OPLS-AA/M force field (`OPLSAAM`).

## Step 1: build the mutant side chains

QresFEP needs one PDB file per mutant residue, placed on the wild-type
backbone. `setup_resFEP` builds them with the mutagenesis wizard of PyMOL and
writes them to `mutants/`. To use your own files instead, pass their directory
with `-mp`; the files are named `<MUTANT><POSITION>.pdb`, such as `ALA39.pdb`.

The PyMOL wizard knows the 20 standard amino acids, not the protonation
variants of Q (`ASH`, `GLH`, `HIP`, `LYN`). Supply those residues with `-mp`.

## Step 2: prepare a sphere for each mutation

Each mutation gets its own water sphere, centred on the mutated residue: on
Cβ, or on a side-chain hydrogen for glycine. `qprep_prot` builds the sphere
and neutralises the ionisable residues in its outer 3 Å; see
[The spherical boundary](../concepts/spherical-boundary.md) and
[Electrostatics](../concepts/electrostatics.md).

The centre matters. A residue near or beyond the boundary is neutralised, and
a neutralised residue is no longer the residue that the mutation names; in
that case `qresfep` stops with an error.

`qprep_prot` also writes `prep.json`, which records the sphere centre and
radius, the charge inside the sphere, the disulfide bonds and the mapping from
the residue numbers of the input PDB to those of the Q topology. `qresfep`
uses it so that mutations can be given in the numbering of the input PDB.

To keep one sphere for a whole series, for example for residues around a
bound ligand, pass `-cog x y z` to `setup_resFEP`. It then checks every
mutation against that sphere before writing any file, and rejects the series
if a mutated residue reaches beyond the radius or into the outer layer.

## Step 3: generate the FEP inputs

`setup_resFEP` does steps 1 to 3 for a list of mutations:

```bash
setup_resFEP -i protein.pdb -M mutations.txt -mc A -FF OPLSAAM -r 25 \
    -c SNELLIUS -w 25 -s exponential -T 298 -R 10
```

It writes:

```text
mutants/                    mutant residue PDBs built by PyMOL
work/<MUTATION>/            the prepared sphere of each mutation
protein/FEP_<MUTATION>/     protein leg
tripeptide/FEP_<MUTATION>/  reference-peptide leg
```

Before it writes anything, `setup_resFEP` checks the residue names against the
force-field library and every mutation against the structure, and reports all
problems together. `-n` runs only these checks.

`qresfep` sets up one leg of one mutation in a directory that already holds
the output of `qprep_prot`. Each `FEP_<MUTATION>` directory contains:

| File | Contents |
|---|---|
| `inputfiles/*.lib` | the hybrid residue: both side chains in one library entry |
| `inputfiles/FEP1.fep`, `FEP2.fep` | the two stages |
| `inputfiles/OPLSAAM_merged.prm` | force field, with zero-force terms for the bonded terms that span the two side chains |
| `inputfiles/eq*.inp` | equilibration |
| `inputfiles/md_{1,2}_*.inp` | production, per stage and λ window |
| `inputfiles/qfep.inp` | free-energy analysis input |
| `inputfiles/run<CLUSTER>.sh` | Slurm script that runs both stages |
| `inputfiles/resfep_config.json` | the settings the directory was generated with |
| `FEP_submit.sh` | submits the replicates |

## Protocol choices

| Choice | Option | Default | Notes |
|---|---|---|---|
| Reference peptide | `-t` | `A` | `A`: Ala–X–Ala. `G`: Gly–X–Gly. `X`: the capped residue alone. `Z`: the native neighbours |
| λ windows | `-w` | 25 | per stage; the total is twice this number |
| λ spacing | `-s` | `exponential` | exponential in stage 1 and reverse exponential in stage 2, which places more windows where each stage changes fastest |
| Production length | `-ps` | 10 000 steps | per window; 20 ps at 2 fs |
| Final equilibration | `-eqs` | 2.5 ns at 2 fs | length of `eq5`, in steps |
| Timestep | `-ts` | `2fs` | |
| Replicates | `-r` (`qresfep`), `-R` (`setup_resFEP`) | 10 | |
| Random seeds | `-rs` or `--seeds` | none | for reproducible runs |
| `separate_scaling` of `qdyn` | `--separate-scaling` | `on` | |

The published QresFEP-2 protocol used 50 windows per stage, 20 ps per window
and 10 replicates.

### Charge-changing mutations

A mutation to or from a charged residue changes the net charge of the system.
For these mutations QresFEP adds sodium or chloride ions to the sphere of the
reference peptide, so that its total charge equals that of the protein sphere.
The ions are placed inside the solvent boundary and restrained away from its
outer layer. This is automatic.
[Electrostatics](../concepts/electrostatics.md) explains why the two spheres
of a charge-changing perturbation should have the same charge.

`--manuscript-settings` selects, in one option, the settings of a
charge-changing protocol that is described in a manuscript in preparation:
a 25 Å sphere, the Ala–X–Ala reference peptide, 50 windows per stage, a 2 fs
timestep, 10 000 steps per window, `separate_scaling off` and ten fixed random
seeds.

## Step 4: submit

Each `FEP_<MUTATION>` directory has its own submission script:

```bash
for d in protein/FEP_* tripeptide/FEP_*; do (cd "$d" && bash FEP_submit.sh); done
```

Every job runs both stages of one replicate. Test one mutation in both legs
before submitting a whole series.

## Step 5: analyse

```bash
qresfep_analyze -p protein -t tripeptide -T 298
```

`qresfep_analyze` finds the mutations from the `FEP_<WT><POS><MUT>` directory
names and writes `resfep_results.csv`, with one row per mutation: the free
energy and standard error of each leg, the number of replicates, and
`ddG_fold` with its propagated standard error. A mutation with a missing leg
gets `NaN` and is listed in a warning.

Compare the forward and reverse estimates (`-e dG_forward`, `-e dG_reverse`).
A large difference between them is the clearest sign that a stage has not
converged. With `-exp`, the results are scored against experimental values.

See the [CLI reference](../cli/qresfep.md) for all options.
