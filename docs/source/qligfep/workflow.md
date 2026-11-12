# QligFEP workflow

QligFEP prepares dual-topology relative ligand free-energy calculations for Q.
For every mapping edge it builds the same alchemical change in a solvated ligand
(**water leg**) and receptor-bound complex (**protein leg**), then subtracts the
leg free energies.

## Inputs

A complete project directory normally contains:

| Input | Purpose |
|---|---|
| multi-molecule ligand SDF | parameterization and network construction |
| one individual SDF per ligand | atom mapping/restraints during edge setup |
| one `.pdb`, `.lib`, `.prm` per ligand | Q coordinates and force-field data |
| `protein.pdb` | prepared receptor coordinates |
| `water.pdb` | prepared spherical solvent coordinates |
| mapping JSON | directed perturbation edges and metadata |
| cluster profile | executable paths and Slurm settings |

Ligand names must be consistent across SDF titles, filenames, generated
parameter files, and mapping nodes/edges.

## Stage 1: curate structures

Before automation:

- choose protonation, tautomeric, stereochemical, and formal-charge states;
- use one chemically consistent ligand series;
- place ligands in a common binding-site frame;
- resolve receptor alternate locations, missing atoms, termini, cofactors,
  disulfides, and protonation;
- retain a written record of every modeling choice.

QligFEP does not turn ambiguous chemistry into an unambiguous thermodynamic
cycle.

## Stage 2: parameterize ligands

```bash
qparams -i ligands.sdf -p 4
```

This uses OpenFF 2.3.0 and the
`openff-gnn-am1bcc-1.0.0.pt` NAGL model by default. `-am1bcc` selects AM1-BCC.
Each molecule produces `<name>.lib`, `<name>.prm`, and `<name>.pdb` beside the
input SDF.

Check atom ordering, geometry, formal/net partial charge, atom types, and any
unusual torsions. Cofactor mode (`-pcof`) writes a combined protein/cofactor
library and parameter set and currently supports a limited number of cofactors.

## Stage 3: design the network

Use either:

```bash
qlomap -i ligands.sdf
```

or:

```bash
qkonnektor -i ligands.sdf -n mst -s combined
```

A minimal spanning tree uses few simulations but has no redundant cycle closure.
Redundant, fixed-connectivity, cyclic, or star designs can improve diagnostics at
additional cost. Inspect transformations rather than selecting solely by score.

Charge-changing edges need an explicit protocol decision. `qkonnektor` blocks
them by default in the combined scorer (`--charge_changes_score 0.0`) and can
split formal-charge groups with `--separate_charges`.

## Stage 4: prepare the sphere

Get a center from the aligned series:

```bash
qcog -i ligands.sdf
```

For AMBER naming, normalize the receptor if needed:

```bash
pdb2amber -i receptor.pdb -o protein.pdb
```

Prepare the sphere:

```bash
qprep_prot -i protein.pdb -FF AMBER14sb \
  -cog X Y Z -r 25
```

The wrapper removes crystal waters, handles selected boundary ionizable groups,
can merge cofactors, writes qprep input, runs qprep, extracts `water.pdb`, and
filters fragments outside the sphere. Review `qprep.inp`, `qprep.out`, the final
PDB, charge, and boundary selection.

## Stage 5: generate edge directories

For all mapping edges and both legs:

```bash
setupFEP -FF AMBER14sb -j mapping.json -c MY_CLUSTER -r 25
```

For a single edge/leg:

```bash
qligfep -l1 ligandA -l2 ligandB -FF AMBER14sb \
  -s protein -c MY_CLUSTER -r 25
```

`setupFEP` invokes `qligfep` once per edge **per selected leg**—twice per edge
when both water and protein are prepared—and moves results under `1.water` and
`2.protein`. It prepares files but does not submit them.

## Stage 6: audit generated inputs

Every edge stores `inputfiles/fep_config.json`, including the setup command,
parameters, timestamp, and package version. Also inspect:

- qprep log and generated topology;
- combined ligand coordinates;
- `FEP1.fep` state definitions;
- mapped-atom restraints in equilibration/production inputs;
- lambda schedule or switching schedule;
- force field, charge method, and softcore method;
- Slurm modules, executable paths, resources, and cleanup.

Validate one edge end to end before expanding to the network.

## Stage 7: run both legs

Each edge has `FEP_submit.sh`. Submit according to local cluster policy. Start
with the protein leg because receptor/ligand clashes often fail there first.
Do not treat Slurm `COMPLETED` as proof of Q normal termination.

## Stage 8: analyze

Equilibrium:

```bash
qligfep_analyze -j mapping.json -t TARGET -m ddGbar
```

Non-equilibrium:

```bash
qligfep_neq_analyze -j mapping.json -t TARGET \
  -p 2.protein -w 1.water -T 298 -u kcal
```

Keep per-leg and per-replicate details. Investigate missing edges, hysteresis,
poor overlap, and outliers before interpreting network statistics.

## Output layout

```text
project/
├── mapping.json
├── protein.pdb
├── water.pdb
├── ligandA.{sdf,pdb,lib,prm}
├── ligandB.{sdf,pdb,lib,prm}
├── 1.water/
│   └── FEP_ligandA_ligandB/
│       ├── inputfiles/
│       └── FEP_submit.sh
└── 2.protein/
    └── FEP_ligandA_ligandB/
        ├── inputfiles/
        └── FEP_submit.sh
```
