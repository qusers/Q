# Quickstart

This quickstart prepares a relative binding free-energy network using the Tyk2
example. It demonstrates file generation; production sampling still requires a
validated cluster configuration and deliberate protocol choices.

## 1. Install and build

Follow [Installation](installation.md), then verify both the Python commands and
the required Q binaries. For a cluster calculation, build `qdynp` as well as the
serial programs.

## 2. Parameterize the ligand series

From the repository root:

```bash
cd tutorials/Tyk2/ligprep
qparams -i tyk2_ligands.sdf -p 4
```

The default partial-charge method is NAGL. Add `-am1bcc` when AM1-BCC is required.
For each molecule, inspect the generated `.pdb`, `.lib`, and `.prm` files and
confirm the formal and partial charges.

## 3. Build the perturbation network

Choose one network generator.

### LoMap

```bash
qlomap -i tyk2_ligands.sdf --exp-key r_exp_dg
```

This splits the input SDF into `tyk2_ligands/` and writes
`tyk2_ligands/lomap.json`.

### Konnektor

```bash
qkonnektor -i tyk2_ligands.sdf --exp-key r_exp_dg
```

This writes `tyk2_ligands/mapping.json`. Inspect every edge for chemical
plausibility, alignment quality, charge changes, and network connectivity.

## 4. Calculate the sphere center

```bash
qcog -i tyk2_ligands.sdf
```

`qcog` prints the average center of geometry. Save the values produced by your
input; do not copy coordinates blindly from an example.

## 5. Prepare the receptor and water sphere

Move to the AMBER preparation example and substitute your `qcog` coordinates:

```bash
cd ../setupFEP/amber
qprep_prot -i protein.pdb -FF AMBER14sb \
  -cog -4.689 26.119 -30.570 -r 25
```

Important outputs include `water.pdb`, `dualtop.top`, `qprep.out`, and a prepared
structure. Read `qprep.out` and visualize the structure before continuing.

:::{warning}
Although the parser currently labels `-cog` optional, the implementation requires
it in practice. Always supply all three coordinates.
:::

## 6. Assemble a clean setup directory

A setup directory needs:

- `protein.pdb` and `water.pdb`;
- each ligand's `.pdb`, `.lib`, `.prm`, and individual `.sdf`;
- one network JSON (`lomap.json` or `mapping.json`).

The checked-in tutorial describes the exact copy commands for its layout. For a
new project, create a new directory rather than writing production output into
the repository tutorial.

## 7. Generate both thermodynamic legs

For a conventional equilibrium-window calculation:

```bash
setupFEP -FF AMBER14sb -r 25 -ts 2fs \
  -j lomap.json -rs 42 -c SNELLIUS
```

Replace `SNELLIUS` with a configured cluster. This creates:

```text
1.water/FEP_<from>_<to>/
2.protein/FEP_<from>_<to>/
```

Each edge contains `inputfiles/`, `FEP_submit.sh`, and the generated Q inputs.
`setupFEP` does **not** submit jobs.

For one edge or one leg, call `qligfep` directly:

```bash
qligfep -l1 ejm_31 -l2 ejm_42 -FF AMBER14sb \
  -s protein -c SNELLIUS -r 25 -rs 42
```

## 8. Inspect before submitting

For at least one edge in each leg:

```bash
less 2.protein/FEP_ejm_31_ejm_42/inputfiles/fep_config.json
less 2.protein/FEP_ejm_31_ejm_42/inputfiles/FEP1.fep
less 2.protein/FEP_ejm_31_ejm_42/inputfiles/eq1.inp
less 2.protein/FEP_ejm_31_ejm_42/FEP_submit.sh
```

Confirm:

- edge direction and ligand names;
- topology and FEP atom indices;
- mapped-atom distance restraints;
- sphere radius, temperature, timestep, windows, and replicates;
- softcore and charge-handling methods;
- executable paths, modules, account, wall time, and task count;
- no unintended cleanup command.

## 9. Submit and monitor

Submit `FEP_submit.sh` from each edge directory according to local Slurm policy.
Start with a single protein edge. Confirm equilibration and normal termination
before submitting the full protein and water networks.

## 10. Analyze equilibrium results

From the directory containing `1.water`, `2.protein`, and the mapping:

```bash
qligfep_analyze -t Tyk2 -j lomap.json \
  -exp ddg_value -m ddGbar
```

Without experimental values, omit `-exp`. The analyzer writes
`Tyk2_FEP_results.json` and `lomap_ddG.json`; with experimental values it also
writes `Tyk2_ddG_plot.png`.

## Non-equilibrium alternative

Use the same preparation and add `--neq` during setup:

```bash
setupFEP -FF AMBER14sb -r 25 -ts 2fs \
  -j lomap.json -rs 42 -c SNELLIUS \
  --neq --neq-reps 5 --neq-steps 50000 \
  --neq-eq-steps 1000 --neq-relax-steps 5000 \
  --neq-steepness 8 --neq-schedule sigmoidal
```

Analyze completed forward/reverse switches with:

```bash
qligfep_neq_analyze -p 2.protein -w 1.water \
  -T 298 -u kcal -o neq_results.csv \
  -j lomap.json -exp ddg_value -t Tyk2
```

Read the dedicated [equilibrium](qligfep/equilibrium.md) and
[non-equilibrium](qligfep/neq.md) pages before choosing a production protocol.
