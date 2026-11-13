# Tyk2 tutorial

The repository includes a Tyk2 ligand series and prepared example files under
`tutorials/Tyk2/`. Run this tutorial in a disposable copy because commands create
or overwrite parameter, mapping, topology, and run directories.

## Prerequisites

From the repository root, install QligFEP and build Q:

```bash
python -m pip install -e .
make -C src/q6 all COMP=gcc
make -C src/q6 mpi COMP=gcc
```

Choose a valid cluster profile before generating production scripts.

## 1. Ligand parameters

```bash
cd tutorials/Tyk2/ligprep
qparams -i tyk2_ligands.sdf -p 4
```

The repository also contains generated ligand files for inspection. Compare new
output rather than assuming it is byte-identical across dependency versions.

## 2. Perturbation network

### LoMap path

```bash
qlomap -i tyk2_ligands.sdf --exp-key r_exp_dg
```

Expected mapping: `tyk2_ligands/lomap.json`.

### Konnektor path

```bash
qkonnektor -i tyk2_ligands.sdf --exp-key r_exp_dg
```

Expected mapping: `tyk2_ligands/mapping.json`.

Try a redundant network only after inspecting the default:

```bash
qkonnektor -i tyk2_ligands.sdf -n cyclic \
  --exp-key r_exp_dg -o tyk2_cyclic
```

## 3. Optional browser visualization

```bash
cd ..
qmapfep -i ligprep/tyk2_ligands.sdf \
  -l ligprep/tyk2_ligands/lomap.json -wd .
```

Open `qmapfep.html`, import the generated JSON, review edge mappings, and export
an edited JSON if needed.

## 4. Sphere center

```bash
qcog -i ligprep/tyk2_ligands.sdf
```

A checked-in coordinate set prints a center close to
`[-4.689, 26.119, -30.570]`, but use the value from your run.

## 5. Receptor and water sphere

```bash
cd setupFEP/amber
qprep_prot -i protein.pdb -FF AMBER14sb \
  -cog -4.689 26.119 -30.570 -r 25
less qprep.out
```

Visualize the prepared protein and `water.pdb`. Confirm the sphere contains the
binding site and inspect boundary-neutralization decisions.

## 6. Assemble setup inputs

The checked-in `tutorials/Tyk2/setupFEP/` contains individual ligand SDF/PDB/LIB/
PRM files, `protein.pdb`, `water.pdb`, and example mappings. For newly generated
files, copy the outputs from `ligprep` and sphere preparation into a fresh setup
directory.

```text
Tyk2-run/
├── protein.pdb
├── water.pdb
├── lomap.json
├── ejm_31.{sdf,pdb,lib,prm}
├── ejm_42.{sdf,pdb,lib,prm}
└── ...
```

## 7. Prepare an equilibrium network

From that setup directory:

```bash
setupFEP -FF AMBER14sb -r 25 -ts 2fs \
  -j lomap.json -rs 42 -c SNELLIUS
```

Use your cluster instead of `SNELLIUS`. Omit `-clean` for the first run so that
trajectories remain available for diagnosis.

Inspect one small edge, for example `ejm_31` → `ejm_42`:

```bash
less 2.protein/FEP_ejm_31_ejm_42/inputfiles/fep_config.json
less 2.protein/FEP_ejm_31_ejm_42/inputfiles/FEP1.fep
less 2.protein/FEP_ejm_31_ejm_42/inputfiles/eq5.inp
less 2.protein/FEP_ejm_31_ejm_42/FEP_submit.sh
```

Submit one protein edge according to local Slurm practice. After it equilibrates
and terminates normally, test its water counterpart, then scale out.

## 8. Analyze equilibrium output

```bash
qligfep_analyze -t Tyk2 -j lomap.json \
  -exp ddg_value -m ddGbar -v
```

Expected principal products are `Tyk2_FEP_results.json`, `lomap_ddG.json`, and
`Tyk2_ddG_plot.png`, plus verbose CSVs.

## 9. Prepare the non-equilibrium alternative

Use a fresh directory containing the same prepared inputs:

```bash
setupFEP -FF AMBER14sb -r 25 -ts 2fs \
  -j lomap.json -rs 42 -c SNELLIUS \
  --neq --neq-reps 5 --neq-steps 50000 \
  --neq-eq-steps 1000 --neq-relax-steps 5000 \
  -L 8 --neq-schedule sigmoidal
```

Confirm switching input:

```bash
tail -n 4 2.protein/FEP_ejm_31_ejm_42/inputfiles/neq_0.inp
```

After completing both legs:

```bash
qligfep_neq_analyze -p 2.protein -w 1.water \
  -T 298 -u kcal -j lomap.json -exp ddg_value \
  -t Tyk2 -o neq_results.csv
```

Inspect switch counts, failed switches, and protein/water overlap before
interpreting the plot.

## Expected outcomes

Before treating the tutorial as complete, you should have:

- reviewed ligand chemical states, generated charges, and every mapping edge;
- inspected the sphere and boundary residues visually;
- run at least one edge successfully through preparation, equilibration, and
  analysis in both legs;
- accounted for every missing or failed network edge; and
- confirmed that replicate and overlap diagnostics support the reported
  uncertainty.
