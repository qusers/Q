# CCR2 Demo

## Align CCR2 ligands in [ligand_sdf_files](./ligand_sdf_files) directory

This is done in a jupyter notebook, please have a look at the [align_CCR2_ligands.ipynb](./align_CCR2_ligands.ipynb).

## Prepare a new water.pdb file

```bash
mkdir qprep && cp protein.pdb qprep && cd qprep
qprep_prot -i protein.pdb -FF OPLS2015 -cog 24.357 51.502 31.997
# COG taken from qcog -i aligned_ligands.sdf
```
## generate the ligand perturbation network

First we create a dedicated directory for all the perturbation files.

``` bash
mkdir ../perturbation && cd ../perturbation
```

Once in the perturbation directory, we move the ligands to it and generate the perturbation network using the `qlomap` command.

```bash
cp ../aligned_ligands.sdf .
qlomap -i aligned_ligands.sdf
```

To generate the QmapFEP html file:
```bash
qmapfep -i aligned_ligands.sdf -wd . -l aligned_ligands/lomap.json
```

## create ligand parameters

```bash
qparams -i aligned_ligands.sdf -p 4 -nagl
```
## move the sdfs & lomap.json to the same directory to run the perturbations

```bash
cp aligned_ligands/*.sdf .
cp aligned_ligands/lomap.json .
cp ../qprep/water.pdb .
cp ../qprep/protein.pdb .
```
## run setup FEP

```bash
setupFEP -c SNELLIUS -FF OPLS2015 -rs 42 -j lomap.json -rest hybridization_p -clean .dcd
```
