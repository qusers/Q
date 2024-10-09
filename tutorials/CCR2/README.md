# CCR2 Demo

## Prepare a new water.pdb file

```bash
mkdir qprep && cp protein.pdb qprep && cd qprep
qprep_prot -i protein.pdb -FF OPLS2015 -cog 24.357 51.502 31.997
# COG taken from qcog -i aligned_ligands.sdf
```
## generate the ligand perturbation network

```bash
qlomap -i aligned_ligands.sdf
```
## create the parameters

```bash
qparams -i aligned_ligands.sdf -p 4 -nagl
```
## move the sdfs & lomap.json to the same directory to run the perturbations

```bash
cp aligned_ligands/*.sdf .
```
## run setup FEP

```bash
setupFEP -c SNELLIUS -FF OPLS2015 -rs 42 -j lomap.json -rest hybridization_p
```