# Scoring

Three scoring functions are implemented in `qcalc`: X-Score [27],
ChemScore [28] and PMF-Score [29]. X-Score and ChemScore are empirical;
PMF-Score is knowledge-based.

All of them need a loaded topology and a mask. The initial topology, with the
coordinates of the `.top` file, can be scored to verify the atom typing.

Trajectory and restart files can both be scored. When you give the file names
in `qcalc`:

- `,frames=every n` after a name restricts the calculation to every *n*-th
  frame;
- `,frames=n-m` restricts it to frames *n* to *m*;
- `mean` instead of a file name gives the mean of all frames processed since
  the start, or since the last `mean`.

The input is similar for the three functions. The functions read data files
(atom definitions, translation tables, parameter files) that are not
distributed in this repository; you must supply them.

## X-Score

### Input

| Prompt | Input |
|---|---|
| Topology file: | `adm.top` |
| Qcalc> | `xscore` |
| Mask: | `residue 1 407` |
| Mask: | `.` |
| Score initial topology? (yes/no) | `yes` |
| Q-atom (FEP) file: | `lig.fep` |
| Cofactor (. or EOL terminates): | `restype=HEM` |
| Cofactor (. or EOL terminates): | `.` |
| FF translation key: | `qoplsaa` |
| Scoring parameters: | `xscore_default.input` |
| Qcalc> | `go` |
| Enter names of coordinate or restart files | `md.dcd,frames=every 5` then `mean` |

### Cofactors

X-Score uses different typing schemes for protein and ligand atoms. Protein
atoms are typed per residue, from the residue library
`RESIDUE_DEF_XTOOL.dat`. Ligand atoms are typed per atom, from
`ATOM_DEF_XTOOL.dat`. A part of the protein that is not in the residue
library, such as the haem of P450, can be typed with the ligand procedure by
defining it as a cofactor. A cofactor is defined on its own line as
`restype=RES`, where `RES` is the residue name. All atoms in residues of that
name are typed as ligand atoms but otherwise count as part of the protein.
Any number of cofactors can be defined.

### Force field

A force-field translation key is needed to translate Q atom types to Sybyl
atom types. The tables are in `ATOM_TYPE_CONVERSIONS.dat`, shared with
PMF-Score. Another file can be named in the parameter input.

### Parameters

Scoring parameters, output options and data files are given in an input
file. Entering `default` at the prompt uses these values:

| Keyword | Default | Meaning |
|---|---|---|
| `SHOW_ABS` | `NO` | Show the binding score of each atom |
| `SHOW_TOTAL` | `YES` | Show the total binding score |
| `SHOW_LIGAND` | `YES` | Show ligand atoms |
| `SHOW_PROTEIN` | `NO` | Show protein atoms |
| `SHOW_COFACTOR` | `YES` | Show cofactor atoms |
| `APPLY_HPSCORE` | `YES` | Use the hydrophobic contact algorithm |
| `APPLY_HMSCORE` | `YES` | Use the hydrophobic matching algorithm |
| `APPLY_HSSCORE` | `YES` | Use the hydrophobic surface algorithm |
| `RESIDUE_DEFINITIONS` | `residue_def_xtool.dat` | |
| `ATOM_DEFINITIONS` | `atom_def_xtool.dat` | |
| `LOGP_DEFINITIONS` | `atom_def_xlogp.dat` | |
| `SURFACE_DEFINITIONS` | `surface_def_xtool.dat` | |
| `ATOM_TRANSLATIONS` | `atom_type_conversions.dat` | |

Applying more than one hydrophobic algorithm gives a consensus score. The
default scoring coefficients are those of [27].

### Output

With `SHOW_LIGAND` or `SHOW_PROTEIN`, the atoms are listed with their type,
residue, properties, neighbours and bonds, followed by the bonds and aromatic
rings that were detected.

For each frame, the contribution of each term is shown with the total: van
der Waals (VDW), hydrogen bonding (HB), hydrophobic contact (HP), matching
(HM) and surface (HS), and rotation (RT). With `SHOW_ABS` the contribution of
every ligand atom is shown. The atomic binding score is always shown for the
initial configuration.

X-Score results are in p*K*<sub>d</sub> units.

## ChemScore

### Input

| Prompt | Input |
|---|---|
| Topology file: | `adm.top` |
| Qcalc> | `chemscore` |
| Mask: | `residue 1 407` |
| Mask: | `.` |
| Score initial topology? (yes/no) | `yes` |
| Q-atom (FEP) file: | `lig.fep` |
| Parameter file: | `chemscore_oplsaa.prm` |
| Qcalc> | `go` |
| Enter names of coordinate or restart files | `md.dcd,frames=every 5` then `mean` |

### Parameter file

ChemScore reads all atom parameters from one parameter file, which defines the
properties of all atom types. There is a file for each force field.

### Output

Before scoring, ChemScore prints the atom type and bond information of all
ligand atoms and the rings it detected. For every frame it prints the
contribution of each term (hydrogen bonds, metal contacts, lipophilic
contacts and frozen rotatable bonds) and the total.

ChemScore results are in kJ/mol. A negative score is a negative energy.

## PMF-Score

### Input

| Prompt | Input |
|---|---|
| Topology file: | `adm.top` |
| Qcalc> | `pmfscore` |
| Mask: | `residue 1 407` |
| Mask: | `.` |
| Score initial topology? (yes/no) | `yes` |
| Q-atom (FEP) file: | `lig.fep` |
| Force field translation key: | `qoplsaa` |
| Scoring parameters: | `pmfscore_default.input` |
| Qcalc> | `go` |
| Enter names of coordinate or restart files | `md.dcd,frames=every 5` then `mean` |

All atoms defined as solvent are ignored. Water molecules that matter must be
defined as part of the protein, or they are left out.

### Force field

As for X-Score, a translation key is needed. The Sybyl type is used only to
determine the hybridisation of carbon and nitrogen atoms.

### Scoring parameters

Output options, data files and the maximum ring size are given in an input
file. The output options are similar to those of X-Score. The maximum ring
size is the number of steps the ring search takes in each direction. Too low
a value misses rings; too high a value makes the search slow.

### Output

With `SHOW_LIGAND` or `SHOW_PROTEIN`, the atoms are listed as for X-Score,
followed by the rings that were detected and, if requested, the bonds. The
atomic binding score is shown only for the initial configuration.

PMF-Score results should be read as rankings, in which a more negative score
means better binding. For the conversion to a binding free energy see [29].
