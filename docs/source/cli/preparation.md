# Preparation commands

## `qcog`

Calculate and print a center of geometry.

```text
qcog -i FILE [-incl RECORDS] [-log LEVEL]
```

| Option | Required | Default | Description |
|---|---:|---|---|
| `-i`, `--input` | yes | — | PDB or SDF input |
| `-incl`, `--include` | no | `ATOM,HETATM` | Comma-separated PDB record prefixes |
| `-log`, `--log-level` | no | `info` | Logging threshold |

For PDB, coordinates from matching record prefixes are averaged. For SDF, RDKit
removes hydrogens, calculates each molecule center, then averages the molecule
centers. Output is logged as `Center of geometry: [x y z]`; no file is written.
Only `.pdb` and `.sdf` suffixes are accepted.

```bash
qcog -i ligands.sdf
qcog -i complex.pdb -incl HETATM
```

## `pdb2amber`

Normalize common PDB residue/atom conventions for the bundled AMBER definitions.

```text
pdb2amber -i INPUT_PDB [-o OUTPUT_PDB]
```

| Option | Required | Default | Description |
|---|---:|---|---|
| `-i`, `--input` | yes | — | Input PDB |
| `-o`, `--output` | no | `<input-stem>_renamed.pdb` | Output path/name |

The command filters alternate locations, checks atom-name uniqueness, detects
ASP/GLU/histidine protonation from present atoms, renames termini, normalizes
amino-acid/ion names, and converts several proton/oxygen conventions. Explicit
`-o` is resolved relative to the current working directory.

```bash
pdb2amber -i receptor.pdb -o protein.pdb
```

Inspect protonation and deleted alternate locations in the result.

## `qparams`

Generate Q ligand parameter, library, and coordinate files from SDF.

```text
qparams -i INPUT_SDF [-lff OpenFF] [-p JOBS] [-am1bcc]
        [-ff OFFXML] [-nagl-model MODEL] [-pcof] [-pff PROTEIN_FF]
        [-log LEVEL]
```

| Option | Required | Default | Description |
|---|---:|---|---|
| `-i`, `--input` | yes | — | One- or multi-molecule SDF |
| `-lff`, `--ligand-forcefield` | no | `OpenFF` | Ligand parameterizer; currently only `OpenFF` |
| `-p`, `--parallel` | no | `1` | Parallel charge-assignment jobs |
| `-am1bcc`, `--use_am1bcc` | no | off | Use AM1-BCC instead of NAGL |
| `-ff`, `--forcefield` | no | `openff-2.3.0.offxml` | OpenFF force-field file |
| `-nagl-model`, `--nagl-model` | no | `openff-gnn-am1bcc-1.0.0.pt` | NAGL model |
| `-pcof`, `--parametrize-cofactors` | no | off | Write protein + cofactor files |
| `-pff`, `--protein-forcefield` | no | `AMBER14sb` | Protein force field for cofactor mode |
| `-log`, `--log-level` | no | `info` | Logging threshold |

Normal mode writes `<molecule>.lib`, `<molecule>.prm`, and `<molecule>.pdb` beside
the input SDF. Partial charges are rounded while preserving integer net charge.
Cofactor mode accepts a limited set (currently at most six) and writes
`<protein-FF>_plus_cofactor.lib`, `<protein-FF>_plus_cofactor.prm`, and
`all_cofactors.pdb`.

```bash
qparams -i ligands.sdf -p 4
qparams -i difficult_ligand.sdf -am1bcc
qparams -i cofactors.sdf -pcof -pff AMBER14sb
```

(qprep-prot-cli)=
## `qprep_prot`

Prepare a spherical protein system and run bundled qprep.

```text
qprep_prot -i PROTEIN -cog X Y Z
  [-FF FORCEFIELD] [-r RADIUS] [-b BONDS] [-sp DISTANCE]
  [-cof [PDB ...]] [-skip-n] [-nbo OFFSET] [-sbc DISTANCE]
  [-skip-ff] [-log LEVEL]
```

| Option | Required | Default | Description |
|---|---:|---|---|
| `-i`, `--input_pdb_file` | yes | — | Protein PDB |
| `-FF`, `--forcefield` | no | `AMBER14sb` | Bundled name or path base |
| `-cog`, `--center_of_geometry` | **yes in practice** | — | Sphere center `X Y Z` |
| `-r`, `--sphereradius` | no | `25` | Sphere radius in Å |
| `-b`, `--cysbond` | no | `auto` | `atom1_atom2,...` or automatic disulfide detection |
| `-sp`, `--solvent_pack` | no | `3.0` | Minimum solvent packing distance in Å |
| `-cof`, `--cofactors` | no | none | Zero or more cofactor PDB files |
| `-skip-n`, `--skip-neutralization` | no | off | Disable boundary neutralization |
| `-nbo`, `--neutralize_boundary_offset` | no | `3.0` | Neutralization offset inside sphere, Å |
| `-sbc`, `--salt_bridge_cutoff` | no | `4.0` | Salt-bridge cutoff, Å |
| `-skip-ff`, `--skip-fragment-filter` | no | off | Keep entirely out-of-sphere fragments |
| `-log`, `--log-level` | no | `info` | `debug`, `info`, `warning`, `error`, or `critical` |

Valid bundled force-field names are `AMBER14sb`, `CHARMM36`, `OPLS2005`, and
`OPLS2015`. Automatic disulfide detection uses 1.8–2.5 Å for ordinary or mixed
cysteine labels and 1.8–4.0 Å only when both residues are already labeled
`CYX`/`CYD`; inspect every detected bond. The neutralization offset is validated
between zero and the sphere radius.

The command removes crystal waters, merges cofactors, handles boundary
ionizable groups, optionally converts POPC naming, reindexes residues, filters
out-of-sphere fragments, writes `qprep.inp`, and captures qprep output in
`qprep.out`. Expected products include `dualtop.top`, `top_p.pdb`,
`complexnotexcluded.pdb`, and `water.pdb`; preprocessing can also write
`protein_processed.pdb`.

```bash
qprep_prot -i protein.pdb -FF AMBER14sb \
  -cog -4.689 26.119 -30.570 -r 25
```
