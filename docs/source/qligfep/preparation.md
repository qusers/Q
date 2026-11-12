# System preparation

## Ligand SDF requirements

Use one molecule per SDF record with a unique title/name. Preserve hydrogens and
stereochemistry appropriate to the chosen chemical state. If experimental values
will be propagated into the network, store them under one consistently named SDF
property.

`qparams` writes outputs using molecule names, so duplicates or filesystem-unsafe
names will cause collisions or awkward edge paths.

## Protein naming

Qprep requires residue and atom names matching the selected `.lib` file.
`pdb2amber` converts common PDB conventions to those expected by the bundled
AMBER library. It also resolves alternate locations and identifies protonation
from present atoms; this makes the input's hydrogen model consequential.

```bash
pdb2amber -i receptor.pdb
# writes receptor_renamed.pdb
```

Always compare the result to the original structure. Pay special attention to
histidines, acidic residues, termini, ions, NME/NMA caps, alternate locations,
and atom-name uniqueness.

## Center of geometry

```bash
qcog -i ligands.sdf
```

For SDF input, RDKit removes hydrogens, computes each molecule's center, then
averages the molecule centers. For PDB input, records beginning with `ATOM` or
`HETATM` are included by default; change this with `-incl`.

A center computed from multiple aligned ligands is useful for a series. Poorly
aligned outliers can displace it, so visualize the sphere.

## Protein sphere preparation

```bash
qprep_prot -i protein.pdb -FF AMBER14sb \
  -cog X Y Z -r 25
```

Useful controls include:

- `-b auto` or explicit disulfide atom-index pairs;
- `-cof cofactor1.pdb cofactor2.pdb`;
- `-sp 3.0` minimum solvent packing distance;
- `-nbo 3.0` boundary neutralization offset;
- `-sbc 4.0` salt-bridge distance cutoff;
- `-skip-n` to disable boundary neutralization;
- `-skip-ff` to retain fragments even if they are entirely outside the sphere.

The neutralization offset must lie between zero and the sphere radius.

## Boundary neutralization

A truncated spherical model can leave ionizable groups near or outside the
boundary in an unfavorable electrostatic environment. The preparation wrapper
neutralizes selected boundary groups by modifying residue forms and considers
salt-bridge pairing. This changes the model and must be recorded.

After preparation, confirm:

- which residues were changed;
- whether a buried or catalytic group was altered unexpectedly;
- whether salt bridges were preserved;
- the final system charge;
- that complete charge groups are treated consistently near the cutoff.

Use `-skip-n` only as a deliberate protocol choice, not as a workaround for an
unreviewed preparation error.

## Disulfides and cofactors

Automatic disulfide detection considers close cysteine sulfur atoms. Confirm
every proposed bond and supply explicit atom pairs when needed. Cofactors require
coordinates plus force-field-compatible library/parameter definitions; the PDB
alone is not enough.

## Water and exclusion checks

QligFEP removes waters too close to the system during edge preparation. The
default threshold is 1.4 Å (`-wath`); `-wath-ligo` limits that check to ligand
atoms. Review waters in the binding site rather than assuming all removals are
benign.

Key preparation products include:

- `qprep.inp` — exact qprep commands;
- `qprep.out` — warnings, charge, and topology-generation output;
- `dualtop.top` — prepared topology;
- `complexnotexcluded.pdb` — active/unexcluded system used for index checks;
- `water.pdb` — extracted water sphere;
- `protein_processed.pdb` when preprocessing modifies the input.

## Before continuing to FEP setup

Review the preparation in four areas:

| Area | What to verify |
|---|---|
| Topology | Every residue and atom matched a library entry; disulfides and other covalent links are correct; no duplicate or coincident atoms remain. |
| Chemical state | Ligand and cofactor net charges are correct; histidines, acidic residues, and termini have the intended protonation. |
| Spherical model | The center and radius contain the intended flexible region; boundary-neutralized residues and retained binding-site waters were reviewed. |
| Structure quality | High bonded-energy diagnostics were investigated and the prepared coordinates were inspected in a molecular viewer. |

Proceed only when the generated structure and `qprep.out` support the intended
physical model.
