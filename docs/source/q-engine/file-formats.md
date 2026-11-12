# Q file formats

## Format overview

| File | Format | Produced/consumed by |
|---|---|---|
| PDB | Text, fixed columns | Input/output of qprep; workflow structures |
| `.lib` | Text, residue sections | qprep, QligFEP |
| `.prm` | Text, parameter sections | qprep, QligFEP |
| `.top` | Q binary topology | qprep → qdyn/qcalc |
| `.inp` | Text, bracketed sections | qprep command stream or qdyn controls |
| `.fep` | Text, state definitions | qdyn/qdynp |
| `.re` | Q binary restart | qdyn/qprep/qcalc |
| `.dcd` | Binary trajectory | qdyn → qcalc/viewers |
| `.en` | Q binary state energies | qdyn → qfep |
| mapping `.json` | JSON nodes and edges | QligFEP network/setup/analysis |

Binary Q files are not a substitute for retaining the text inputs and version
that define how they were generated.

## Qdyn input

The main sections are:

| Section | Purpose |
|---|---|
| `[MD]` | Steps, timestep, temperature, thermostat coupling, SHAKE, minimization |
| `[PBC]` | Periodic box and pressure controls |
| `[cut-offs]` | Ordinary, Q-atom, and LRF cutoffs |
| `[sphere]` | Solute shell boundary controls |
| `[solvent]` | Spherical solvent boundary/polarization controls |
| `[intervals]` | Logging, trajectory, energy, pair-list, and volume frequencies |
| `[files]` | Topology, restart, final, trajectory, energy, FEP, restraint paths |
| `[trajectory_atoms]` | Atom mask written to the trajectory |
| `[lambdas]` | One weight per state |
| `[sequence_restraints]` | Range restraints |
| `[atom_restraints]` | Cartesian restraints |
| `[distance_restraints]` | Pair flat-bottom/harmonic restraints |
| `[angle_restraints]` | Angular restraints |
| `[wall_restraints]` | Radial wall restraints |
| `[lambda_scaling]` | Continuous two-state switching schedule |

Required core values are `steps`, `stepsize`, `temperature`, `topology`, and
`final`. Additional files become required when their interval or feature is
enabled.

## FEP file

The minimum conceptual structure is:

```text
[FEP]
states 2

[atoms]
! q-index topology-index
1 4687
2 4688

[change_charges]
! q-index state-1 state-2
1 -0.12 0.00
2  0.00 -0.12
```

Common optional sections define atom types, charge changes, softcore values,
state-dependent bonded terms, exclusions, restraints, and EVB couplings. Two
numbering schemes coexist:

- `[atoms]`, `[change_atoms]`, `[change_charges]`, `[softcore]`, and selected
  Q–Q controls use the compact Q-atom index;
- state-dependent bonds, angles, torsions, impropers, and general exclusions use
  topology atom numbers.

Read generated FEP files alongside the prepared PDB to avoid confusing them.

## Fragment library

```text
{RES}
[info]
SYBYLtype GROUP

[atoms]
1 C1 c3 -0.1200
2 H1 hc  0.0400

[bonds]
C1 H1

[charge_groups]
C1 H1
```

Entry names and atom names must match the PDB. `[connections] head/tail` links
polymer residues. `[impropers]` supplies explicitly enumerated improper torsions
where the force field requires them.

## Parameter file

```text
[options]
vdw_rule arithmetic
scale_14 0.833333
switch_atoms off

[atom_types]
! force-field-specific LJ columns followed by mass
...

[bonds]
...
[angles]
...
[torsions]
...
[impropers]
...
```

`vdw_rule` determines whether atom-type LJ fields represent geometric `A/B` or
arithmetic `R*/epsilon` combination. Follow an existing parameter file from the
same force field rather than translating columns by appearance.

## Atom masks

A mask line combines properties such as `solute`, `heavy`, `excluded`, and
`restrained`; prefix with `not` to negate. Optional atom/residue ranges narrow
the selection. Multiple lines are combined as a union.

## Comments and names

Q text parsers conventionally accept `!`, `#`, and `*` comments. Bracketed
section matching is case-insensitive. In the current parameter parser,
numeric/logical keys are case-sensitive while string keys are normalized; use
the exact canonical key spelling from this reference. File paths, PDB names, and
force-field atom/residue names may also be case-sensitive. Preserve canonical
spelling and verify the initialization log.

The [Q manual](../manual/index.md) contains the fuller section and column
reference.
