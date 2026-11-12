# Q 6 and QligFEP manual

**Molecular dynamics and free-energy calculations in biomolecular systems**

Updated for the current Q/QligFEP source tree (Q engine 6.0.1; QligFEP 2.x).

:::{important}
This is the maintained manual. It replaces the OCR conversion of the 17 August
2017 Q 6 manual. Command names, build paths, supported features, and defaults in
this edition have been checked against the current source tree. The original PDF
is retained in `q-manual/qman.pdf` for historical reference.
:::

## How to use this manual

Q is a Fortran molecular-dynamics suite designed around alchemical free-energy
calculations. QligFEP is the Python command-line workflow in this repository that
prepares, runs, and analyzes relative ligand binding free-energy calculations
with Q.

Use this manual when you need the engine-level concepts and file formats. Use the
separate QligFEP workflow and CLI pages for automation details. Input files are
plain text unless noted otherwise; paths are interpreted relative to the working
directory in which the executable is run.

## Program overview

The current build produces these programs:

| Program | Purpose | Typical invocation |
|---|---|---|
| `qprep` | Build a topology and prepared coordinates | `qprep prep.inp > prep.log` |
| `qdyn` | Serial minimization, MD, equilibrium FEP/EVB, and non-equilibrium switching | `qdyn md.inp > md.log` |
| `qdynp` | MPI version of the MD engine | `mpirun -np 8 qdynp md.inp > md.log` |
| `qdum` | Initialize and validate a `qdyn` input without running the dynamics loop | `qdum md.inp > check.log` |
| `qfep` | Analyze Q energy files | `qfep < qfep.inp > qfep.out` |
| `qcalc` | Analyze restart and trajectory files | `qcalc < qcalc.inp > qcalc.out` |

`qprep`, `qdyn`, `qfep`, and `qcalc` form the usual data flow:

```{mermaid}
flowchart TD
    inputs["PDB coordinates + fragment libraries + parameters"] --> qprep["qprep"]
    qprep --> system["Topology .top + prepared coordinates"]
    system --> dynamics["qdyn serial / qdynp MPI"]
    dynamics --> restart["Restart .re"]
    dynamics --> trajectory["Trajectory .dcd"]
    dynamics --> energy["Q-state energies .en"]
    energy --> qfep["qfep"]
    qfep --> results["Free-energy estimates + profiles"]
```

QligFEP adds ligand parameterization, network generation, dual-topology setup,
Slurm run scripts, and result aggregation around this engine workflow.

## Scientific scope

Q supports the following main calculation families:

- conventional molecular dynamics;
- free-energy perturbation (FEP), including multi-state descriptions of a
  designated perturbation region;
- empirical valence bond (EVB) calculations through state-dependent valence
  terms and off-diagonal couplings;
- linear interaction energy (LIE) workflows;
- spherical boundary conditions for a local region such as a binding site;
- periodic boundary conditions;
- serial and MPI dynamics;
- equilibrium lambda windows and two-state non-equilibrium lambda switching.

Q is force-field agnostic in the sense that atom, bonded, and nonbonded
parameters are read from external library and parameter files. This repository
ships protein definitions for `AMBER14sb`, `CHARMM36`, `OPLS2005`, and
`OPLS2015`. Ligands can be parameterized through the QligFEP `qparams` command.
Do not mix naming conventions or parameter files from different force fields.

## Installation and compilation

### Requirements

The Python workflow currently requires Python 3.13. Building the serial Fortran
programs requires `gfortran`; building `qdynp` additionally requires an MPI
Fortran wrapper such as `mpif90` from an MPI installation built with a compatible
compiler.

The complete Python environment includes scientific Python dependencies plus
OpenFF Toolkit/NAGL, LoMap, Kartograf, FKComBU, and Konnektor. Follow the
installation page in the online documentation for the current environment
command.

### Install QligFEP

From the repository root:

```bash
python -m pip install -e .
```

Installation attempts to compile `qprep`, because the Python preparation
commands require it. A successful Python installation does **not** imply that all
simulation executables have been built.

### Build all serial programs

```bash
cd src/q6
make all COMP=gcc
```

The makefile places executables in `src/q6/bin/q6/`. Verify them with:

```bash
ls -l bin/q6/{qprep,qdyn,qdum,qfep,qcalc}
```

### Build MPI dynamics

Load the compiler and MPI modules used on your cluster, then run:

```bash
cd src/q6
make mpi COMP=gcc
ls -l bin/q6/qdynp
```

The compiler used by `mpif90` must be compatible with the compiler and runtime
used for the serial objects. Cluster module names are site-specific.

### Build checks

```bash
cd src/q6
make test COMP=gcc
```

For a quick input initialization check, use `qdum` on a representative `qdyn`
input. A production installation should also run a small serial and MPI
simulation before launching a network.

## Units and conventions

| Quantity | Unit |
|---|---|
| Time and timestep | fs |
| Temperature | K |
| Distance | Å |
| Energy | kcal mol$^{-1}$ |
| Harmonic distance force constant | kcal mol$^{-1}$ Å$^{-2}$ |
| Charge | elementary charge, e |
| Atomic mass | unified atomic mass unit |
| Angles in input parameter files | degrees unless stated otherwise |

Atom and residue numbering in Q topology-dependent files is one-based. Qprep may
renumber both atoms and residues, so always take indices used in MD restraints or
FEP files from Qprep's prepared structure, not from the original PDB.

## Preparing a topology with `qprep`

### Inputs

Qprep normally reads:

- a fixed-column PDB coordinate file;
- one or more fragment-library (`.lib`) files;
- one force-field parameter (`.prm`) file;
- optionally a solvent coordinate file, topology, restart, or trajectory.

Only `ATOM` and `HETATM` records supply atoms. Residue and atom names must match a
loaded library entry. Separate disconnected molecules with a line containing
`GAP`. Inspect alternate locations, insertion codes, missing heavy atoms, and
protonation states before topology generation.

### Interactive and scripted use

Start an interactive session with:

```bash
qprep
```

For reproducible work, put the same commands in a file:

```bash
qprep prep.inp > prep.log
```

A minimal spherical-solvation example is:

```text
! Ligand-specific and protein/solvent libraries
readlib ligand.lib
readlib AMBER14sb.lib

! One parameter set containing all required types
readprm AMBER14sb.prm

readpdb complex.pdb

! Boundary type 1: center x/y/z and radius
boundary 1 -4.689 26.119 -30.570 25.0
solvate -4.689 26.119 -30.570 25.0 1 HOH

maketop complex
writetop complex.top

checkbonds 6
checkangs 6
checktors 6
checkimps 6

writepdb complex_prepared.pdb y
quit
```

Aliases such as `rl`, `rprm`, `rp`, `mt`, `wt`, and `wp` are accepted, but full
command names are clearer in archived preparation inputs.

### Frequently used Qprep commands

| Command | Common alias | Function |
|---|---|---|
| `readlib FILE` | `rl` | Load a fragment library; repeat for additional libraries |
| `readprm FILE` | `rprm`, `rff` | Load force-field parameters |
| `readpdb FILE` | `rp` | Load PDB coordinates |
| `readtop FILE` | `rt` | Load an existing Q topology |
| `readx FILE` | `rx` | Load coordinates from a restart |
| `trajectory FILE` | `tr` | Open a trajectory |
| `readframe FRAME` | `rf` | Load a numbered frame from the trajectory opened by `trajectory FILE` |
| `readnext` | `rn` | Load the next frame from the open trajectory |
| `boundary ...` | `bc` | Define sphere or box boundary conditions |
| `solvate ...` | `so` | Fill a sphere or box with solvent |
| `addbond I J` | `ab` | Add a cross-link, for example a disulfide bond |
| `xlink` | `xl` | Search interactively for close heavy-atom cross-links |
| `maketop [TITLE]` | `mt` | Generate the topology |
| `checkbonds E` | `cb` | Report bond terms above threshold `E` |
| `checkangs E` | `ca` | Report angle terms above threshold `E` |
| `checktors E` | `ct` | Report torsion terms above threshold `E` |
| `checkimps E` | `ci` | Report improper terms above threshold `E` |
| `writepdb FILE [y\|n]` | `wp` | Write selected coordinates |
| `writemol2 ...` | `wm` | Write a MOL2 structure |
| `writetop FILE [TITLE]` | `wt` | Write a binary Q topology |
| `listseq` | `ls` | List the residue sequence |
| `listres N` | `lr` | List atoms in residue `N` |
| `mask ...` | `ma` | Define the output atom selection |
| `preferences` | `prefs` | Show Qprep preferences |
| `set ...` | `s` | Set a preference |
| `help` | `h`, `?` | Show available commands |
| `quit` | `q` | Exit |

Run `help` in the installed executable for the definitive command set of that
build.

### Atom masks

Masks select subsets for Qprep output, Qdyn trajectories, and several Qcalc
operations. Selectors include `solute`, `heavy`, `excluded`, and `restrained`;
`not` negates a selector. Atom or residue ranges can be appended. Multiple mask
lines are additive.

Examples:

```text
heavy
not excluded
solute residue 10 25
```

An omitted mask means all atoms. Verify any mask by writing a PDB before relying
on it for a production trajectory.

### Validation

Do not continue merely because Qprep wrote a topology. Read its complete log and:

1. confirm every residue and atom matched a library definition;
2. inspect all added hydrogens and cross-links;
3. inspect the prepared PDB in a molecular viewer;
4. investigate high bond, angle, torsion, and improper energies;
5. confirm the sphere center, radius, solvent count, and excluded region;
6. confirm the system and Q-region charges.

QligFEP's `pdb2amber` and `qprep_prot` automate common protein naming,
neutralization, sphere construction, and water extraction tasks, but their output
still requires these checks.

## Force-field files

### Fragment libraries (`.lib`)

A library defines molecular building blocks. Each entry begins with a residue
name in braces and contains sections such as:

```text
{LIG}
[info]
SYBYLtype GROUP

[atoms]
! index name type charge
1 C1 c3 -0.1200
2 H1 hc  0.0400

[bonds]
C1 H1

[charge_groups]
C1 H1
```

Important sections are:

| Section | Meaning |
|---|---|
| `[info]` | Optional metadata, including solvent and output record type |
| `[atoms]` | Sequential index, atom name, atom type, and partial charge |
| `[bonds]` | Intrafragment bonds by atom name |
| `[connections]` | `head` and `tail` atoms for polymer links |
| `[impropers]` | Explicit improper definitions |
| `[charge_groups]` | Nonbonded cutoff groups; switching atom appears first |
| `[build_rules]` | Rules for coordinates of generated atoms |

PDB atom names must match `[atoms]`. The charge of each charge group should be
consistent with the intended force-field convention; incomplete groups near a
cutoff can produce electrostatic artifacts.

### Parameter files (`.prm`)

The parameter file contains these principal sections:

| Section | Contents |
|---|---|
| `[options]` | van der Waals combination rule, 1–4 electrostatic scale, cutoff convention |
| `[atom_types]` | LJ parameters, mass, and optional output type |
| `[atom_aliases]` | Alternate atom-type names |
| `[LJ_type2_pairs]` | Pairs using the alternate LJ set |
| `[bonds]` | atom types, harmonic force constant, equilibrium length |
| `[angles]` | three atom types, force constant, equilibrium angle |
| `[torsions]` | four types/wildcards, force constant, periodicity, phase, paths |
| `[impropers]` | central type matching and improper parameters |

Q supports geometric $A/B$ and arithmetic $R^*/\epsilon$ Lennard-Jones
representations. The `vdw_rule` in `[options]` controls how unlike atom types are
combined. Never infer the rule from a force-field name; inspect the parameter
file.

When extending a force field, check for duplicate definitions, wildcard
precedence, 1–4 values, net residue charge, and every bonded type required by a
Qprep test structure.

## Running molecular dynamics

### Minimal `qdyn` input

Qdyn reads an INI-like file. Section names are case-insensitive. Use the exact
canonical key spelling shown here: numeric and logical key lookup in the current
parser is case-sensitive, while string-key lookup is normalized. `!`, `#`, and
`*` introduce comments. Unknown optional keys can be ignored, so always inspect
the initialization log to confirm every intended setting.

```text
[MD]
steps                 10000
stepsize              2.0
temperature           298
initial_temperature   10
random_seed           12345
bath_coupling         100
shake_solvent         on
shake_hydrogens       on
lrf                    on

[cut-offs]
solute_solute         10
solute_solvent        10
solvent_solvent       10
q_atom                99
lrf                   99

[sphere]
shell_radius          22
shell_force           20

[intervals]
non_bond              25
output                50
energy                10
trajectory            100

[files]
topology              complex.top
final                 md.re
trajectory            md.dcd
energy                md.en
fep                    transformation.fep

[lambdas]
1.000 0.000
```

Without a restart, provide `initial_temperature` so velocities can be generated.
For continuation runs, replace it with `restart previous.re`. `[files] topology`
and `[files] final` are required.

### Core `[MD]` keys

| Key | Meaning | Current engine behavior |
|---|---|---|
| `steps` | Number of integration steps | Required |
| `stepsize` | Timestep in fs | Required |
| `temperature` | Target temperature in K | Required |
| `bath_coupling` | Weak-coupling relaxation time in fs | Default 10; must be at least the timestep |
| `initial_temperature` | Maxwell velocity temperature | Needed when no restart supplies velocities |
| `random_seed` | Velocity/random-operation seed | Default 1 when generating velocities |
| `separate_scaling` | Scale solute and solvent temperatures separately | Default off |
| `shake_solvent` | Constrain solvent geometry | Default on |
| `shake_solute` | Constrain solute bonds | Default off |
| `shake_hydrogens` | Constrain bonds involving hydrogen | Default off |
| `lrf` | Enable local reaction field | Default off |
| `force_rms` | Calculate and print the root-mean-square force | Default off |
| `minimize` | Run steepest-descent minimization before MD | Default off |
| `max_minimize_steps` | Minimizer iteration limit | Default 1000 |
| `minimize_tolerance` | Minimizer force tolerance | Default 0.1 kcal mol⁻¹ Å⁻¹ |
| `minimize_step_size` | Minimizer coordinate step | Default 0.001 Å |

The current source implements the established leapfrog-style MD update,
Berendsen weak temperature coupling, and SHAKE constraints. Do not use input keys
for alternative integrators, thermostats, or constraint methods from older
manual drafts unless they are present in the source version you compiled.
Minimization cannot be combined with non-equilibrium `[lambda_scaling]` mode.

### Nonbonded cutoffs

The `[cut-offs]` section accepts:

- `solute_solute`;
- `solute_solvent`;
- `solvent_solvent`;
- `q_atom`;
- `lrf`.

The conventional non-Q defaults are 10 Å. A value of 99 Å for `q_atom` is used to
make the perturbation region interact with essentially the full simulated
system. The LRF radius cannot be shorter than the ordinary interaction cutoffs.
Cutoff behavior also depends on charge-group and switch-atom settings in the
force field.

### Intervals and output

`[intervals]` controls the following step intervals:

| Key | Default when omitted |
|---|---:|
| `non_bond` | 10 |
| `output` | 10 |
| `temperature` | the `output` interval |
| `energy` | 0 (disabled) |
| `trajectory` | 0 (disabled) |
| `volume_change` | 10 when constant pressure is active |

An interval of zero disables the corresponding saved data. Energy files are
necessary for FEP analysis. Trajectories can dominate disk usage; choose a
scientifically useful interval rather than deleting data before basic quality
control.

The `[files]` section accepts:

| Key | Requirement and purpose |
|---|---|
| `topology` | Required input `.top` |
| `restart` | Initial coordinates/velocities; required if `initial_temperature` is absent |
| `final` | Required final restart path |
| `trajectory` | Required when the trajectory interval is positive; invalid when it is disabled |
| `energy` | Required when the energy interval is positive; invalid when it is disabled |
| `fep` | FEP strategy; required when `[lambdas]` is supplied |
| `restraint` | Optional restart whose coordinates define restraint references |

`[trajectory_atoms]` contains an atom mask and defaults to all atoms when
trajectory output is enabled. Common output extensions are:

| Extension | Contents |
|---|---|
| `.log` | Text written to standard output |
| `.re` | Binary restart coordinates and velocities |
| `.dcd` | Binary trajectory |
| `.en` | Binary state-energy records for qfep |

### Spherical boundary conditions

A spherical system focuses computational effort on the active region. Qprep
stores the sphere center and outer boundary in the topology. During dynamics:

- solute atoms outside the sphere are excluded from nonbonded interactions and
  strongly restrained;
- a shell near the boundary restrains active/inactive covalent interfaces;
- solvent receives radial boundary restraints to maintain the sphere;
- polarization restraints model the orientational behavior of boundary water;
- LRF can approximate electrostatics beyond the explicit cutoff.

The `[sphere]` section controls:

| Key | Default |
|---|---|
| `shell_radius` | 85% of the exclusion-shell radius |
| `shell_force` | 10 kcal mol⁻¹ Å⁻² |
| `exclude_bonded` | off |

The optional `[solvent]` section controls:

| Key | Default/meaning |
|---|---|
| `radius` | Override target solvent radius stored in the topology |
| `radial_force` | 60 kcal mol⁻¹ Å⁻² after sphere initialization |
| `polarization` | on |
| `charge_correction` | follows `polarization`; cannot be on if polarization is off |
| `polarization_force` | 20 kcal mol⁻¹ rad⁻² after sphere initialization |
| `morse_depth` | Radius-dependent engine value |
| `morse_width` | Radius-dependent engine value |

The legacy `model` key is recognized only to report that it is obsolete; it does
not select a solvent model. The topology and loaded solvent library define the
model.

Sphere placement is part of the model, not only a performance option. Confirm
that the ligand and relevant flexible residues are comfortably inside the
unrestrained region. QligFEP commonly uses a 25 Å sphere centered on the ligand
series and neutralizes ionizable groups near the boundary.

### Periodic boundary conditions

A `[PBC]` section activates periodic handling. The box can be rectangular.
Nonbonded cutoffs must not exceed half the shortest box side. Relevant keys are:

| Key | Meaning |
|---|---|
| `rigid_box_center` | Keep the topology box center fixed |
| `constant_pressure` | Enable Monte Carlo volume moves |
| `max_volume_displ` | Required maximum volume displacement when pressure control is on |
| `pressure_seed` | Volume-move random seed; default 3781 |
| `pressure` | Target pressure in bar; default 1.0 |
| `control_box` | Three target box lengths for an NVT box-size change |
| `put_solute_back_in_box` | Wrap solute molecules back into the primary box |
| `put_solvent_back_in_box` | Wrap solvent molecules back into the primary box |

Specify both `put_*_back_in_box` booleans explicitly in a PBC input. The
`control_box` key is used only when `constant_pressure` is off. Spherical solvent
restraints and periodic settings are mutually exclusive. Test periodic
free-energy protocols independently rather than assuming spherical workflow
defaults transfer unchanged.

### Restraints

Qdyn supports these tabular restraint sections and column orders:

| Section | Columns |
|---|---|
| `[sequence_restraints]` | `first_atom last_atom force H_flag to_center` |
| `[atom_restraints]` | `atom x y z force_x force_y force_z state` |
| `[distance_restraints]` | `atom1 atom2 lower upper force state` |
| `[angle_restraints]` | `atom1 atom2 atom3 target_angle force state` |
| `[wall_restraints]` | `first_atom last_atom distance force morse_width morse_depth H_flag` |

For sequence restraints, `H_flag` controls hydrogen inclusion and `to_center`
selects per-atom coordinates (`0`), mass-weighted group-center behavior (`1`), or
center-of-mass behavior (`2`). Distance restraints also accept `residue:atom`
references for the first two columns. `lower` and `upper` define the unrestrained
interval; equal limits make the ordinary harmonic form. `state = 0` applies in
all states, while a positive state number makes the term state-specific and
lambda-scaled. Wall restraints are available only for spherical systems.

QligFEP writes dual-topology interligand distance restraints automatically. Their
indices refer to `complexnotexcluded.pdb`/the prepared topology.

## Free-energy perturbation and EVB

### Q-atoms and states

The FEP file marks a subset of topology atoms as Q-atoms. Charges, atom types,
bonds, angles, torsions, impropers, exclusions, constraints, and restraints can
then differ among states. Qdyn evaluates the state energies separately while
using the lambda-weighted force mixture to propagate coordinates.

Lambda weights are specified in the Qdyn input:

```text
[lambdas]
0.750 0.250
```

The weights must lie in `[0, 1]` and sum to one. In two-state calculations,
state 1 is fully active at `(1, 0)` and state 2 at `(0, 1)`. In `[FEP]`, `offset`
adds an atom-number offset; alternatively, `offset_residue` or `offset_name`
locates the first atom of a residue. Under PBC, the FEP file also requires
`[PBC] switching_atom` for Q-atom pair-list construction.

### FEP file syntax

FEP files use bracketed sections. The compact example below changes Q-atom types
and charges between two states:

```text
[FEP]
states 2
softcore_use_max_potential on
softcore_method standard

[atoms]
1 4687
2 4688

[atom_types]
! name  LJ fields ... mass
DUM 0.0 0.0 0.0 0.0 0.0 0.0 0.0
C1  1.9 0.1 0.0 0.0 1.9 0.1 12.011

[change_atoms]
1 C1  DUM
2 DUM C1

[change_charges]
1 -0.12 0.00
2  0.00 -0.12

[softcore]
1 0.0 50.0
2 50.0 0.0
```

The exact atom-type columns depend on the chosen LJ representation. Generated
files are the safest starting point for ligand work.

### Principal FEP sections

| Section | Purpose | Atom numbering |
|---|---|---|
| `[FEP]` | Number of states, offsets, library-charge and softcore controls | — |
| `[atoms]` | Map Q-atom index to topology atom | Q index + topology index |
| `[atom_types]` | Define Q-region atom types | names |
| `[change_atoms]` | Atom type in each state | Q-atom index |
| `[change_charges]` | Partial charge in each state | Q-atom index |
| `[softcore]` | State-specific softcore value | Q-atom index |
| `[soft_pairs]` | Exponential repulsion pair | Q-atom indices |
| `[excluded_pairs]` | State-dependent nonbonded exclusion | topology indices |
| `[el_scale]` | Scale selected Q–Q electrostatics | Q-atom indices |
| `[bond_types]` / `[change_bonds]` | Define and assign state bonds | topology indices in assignments |
| `[angle_types]` / `[change_angles]` | Define and assign state angles | topology indices |
| `[torsion_types]` / `[change_torsions]` | Define and assign state torsions | topology indices |
| `[improper_types]` / `[change_impropers]` | Define and assign state impropers | topology indices |
| `[shake_constraints]` | Lambda-weighted constraint distances | topology indices |
| `[off_diagonals]` | EVB off-diagonal couplings | state and Q-atom indices |
| `[monitor_groups]` | Define nonbonded reporting groups | topology indices |
| `[monitor_group_pairs]` | Select pairs to report | group indices |

A type code of zero in a state removes the corresponding bonded term. Coupling
sections can scale angles, torsions, or impropers with a changing Morse bond.

### Softcore methods

The current QligFEP CLI exposes three engine methods:

| CLI value | Engine concept |
|---|---|
| `standard` | Established Lennard-Jones softcore |
| `beutler_coul` | Softened LJ and Coulomb effective distance |
| `gapsys` | Short-range force linearization for LJ and Coulomb |

`[FEP]` recognizes `qq_use_library_charges`,
`softcore_use_max_potential`, and `softcore_method`. With
`softcore_method gapsys`, the optional controls and defaults are:

| Key | Default |
|---|---:|
| `gapsys_scale_linpoint_lj` | 0.85 |
| `gapsys_scale_linpoint_q` | 0.3 |
| `gapsys_sigma_lj` | 0.3 |
| `gapsys_sigma_q` | 1.0 |

Softcore potentials prevent singular energies when atoms appear, disappear, or
overlap near an endpoint. Method parameters and endpoint behavior are
scientifically significant; preserve the method in provenance and do not combine
results from different methods as if they were identical protocols.

### Equilibrium window protocol

A two-state equilibrium FEP runs many fixed lambda values. Each window must save
`.en` records at a suitable interval. QligFEP supports linear, sigmoidal,
exponential, and reverse-exponential spacing and normally launches independent
replicates. Adequate overlap between neighboring windows is essential.

### Non-equilibrium switching

A `[lambda_scaling]` section activates continuous two-state switching in serial
`qdyn`:

```text
[lambda_scaling]
scaling_parameter sigmoidal
L_sigmoid 8.0
```

The mode requires two FEP states. QligFEP equilibrates endpoints, launches forward
and reverse switches, extracts work, and applies Bennett acceptance ratio (BAR).
Fixed-lambda equilibration may use `qdynp`, but switching itself uses serial
`qdyn`; QligFEP distributes independent switches across allocated resources.

## Free-energy analysis with `qfep`

Qfep reads standard input:

```bash
qfep < qfep.inp > qfep.out
```

The input gives the number of energy files and states, thermal energy `kT`, the
number of initial samples to skip, energy-gap binning, state energy shifts,
optional off-diagonal terms, reaction-coordinate coefficients, and the ordered
list of `.en` files.

Qfep reports several estimators and diagnostics, including:

- forward and reverse Zwanzig exponential averaging;
- thermodynamic integration;
- overlap sampling;
- Bennett acceptance ratio;
- EVB/reaction free-energy profiles binned by energy gap.

The current qfep analysis arrays support up to four states and 600 bins. Read the
entire output: a final number without adequate neighboring-window overlap,
replicate consistency, and stable forward/reverse estimates is not evidence of a
converged calculation.

QligFEP's `qligfep_analyze` runs/reads qfep for the water and protein legs,
calculates

$$
\Delta\Delta G = \Delta G_{\mathrm{protein}} - \Delta G_{\mathrm{water}},
$$

and writes a result JSON plus an annotated copy of the mapping network.

## Structure and trajectory analysis with `qcalc`

Qcalc is stdin-driven:

```bash
qcalc < qcalc.inp > qcalc.out
```

A session first selects a topology, then one or more calculations, then the
trajectory or restart files to process. Implemented analyses include RMSD, RMSF,
least-squares fitting, distance, angle, torsion, Schlitter entropy, radial
distribution functions, center-of-mass quantities, nonbonded interactions and
residue/Q-atom contributions, plus the bundled X-Score, ChemScore, and PMF-Score
calculations.

Because Q trajectories can contain a mask rather than every topology atom, make
sure the selected trajectory atoms are compatible with the requested analysis.
Use Qprep's `readx`, `trajectory`, and `writepdb` commands when a restart or frame
must be converted for visualization.

## QligFEP workflow summary

A relative binding free-energy network normally follows these steps:

1. curate and align a congeneric ligand series in SDF format;
2. use `qparams` to write one `.lib`, `.prm`, and `.pdb` per ligand;
3. use `qlomap` or `qkonnektor` to generate a mapping JSON;
4. use `qcog` to obtain a sphere center;
5. prepare the receptor and water sphere with `pdb2amber` and/or `qprep_prot`;
6. collect `protein.pdb`, `water.pdb`, ligand files, SDFs, and mapping JSON in a
   setup directory;
7. use `setupFEP` for every network edge, or `qligfep` for one edge/leg;
8. inspect generated topologies, logs, FEP files, restraints, and job scripts;
9. submit both protein and water legs;
10. analyze with `qligfep_analyze` or `qligfep_neq_analyze`.

`setupFEP` prepares jobs but does not submit them. The normal output layout is:

```text
1.water/
  FEP_ligandA_ligandB/
    inputfiles/
    FEP_submit.sh
2.protein/
  FEP_ligandA_ligandB/
    inputfiles/
    FEP_submit.sh
```

## Reproducibility record

For every calculation retain:

- Q and QligFEP version/commit;
- compiler, compiler flags, MPI implementation, and cluster profile;
- input SDF/PDB files and protonation decisions;
- all `.lib`, `.prm`, `.top`, `.fep`, and `.inp` files;
- mapping JSON and edge direction;
- sphere center/radius and prepared PDB;
- random seeds, temperatures, replicates, windows/switch counts, and softcore
  method;
- complete preparation, simulation, scheduler, and analysis logs;
- failed-run records and all exclusions from the final analysis.

QligFEP writes `inputfiles/fep_config.json` for each edge. Preserve it with the
run, but also archive the actual generated files because they are the executable
protocol.

## Troubleshooting

### Qprep cannot find a residue or atom

Confirm that the correct force-field library was loaded before `readpdb`, names
are uppercase/fixed-column where required, alternate locations were removed, and
each disconnected molecule is separated. Do not rename atoms by guesswork.

### High initial bonded energy or SHAKE failure

Inspect the prepared PDB for missing atoms, duplicated coordinates, bad
cross-links, ligand/protein clashes, and waters overlapping the ligand. Begin
with conservative equilibration or minimization; do not hide a bad topology by
loosening convergence indefinitely.

### Simulation heats rapidly

Check units, timestep, constraints, initial coordinates, thermostat coupling,
and restart compatibility. Find the first abnormal energy component and first
hot atom in the log.

### Sphere-boundary artifacts

Check that the region of interest is inside the unrestrained sphere, boundary
ionizable groups were handled consistently, charge groups are complete, and the
water sphere has no voids. Compare sensitivity to a larger sphere when boundary
placement may affect the result.

### qfep gives unstable or inconsistent estimates

Check completed energy files, ordering of lambda windows, equilibration discard,
neighbor overlap, forward/reverse agreement, replicate spread, and endpoint
energies. Do not report a single pooled estimate without replicate diagnostics.

### A charge-changing edge behaves anomalously

Record the charge-handling method, inspect counter-ion/co-alchemical definitions,
confirm both legs use the intended net charge, and compare forward/reverse edge
behavior. Charge-changing transformations require additional physical and
statistical care.

## Citation

Please cite the methods you use, not only the software repository.

- **Q:** Marelius, Kolmodin, Feierberg, and Åqvist, *Q: A Molecular Dynamics
  Program for Free Energy Calculations and Empirical Valence Bond Simulations in
  Biomolecular Systems*, Journal of Molecular Graphics and Modelling (1998),
  <https://doi.org/10.1016/S1093-3263(98)80006-5>.
- **QligFEP:** Jespers *et al.*, *QligFEP: an automated workflow for small
  molecule free energy calculations in Q*, Journal of Cheminformatics (2019),
  <https://doi.org/10.1186/s13321-019-0348-5>.
- **Current QligFEP protocol:** Araripe *et al.*, *Doing More with Less: Accurate
  and Scalable Ligand Free Energy Calculations by Focusing on the Binding Site*,
  Journal of Chemical Information and Modeling (2026),
  <https://doi.org/10.1021/acs.jcim.5c02932>.
- **QresFEP:** Jespers *et al.*, Journal of Chemical Theory and Computation
  (2019), <https://doi.org/10.1021/acs.jctc.9b00538>.

Consult the repository `LICENSE` file for distribution terms.
