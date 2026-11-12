# Command-line reference

Installing QligFEP declares eleven console scripts.

## Preparation

| Command | Purpose |
|---|---|
| [`qcog`](preparation.md#qcog) | Calculate a ligand/PDB center of geometry |
| [`pdb2amber`](preparation.md#pdb2amber) | Normalize PDB residue and atom names for AMBER |
| [`qparams`](preparation.md#qparams) | Generate Q ligand/cofactor parameters with OpenFF |
| {ref}`qprep_prot <qprep-prot-cli>` | Prepare and solvate a spherical protein system |

## Network generation and visualization

| Command | Purpose |
|---|---|
| [`qlomap`](networks.md#qlomap) | Generate a LoMap perturbation network |
| [`qkonnektor`](networks.md#qkonnektor) | Generate a configurable Konnektor network |
| [`qmapfep`](networks.md#qmapfep) | Generate or render a mapping for browser editing |

## FEP setup

| Command | Purpose |
|---|---|
| [`qligfep`](setup.md#qligfep) | Prepare one directed edge and one thermodynamic leg |
| [`setupFEP`](setup.md#setupfep) | Prepare all mapped edges for water/protein legs |

## Analysis

| Command | Purpose |
|---|---|
| {ref}`qligfep_analyze <qligfep-analyze-cli>` | Analyze equilibrium-window simulations |
| {ref}`qligfep_neq_analyze <qligfep-neq-analyze-cli>` | Analyze non-equilibrium switching work |

Every command accepts `-h`/`--help`. Most commands use `-log`/`--log-level` with
`trace`, `debug`, `info`, `warning`, `error`, or `critical`; `qprep_prot` omits
`trace`. The reference records implementation behavior even where an argument's
current help text is less strict.
