# Command-line reference

Installing the `QligFEP` Python package provides fourteen commands, for both
the QligFEP and the QresFEP workflow.

## Preparation

| Command | Purpose |
|---|---|
| [`qcog`](preparation.md#qcog) | Calculate a ligand/PDB center of geometry |
| [`pdb2amber`](preparation.md#pdb2amber) | Normalize PDB residue and atom names for AMBER |
| [`qparams`](preparation.md#qparams) | Generate Q ligand/cofactor parameters with OpenFF |
| [`qprep_prot`](#qprep-prot-cli) | Prepare and solvate a spherical protein system |

## Network generation and visualization

| Command | Purpose |
|---|---|
| [`qlomap`](networks.md#qlomap) | Generate a LoMap perturbation network |
| [`qkonnektor`](networks.md#qkonnektor) | Generate a configurable Konnektor network |
| [`qmapfep`](networks.md#qmapfep) | Generate or render a mapping for browser editing |

## Ligand FEP setup

| Command | Purpose |
|---|---|
| [`qligfep`](setup.md#qligfep) | Prepare one directed edge and one thermodynamic leg |
| [`setupFEP`](setup.md#setupfep) | Prepare all mapped edges for water/protein legs |

## Ligand FEP analysis

| Command | Purpose |
|---|---|
| [`qligfep_analyze`](#qligfep-analyze-cli) | Analyze equilibrium-window simulations |
| [`qligfep_neq_analyze`](#qligfep-neq-analyze-cli) | Analyze non-equilibrium switching work |

## Residue FEP

| Command | Purpose |
|---|---|
| [`qresfep`](#qresfep-cli) | Prepare one leg of one amino-acid mutation |
| [`setup_resFEP`](#setup-resfep-cli) | Prepare both legs for a list of mutations |
| [`qresfep_analyze`](#qresfep-analyze-cli) | Combine both legs into folding free-energy changes |

Every command accepts `-h`/`--help`. Most commands use `-log`/`--log-level` with
`trace`, `debug`, `info`, `warning`, `error`, or `critical`; `qprep_prot` omits
`trace`.
