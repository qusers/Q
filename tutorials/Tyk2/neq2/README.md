# Non-equilibrium FEP with Tyk2

This tutorial walks through relative binding free energy (RBFE) calculation on the Tyk2 benchmark set using the **non-equilibrium (NEQ)** protocol in Q ($NEQ^{2}$), from edge preparation through submission to analysis. It reuses the same inputs as the [equilibrium Tyk2 tutorial](../README.md); only the `setupFEP` and analysis steps differ.

## What NEQ FEP is

The equilibrium protocol samples ~100 fixed $\lambda$ windows per edge and combines them sequentially with `qfep`. NEQ instead drives $\lambda$ continuously from one end state to the other over a single short simulation (a *switch*) with the `qdyn_neq` binary, accumulating the work done on the system. Running many forward ($\lambda: 1 \to 0$) and reverse ($\lambda: 0 \to 1$) switches yields two work distributions, and the free energy is recovered from the Bennett Acceptance Ratio (BAR) over them. The relative binding free energy is $\Delta\Delta G = \Delta F_{\text{protein}} - \Delta F_{\text{water}}$, exactly as in the equilibrium case.

Each leg is written $\Delta F$, a *Helmholtz* free energy, rather than $\Delta G$ because the work theorems BAR is built on ([Jarzynski](https://doi.org/10.1103/PhysRevLett.78.2690), [Crooks](https://doi.org/10.1103/PhysRevE.60.2721)) are derived at constant volume, where $F = U - TS$ (not the Gibbs $G = F + pV$) is the natural free energy. Q simulates inside a fixed-radius solvent sphere with no pressure coupling, so the volume is effectively constant and $F$ and $G$ differ only by a negligible $p\Delta V$ term; the combined binding estimate is therefore reported as the familiar $\Delta\Delta G$.

The switch advances $\lambda$ on a sigmoidal schedule of the progress fraction $p$ (the share of switch steps completed). The raw sigmoid $s(p) = 1 / \left[1 + e^{-(2p-1)L}\right]$ is rescaled to hit the endpoints exactly, $S(p,L) = \left[s(p) - s(0)\right] / \left[s(1) - s(0)\right]$, and $\lambda$ is the linear interpolation between the starting endpoint and its complement by $S(p,L)$. A larger steepness $L$ spends more time near the end states. Because each edge needs only a handful of short switches rather than ~100 fully-equilibrated windows, NEQ reaches comparable *ranking* power at roughly 20-40% less compute.

## Prerequisites

- A working QligFEP install with the Q binaries built (`make all` in `src/q6` builds `qdyn_neq` alongside the others).
- The analysis CLI `qligfep_neq_analyze` (installed with the package).

## 1. Preparation (reuse the Tyk2 inputs)

NEQ uses the **exact same preparation** as the equilibrium workflow: ligand parameters, the perturbation network, and the water sphere. Rather than repeat it here, follow the main tutorial through the steps below, which start from the inputs already in this repository:

| Input                | File | Step in the main tutorial |
|----------------------|------|---------------------------|
| Ligands              | [`../ligprep/tyk2_ligands.sdf`](../ligprep/tyk2_ligands.sdf) | [Ligand parameters](../README.md#ligand-parameters) |
| Perturbation network | `lomap.json` (generated) | [Visualize the perturbation network](../README.md#visualize-the-perturbation-network) |
| Protein / water      | [`../setupFEP/amber/protein.pdb`](../setupFEP/amber/protein.pdb), [`../setupFEP/water.pdb`](../setupFEP/water.pdb) | [Water sphere](../README.md#water-sphere) |

Work in a fresh directory so the NEQ run does not collide with the equilibrium `1.water`/`2.protein` folders. The equilibrium tutorial already gathered every input `setupFEP` needs into `../setupFEP`, so pull them in from there (ligand `.pdb`/`.prm`/`.lib`, the `.sdf` files, `lomap.json`, and the prepared `protein.pdb`/`water.pdb`):

```bash
cd tutorials/Tyk2/neq2
cp ../setupFEP/*.pdb ../setupFEP/*.prm ../setupFEP/*.lib ../setupFEP/*.sdf ../setupFEP/lomap.json .
```

Stop once you have `lomap.json` and the prepared structures in place, then continue here.

## 2. Setup NEQ FEP

From the directory holding `lomap.json`, run `setupFEP` with `--neq`:

```bash
setupFEP -FF AMBER14sb -r 25 -ts 2fs -j lomap.json -rs 42 -c SNELLIUS \
    --neq --neq-reps 5 --neq-steps 50000 --neq-eq-steps 1000 --neq-relax-steps 5000 \
    -L 8 --neq-schedule sigmoidal
```

The shared flags (`-FF`, `-r`, `-ts`, `-j`, `-rs`, `-c`) behave exactly as in the
[equilibrium setup](../README.md#setup-fep). The NEQ-specific flags are:

- `--neq`: switch QligFEP into non-equilibrium mode. The windowed flags `-w/--windows` and `-S/--sampling` are ignored in this mode;
- `--neq-reps 5`: number of forward/reverse switching pairs per replicate;
- `--neq-steps 50000`: length of each switch in MD steps (100 ps at 2 fs; recommended > 16000);
- `--neq-eq-steps 1000`: endpoint equilibration steps between successive switches (the tEQ spacing; recommended > 250);
- `--neq-relax-steps 5000`: length of the one-time endpoint relaxation at $\lambda = 0$ and $\lambda = 1$ before the first switch (~10 ps at 2 fs), which settles the nearly-decoupled ligand at each endpoint. The first switching iteration uses this longer relaxation; later iterations use the shorter `--neq-eq-steps` spacing;
- `-L 8` (`--neq-steepness`): steepness of the sigmoidal schedule;
- `--neq-schedule sigmoidal`: switching schedule (`sigmoidal` or `linear`).

### Recommended defaults: the paper's production protocol

The flag values above are also the CLI **defaults** — they are spelled out only for clarity. They reproduce the final $NEQ^{2}$ manuscript protocol, benchmarked across the full ~200-edge JACS FEP set, so `setupFEP --neq` with no overrides already gives the recommended protocol.

| Flag | Default | Why this value |
|------|---------|----------------|
| `--neq-steps 50000` | 100 ps @ 2 fs | Switch length (tNEQ). Longer switches lower the irreversible work and increase forward/reverse overlap, which BAR needs. The grid search started at 32 ps, but 100 ps was required to recover equilibrium-level ranking (Kendall $\tau$ = 0.46 vs QligFEP 0.47). |
| `-L 8` | 8 | Sigmoidal steepness. L=8 balances the work-accumulation gradient between the endpoints (too shallow) and the midpoint (too steep), giving the best forward/reverse overlap. |
| `--neq-schedule sigmoidal` | sigmoidal | Spends more sampling near the end states, where forces change fastest. |
| `--neq-eq-steps 1000` | 1000 steps | Decorrelation spacing (tEQ) between successive switches; shorter intervals gave irreproducible correlations. |
| `--neq-reps 5` | 5 | Iterations (N), set to cost-match equilibrium (see below). |
| `--replicates 10` | 10 | Independent seeds for the uncertainty (SEM), as in the QligFEP benchmark. |

`--neq-relax-steps 5000` is the one exception — a one-time longer endpoint relaxation before the first switch, an implementation detail not part of the manuscript's parameter set.

**Cost-matching N to equilibrium.** N is chosen so the switching work roughly equals the equilibrium production cost: 5 iterations × 2 switches × 50,000 steps = **500,000 steps**, against ~505,000 for the ~100 equilibrium windows. If you change `--neq-steps`, rescale `--neq-reps` to preserve the match (halving the switch length → roughly double the reps). Cutting compute below the cost-matched point is a deliberate accuracy/cost trade.

This creates `1.water` and `2.protein`, each with one `FEP_<lig1>_<lig2>` folder per edge. Every `inputfiles/` directory now holds the equilibration files (`eq1`-`eq5`), the one-time endpoint-relaxation templates (`relax_0`, `relax_1`), the endpoint-spacing templates (`eq6_0`, `eq6_1`), and the switching templates (`neq_0`, `neq_1`) — instead of the ~100 `md_xxxx_xxxx.inp` window files of the equilibrium run.

The dual-topology distance restraints are applied exactly as in the equilibrium protocol (see [Restraint setting](../README.md#restraint-setting)): 1.5 kcal/mol/Å² during eq1-eq4 and 0.5 kcal/mol/Å² for eq5 and every NEQ step. Confirm the switching schedule made it into the inputs:

```bash
tail -n 4 2.protein/FEP_ejm_31_ejm_42/inputfiles/neq_0.inp
```

```text
[lambda_scaling]
scaling_parameter          sigmoidal
L_sigmoid        8.0
```

## 3. Job submission

Submission is identical to the equilibrium workflow — each `FEP_<lig1>_<lig2>` directory has a `FEP_submit.sh` script, submitted per leg with the [`submitFEPjobs`](../README.md#job-submission) helper:

```bash
cd 2.protein
submitFEPjobs
# once the protein leg is running cleanly, do the water leg
cd ../1.water
submitFEPjobs
```

Each replicate (a SLURM array task) runs `eq1`-`eq5`, performs the one-time endpoint relaxation, then loops the requested forward/reverse switches. The switching work is written to **`neq_1_*.log`** (forward, $\lambda: 1 \to 0$) and **`neq_0_*.log`** (reverse, $\lambda: 0 \to 1$) under each replicate's run directory.

**How the switches map to cores.** Each iteration fires a forward *and* a reverse switch, so `--neq-reps 5` gives `2 × 5 = 10` independent switches per replicate. Because the switches are serial and independent, `setupFEP` sizes the allocation to run them all in **one concurrent wave**, setting `--ntasks-per-node` to `max(2 × neq_reps, cluster default)`. The cluster default (16 on SNELLIUS) is the *billed* core count from [`settings.py`](../../../src/QligFEP/settings/settings.py): flooring there means a cluster that bills a fixed minimum never runs fewer cores than you pay for, and the equilibration (`eq1`-`eq5`, on the parallel `qdynp` engine) uses all of them. The count grows to `2 × neq_reps` only when the protocol needs more cores than the default (e.g. 10 switches where the default is 8), keeping the switches in one wave rather than a slow tail. A single `mpirun` then streams the switches onto the cores with no cross-rank communication — the independence behind NEQ's ~20–40% compute saving.

**Adjusting for your cluster.** If a node cannot grant `2 × neq_reps` cores, lower `--ntasks-per-node` in the generated `run*.sh`; the switches then run in successive waves — correct, just slower. Use a **divisor of `2 × neq_reps`** to keep the waves even; the script prints a warning before launching if your count leaves a partial tail wave with idle billed cores. For a new cluster, set its `NTASKS` (the billed floor) in [`settings.py`](../../../src/QligFEP/settings/settings.py) so `setupFEP` sizes the allocation automatically.

❗**Tip**❗ As with equilibrium FEP, protein legs are more failure-prone than water legs (residual clashes, a nearly-decoupled ligand at switch start). We recommend first submitting the protein leg.

## 4. Analysis

Once the switches finish, estimate the free energies with `qligfep_neq_analyze`. It reads the work from the switching logs, runs BAR with a bootstrap uncertainty, combines the legs into $\Delta\Delta G = \Delta F_{\text{protein}} - \Delta F_{\text{water}}$, and — when given the mapping JSON — compares to experiment and saves the correlation plot.

```bash
qligfep_neq_analyze -p 2.protein -w 1.water -T 300 -u kcal \
    -j lomap.json -exp ddg_value -t Tyk2 -o neq_results.csv
```

Options:

- `-p 2.protein` / `-w 1.water`: the protein- and water-leg directories holding the
  `FEP_*` edges (same flags as `qligfep_analyze`);
- `-T 300`: temperature (K) used for the kcal/mol conversion (and for $\beta$ when not `kT`);
- `-u kcal`: work units for BAR — the default ($\beta = 1/k_B T$, the physically consistent factor); `kT` uses $\beta = 1$;
- `-j lomap.json` / `-exp ddg_value`: the mapping JSON and the edge key holding the
  experimental ddG. With both, the analyzer attaches the experimental values, logs the
  correlation metrics, and writes `Tyk2_neq_ddG_plot.png`. Omit them to skip the comparison;
- `-t Tyk2`: target name used in the plot title and output file names;
- `-o neq_results.csv`: the per-edge results table.

### Outputs

**`neq_results.csv`** — one row per edge, with the estimate **and** run diagnostics:

| Column                               | Meaning                                                          |
|--------------------------------------|------------------------------------------------------------------|
| `ddG_kcal`, `ddG_err_kcal`           | BAR $\Delta\Delta G$ and its bootstrap standard error            |
| `dF_protein_kcal`, `dF_water_kcal`   | per-leg BAR free energies                                        |
| `overlap_protein`, `overlap_water`   | forward/reverse work-distribution overlap (0-1)                  |
| `n_forward_*`, `n_reverse_*`         | completed switches per leg                                       |
| `n_failed_protein`, `n_failed_water` | switches that started but did not finish (SHAKE/overlap crashes) |
| `ddG_exp`                            | experimental value (only with `-j`/`-exp`)                       |

**`neq_results_run_data.csv`** — per-replicate `runtime`, `seed`, and `status`
(`SUCCESS`/`TIMEOUT`/`OOM`/`CANCELLED`/`CRASHED`) parsed from the `slurm*.out` files. Pass
`-norun` to skip this if the slurm logs are not present.

**`Tyk2_neq_ddG_plot.png`** — experimental vs calculated $\Delta\Delta G$, with the correlation statistics
(Kendall $\tau$, RMSE, MUE) when experimental data is supplied. The correlation metrics
(R², Pearson, Spearman, Kendall, RMSE, MAE) are also logged to the console.

## Interpreting the results

- **Overlap and failed switches are the reliability check.** A low `overlap_*` (well under
  ~0.3) or a non-zero `n_failed_*` means the forward and reverse work distributions barely
  meet, so the BAR estimate for that leg is poorly determined. These flag the
  switch-initiation instability — SHAKE errors and high-energy nearly-decoupled-ligand
  configurations — that the relaxation step is meant to reduce. Lengthen `--neq-steps` if
  forward/reverse overlap is poor.
