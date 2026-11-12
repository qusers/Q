# FEP and analysis

## Equilibrium FEP

For states with potential energies $U_i(x)$, Qdyn samples a lambda-weighted
state mixture and records the state energies required for free-energy analysis.
A two-state calculation traverses values such as

```text
1.000 0.000
0.950 0.050
...
0.000 1.000
```

in separate windows. QligFEP can use linear, sigmoidal, exponential, or reverse
exponential spacing.

### Run qfep

Qfep reads its control file from standard input:

```bash
src/q6/bin/q6/qfep < qfep.inp > qfep.out
```

The control input supplies:

1. number of energy files and FEP states;
2. thermal energy `kT` and samples to skip;
3. energy-gap binning controls;
4. state energy shifts;
5. optional off-diagonal EVB parameters;
6. reaction-coordinate coefficients;
7. the ordered energy-file list.

Use a generated QligFEP `qfep.inp` as a known-current example. File order must
follow the lambda path used in the simulations.

### Estimators

Qfep reports:

- forward and reverse Zwanzig exponential averages;
- thermodynamic integration (TI);
- overlap sampling (Lu–Kofke–Woolf);
- Bennett acceptance ratio (BAR);
- energy-gap free-energy profiles for EVB/state-mixing calculations.

Do not choose an estimator only because it gives the expected answer. Report the
specified primary method and inspect agreement, overlap, hysteresis, and
replicate dispersion.

## Relative binding free energy

For the same alchemical edge in the bound and unbound systems:

$$
\Delta\Delta G = \Delta G_\mathrm{protein} - \Delta G_\mathrm{water}.
$$

Edge direction matters. QligFEP mapping edges use `from` → `to`, and experimental
`ddg_value` follows $G_\mathrm{to} - G_\mathrm{from}$.

## Non-equilibrium FEP

Non-equilibrium calculations begin from equilibrium endpoint configurations and
drive lambda through many independent forward and reverse switches. The work
distributions are combined by BAR. QligFEP reports per-leg free energies,
protein-minus-water `ddG_kcal`, uncertainty, sample counts, and overlap.

Low overlap, few successful switches, or one-direction failures are convergence
problems, not merely plotting issues.

## Qcalc

Qcalc reads a command sequence from standard input:

```bash
src/q6/bin/q6/qcalc < qcalc.inp > qcalc.out
```

It supports:

- RMSD and RMSF;
- least-squares coordinate fitting;
- distances, angles, and torsions;
- Schlitter entropy;
- radial distribution functions;
- center-of-mass position and kinetic energy;
- nonbonded energies and residue/Q-atom interaction decompositions;
- X-Score, ChemScore, and PMF-Score.

Qcalc is interactive when no redirected control file is supplied. Capture the
interactive prompts once, then archive a scripted input for repeatable analysis.

## Minimum quality-control record

For every reported edge retain:

- all per-replicate estimates;
- number and order of lambda windows, or forward/reverse switch counts;
- discarded equilibration samples;
- failed/missing windows or switches;
- forward/reverse agreement;
- overlap diagnostics;
- water- and protein-leg values before subtraction;
- mapping direction and experimental-value convention;
- exact analysis command and software version.
