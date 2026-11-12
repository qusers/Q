# How Q does FEP

Q describes an alchemical transformation with **states**, not with a single
set of λ-dependent parameters. This page defines the terms, shows what `qdyn`
simulates and stores, and shows how `qfep` turns the stored numbers into a free
energy. The [soft-core page](softcore.md) builds on it.

## Q-atoms and states

The FEP file designates some atoms of the topology as **Q-atoms**. For each
Q-atom it can give a different charge, atom type (Lennard-Jones parameters) and
bonded terms in each of $n$ **states**. Everything else in the system has one
description, taken from the topology.

For a configuration $x$, `qdyn` evaluates one **state energy** per state:

$$
\varepsilon_s(x) = E^{\mathrm{bonded}}_{s} + E^{qq}_{s} + E^{qp}_{s} + E^{qw}_{s} + E^{\mathrm{restr}}_{s},
\qquad s = 1,\dots,n
$$

where $qq$, $qp$ and $qw$ are the nonbonded (electrostatic plus Lennard-Jones)
energies of Q-atoms with Q-atoms, with the rest of the solute, and with the
solvent. Each $\varepsilon_s$ is computed with the *fixed* parameters of state
$s$. It contains no λ.

## The mapping potential

The simulation is run on the **mapping potential**, a linear mixture of the
state energies, with the weights given in the `[lambdas]` section of the `qdyn`
input:

$$
V_m(x) = \sum_{s} \lambda^{(m)}_{s}\,\varepsilon_s(x) + V_{\mathrm{env}}(x),
\qquad \sum_s \lambda^{(m)}_{s} = 1 .
$$

$V_{\mathrm{env}}$ collects every interaction that does not involve a Q-atom
and is the same in all states. The index $m$ labels a **window**: one
simulation at one fixed weight vector. A two-state calculation runs a series
of windows from $(1, 0)$ to $(0, 1)$.

The forces follow the same rule: the force on an atom is the λ-weighted sum of
the forces from each state.

## What is written to the energy file

At every energy-output step `qdyn` writes the state energies
$\varepsilon_1(x), \dots, \varepsilon_n(x)$ and the weights of the window to the
binary `.en` file. It does not write $V_m$ for other windows.

That is sufficient, because $V_m$ is linear in the weights. For any other
window $m'$,

$$
V_{m'}(x) - V_m(x) = \sum_s \bigl(\lambda^{(m')}_{s} - \lambda^{(m)}_{s}\bigr)\,\varepsilon_s(x),
$$

so the energy of a stored configuration under a neighbouring window can be
reconstructed afterwards without re-running anything. `qfep` relies on this.

## Dynamics and bookkeeping

It helps to keep two things apart.

- The **dynamics**: which potential moves the atoms. In Q this is the
  mapping potential, a weighted sum of state energies.
- The **bookkeeping**: which numbers are stored, and how the analysis turns
  them into the energy of a configuration under another window. In Q the
  state energies are stored and the other window is reconstructed by
  changing the weights.

The bookkeeping is the stricter of the two. The dynamics only needs each
state to have an energy and a force. The bookkeeping needs each state energy
to be the same function of the coordinates in every window, so that a value
stored in one window can be reused for another. Anything that makes a state
energy depend on the weights is fine for the dynamics and wrong for the
bookkeeping.

<!-- FIGURE TODO: dynamics-vs-bookkeeping
What: a two-row flow diagram for one stored configuration x of window m.
  Row 1 (dynamics, inside qdyn): x -> eps_1(x), eps_2(x) -> V_m = l1*eps_1 +
  l2*eps_2 -> forces -> next x.
  Row 2 (bookkeeping): eps_1(x), eps_2(x) written to the .en file -> qfep
  re-weights with the lambdas of window m+1 -> V_(m+1)(x) -> exp average / BAR.
  Mark with a warning symbol the arrow "re-weights", annotated "valid only if
  eps_s does not depend on lambda".
Why: this distinction is the key to understanding both the energy files and
  why a lambda-dependent soft-core cannot simply be switched on.
How: a mermaid flowchart in the page, or a drawn schematic.
-->

## What qfep computes

For neighbouring windows $m$ and $m+1$, `qfep` forms the energy difference above
for every stored configuration and applies the usual estimators:

- **Zwanzig (exponential averaging)**, forward and reverse,

  $$
  \Delta G_{m \to m+1} = -RT \ln \left\langle e^{-(V_{m+1} - V_m)/RT} \right\rangle_m ;
  $$

- **thermodynamic integration**, using $\partial V/\partial\lambda$, which for
  two states is simply $\varepsilon_2 - \varepsilon_1$;
- **overlap sampling** and **Bennett acceptance ratio (BAR)**, which combine
  the forward difference sampled in window $m$ with the reverse difference
  sampled in window $m+1$.

The total free energy is the sum over neighbouring pairs. Every one of these
estimators uses the linear reconstruction, so every one of them assumes that
$\varepsilon_s(x)$ does not depend on the window it was recorded in.

## How this compares with other codes

Most other packages present alchemical FEP as one Hamiltonian $H(\lambda)$
whose parameters depend on a scalar λ. For interactions that are linear in the
perturbed parameters the two pictures are the same thing. Take an atom that is
present in state 1 and absent in state 2, and write $\lambda_1 = 1 - \lambda$:

$$
\lambda_1\,V_{\mathrm{LJ}}(r) + \lambda_2 \cdot 0 = (1-\lambda)\,V_{\mathrm{LJ}}(r).
$$

The practical differences are:

| | Q | Typical single-λ code |
|---|---|---|
| Path variable | a weight per state, summing to 1 | one scalar λ (sometimes several) |
| What is stored | state energies $\varepsilon_s$ | energy differences to other λ values, or $\partial H/\partial\lambda$ |
| Cross-window energies | reconstructed in analysis from linearity | evaluated during the run |
| Requirement | $V_m$ must be linear in the weights | none |

The last row is the important one. Any term that makes a state energy depend
on λ, such as a λ-dependent soft-core, breaks the reconstruction; see
[Soft-core](softcore.md).

The state picture is not an FEP-specific trick. It comes from the empirical
valence bond (EVB) method, where the states are valence-bond structures of
reactants and products, and it is what lets Q handle more than two states and
off-diagonal coupling between them.

## Dual topology in QligFEP

QligFEP builds a relative binding free-energy calculation as a **dual
topology**: both ligands are in the system at once, on top of each other.

- Atoms of ligand A have their real atom types and charges in state 1. In
  state 2 they have the atom type `DUM` and charge 0.
- Atoms of ligand B are the mirror image: `DUM` in state 1, real in state 2.
- `DUM` has zero Lennard-Jones parameters, so a dummy atom has no nonbonded
  interaction with anything. It keeps its mass and its bonded terms.
- Distance restraints between matched atoms keep the two ligands superimposed.

```{figure} ../_static/figures/dual-topology-states-light.svg
:alt: Table of two ligands by two states. Ligand A is drawn in colour in state 1 and as a grey shade in state 2. Ligand B is drawn as a grey shade in state 1 and in colour in state 2.
:width: 85%

The atoms of both ligands are in the topology and are simulated in every
window. Each ligand has two descriptions, one per state. Colour: the atoms
have their real atom types and charges. Grey: the atoms have type `DUM` and
charge 0, and keep their mass and bonded terms. State 1 is ligand A with a
non-interacting copy of ligand B; state 2 is ligand B with a non-interacting
copy of ligand A. The ligands are drawn apart here; in the simulation they are
superimposed. Generated by `docs/scripts/plot_dual_topology_states.py`.
```

At $(1, 0)$ the system is ligand A in its environment plus a non-interacting
copy of ligand B; at $(0, 1)$ it is the reverse. In between, each ligand
interacts with its surroundings with its real potential scaled by its state
weight. The two ligands never interact with each other through nonbonded
terms, because in every state one of them is `DUM`.

The same calculation is run for the ligand pair in the protein and in water,
and

$$
\Delta\Delta G_{\mathrm{bind}} = \Delta G_{\mathrm{protein}} - \Delta G_{\mathrm{water}} .
$$
