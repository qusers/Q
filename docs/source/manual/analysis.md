# Analysis of results

## Structures from the simulation

`qprep` converts the binary restart and trajectory files of `qdyn` to PDB or
mol2 files for molecular graphics programs.

### PDB or mol2 files from restart or trajectory files

1. Load the topology with `readtop`. The fragment libraries used to make it
   are loaded automatically if they are available; otherwise load them with
   `readlib`.
2. Load coordinates:
   - from a restart file with `readx`, or
   - open a trajectory with `trajectory` and read a frame with `readframe`.
     You are asked whether to use the atom mask of the trajectory, so that
     only atoms in the trajectory appear in the files you write.
3. To select the atoms to write, use `mask`. Enter `mask none` to clear the
   current mask, then add atoms with the
   [atom mask](reference/qprep.md#atom-masks) syntax.
4. Write the file with `writepdb` (with or without `TER` cards) or
   `writemol2`.
5. Repeat steps 2 and 4 for more files. `readnext` reads the next frame of a
   trajectory.

`CONECT` records are written only for fragments whose library entry has
`PDBtype HETATM` in its `[info]` section.

### Trajectory animation

Programs such as VMD [21, 22] and PyMOL read the DCD trajectories of `qdyn`. They need a
PDB file with the same atoms as the trajectory:

1. Load the topology with `readtop`.
2. Open the trajectory with `trajectory` and answer `y` to use its atom mask.
3. To use coordinates other than those of the topology, use `readframe` or
   `readx`.
4. Write the file with `writepdb`, without `TER` cards.

## Free energies with qfep

An FEP calculation is a series of `qdyn` runs in which the weight vector λ
goes from one end state to the other, usually from (1, 0) to (0, 1) for two
states. `qfep` reads the energy files and calculates the free-energy change
for the complete perturbation from state 1, with energy $\varepsilon_1$, to
state 2, with energy $\varepsilon_2$.

### The perturbation formula

The free-energy difference is the sum of the differences between adjacent
windows, each from Zwanzig's formula:

$$
\Delta G = \sum \Delta g = \sum -RT \ln \left\langle e^{-\Delta V_{\mathrm{eff}}/RT} \right\rangle_A
$$ (qfep-zwanzig)

where

$$
V_{\mathrm{eff}} = \lambda_1 \varepsilon_1 + \lambda_2 \varepsilon_2, \qquad \sum_n \lambda_n = 1,
$$

and $\Delta V_{\mathrm{eff}}$ is the difference in $V_{\mathrm{eff}}$ between
two adjacent windows.

`qfep` first lists the average energies and λ values of each file. Then it
summarises the free-energy change of each step, relative to both the previous
and the following window (`dGf` and `dGr`, forward and reverse). It gives the
accumulated sums (`sum(dGf)` and `sum(dGr)`) and the average of the forward
and reverse accumulated change, `<dG>`.

The output is divided into parts:

| Part | Content |
|---|---|
| 0 | Average energies of all states in all files |
| 1 | Free energy perturbation summary (`dGf`, `dGr`, `<dG>`) |
| 2 | Reaction free energy summary |
| 3 | Bin-averaged summary |
| 4 | Thermodynamic integration |
| 5 | Overlap sampling |
| 6 | Bennett acceptance ratio (BAR) |

See [How Q does FEP](../concepts/fep-in-q.md#what-qfep-computes) for the
estimators.

### Free-energy profiles

`qfep` also calculates free-energy functions, or potentials of mean force:

$$
\Delta G(X_m) = \Delta G(\lambda_i) - RT \ln \left\langle e^{-\left(E_g(X_m) - V_i(X_m)\right)/RT} \right\rangle_i
$$ (qfep-pmf)

The reaction coordinate $X$ is the energy gap between the states,
$X = \Delta V = \varepsilon_1 - \varepsilon_2$, divided into intervals $X_m$
(bins). The first term is the free-energy difference between the initial
state $\varepsilon_1$ and the mapping potential $V_i$:

$$
\Delta G(\lambda_i) = -RT \sum_{n=0}^{i-1} \ln \left\langle e^{-(V_{n+1} - V_n)/RT} \right\rangle_n
$$ (qfep-dglambda)

The second term is the free-energy difference between the mapping potential
$V_i$ and the ground-state potential $E_g$. Its average is taken only over
configurations whose $X$ falls in $X_m$. $E_g$ is the solution of the secular
determinant

$$
\begin{vmatrix}
\varepsilon_1 - E_g & \cdots & H_{1n} \\
\vdots & \ddots & \vdots \\
H_{n1} & \cdots & \varepsilon_n - E_g
\end{vmatrix} = 0
$$ (qfep-secular)

which for two states is

$$
E_g = \tfrac{1}{2}(\varepsilon_1 + \varepsilon_2) - \tfrac{1}{2}\sqrt{(\varepsilon_1 - \varepsilon_2)^2 + 4 H_{12}^2} .
$$ (qfep-eg)

$H_{ij}$ is the off-diagonal matrix element that represents the
quantum-mechanical coupling of the states. It is zero in a normal FEP
calculation. $H_{ij} \neq 0$ mixes states $i$ and $j$, which is what is wanted
for reaction free-energy profiles in the EVB model. In `qfep` the
off-diagonal element has the form

$$
H_{ij} = A_{ij} \cdot e^{-\left(\mu (r_{ij} - r_0) + \eta (r_{ij} - r_0)^2\right)}
$$ (qfep-hij)

where $r_{ij}$ is the distance between the reacting atoms. Depending on $\mu$
and $\eta$, $H_{ij}$ is a constant, an exponential or a Gaussian.

EVB allows a simulated reference reaction to be calibrated against
experimental data from the gas phase or solution. The two EVB parameters
$H_{ij}$ (mostly $A_{ij}$) and $\Delta\alpha_{ij}$ are varied until the
calculated profile matches the data. $\Delta\alpha_{ij}$ is a constant energy
shift between the states, which represents their difference in heat of
formation and is not part of the force field. Broadly, $\Delta\alpha_{ij}$
sets the level of $\Delta G^{\circ}$ and $H_{ij}$ sets the degree of mixing at
the transition state, that is, the level of $\Delta G^{\ddagger}$.

The last table of the `qfep` output summarises the FEP functions and the
reaction free-energy profile. Each bin has contributions from several values
of λ, and each value of λ contributes to several bins.

More than two valence-bond states can be handled in `qdyn` and `qfep`, but
sampling and calibration are then difficult.

### Input

`qfep` reads its input from standard input, one item per line:

```bash
qfep < qfep.inp > qfep.out
```

| Line | Example | Meaning |
|---|---|---|
| 1 | `11` | Number of energy files |
| 2 | `2 0` | Number of states; number of predefined off-diagonal elements |
| 3 | `0.596 100` | $kT$; number of points to skip in each energy file. An optional third value of 1 restricts the analysis to the Q–Q energies. |
| 4 | `40` | Number of gap bins |
| 5 | `20` | Minimum number of points per bin |
| 6 | `12.3` | Energy shift $\Delta\alpha$ of state 2. One line for each state after the first. |
| 7 | `1` | Number of off-diagonal elements that are not zero |
| 7.1 | `1 2 18.1 0.32 0.0 2.0` | For each of them: $i$, $j$, $A_{ij}$, $\mu_{ij}$, $\eta_{ij}$, $r_{0,ij}$ |
| 8 | `1 -1` | Linear combination of states that defines the reaction coordinate ($\varepsilon_1 - \varepsilon_2$) |
| 9… | `md_00.en` … | The energy files, one per line, in the order of the λ path |

If the second value of line 2 is not zero, the off-diagonal elements are
taken from the energy files. Lines 7 and 7.1 are then replaced by one line
with a scaling factor for $H_{ij}$.

For a plain two-state FEP calculation the input is therefore:

```text
51
2 0
0.596 100
40
20
0
0
1 -1
md_1000_0000.en
md_0998_0002.en
...
```

### Output

```text
--> Name of file number   1:
...file opened, and lambdas: md_00.ene  0.00  1.00
Number of points: 1000

state    EQtot   EQbond   EQang  EQtor  EQimp     EQel   EQvdW
    1 -1087.71  -578.19    7.68   0.63   0.00  -541.46   22.70
    2 -1145.21  -586.25   31.48   1.16   0.00  -603.21   10.68

state   Eel_qq  EvdW_qq  Eel_qp  EvdW_qp   Eel_qw  EvdW_qw  Eqrstr
    1   112.02     0.57  -18.78    -1.97  -634.70    24.10    0.00
    2     0.00    -0.30   -5.32    -2.51  -597.89    13.50    0.00
...

lambda(1)    dGf  sum(dGf)     dGr  sum(dGr)    <dG>
     0.00   0.00      0.00   -5.71   -284.20    0.00
     0.01   5.64      5.64   -5.54   -278.49    5.67
     0.02   5.63     11.27   -5.45   -272.95   11.26
     0.03   5.48     16.74   -5.28   -267.50   16.72
     0.04   5.38     22.12  -10.41   -262.22   22.05
     0.06  10.37     32.49  -10.07   -251.80   32.44
...

Min energy-gap is: -292.238384870069
Max energy-gap is:  330.166319351831

Lambda1  Ibin  Energy gap     dGa     dGb     dGg  # pts
   0.00     1     -284.46  -31.28  247.20  -47.72     65
   0.00     2     -268.90  -31.28  231.53  -48.50    200
   0.00     3     -253.34  -31.28  216.31  -49.35    231
   0.00     4     -237.78  -31.28  201.25  -50.31     88
   0.01     2     -268.90  -31.40  231.23  -48.73     79
   0.01     3     -253.34  -31.26  216.16  -49.42    303
   0.01     4     -237.78  -31.13  201.59  -50.14    181
   0.01     5     -222.22  -30.99  189.55  -51.00     25
...
```
