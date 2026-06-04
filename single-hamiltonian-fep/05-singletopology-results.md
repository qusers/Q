# Single-topology single-Hamiltonian FEP — it converges

Date: 2026-06-04. TYK2 `ejm_31 -> ejm_47`, protein leg, 15 Å sphere, AMBER14sb.
Builder: `Q-softcore-verify/experiments/single_h/build_singletop.py` +
`fix_boundary_params.py` + `run_singletop.py`. Engine: worktree `src/q6/qdyn`.

## Build

A hand-merge along the rdFMCS atom map produces a 43-atom hybrid: 28 shared core
atoms that morph A->B, 4 vanishing (ejm_31's methyl C20+H30/31/32), 11 appearing
(ejm_47's R-group). ejm_47 types are X-prefixed to avoid the native name clash;
the existing merged prm is reused. `qprep` builds the topology after one
single-topology wrinkle is resolved: bonded terms at the core/R-group junction mix
shared (ejm_31) and appearing (ejm_47) types, which the prm lacks. The 28 such
terms split into shared<->appearing (copied from their all-X equivalent) and 7
vanishing<->appearing cross-terms (zeroed -- they couple the methyl and R-group,
which never coexist). Single topology needs no ligand tether (one molecule), only
a position restraint.

## Result: 21 uniform windows, single_hamiltonian + gapsys

| | dual-topology | single-topology |
|---|---:|---:|
| BAR dG | 30.3 | **28.8** |
| MBAR dG | 30.5 | **28.5** |
| MBAR iterations | 2000 (capped, not converged) | **128 (converged)** |
| windows flagged POOR | 19 of 20 | **0 of 20** |
| overlap gap range | 5–13 kT | **0.1–1.5 kT** |
| free-energy profile | +79 kcal barrier | **monotonic, no barrier** |

Single topology converges: BAR and MBAR agree to 0.3 kcal, every adjacent window
overlaps well, and the profile climbs smoothly from 0 to 28.8. This is the
multistate-quality convergence the whole investigation pointed to, now obtained
*with* a load-bearing soft-core on the appearing R-group -- the combination
dual-topology cannot reach, because it ramps each whole ligand's interactions
(~217 kcal, ~18 kT/window) regardless of path.

## The absolute value differs from dual-topology, by design

Single-topology dG (28.6) is not the dual/multistate 53.7, and should not be:
the two schemes perturb different things. Dual-topology annihilates the entire
ejm_31 and creates the entire ejm_47 -- including destroying and rebuilding the
shared core, the expensive round-trip that is the +79 barrier. Single-topology
morphs only the difference (methyl->R-group plus a ~0.008-per-atom core charge
shift); the core's interactions persist. The endpoint q-energies nearly match
(-275/-250 single vs -280/-252 dual -- the real atoms are identical, only the
never-interacting dummy atoms differ in count), but the work along the path is
far smaller for single topology. Different dummy content gives different absolute
leg dG; this is standard single-vs-dual topology behaviour. Only
ddG = dG_prot - dG_wat is physically comparable, because the convention-dependent
parts cancel between the two legs.

## Next

Run the single-topology water leg at the same protocol, form
ddG = dG_prot - dG_wat, and compare to the dual-topology ddG and to experiment.
The protein-leg crash that blocked dual-topology ddG (dense-solvent clash of a
co-present whole ligand) should not recur here: single topology grows only the
small R-group, with soft-core, into solvent.
