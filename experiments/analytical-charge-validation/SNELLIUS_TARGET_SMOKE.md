# Eg5 and c-Met: bounded cluster compatibility gate

The user authorized staging and submission under
`/projects/prjs2157/astra-charge-change-perturbation` on `mysnellius`.
Existing runs under `/projects/prjs2157/charge-change/runs` remain read-only.
No build or result is placed in the user's home directory. The existing Python
environment there is read, not modified.

This is a high-performance computing (HPC) compatibility test of Q molecular
dynamics (MD), not a binding free-energy experiment. It must not be reported as
showing that the physical correction works on either target.

## Selection and paired controls

Select the first listed charge-changing edge in each forward pilot mapping:

- c-Met: `CHEMBL3402742_23 -> CHEMBL3402744_300`.
- Eg5: `CHEMBL1085666 -> CHEMBL1089056`.

Selection does not use experimental affinity, historical agreement or test
outcome. Use both forward and reverse setup directions, both protein and water
legs, and reference replica 1. This gives eight input cases. Within each case,
run an integrated-Born and a post-hoc-Born control from identical restarts at
lambda 0.5. Forward/reverse archived restarts are not claimed to be matched
equilibrium samples or a thermodynamic closure test.

Retain original topology, free-energy perturbation (FEP) definition, and
nonempty sequence/distance/wall restraints. These archived FEP files use dual
topologies and Gapsys softcore. They are not the no-softcore ghost-endpoint case
reported by the user. No exact endpoint is sampled or trimmed in this gate.
Do not use these midpoint records to claim a free-energy estimate or to justify
post-hoc lambda scaling across softcore-dependent Hamiltonians.

## Explicit adaptations

Each run uses twenty 1-femtosecond steps at 298 kelvin with repaired SHAKE,
the original thermostat settings, retained restart velocities, and pair-list
updates every step. Disable the local reaction field (LRF), set all pair
cutoffs to 99 angstrom, and check actual whole-trajectory distance coverage.
Retain existing radial/angular force constants and the original sphere/solute
restraints. Record Q's actual solvent radius, not the nominal 20-angstrom label.

Per-state polarization is enabled and its adaptation disabled. In a new restart
copy only, replace the offset record with zeros, preserving coordinate and
velocity records byte for byte. Original offsets and all source hashes are
retained. This deliberately defines a new fixed-offset Hamiltonian
(potential-energy model), not an equilibrated continuation of the old one.
The angular target remains Q-region-only; no total-enclosed-charge angular
proposal is implemented here.

Integrated Born is the only difference between paired runs. Check native
included/excluded non-Q solute charges separately, pure-state charges, topology
Coulomb constant, actual radius, dielectric 80, and
`B_s = -k_e (1-1/80) (Q_included+q_s)^2/(2 R)`.
Require identical final restart bytes between the pair and the corresponding
constant difference in every saved pure-state total. All other saved fields
must match within the existing accounting tolerance. This is not a comparison
between separate equilibrium distributions.

## Budget, provenance and failure handling

The single Slurm scheduler job requests one node, four CPU cores, 4 gigabytes
of memory, and twenty minutes on `rome`, using the verified user association
`ugsei19097`. This is a resource request, not a measured runtime forecast.
The job does not install or modify software environments.

It first builds serial Qdyn/Qprep from a transferred Git archive using the
cluster compiler, runs existing LINCS (Linear Constraint Solver) and SETTLE
(rigid-water solver) tests, then four fresh 10/14-angstrom charge-only probes
with both charge signs. The eight real-target cases run with at most four
concurrent native processes. Each native target invocation has a 120-second
limit; there is no automatic retry.

Total maximum MD budget: `4*20 + 8*2*20 = 400` steps, or 0.4 picoseconds
aggregate. No charge ladder, long equilibration, or production experiment is
submitted by this package. Failures are preserved and block progression, not
hidden by dropping a target. Other cases still report their own outcomes.

The release includes exact Git source/native archives and SHA-256
(cryptographic hash) fingerprints. The job verifies the unpacked sources,
records compiler/build/binary identity, and retains unmodified references,
adapted inputs, initialization audits, all saved energies, final restarts and
per-case outcome reports. Existing releases/attempts cannot be overwritten.
The build archive is validated without needing a remote Git checkout.

The cluster Python environment lacks the pytest test runner, so this job
does not claim to repeat the entire local Python regression suite. It executes
the audited smoke driver plus the native solver tests.

## Next decision

A passing result permits preparing a real-target sampling pilot; it does not
qualify the correction on oppositely charged proteins. Before that pilot,
define interior-endpoint handling for the actual selected FEP dialect,
matched replica/direction controls, adequate equilibration and uncertainty,
and charge-background/radius comparisons. Do not tune boundary parameters to
make one target agree. A result that works for only one target is not a pass.
