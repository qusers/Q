# Approved goal — persistent tracker activation pending

The user approved this goal and its scope limits. The previous goal remains
paused in the thread's tracker, which refuses to register a replacement while
an unfinished goal exists. Do not mark the previous goal achieved to bypass
this restriction. Register this approved objective once the previous tracker
entry is cleared; approval is not a claim that implementation is complete.

Implement and validate analytical corrections for charged perturbations in Q's
existing molecular dynamics (MD) workflow with finite spherical boundaries,
building on `charged-boundary-hamiltonian`. Establish whether the correction and
its energy bookkeeping explain the remaining positive/negative charge asymmetry,
make only evidence-supported minimal changes, and prepare a reproducible
high-performance computing (HPC) validation protocol.

## Scientific questions

1. Check the physical assumptions of the correction against the Surface Constraint
   All-Atom Solvent (SCAAS) model and finite spherical boundary. Treat the supplied
   slide equations as hypotheses, not ground truth.
2. Audit the simulated Hamiltonian (the energy function), forces, pure-state saved
   energies and their lambda weighting (the interpolation between charge states).
   Identify missing terms, double counting, inconsistent constants or state
   definitions. Do not assume every remaining error is a bookkeeping problem.
3. Investigate positive/negative asymmetry using existing results and matched
   charge-perturbation checks. Do not fit charged, protein, ligand, binding or
   experimental results to manufacture agreement.
4. Keep exact ghost-endpoint overlap separate. Historical interior-only Bennett
   acceptance ratio (BAR) analysis estimates a truncated free-energy interval,
   not the full transformation; missing endpoint contributions need not cancel.
   Use charge-only tests with identical real Lennard-Jones (LJ) interactions in
   both states to isolate this issue from new correction validation.

## Scope limits and completion evidence

Use the existing Q dynamics. No Monte Carlo sampler, thermostat replacement,
event-handling integrator, solvent-model redesign or fixed-radial-wall change is
included. If evidence requires a change beyond bookkeeping or the analytical
correction, explain the need and request separate approval before implementation.
Do not resume large neutral-shell calibration to answer a bookkeeping question.

The intended deliverables are an evidence-backed correction/limitations account,
minimal code changes with relevant native and analysis tests, and a reproducible
HPC experiment package: both charge signs, appropriate controls, consistent
boundary/state definitions, endpoint handling, source/input provenance, analysis,
uncertainty/convergence/failure criteria and a justified compute budget. Passing
software tests alone does not validate the physical correction.

Request approval before substantial HPC computation or submission. Introduce
new technical terms and expand acronyms on first use in each documentation file.
The scope above is approved. Any expansion still requires separate approval.

## Clean branch and carried work

Branch: `fix/analytical-charge-corrections`.
Base: `charged-boundary-hamiltonian`, commit `55aa555c`.

Selected changes are copied from the state at `b0e5544d`, before any sampler
implementation, and will have a new scoped commit rather than importing the
exploratory branch history:

- Provenance-checked interior-only BAR utility, documentation and tests.
- Historical charged-perturbation reproduction scripts, frozen records and tests.
- Topology-consistent Coulomb constant in the Born correction and its native tests.

No sampler commits, sampler code, crossing-propagation machinery, fixed-radial
boundary changes or obsolete sampler readiness plans are carried. The previous
worktree, including its uncommitted changes, is preserved separately. Historical
raw data/build artifacts are not part of this initial tracked-file transfer.

Transfer validation: serial Qdyn and Qprep build successfully with gfortran-11;
50 focused tests pass in 6.64 seconds. Two optional integration tests are skipped
because the QFEP native endpoint-test binary and historical reproduction runtime
assets have not been installed in this clean worktree. This checks the transfer,
not completion of the proposed scientific goal or a fresh reproduction of the
archived results.
