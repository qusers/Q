# Soft-core in qgpu: what is left

Branch `feature/qgpu_support_softcore`. Follows the Fortran implementation in
[qusers/Q#140](https://github.com/qusers/Q/pull/140), which is still open: its
input keys and the `.fhe` format may change before it is merged.

## Done

- Gapsys soft-core for Q atoms (`softcore_method gapsys` in the FEP file): GPU
  pair kernel, LRF exact tiles, CPU pair routine, 1-4 pairs.
- Foreign work `U(lambda') - U(lambda)` for the windows in `[foreign_lambdas]`,
  written to `<energy file>.fhe` and returned in `RunResult`.
- `hrex.py` and `fep.py` use the foreign work.
- A run without soft-core has the same energies and speed as before.

## Next: fixed lambda windows in `fep.py`

Replace the adaptive walk by a fixed schedule. Every window then knows its
neighbours, has exact foreign work to them, and its `.fhe` file can go straight
into BAR. Three things to decide first:

- [ ] Default number of windows and their spacing (proposal: 12, evenly spaced).
- [ ] Starting configuration of each window (proposal: from 0.5 outwards, each
      window from the end of the one before).
- [ ] End windows at exactly 0 and 1 (proposal: yes).

Then:

- [ ] Add `schedule(n)` and `run_schedule(engine, lambdas, out, ...)`.
- [ ] Remove `walk`, `next_lambda`, `_pair`, and the interpolation in `work()`
      with `foreign_targets` and `_soft_coordinate`. `overlap` takes two windows
      that are each other's foreign window.
- [ ] `run_hrex.py`: call `run_schedule` instead of `walk`.
- [ ] Move `read_config` and `write_qfep_inp` out of `run_walk.py`, then delete it.
- [ ] `test/test_fep.py`: drop the walk test.

## Clean up

- [ ] `energy.cpp`: the zeroing of states with lambda = 0 is commented out;
      delete it, and fix the comment above `unpack()` and the one in
      `Handler::update_energy_totals()` that still mention it.
- [ ] `fep.py`, `hrex.py`: remove the "lambda of exactly 0 cannot be used" checks.
- [ ] `test/soft_pairs/README.md` mentions the endpoint zeroing.
- [ ] Rebuild the Python extension: the `.so` files next to the package predate
      `RunResult.foreign_work`.

## Not verified yet

- [ ] A complete free energy with Gapsys: many windows, long enough to compare
      with a reference. So far `analyze_fhe` was only shown to read the files.
- [ ] `bar`, `neighbour_work` and `free_energy` in `hrex.py`, and the
      `hrex.tsv` / `hrex_dg.txt` output of `run_hrex.py`.
- [ ] The existing tests in `python/test`: they need the benchmark folder
      outside the repository.
- [ ] `walk` end to end (only its parts were run); moot once it is removed.
- [ ] Forces of the whole system against finite differences. The pair functions
      are checked that way, and CPU and GPU energies agree.
- [ ] 1-4 soft pairs on their own; they are covered only as part of the totals.
- [ ] Timing in SPFP, and of Gapsys together with LRF.
- [ ] Whether Fortran Q reports the energies of a state with lambda = 0;
      `test/runTEST.py` compares the logs.

## Tests to add to the repository

- [ ] Python test of `foreign_work` against a direct evaluation. The script
      used so far needs `1_eq5/`, which is not tracked; pick a tracked system.
- [ ] End-to-end test of the `.fhe` file.

## Known limits

- `softcore_method standard` is not implemented: `[softcore]` is read and
  ignored, without a message. Beutler is rejected.
- `[softcore]` rows are matched to Q atoms by their order; the atom number in
  the first column is not used. The section must list every Q atom once.
- A flag on a state where the atom is a dummy does nothing. FEP files written
  for the standard method have the flags there; adding `softcore_method gapsys`
  to such a file runs the hard potential. This is not an error on purpose: not
  flagging an atom is a valid choice.
- Only two FEP states, no `[off_diagonals]`, no `[soft_pairs]`.
- At most 8 windows in `[foreign_lambdas]`.
- `foreign_hamiltonian_mode` in the md input is ignored.
- The foreign work leaves out the lambda dependence of the water boundary
  terms: about 1e-4 kcal/mol on cdk2.
- Energies above about 8.6e9 kcal/mol overflow the fixed-point accumulators.
  A state with lambda = 0 exerts no force, so atoms can overlap its hard
  (unflagged) atoms; its reported energy and the foreign work are then garbage.
  The total potential is not affected.
- `work()` in `fep.py` interpolates between foreign windows: up to 0.08 kT off
  on cdk2. Fine for placing windows, not for free energies.

## Later

- [ ] Compute energies only on the steps that use them. Dropping the energy
      accumulation from the pair kernel alone saves 3-5% on cdk2;
      `Handler::output_requirements(iteration).energy` is the switch.
- [ ] Beutler soft-core. Its Coulomb softening reaches about 10 A, so it needs
      a decision on the pairs that go through LRF.
- [ ] Keep the atom number of `[softcore]` rows in the parser.
- [ ] QligFEP: generate qgpu runs with `-sc gapsys` once #140 is merged.

## Seen on the way, not from this work

- CPU and GPU differ in the kinetic energy at step 0 by 18.8 kcal/mol (0.8 K)
  on cdk2, with and without soft-core; the trajectories then diverge.
- Two runs of the same input on the GPU give `.en` files that differ in a few
  bytes from about frame 13 on.
