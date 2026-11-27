"""Replica exchange FEP for one leg and one replicate.

    python -m fep.run_hrex <input_dir> --windows 16 --replicate 1 --out <folder>

The run heats both end states, finds the lambda windows with a sweep between them, and runs
replica exchange on them. Written to <folder>/rep_<replicate>:

    run.log                 what the run did and its status, one line per event
    engine.log              what the engine of this process printed (heating and sweep)
    end0.re, end1.re        the heated end states; a rerun loads them
    windows.txt             the lambda windows and the swap acceptance expected for them
    production.npz          the work between neighbours, the swap history and the acceptance;
                            rewritten during the run (see replicate_exchange.save_result)
    result.txt              the free energy and a line for every pair of neighbours
"""
import os
import shutil
import socket
import time
from dataclasses import replace
from pathlib import Path

import numpy as np

from qgpu import Config, Engine

from fep.async_engine import AsyncEngine
from fep.equilibration import equilibrate
from fep.log import log, setup
from fep.parser import parse_args, q_atoms
from fep.restraints import (release_past_joints, release_rotors, restrain_core, topology_bonds,
                            topology_masses)
from fep.replicate_exchange import Window, replica_exchange, free_energy, save_result, summary
from fep.windows import get_windows, place_windows, refine_table

SEEDS = (2924, 25360, 21448, 9187, 31337, 17761, 4421, 12903, 27077, 6359)    # one per replicate


def read_config(template, seed):
    """The physics of the md template, without its restart and without its own output files."""
    # restart="": start from the coordinates in the topology; the workflow equilibrates itself,
    # from 1 K (both are fixed once the engine exists)
    config = Config.from_inp(template, restart="", seed=seed, initial_temperature=1.0)
    return replace(config, output=1000, energy=0, energy_file="",
                   trajectory=0, trajectory_file="", final_file="")


# REST2 windows sit this far from their end state in lambda: the engine needs every window it
# computes work for to have its own lambda, and 1e-6 in lambda changes the energy by ~1e-5 kcal/mol.
REST_OFFSET = 1e-6


def chain_with_rest(lambdas, states, rest):
    """The windows and their starting states: the lambda windows, with a REST2 branch at both ends.

    rest: the scales of each branch, from the end state outwards (e.g. (0.8, 0.64)).
    """
    chain = [Window(float(lam)) for lam in lambdas]
    states = list(states)
    for i, scale in enumerate(rest, 1):
        chain.insert(0, Window(i * REST_OFFSET, scale))
        states.insert(0, states[0])
        chain.append(Window(1.0 - i * REST_OFFSET, scale))
        states.append(states[-1])
    return chain, states


def build_engine(config, device, precision, tries=4):
    """The engine of this process; a GPU that is not ready yet (MPS starting) is tried again."""
    from qgpu import CudaError
    for attempt in range(1, tries + 1):
        try:
            return Engine(config, device=device, precision=precision)
        except CudaError as error:
            if attempt == tries:
                raise
            log.info("no GPU yet (%s), trying again in 15 s", error)
            time.sleep(15)


def run_one_replicate(inputs, windows, device, precision, seed, out, settle_cycles, cycles, block,
                      rest=(), save_atoms=(), refine=False, rotors=False, core=False, joints=False):
    """Heating, windows and replica exchange of one replicate; everything is written into `out`."""
    started = time.time()
    config = read_config(inputs.template, seed)
    if joints:
        kept, dropped = release_past_joints(config.dist_restraints, topology_bonds(inputs.topology),
                                            q_atoms(inputs.fep), topology_masses(inputs.topology))
        config = replace(config, dist_restraints=tuple(kept))
        log.info("pair restraints between the ligand copies: kept %d, released %d past an sp3 joint %s",
                 len(kept), len(dropped), [(r.ai, r.aj) for r in dropped])
    elif core:
        kept, dropped = restrain_core(config.dist_restraints, topology_bonds(inputs.topology),
                                      q_atoms(inputs.fep), topology_masses(inputs.topology))
        config = replace(config, dist_restraints=tuple(kept))
        log.info("pair restraints between the ligand copies: kept %d in the rigid core %s, released %d",
                 len(kept), [(r.ai, r.aj) for r in kept], len(dropped))
    elif rotors:
        kept, dropped = release_rotors(config.dist_restraints, topology_bonds(inputs.topology),
                                       q_atoms(inputs.fep))
        config = replace(config, dist_restraints=tuple(kept))
        log.info("pair restraints between the ligand copies: kept %d, released %d past a rotatable bond %s",
                 len(kept), len(dropped), [(r.ai, r.aj) for r in dropped])
    # The engine of this process prints its energies to standard output: keep them in a file.
    os.dup2(os.open(out / "engine.log", os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644), 1)
    engine = build_engine(config, device, precision)
    log.info("system: %d atoms, %g K, %g fs steps", engine.n_atoms, config.temperature, config.stepsize)

    log.info("== 1. heating the two end states")
    start = engine.get_state()
    end0 = equilibrate(engine, config, (1.0, 0.0), out / "end0.re")
    engine.set_state(start)
    end1 = equilibrate(engine, config, (0.0, 1.0), out / "end1.re")

    log.info("== 2. finding %d windows", windows)
    lambdas, states, acceptance, table = get_windows(engine, config, end0, end1, out / "sweep.en", windows,
                                                     with_table=True)
    (out / "windows.txt").write_text(
        f"# expected swap acceptance between neighbours: {acceptance:.3f}\n"
        + "".join(f"{lam:.6f}\n" for lam in lambdas))
    for scratch in (out / "sweep.en", out / "sweep.en.fhe"):
        scratch.unlink(missing_ok=True)

    log.info("== 3. replica exchange")
    chain, states = chain_with_rest(lambdas, states, rest)
    if rest:
        log.info("REST2 branch at each end state, scales %s: %d windows in all", rest, len(chain))
    atoms = sorted(set(q_atoms(inputs.fep)) | set(save_atoms))
    log.info("saving the coordinates of %d atoms (the Q atoms and %d more) after every cycle",
             len(atoms), len(set(save_atoms) - set(q_atoms(inputs.fep))))
    (out / "saved_atoms.txt").write_text("".join(f"{a + 1}\n" for a in atoms))
    began = time.time()
    # one engine per window, each in its own process
    engines = [AsyncEngine(config, device=device, precision=precision) for _ in chain]
    try:
        log.info("started %d engine processes (%d atoms each) in %.0f s", len(engines),
                 engines[-1].n_atoms, time.time() - began)
        settled = replica_exchange(engines, config, chain, states, out / "scratch", cycles=settle_cycles,
                                   block=block, seed=seed + 1, name="settling")
        log.info("settling: acceptance per pair %s", " ".join(f"{a:.2f}" for a in settled.acceptance))
        if refine:
            settled = _refine(engines, config, chain, settled, table, lambdas, windows, rest, out,
                              settle_cycles, block, seed)
            chain = settled.windows
        result = replica_exchange(engines, config, chain, settled.states, out / "scratch", cycles=cycles,
                                  block=block, seed=seed, name="production", save=out / "production.npz",
                                  save_atoms=atoms)
    finally:
        for async_engine in engines:
            async_engine.close()

    shutil.rmtree(out / "scratch", ignore_errors=True)   # the engines' own energy files
    save_result(result, out / "production.npz")
    (out / "result.txt").write_text(summary(result, config.temperature))
    log.info("production: acceptance per pair %s", " ".join(f"{a:.2f}" for a in result.acceptance))
    log.info("== done in %.1f min: dG %.3f kcal/mol (first half %.3f, second half %.3f), round trips %d",
             (time.time() - started) / 60, free_energy(result, config.temperature),
             free_energy(result, config.temperature, 0.0, 0.5), free_energy(result, config.temperature, 0.5, 1.0),
             int(result.round_trips.sum()))
    log.info("results in %s", out.resolve())
    return result


def _refine(engines, config, chain, settled, table, lambdas, n_windows, rest, out, settle_cycles, block, seed):
    """Place the windows again with the work measured during settling, and settle the new ones."""
    first = len(rest)                       # the lambda windows sit between the two REST branches
    up = settled.work_up[first:first + n_windows - 1]
    down = settled.work_down[first:first + n_windows - 1]
    refined, measured, predicted = refine_table(table, lambdas, up, down, config.temperature)
    log.info("refining the windows: acceptance from the settling work per pair %s",
             " ".join(f"{a:.2f}" for a in measured))
    log.info("    (the sweep had predicted %s)", " ".join(f"{a:.2f}" for a in predicted))
    new, expected = place_windows(refined, n_windows)
    log.info("placed %d windows again, expected swap acceptance %.2f", n_windows, expected)
    log.info("    lambdas: %s", " ".join(f"{lam:.3f}" for lam in new))
    # every new window starts from the settled state of the old lambda window nearest to it
    old_states = settled.states[first:first + n_windows]
    states = [old_states[int(np.argmin(np.abs(np.asarray(lambdas) - lam)))] for lam in new]
    chain, states = chain_with_rest(new, states, rest)
    # the REST branches keep their settled states
    for i in range(first):
        states[i] = settled.states[i]
        states[-1 - i] = settled.states[-1 - i]
    (out / "windows.txt").write_text(
        f"# placed again after settling; expected swap acceptance between neighbours: {expected:.3f}\n"
        + "".join(f"{lam:.6f}\n" for lam in new))
    return replica_exchange(engines, config, chain, states, out / "scratch", cycles=settle_cycles,
                            block=block, seed=seed + 2, name="settling again")


def main(argv=None):
    args = parse_args(argv)
    if args.replicate > len(SEEDS):
        raise SystemExit(f"--replicate {args.replicate}: there are {len(SEEDS)} seeds")
    seed = SEEDS[args.replicate - 1]
    out = args.out / f"rep_{args.replicate}"
    out.mkdir(parents=True, exist_ok=True)
    setup(out / "run.log")
    log.info("fep.run_hrex on %s, GPU device %d, %s precision", socket.gethostname(), args.device, args.precision)
    log.info("input %s", args.input_dir.directory)
    log.info("replicate %d (seed %d), %d windows, settling %d and production %d cycles of %d steps, "
             "REST2 end branches %s", args.replicate, seed, args.windows, args.settle_cycles, args.cycles,
             args.block, args.rest or "none")
    try:
        run_one_replicate(args.input_dir, args.windows, args.device, args.precision, seed, out,
                          args.settle_cycles, args.cycles, args.block, args.rest, args.save_atoms,
                          args.refine, args.release_rotors, args.restrain_core, args.release_joints)
    except BaseException:
        log.exception("the run failed")
        raise


if __name__ == "__main__":
    main()
