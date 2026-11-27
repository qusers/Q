"""Replica exchange along a chain of windows.

    windows = [Window(lam) for lam in lambdas]
    engines = [AsyncEngine(config) for _ in windows]
    result = replica_exchange(engines, config, windows, states, scratch, cycles=1000)
    dg = free_energy(result, config.temperature)

A window is a lambda and a REST2 scale. Here lambda is the lambda of state 2: a window has
lambdas = (1 - lambda, lambda). The windows are given in the order of the chain: neighbours in
the list are the ones that swap. A chain is usually the lambda windows from 0 to 1, all with
scale 1. REST2 windows can be put in front of the first and behind the last one, with the same
lambda and falling scale: they heat the Q atoms of an end state and are not part of the free
energy (see physical_ends).

Every window has an engine of its own and keeps its lambda and scale; the configurations move.
Swapping the configurations x_a, x_b of windows a and b changes the total energy by

    [U_b(x_a) - U_a(x_a)] + [U_a(x_b) - U_b(x_b)]

Each bracket is the work of one window to the other at its coordinates, which the engine
computes (foreign work); it is right with soft-core and with REST2.
"""
import os
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from fep.log import log
from fep.metrics import KB, bar


@dataclass(frozen=True)
class Window:
    lam: float              # lambda of state 2
    scale: float = 1.0      # REST2 scale of the Q atoms; 1 = not scaled


@dataclass(frozen=True)
class Result:
    """What replica_exchange() returns. n = number of windows."""
    windows: list
    work_up: np.ndarray      # (n - 1, frames): U(window k + 1) - U(window k) at the frames of k, kcal/mol
    work_down: np.ndarray    # (n - 1, frames): U(window k) - U(window k + 1) at the frames of k + 1
    replica: np.ndarray      # (cycles + 1, n): the configuration that each window holds
    tried: np.ndarray        # (n - 1,) swap attempts of window k with k + 1
    accepted: np.ndarray     # (n - 1,) accepted ones
    states: list             # the state of every window after the last cycle
    coords: np.ndarray = None    # (saved cycles, n, atoms, 3) float32, see replica_exchange(save_atoms=)

    @property
    def acceptance(self):
        return self.accepted / np.maximum(self.tried, 1)

    @property
    def round_trips(self):
        """How often each configuration went from one end of the chain to the other and back."""
        return round_trips(self.replica)


def round_trips(replica):
    """How often each configuration went from one end of the chain to the other and back.

    replica: (cycles + 1, n), the configuration that each window holds after every cycle.
    Returns one count per configuration.
    """
    replica = np.asarray(replica)
    n = replica.shape[1]
    last_end = np.full(n, -1)           # the end each configuration was at last: 0, 1, or -1 for none yet
    crossings = np.zeros(n, dtype=int)  # times it arrived at the other end
    for row in replica:
        for end, window in ((0, 0), (1, n - 1)):
            configuration = row[window]
            if last_end[configuration] not in (-1, end):
                crossings[configuration] += 1
            last_end[configuration] = end
    return crossings // 2


def replica_exchange(engines, config, windows, states, scratch, cycles, block=500, energy=10,
                     seed=0, exchange=True, name="replica exchange", report=50, save=None,
                     save_atoms=None, save_every=1):
    """Run `cycles` cycles of replica exchange over `windows`, from their `states`.

    engines: one AsyncEngine per window (fep.async_engine); window k stays in engines[k] with
    its lambda and scale for the whole run. A cycle runs `block` steps in every window, all at
    the same time as far as AsyncEngine allows, and then tries to swap neighbours: windows
    (0, 1), (2, 3), ... after even cycles and (1, 2), (3, 4), ... after odd ones. An accepted
    swap moves the two configurations between their engines. The work to both neighbours is
    recorded every `energy` steps. exchange=False runs the same without swaps.

    `config` is the template: the run uses its step size, temperature and restraints, and no
    solute constraints other than bonds to hydrogen. `scratch` is a folder for the energy file
    of every engine, which the engine overwrites in every block. `seed` is for the swap
    decisions only. To continue a run (for example after a part that is thrown away), pass the
    states of its result.

    Every `report` cycles the status goes to the log under `name`, and with `save` (a file name)
    the result so far is written there (see save_result), so that a run that is cut off keeps
    what it has done. save_atoms: 0-based atom indices whose coordinates every window keeps after
    every `save_every`-th cycle (Result.coords; row i is the end of cycle (i + 1) * save_every).
    """
    windows = list(windows)
    n = len(windows)
    if n < 2:
        raise ValueError("replica exchange needs at least two windows")
    if len(states) != n or len(engines) != n:
        raise ValueError("one state and one engine per window")
    if cycles < 1 or energy < 1 or block % energy:
        raise ValueError("need cycles >= 1 and a block that is a multiple of the energy interval")

    kT = KB * config.temperature
    rng = np.random.default_rng(seed)
    frames = block // energy
    work_up = np.empty((n - 1, cycles * frames))
    work_down = np.empty((n - 1, cycles * frames))
    replica = np.empty((cycles + 1, n), dtype=np.int32)
    holder = np.arange(n)      # the configuration held by each window
    replica[0] = holder
    tried = np.zeros(n - 1, dtype=np.int64)
    accepted = np.zeros(n - 1, dtype=np.int64)
    atoms = None if save_atoms is None else np.asarray(save_atoms, dtype=np.int64)
    coords = np.zeros((cycles // save_every, n, 0 if atoms is None else len(atoms), 3), dtype=np.float32)

    scratch = Path(scratch)
    scratch.mkdir(parents=True, exist_ok=True)
    neighbours = [[j for j in (k - 1, k + 1) if 0 <= j < n] for k in range(n)]
    for k, (engine, window) in enumerate(zip(engines, windows)):
        engine.set_state(states[k])
        engine.update_config(
            steps=block, stepsize=config.stepsize, temperature=config.temperature,
            bath_coupling=config.bath_coupling, constrain_hydrogens=True, constrain_solute=False,
            seq_restraints=config.seq_restraints, dist_restraints=config.dist_restraints,
            lambdas=(1.0 - window.lam, window.lam), rest_scale=window.scale,
            foreign_lambdas=[windows[j].lam for j in neighbours[k]],
            foreign_rest=[windows[j].scale for j in neighbours[k]],
            energy=energy, energy_file=str(scratch / f"window_{k}.en"))

    log.info("%s: %d windows, %d cycles of %d steps (%.0f ps per window)", name, n, cycles, block,
             cycles * block * config.stepsize / 1000)
    started = time.time()
    for cycle in range(cycles):
        at = slice(cycle * frames, (cycle + 1) * frames)
        # work of every window to the next / previous one at its last frame
        final_up, final_down = np.zeros(n), np.zeros(n)
        handles = [engine.run() for engine in engines]
        for k, handle in enumerate(handles):
            result = handle.result()
            # the engine gives the frames inside the block and, apart from them, the last one
            work = np.vstack([result.foreign_work, result.final_foreign_work])
            for column, j in enumerate(neighbours[k]):
                if j > k:
                    work_up[k, at] = work[:, column]
                    final_up[k] = work[-1, column]
                else:
                    work_down[j, at] = work[:, column]
                    final_down[k] = work[-1, column]

        if exchange:
            for k in range(cycle % 2, n - 1, 2):
                change = final_up[k] + final_down[k + 1]
                tried[k] += 1
                if change <= 0 or rng.random() < np.exp(-change / kT):
                    accepted[k] += 1
                    lower, upper = engines[k].get_state(), engines[k + 1].get_state()
                    engines[k].set_state(upper)
                    engines[k + 1].set_state(lower)
                    holder[[k, k + 1]] = holder[[k + 1, k]]
        replica[cycle + 1] = holder
        if atoms is not None and (cycle + 1) % save_every == 0:
            for k, engine in enumerate(engines):
                coords[(cycle + 1) // save_every - 1, k] = engine.coordinates(atoms)

        done = cycle + 1
        if done % report == 0 or done == cycles:
            so_far = Result(windows, work_up[:, :done * frames], work_down[:, :done * frames],
                            replica[:done + 1], tried.copy(), accepted.copy(), [],
                            coords[:done // save_every])
            _report(name, so_far, config, done, cycles, block, time.time() - started)
            if save is not None:
                save_result(so_far, save)

    final = [engine.get_state() for engine in engines]
    for engine in engines:
        engine.update_config(energy=0, energy_file="", foreign_lambdas=(), rest_scale=1.0)
    return Result(windows, work_up, work_down, replica, tried, accepted, final, coords)


def _report(name, result, config, done, cycles, block, seconds):
    """One status line of a run in the log."""
    n = len(result.windows)
    rate = n * done * block / seconds
    line = (f"{name}: cycle {done} of {cycles}, {seconds / 60:.1f} min, {rate:.0f} steps/s in all "
            f"({rate * config.stepsize * 86400 / 1e6:.0f} ns/day), about {seconds / done * (cycles - done) / 60:.0f} min left; "
            f"acceptance mean {result.acceptance.mean():.2f}, lowest {result.acceptance.min():.2f}; "
            f"round trips {int(result.round_trips.sum())}")
    try:
        line += f"; dG so far {free_energy(result, config.temperature):.2f} kcal/mol"
    except ValueError:      # a chain without both physical end windows
        pass
    log.info(line)


def save_result(result, path):
    """Write a result (without the states) to an .npz file: lambdas, scales, work_up, work_down,
    replica, tried, accepted, as in Result. The file is replaced in one step."""
    path = Path(path)
    partial = path.with_name(path.name + ".partial.npz")
    np.savez(partial,
             lambdas=np.array([w.lam for w in result.windows]),
             scales=np.array([w.scale for w in result.windows]),
             work_up=result.work_up.astype(np.float32), work_down=result.work_down.astype(np.float32),
             replica=result.replica, tried=result.tried, accepted=result.accepted,
             coords=np.zeros((0, 0, 0, 3), np.float32) if result.coords is None else result.coords)
    os.replace(partial, path)


def load_result(path):
    """The result that save_result wrote (its states are empty)."""
    with np.load(path) as data:
        windows = [Window(float(lam), float(scale)) for lam, scale in zip(data["lambdas"], data["scales"])]
        return Result(windows, data["work_up"].astype(float), data["work_down"].astype(float),
                      data["replica"], data["tried"], data["accepted"], [],
                      data["coords"] if "coords" in data else None)


def summary(result, temperature):
    """The outcome of a run as text: the free energy, and a line for every pair of neighbours."""
    steps = free_energy_steps(result, temperature)
    trips = result.round_trips
    lines = [f"windows                {len(result.windows)}",
             f"frames per window      {result.work_up.shape[1]}",
             f"swap cycles            {len(result.replica) - 1}",
             f"round trips            {int(trips.sum())} (per configuration: {' '.join(str(t) for t in trips)})",
             f"acceptance             mean {result.acceptance.mean():.3f}, lowest {result.acceptance.min():.3f}"]
    try:
        lines += [f"dG state 1 -> state 2  {free_energy(result, temperature):.4f} kcal/mol",
                  f"dG first half          {free_energy(result, temperature, 0.0, 0.5):.4f}",
                  f"dG second half         {free_energy(result, temperature, 0.5, 1.0):.4f}"]
    except ValueError:
        pass
    lines += ["", "lambda   scale   next_lambda  acceptance  dF_to_next"]
    for k, step in enumerate(steps):
        a, b = result.windows[k], result.windows[k + 1]
        lines.append(f"{a.lam:.4f}   {a.scale:.3f}   {b.lam:.4f}       {result.acceptance[k]:.3f}       {step:9.4f}")
    return "\n".join(lines) + "\n"


def free_energy_steps(result, temperature, first=0.0, last=1.0):
    """F(window k + 1) - F(window k) for every neighbouring pair, kcal/mol (BAR on the work).

    first, last: the part of the frames that is used, as fractions of the run; (0, 0.5) is the
    first half.
    """
    frames = result.work_up.shape[1]
    part = slice(int(first * frames), int(last * frames))
    return np.array([bar(up[part], down[part], temperature)
                     for up, down in zip(result.work_up, result.work_down)])


def physical_ends(windows):
    """(i, j): the windows that are state 1 and state 2, that is lambda 0 and 1 with scale 1."""
    ends = []
    for lam in (0.0, 1.0):
        found = [k for k, w in enumerate(windows) if w.lam == lam and w.scale == 1.0]
        if len(found) != 1:
            raise ValueError(f"need exactly one window with lambda {lam:g} and scale 1")
        ends.append(found[0])
    if ends[0] > ends[1]:
        raise ValueError("the window of lambda 0 must come before that of lambda 1")
    return tuple(ends)


def free_energy(result, temperature, first=0.0, last=1.0):
    """Free energy of state 1 -> state 2, kcal/mol: the steps between the two physical windows."""
    i, j = physical_ends(result.windows)
    return float(free_energy_steps(result, temperature, first, last)[i:j].sum())
