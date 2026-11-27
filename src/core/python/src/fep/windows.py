"""Choose the lambda windows of one leg from a short sweep between the two end states.

    swept = sweep(engine, config, end0, end1, energy_file)
    table = acceptance_table(swept, config.temperature)
    lambdas, acceptance = place_windows(table, n_windows)
    states = [swept.states[nearest(l)] for l in lambdas]

Here lambda is the lambda of state 2: a window has lambdas = (1 - lambda, lambda), lambda = 0 is
state 1 and lambda = 1 is state 2.

The sweep visits the grid 0, 0.01, ..., 1 in both directions, a short run at every point, and
records the work to grid points nearby. From that work comes the swap
acceptance to expect between two windows (acceptance_table), and the windows are placed so that
all neighbours have the same expected acceptance (place_windows). The placement follows
Opt-PSO (Zhang et al., JCTC 2024, 20, 3935) with the swap acceptance in place of the overlap
index and a sweep in place of the burn-in simulation of every window.
"""
import time
from dataclasses import dataclass

import numpy as np

from fep.log import log
from fep.metrics import KB

GRID = 100                 # grid points are k / GRID, k = 0 ... GRID
OFFSETS = (1, 2, 4, 8)     # work is recorded to the points this many grid steps away, both sides


@dataclass(frozen=True)
class Sweep:
    """What sweep() returns. Grid points are numbered k = 0 ... GRID, lambda = k / GRID."""
    work: dict      # (a, b) -> U(b) - U(a) at the frames of point a, kcal/mol
    states: list    # the state at the end of the run at every point


def grid_lambda(k):
    return round(k / GRID, 6)


def nearest(lam):
    """The grid point nearest to lambda."""
    return int(round(lam * GRID))


def _one_pass(engine, config, start, points, energy_file, steps, energy, discard):
    """Run at every grid point of `points` in turn, from the state `start`.

    Returns (work, states) of the points, as in Sweep.
    """
    skip = int(discard * (steps // energy))
    work, states = {}, {}
    engine.set_state(start)
    points = list(points)
    started = time.time()
    log.info("sweep from lambda %.2f to %.2f: %d points of %d steps", grid_lambda(points[0]),
             grid_lambda(points[-1]), len(points), steps)
    for count, k in enumerate(points, 1):
        others = [k + sign * d for d in OFFSETS for sign in (1, -1) if 0 <= k + sign * d <= GRID]
        lam = grid_lambda(k)
        engine.update_config(
            steps=steps, stepsize=config.stepsize, temperature=config.temperature,
            bath_coupling=config.bath_coupling, constrain_hydrogens=True, constrain_solute=False,
            seq_restraints=config.seq_restraints, dist_restraints=config.dist_restraints,
            lambdas=(1.0 - lam, lam), foreign_lambdas=[grid_lambda(o) for o in others],
            energy=energy, energy_file=str(energy_file))
        result = engine.run()
        for column, other in enumerate(others):
            work[k, other] = result.foreign_work[skip:, column].copy()
        states[k] = engine.get_state()
        if count % 25 == 0 or count == len(points):
            log.info("    point %d of %d (lambda %.2f), %.0f s", count, len(points), lam, time.time() - started)
    engine.update_config(energy=0, energy_file="", foreign_lambdas=())
    return work, states


def sweep(engine, config, end0, end1, energy_file, steps=500, energy=10, discard=0.2):
    """Visit every grid point twice: going up from end0 and going down from end1.

    end0, end1: the equilibrated states at lambda 0 and 1. Every point runs `steps` steps from
    the final state of the point before, with the work to its neighbours every `energy` steps;
    the first `discard` (a fraction) of the frames of a point is left out. A point lags behind
    its lambda towards the side it was reached from; with both directions the two lags are on
    opposite sides. The state kept for a point is that of the pass that started nearer to it.

    `config` is the template: the sweep runs with its step size, temperature and restraints,
    and without solute constraints other than bonds to hydrogen. `energy_file` is a scratch
    file that the engine overwrites at every point.
    """
    up, up_states = _one_pass(engine, config, end0, range(0, GRID + 1), energy_file, steps, energy, discard)
    down, down_states = _one_pass(engine, config, end1, range(GRID, -1, -1), energy_file, steps, energy, discard)
    work = {pair: np.concatenate([up[pair], down[pair]]) for pair in up}
    middle = GRID // 2
    states = [up_states[k] if k <= middle else down_states[k] for k in range(GRID + 1)]
    return Sweep(work, states)


def swap_acceptance(work_ab, work_ba, kT):
    """Expected acceptance of a swap between windows a and b.

    work_ab: U(b) - U(a) at frames of a; work_ba: U(a) - U(b) at frames of b. A swap of one
    frame of each is accepted with min(1, exp(-(w_ab + w_ba) / kT)); this is the mean over all
    combinations of the frames.
    """
    total = work_ab[:, None] + work_ba[None, :]
    return float(np.exp(-np.maximum(total, 0.0) / kT).mean())


@dataclass(frozen=True)
class AcceptanceTable:
    """Expected swap acceptance between any two lambdas.

    Between two lambdas a distance d apart and with midpoint m, the acceptance is written as
    exp(-z(m, d) * d**2), which is 1 at d = 0. z[i, j] is z at distance distances[i] and
    midpoint midpoints[j], smoothed along the midpoint.
    """
    distances: np.ndarray   # (n_offsets,)
    midpoints: np.ndarray   # (n_midpoints,)
    z: np.ndarray           # (n_offsets, n_midpoints)

    def __call__(self, lam_a, lam_b):
        d = abs(lam_b - lam_a)
        if d == 0.0:
            return 1.0
        m = 0.5 * (lam_a + lam_b)
        # z at this midpoint for every measured distance, then between the distances
        at_m = [np.interp(m, self.midpoints, row) for row in self.z]
        z = np.interp(d, self.distances, at_m)  # beyond the measured distances: the nearest one
        return float(np.exp(-z * d * d))


def acceptance_table(swept, temperature, width=0.02, floor=1e-4):
    """The acceptance table of a sweep.

    width: the z values of one distance are averaged along the midpoint with a Gaussian of this
    standard deviation (in lambda), because one grid point has few frames. floor: lowest
    acceptance that is used as measured (lower ones are set to it: they are too noisy).
    """
    kT = KB * temperature
    midpoints = np.arange(0.0, 1.0 + 1e-9, 0.5 / GRID)
    distances, rows = [], []
    for offset in OFFSETS:
        d = offset / GRID
        at, z = [], []
        for a in range(0, GRID + 1 - offset):
            b = a + offset
            acceptance = swap_acceptance(swept.work[a, b], swept.work[b, a], kT)
            at.append(0.5 * (a + b) / GRID)
            z.append(-np.log(max(acceptance, floor)) / d ** 2)
        at, z = np.array(at), np.array(z)
        weights = np.exp(-0.5 * ((midpoints[:, None] - at[None, :]) / width) ** 2)
        distances.append(d)
        rows.append((weights * z).sum(axis=1) / weights.sum(axis=1))
    return AcceptanceTable(np.array(distances), midpoints, np.array(rows))


def refine_table(table, lambdas, work_up, work_down, temperature, samples=1000):
    """The table corrected by work measured between neighbouring windows at `lambdas`.

    work_up[k], work_down[k]: work from window k to k + 1 at the frames of k, and back at the
    frames of k + 1 (Result.work_up and work_down of a replica exchange run on these windows).
    Between two neighbours the table is scaled so that it gives the acceptance computed from the
    measured work; the factor applies to all distances at the midpoints between the two windows.
    Returns (table, measured, predicted): the new table and the acceptance per pair before.
    """
    kT = KB * temperature
    factor = np.ones_like(table.midpoints)
    measured, predicted = [], []
    for k in range(len(lambdas) - 1):
        a, b = float(lambdas[k]), float(lambdas[k + 1])
        up = work_up[k][:: max(1, len(work_up[k]) // samples)]
        down = work_down[k][:: max(1, len(work_down[k]) // samples)]
        measured.append(swap_acceptance(up, down, kT))
        predicted.append(table(a, b))
        z_now = -np.log(max(predicted[-1], 1e-6))
        z_new = -np.log(max(measured[-1], 1e-4))
        inside = (table.midpoints >= a) & (table.midpoints <= b)
        factor[inside] = z_new / max(z_now, 1e-6)
    refined = AcceptanceTable(table.distances, table.midpoints, table.z * factor[None, :])
    return refined, np.array(measured), np.array(predicted)


def _next_lambda(table, lam, acceptance):
    """The lambda above `lam` whose expected acceptance with `lam` is `acceptance` (at most 1)."""
    if table(lam, 1.0) >= acceptance:
        return 1.0
    low, high = lam, 1.0  # acceptance(low) >= target > acceptance(high)
    for _ in range(40):
        mid = 0.5 * (low + high)
        low, high = (mid, high) if table(lam, mid) >= acceptance else (low, mid)
    return low


def _walk(table, n_windows, acceptance):
    lambdas = [0.0]
    for _ in range(n_windows - 1):
        lambdas.append(_next_lambda(table, lambdas[-1], acceptance))
    return lambdas


def place_windows(table, n_windows):
    """n_windows lambdas from 0 to 1 whose neighbours all have the same expected acceptance.

    Returns (lambdas, acceptance): the windows and that common acceptance. A low value means
    that n_windows is too small for this leg.
    """
    if n_windows < 2:
        raise ValueError("n_windows must be at least 2")
    # The highest acceptance with which n_windows steps of equal acceptance still reach 1.
    low, high = 0.0, 1.0  # reaches 1 with low, not with high
    for _ in range(40):
        mid = 0.5 * (low + high)
        low, high = (mid, high) if _walk(table, n_windows, mid)[-1] >= 1.0 else (low, mid)
    lambdas = _walk(table, n_windows, low)
    lambdas[-1] = 1.0
    return np.array(lambdas), low


def get_windows(engine, config, end0, end1, energy_file, n_windows, steps=500, energy=10, discard=0.2,
                with_table=False):
    """The windows of this leg: (lambdas, states, acceptance).

    lambdas: n_windows values from 0 to 1; states: a starting state for each; acceptance: the
    swap acceptance to expect between neighbours. with_table=True adds the acceptance table.
    """
    swept = sweep(engine, config, end0, end1, energy_file, steps, energy, discard)
    table = acceptance_table(swept, config.temperature)
    lambdas, acceptance = place_windows(table, n_windows)
    log.info("placed %d windows, expected swap acceptance %.2f", n_windows, acceptance)
    log.info("    lambdas: %s", " ".join(f"{lam:.3f}" for lam in lambdas))
    even = np.linspace(0.0, 1.0, n_windows)
    log.info("    (evenly spaced windows would have: %s)",
             " ".join(f"{table(a, b):.2f}" for a, b in zip(even[:-1], even[1:])))
    states = [swept.states[nearest(lam)] for lam in lambdas]
    if with_table:
        return lambdas, states, acceptance, table
    return lambdas, states, acceptance
