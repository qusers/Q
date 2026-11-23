"""Hamiltonian replica exchange between lambda windows (prototype).

Every window keeps its lambda; the configurations move. One engine runs the windows in
turn: set_state, run a block, get_state. After each cycle, neighbouring windows try to
swap configurations.

Swapping the configurations x_a, x_b of windows a and b changes the total energy by

    [U_b(x_a) - U_a(x_a)] + [U_a(x_b) - U_b(x_b)]

Each bracket is the foreign work of one window to the other at its final coordinates
(RunResult.final_foreign_work), so every window is run with its neighbours as foreign
windows. Without soft-core this is (L_a - L_b) @ (q(x_b) - q(x_a)); with soft-core the Q
energies depend on lambda and only the foreign work is right.
"""

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .errors import InputError
from .fep import K_B, _energy_file, append_fhe


@dataclass(frozen=True)
class HrexResult:
    lambdas: list            # [(lambda_1, lambda_2)] per window
    q_energy: list           # per window: (n_frames, 2), all cycles
    replica: np.ndarray      # (cycles + 1, n_windows): which configuration each window holds
    tried: np.ndarray        # (n_windows - 1,) swap attempts per neighbouring pair
    accepted: np.ndarray     # (n_windows - 1,) accepted swaps
    coords: np.ndarray       # (cycles, n_windows, n_saved, 3) or empty; see save_atoms
    states: list             # Engine.get_state() of every window after the last cycle
    # per window: lambda_2 of the windows its foreign work is for (itself and its neighbours,
    # sorted), and the work U(that window) - U(this window) of every frame, (n_frames, n_targets)
    foreign_lambdas: list
    foreign_work: list

    @property
    def acceptance(self):
        return self.accepted / np.maximum(self.tried, 1)


def swap_energy(work_a_to_b, work_b_to_a):
    """Energy change (kcal/mol) when windows a and b swap their configurations.

    work_a_to_b: U_b - U_a at the configuration of a; work_b_to_a: U_a - U_b at that of b.
    """
    return float(work_a_to_b) + float(work_b_to_a)


def bar(work_forward, work_backward, temperature):
    """Free energy F_b - F_a (kcal/mol) from the work a -> b of the frames of a and
    b -> a of the frames of b (Bennett acceptance ratio)."""
    kt = K_B * temperature
    wf = np.asarray(work_forward, dtype=float) / kt
    wr = np.asarray(work_backward, dtype=float) / kt
    m = np.log(len(wf) / len(wr))

    def excess(df):     # rises with df; zero at the BAR estimate
        forward = np.logaddexp(0.0, m + wf - df)
        backward = np.logaddexp(0.0, -m + wr + df)
        return np.logaddexp.reduce(-forward) - np.logaddexp.reduce(-backward)

    lo, hi = min(wf.min(), -wr.max()) - 1.0, max(wf.max(), -wr.min()) + 1.0
    for _ in range(100):
        middle = 0.5 * (lo + hi)
        lo, hi = (middle, hi) if excess(middle) < 0 else (lo, middle)
    return 0.5 * (lo + hi) * kt


def neighbour_work(result):
    """(up, down): up[k] is the work of every frame of window k to window k + 1, down[k] that
    of window k + 1 to window k; each (n_windows - 1, n_frames), kcal/mol."""
    up, down = [], []
    for k in range(len(result.lambdas) - 1):
        i = np.flatnonzero(np.isclose(result.foreign_lambdas[k], result.lambdas[k + 1][1], rtol=0.0, atol=1e-12))[0]
        j = np.flatnonzero(np.isclose(result.foreign_lambdas[k + 1], result.lambdas[k][1], rtol=0.0, atol=1e-12))[0]
        up.append(result.foreign_work[k][:, i])
        down.append(result.foreign_work[k + 1][:, j])
    return np.array(up), np.array(down)


def free_energy(result, temperature, skip=0):
    """F(next window) - F(window) for every neighbouring pair, from the foreign work.

    result: HrexResult; skip: frames left out at the start of every window. Windows are in
    order of rising lambda_1, so the free energy of state 1 -> state 2 is minus the sum.
    """
    up, down = neighbour_work(result)
    return np.array([bar(u[skip:], d[skip:], temperature) for u, d in zip(up, down)])


def _column(result, lambda2):
    """Index of the foreign window with this lambda_2 in a RunResult."""
    match = np.flatnonzero(np.isclose(result.foreign_lambdas, lambda2, rtol=0.0, atol=1e-12))
    if match.size != 1:
        raise InputError(f"no foreign work for lambda_2 = {lambda2:g} in this run")
    return int(match[0])


def round_trips(replica):
    """How often each configuration went from one end window to the other and back."""
    replica = np.asarray(replica)
    n = replica.shape[1]
    trips = np.zeros(n, dtype=int)
    last_end = np.full(n, -1)           # 0: was last at the first window, 1: at the last
    half = np.zeros(n, dtype=int)       # end-to-end crossings so far
    for row in replica:
        for end, window in ((0, 0), (1, n - 1)):
            r = row[window]
            if last_end[r] not in (-1, end):
                half[r] += 1
            last_end[r] = end
    trips[:] = half // 2
    return trips


def hrex(engine, windows, out, cycles, block=500, energy=10, exchange=True,
         save_atoms=None, seed=0):
    """Run `cycles` cycles of replica exchange over `windows`, starting from their states.

    windows: objects with .lambdas and .state (e.g. the result of fep.walk), any order.
    Each cycle runs `block` steps in every window and then tries to swap neighbours:
    pairs (0,1), (2,3), ... in even cycles and (1,2), (3,4), ... in odd ones.
    exchange=False runs the same schedule without swaps (the control).
    save_atoms: 0-based atom indices whose coordinates are kept after every block.

    Writes md_XXXX_YYYY.en per window into `out`, and next to each its foreign-work file
    (.en.fhe) for the window itself and its neighbours. Afterwards the engine has energy
    output and foreign work switched off.
    """
    if cycles < 1 or block < 1:
        raise InputError("need cycles >= 1 and block >= 1")
    if energy < 1 or block % energy:
        raise InputError("block must be a multiple of the energy interval")
    windows = sorted(windows, key=lambda w: w.lambdas[0])
    if len(windows) < 2:
        raise InputError("replica exchange needs at least two windows")
    for w in windows:
        if 0.0 in w.lambdas:
            raise InputError("a lambda of exactly 0 cannot be used: "
                             "the engine reports no energy for that state")

    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    n = len(windows)
    lambdas = [tuple(w.lambdas) for w in windows]
    states = [w.state for w in windows]
    beta = 1.0 / (K_B * engine.config.temperature)
    rng = np.random.default_rng(seed)

    if len({l[1] for l in lambdas}) != n:
        raise InputError("replica exchange needs windows with different lambdas")
    # each window computes the work to itself (zero; the .fhe file needs the column) and its neighbours
    targets = [tuple(sorted({lambdas[j][1] for j in (k - 1, k, k + 1) if 0 <= j < n})) for k in range(n)]

    files = [_energy_file(out, l) for l in lambdas]
    chunk = out / "hrex_block.en"
    for f in files:
        f.write_bytes(b"")
        Path(str(f) + ".fhe").unlink(missing_ok=True)
    fhe_frames = [0] * n

    holder = np.arange(n)                       # configuration held by each window
    replica = np.empty((cycles + 1, n), dtype=np.int32)
    replica[0] = holder
    tried = np.zeros(n - 1, dtype=np.int64)
    accepted = np.zeros(n - 1, dtype=np.int64)
    parts = [[] for _ in range(n)]
    works = [[] for _ in range(n)]
    save_atoms = None if save_atoms is None else np.asarray(save_atoms, dtype=int)
    coords = np.empty((cycles, n, 0 if save_atoms is None else len(save_atoms), 3), dtype=np.float32)

    for cycle in range(cycles):
        results = []
        for k in range(n):
            engine.set_state(states[k])
            engine.update_config(steps=block, lambdas=lambdas[k], energy=energy,
                                 energy_file=str(chunk), foreign_lambdas=targets[k])
            result = engine.run()
            states[k] = engine.get_state()
            parts[k].append(result.q_energy)
            works[k].append(result.foreign_work)
            results.append(result)
            with open(files[k], "ab") as joined:    # .en files have no header: append
                joined.write(chunk.read_bytes())
            fhe_frames[k] = append_fhe(chunk, files[k], fhe_frames[k])
            if save_atoms is not None:
                coords[cycle, k] = np.asarray(states[k]["coords"]).reshape(-1, 3)[save_atoms]

        if exchange:
            for k in range(cycle % 2, n - 1, 2):
                a, b = results[k], results[k + 1]
                du = swap_energy(a.final_foreign_work[_column(a, lambdas[k + 1][1])],
                                 b.final_foreign_work[_column(b, lambdas[k][1])])
                tried[k] += 1
                if du <= 0 or rng.random() < np.exp(-beta * du):
                    accepted[k] += 1
                    states[k], states[k + 1] = states[k + 1], states[k]
                    holder[[k, k + 1]] = holder[[k + 1, k]]
        replica[cycle + 1] = holder

    engine.update_config(energy=0, energy_file="", foreign_lambdas=())
    chunk.unlink(missing_ok=True)
    return HrexResult(lambdas, [np.concatenate(p) for p in parts], replica, tried, accepted,
                      coords, states, targets, [np.concatenate(w) for w in works])
