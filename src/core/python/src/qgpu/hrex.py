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
import multiprocessing
from pathlib import Path
import shutil
import tempfile
import traceback

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
    # per window: (lambda_1, lambda_2) of the windows its foreign work is for (the one before,
    # itself, the one after), and the work U(that window) - U(this window) of every frame,
    # (n_frames, n_targets)
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
        i = result.foreign_lambdas[k].index(tuple(result.lambdas[k + 1]))
        j = result.foreign_lambdas[k + 1].index(tuple(result.lambdas[k]))
        up.append(result.foreign_work[k][:, i])
        down.append(result.foreign_work[k + 1][:, j])
    return np.array(up), np.array(down)


def free_energy(result, temperature, skip=0):
    """F(next window) - F(window) for every neighbouring pair, from the foreign work.

    result: HrexResult; skip: frames left out at the start of every window. Windows are in
    the order given to hrex (from state 2 to state 1), so the free energy of state 1 ->
    state 2 is minus the sum.
    """
    up, down = neighbour_work(result)
    return np.array([bar(u[skip:], d[skip:], temperature) for u, d in zip(up, down)])


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


@dataclass(frozen=True)
class Start:
    lambdas: tuple      # (lambda_1, lambda_2)
    state: dict         # Engine.get_state() to start the window from


def even_windows(n):
    """lambda_1 of n evenly spaced windows, 0 and 1 included."""
    if n < 2:
        raise InputError("need at least two windows")
    return [round(k / (n - 1), 3) for k in range(n)]


def window_pairs(lambda1):
    """(lambda_1, lambda_2) of the windows at `lambda1`, in order of rising lambda_1."""
    values = sorted({round(float(l), 3) for l in lambda1})
    if len(values) < 2 or values[0] < 0.0 or values[-1] > 1.0:
        raise InputError("need at least two different lambdas in [0, 1]")
    return [(t, round(1.0 - t, 3)) for t in values]


def prepare(engine, windows, steps=5000):
    """Give every window a starting configuration.

    windows: (lambda_1, lambda_2) of every window, in the order replica exchange will use
    (see window_pairs). The window where the two lambdas are closest starts from the
    engine's current state; going outwards to both ends, each window starts from the end of
    the one before it. Every window runs `steps` steps. Returns one Start per window, in the
    same order.
    """
    windows = [(float(a), float(b)) for a, b in windows]
    if len(windows) < 2:
        raise InputError("need at least two windows")
    middle = min(range(len(windows)), key=lambda k: abs(windows[k][0] - windows[k][1]))
    origin = engine.get_state()
    engine.update_config(energy=0, energy_file="", foreign_lambdas=())
    starts = {}
    for order in (range(middle, len(windows)), range(middle - 1, -1, -1)):
        engine.set_state(starts[middle].state if middle in starts else origin)
        for k in order:
            engine.update_config(steps=steps, lambdas=windows[k])
            engine.run()
            starts[k] = Start(windows[k], engine.get_state())
    return [starts[k] for k in range(len(windows))]


@dataclass(frozen=True)
class Block:
    """One window's block: what hrex needs from RunResult, and the state afterwards."""
    q_energy: np.ndarray
    foreign_lambdas: np.ndarray
    foreign_work: np.ndarray
    final_foreign_work: np.ndarray
    state: dict
    fhe_frames: int     # frames in the window's .fhe file so far


def _run_block(engine, chunk, state, lambdas, targets, steps, energy, energy_file, fhe_frames):
    """Run one block of one window and append its energies to the window's files."""
    engine.set_state(state)
    engine.update_config(steps=steps, lambdas=lambdas, energy=energy, energy_file=str(chunk),
                         foreign_lambdas=[t[1] for t in targets])
    result = engine.run()
    with open(energy_file, "ab") as joined:     # .en files have no header: append
        joined.write(Path(chunk).read_bytes())
    fhe_frames = append_fhe(chunk, energy_file, fhe_frames)
    return Block(result.q_energy, result.foreign_lambdas, result.foreign_work,
                 result.final_foreign_work, engine.get_state(), fhe_frames)


class _Serial:
    """One engine runs the windows in turn."""

    def __init__(self, engine):
        self.engine = engine
        self.config = engine.config

    def run_blocks(self, out, jobs):
        chunk = Path(out) / "hrex_block.en"
        blocks = [_run_block(self.engine, chunk, *job) for job in jobs]
        chunk.unlink(missing_ok=True)
        return blocks

    def finish(self):
        self.engine.update_config(energy=0, energy_file="", foreign_lambdas=())


def _serve(pipe, index, config, device, precision):
    """A worker process: one engine, runs what the parent sends until told to stop."""
    try:
        from .engine import Engine
        engine = Engine(config, device=device, precision=precision)
        pipe.send(("ok", None))
        while True:
            task, *args = pipe.recv()
            if task == "stop":
                return
            if task == "prepare":
                pipe.send(("ok", prepare(engine, *args)))
            else:   # "blocks"
                out, jobs = args
                chunk = Path(out) / f"hrex_block_{index}.en"
                blocks = [_run_block(engine, chunk, *job) for job in jobs]
                chunk.unlink(missing_ok=True)
                pipe.send(("ok", blocks))
    except EOFError:        # the parent is gone
        return
    except Exception:
        pipe.send(("error", traceback.format_exc()))


class Workers:
    """Several engines of the same system, each in its own process on the same GPU.

    hrex() takes this in place of an Engine and then runs the windows of a cycle at the
    same time: window k goes to worker k % n. With NVIDIA MPS the processes share the GPU.
    Use as a context manager, or call close().
    """

    def __init__(self, config, n, device=0, precision="spfp"):
        if n < 1:
            raise InputError("need at least one worker")
        self.config = config
        spawn = multiprocessing.get_context("spawn")    # CUDA does not survive a fork
        self._pipes, self._processes = [], []
        for index in range(n):
            ours, theirs = spawn.Pipe()
            process = spawn.Process(target=_serve, args=(theirs, index, config, device, precision), daemon=True)
            process.start()
            theirs.close()
            self._pipes.append(ours)
            self._processes.append(process)
        for pipe in self._pipes:
            self._receive(pipe)

    def __len__(self):
        return len(self._pipes)

    def __enter__(self):
        return self

    def __exit__(self, *error):
        self.close()

    def _receive(self, pipe):
        try:
            status, value = pipe.recv()
        except EOFError:
            raise InputError("a replica-exchange worker died") from None
        if status == "error":
            raise InputError("a replica-exchange worker failed:\n" + value)
        return value

    def prepare(self, lambda1, steps=5000):
        """prepare() in the first worker, from the state of its new engine."""
        self._pipes[0].send(("prepare", lambda1, steps))
        return self._receive(self._pipes[0])

    def run_blocks(self, out, jobs):
        n = len(self._pipes)
        for index, pipe in enumerate(self._pipes):
            pipe.send(("blocks", str(out), jobs[index::n]))
        blocks = [None] * len(jobs)
        for index, pipe in enumerate(self._pipes):
            blocks[index::n] = self._receive(pipe)
        return blocks

    def finish(self):
        pass

    def close(self):
        for pipe in self._pipes:
            try:
                pipe.send(("stop",))
            except (BrokenPipeError, OSError):
                pass
        for process in self._processes:
            process.join(timeout=30)
            if process.is_alive():
                process.terminate()
        self._pipes, self._processes = [], []


def hrex(engine, windows, out, cycles, block=500, energy=10, exchange=True,
         save_atoms=None, seed=0, work_dir=None):
    """Run `cycles` cycles of replica exchange over `windows`, starting from their states.

    engine: an Engine, which runs the windows in turn, or Workers, which run them at the
    same time.
    windows: objects with .lambdas and .state (e.g. the result of prepare), in the order
    of the ladder: neighbours in the list are the ones that swap.
    Each cycle runs `block` steps in every window and then tries to swap neighbours:
    pairs (0,1), (2,3), ... in even cycles and (1,2), (3,4), ... in odd ones.
    exchange=False runs the same schedule without swaps (the control).
    save_atoms: 0-based atom indices whose coordinates are kept after every block.

    Writes md_XXXX_YYYY.en per window into `out`, and next to each its foreign-work file
    (.en.fhe) for the window itself and its neighbours. Every block appends to these files;
    with work_dir they grow in a new folder below it and move to `out` at the end, which
    keeps the many small writes off a shared file system. Afterwards the engine has energy
    output and foreign work switched off.
    """
    if cycles < 1 or block < 1:
        raise InputError("need cycles >= 1 and block >= 1")
    if energy < 1 or block % energy:
        raise InputError("block must be a multiple of the energy interval")
    windows = list(windows)
    if len(windows) < 2:
        raise InputError("replica exchange needs at least two windows")

    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    n = len(windows)
    lambdas = [tuple(w.lambdas) for w in windows]
    states = [w.state for w in windows]
    runner = engine if isinstance(engine, Workers) else _Serial(engine)
    beta = 1.0 / (K_B * runner.config.temperature)
    rng = np.random.default_rng(seed)

    if len(set(lambdas)) != n:
        raise InputError("replica exchange needs windows with different lambdas")
    # each window computes the work to itself (zero; the .fhe file needs the column) and its neighbours
    targets = [tuple(lambdas[j] for j in (k - 1, k, k + 1) if 0 <= j < n) for k in range(n)]

    work = Path(tempfile.mkdtemp(prefix="hrex_", dir=work_dir)) if work_dir is not None else out
    files = [_energy_file(work, l) for l in lambdas]
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
        results = runner.run_blocks(work, [(states[k], lambdas[k], targets[k], block, energy, files[k], fhe_frames[k])
                                          for k in range(n)])
        for k, result in enumerate(results):
            states[k] = result.state
            fhe_frames[k] = result.fhe_frames
            parts[k].append(result.q_energy)
            works[k].append(result.foreign_work)
            if save_atoms is not None:
                coords[cycle, k] = np.asarray(states[k]["coords"]).reshape(-1, 3)[save_atoms]

        if exchange:
            for k in range(cycle % 2, n - 1, 2):
                a, b = results[k], results[k + 1]
                du = swap_energy(a.final_foreign_work[targets[k].index(lambdas[k + 1])],
                                 b.final_foreign_work[targets[k + 1].index(lambdas[k])])
                tried[k] += 1
                if du <= 0 or rng.random() < np.exp(-beta * du):
                    accepted[k] += 1
                    states[k], states[k + 1] = states[k + 1], states[k]
                    holder[[k, k + 1]] = holder[[k + 1, k]]
        replica[cycle + 1] = holder

    runner.finish()
    if work != out:
        for f in files:
            for path in (f, Path(str(f) + ".fhe")):
                if path.exists():
                    shutil.move(str(path), str(out / path.name))
        shutil.rmtree(work, ignore_errors=True)
    return HrexResult(lambdas, [np.concatenate(p) for p in parts], replica, tried, accepted,
                      coords, states, targets, [np.concatenate(w) for w in works])
