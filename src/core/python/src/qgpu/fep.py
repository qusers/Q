"""FEP scheduler driver: scan -> place -> production -> verify.

Design: src/core/adaptive_lambda_hrex_design.md. Only the scan exists so far.
"""

from dataclasses import dataclass
import numpy as np
import tempfile
from .errors import InputError
from pathlib import Path


@dataclass(frozen=True)
class ScanWindow:
    lambdas: tuple  # (lambda_1, lambda_2)
    q_energy: np.ndarray  # (n_frames, 2) kcal/mol; See RunResult.q_energy


@dataclass(frozen=True)
class Stability:
    stable: bool
    mean: float   # of g, kcal/mol
    error: float  # standard error of the mean
    drift: float  # mean of the second half minus mean of the first half


def block_error(x, segments=5):
    x = np.asarray(x, dtype=float)
    n = len(x) // segments
    means = x[len(x) - n * segments:].reshape(segments, n).mean(axis=1)
    return float(means.std(ddof=1) / np.sqrt(segments))


def stability(g, error_limit=1.0, segments=5):
    """Is the mean of g known well enough, or does the window need more steps?

    g: the energy gap q1 - q2 of every frame of the window so far (kcal/mol).
    Stable when the standard error of the mean is below error_limit, and the two halves
    differ by less than 2 * error_limit: a drift then shifts the mean by less than the limit.
    """
    g = np.asarray(g, dtype=float)
    if len(g) < 2 * segments:
        raise InputError(
            f"stability needs at least {2 * segments} frames, got {len(g)}")

    half = len(g) // 2
    error = block_error(g, segments)
    drift = float(g[half:].mean() - g[:half].mean())
    stable = error < error_limit and abs(drift) < 2 * error_limit
    return Stability(stable, float(g.mean()), error, drift)


@dataclass(frozen=True)
class Window:
    lambdas: tuple          # as placed, e.g. (0.0, 1.0)
    q_energy: np.ndarray    # (n_frames, 2) kcal/mol, all blocks
    state: dict             # Engine.get_state() at the end of the window
    steps: int
    stable: bool            # False: stopped at max_blocks
    check: Stability        # the last check


def run_window(engine, lambdas, energy_file, block=1000, energy=10,
               min_blocks=10, max_blocks=50, error_limit=1.0):
    """Run one lambda window from the engine's current state, block by block, until it is stable.

    Stable: at least min_blocks blocks, and stability() passed after each of the last two.
    The window ends at max_blocks in any case. Writes energy_file, the .en file of the
    whole window. Afterwards the engine has energy output switched off.
    """
    if min_blocks < 2 or max_blocks < min_blocks:
        raise InputError("need 2 <= min_blocks <= max_blocks")
    energy_file = Path(energy_file)
    if 0.0 in lambdas:
        raise InputError("run_window cannot use a lambda of exactly 0: "
                         "the engine reports no energy for that state")

    chunks, parts, passed = [], [], 0
    while True:
        # every run() overwrites its energy file: one chunk per block, joined below
        chunks.append(energy_file.with_name(
            f"{energy_file.name}.{len(chunks) + 1:03d}"))
        engine.update_config(steps=block, lambdas=tuple(lambdas), energy=energy,
                             energy_file=str(chunks[-1]))
        parts.append(engine.run().q_energy)
        q_energy = np.concatenate(parts)
        check = stability(q_energy[:, 0] - q_energy[:, 1], error_limit)
        passed = passed + 1 if check.stable else 0
        stable = len(chunks) >= min_blocks and passed >= 2
        if stable or len(chunks) == max_blocks:
            break

    state = engine.get_state()
    engine.update_config(energy=0, energy_file="")

    # .en files have no header and no frame at step 0: appending the chunks gives one file
    with open(energy_file, "wb") as joined:
        for chunk in chunks:
            joined.write(chunk.read_bytes())
            chunk.unlink()

    return Window(tuple(lambdas), q_energy, state, len(chunks) * block, stable, check)


K_B = 0.0019872041  # Boltzmann constant, kcal/(mol K)


@dataclass(frozen=True)
class Overlap:
    fits: bool
    forward: float  # ESS fraction of the frames of a, weighted for b (0..1)
    backward: float  # ESS fraction of the frames of b, weighted for a
    mismatch: float  # kcal/mol, see overlap()


def _reweight(q_energy, lambdas_from, lambdas_to, temperature):
    """ESS fraction and mean gap of the frames sampled at lambdas_from, weighted for lambdas_to."""
    q_energy = np.asarray(q_energy, dtype=float)
    step = np.asarray(lambdas_to, dtype=float) - \
        np.asarray(lambdas_from, dtype=float)
    # energy change of each frame, in kT
    du = q_energy @ step / (K_B * temperature)
    # shifted: the largest weight is 1
    w = np.exp(-(du - du.min()))
    ess = w.sum() ** 2 / (len(w) * (w ** 2).sum())
    gap = q_energy[:, 0] - q_energy[:, 1]
    return float(ess), float((w * gap).sum() / w.sum())


def overlap(a, b, temperature, target=0.5):
    """Do two windows overlap enough to connect them? a and b have .lambdas and .q_energy.

    The frames of a are weighted for the lambda of b, and the frames of b for a. The ESS
    fraction says how many of them still count: 1 = all, near 0 = a few frames carry the
    weight. The pair fits when both fractions reach the target.

    mismatch is a diagnostic, not part of the verdict: the mean gap that one window predicts
    for the other, minus the mean measured there (the larger of the two directions). It should
    be within the errors of the means; a large value means the windows do not sample
    what the weighting assumes.
    """
    if len(a.q_energy) == 0 or len(b.q_energy) == 0:
        raise InputError("overlap needs frames in both windows")
    forward, gap_at_b = _reweight(
        a.q_energy, a.lambdas, b.lambdas, temperature)
    backward, gap_at_a = _reweight(
        b.q_energy, b.lambdas, a.lambdas, temperature)
    miss_b = gap_at_b - float(np.mean(b.q_energy[:, 0] - b.q_energy[:, 1]))
    miss_a = gap_at_a - float(np.mean(a.q_energy[:, 0] - a.q_energy[:, 1]))
    mismatch = miss_b if abs(miss_b) >= abs(miss_a) else miss_a
    return Overlap(forward >= target and backward >= target, forward, backward, mismatch)


def _pair(l1):
    """(lambda_1, lambda_2), each rounded to 3 decimals as in qligfep."""
    return (round(l1, 3), round(1 - l1, 3))


def _energy_file(out, lambdas):
    return Path(out) / f"md_{round(lambdas[0] * 1000):04d}_{round(lambdas[1] * 1000):04d}.en"


def next_lambda(window, end, temperature, aim=0.6, max_step=0.1):
    """lambda_1 of the next window on the way from `window` to `end` (0.0 or 1.0).

    The largest step, in units of 0.001 and at most max_step, for which the frames of
    `window` still have an ESS fraction of `aim` at the new lambda. That is the forward half
    of overlap(), known before the next window is run; `aim` is above the target of
    overlap() to leave room for the backward half.
    """
    l1 = window.lambdas[0]
    direction = 1 if end > l1 else -1
    remaining = round(abs(end - l1), 3)
    best = 0.001
    for k in range(2, round(min(max_step, remaining) * 1000) + 1):
        ess, _ = _reweight(window.q_energy, window.lambdas,
                           _pair(l1 + direction * k / 1000), temperature)
        if ess < aim:
            break
        best = k / 1000
    if best >= remaining:
        return end
    if remaining - best < best / 2:     # no sliver before the end: two equal steps instead
        best = remaining / 2
    return round(l1 + direction * best, 3)


def walk(engine, out, ends=(0.999, 0.001), target=0.5, aim=0.6, max_step=0.1, max_windows=60, dg_limit=0.07, **window_settings):
    """Choose and run the lambda windows of one leg, from the engine's current state (after eq5).

    Runs lambda = 0.5, then walks to 1, then returns to the 0.5 end state and walks to 0.
    Each window runs until stable (run_window; window_settings go to it). Each new window is
    compared with the one before (overlap, with `target`); where a pair does not fit, a window
    is inserted at its midpoint, started from the inner one, until every pair fits. Windows
    that did not fit are kept. Once max_windows exist, no more are inserted.

    Writes md_XXXX_YYYY.en per window into `out`. Returns the windows sorted by lambda_1.
    """
    temperature = engine.config.temperature
    windows = {}

    def run(l1, start):
        step = max_step if start is None else abs(l1 - start.lambdas[0])
        if start is not None:
            engine.set_state(start.state)
        lambdas = _pair(l1)
        windows[lambdas[0]] = run_window(engine, lambdas, _energy_file(out, lambdas),
                                         error_limit=dg_limit / step, **window_settings)
        return windows[lambdas[0]]

    def connect(inner, outer):
        if overlap(inner, outer, temperature, target).fits or len(windows) >= max_windows:
            return
        middle = round((inner.lambdas[0] + outer.lambdas[0]) / 2, 3)
        if middle in (inner.lambdas[0], outer.lambdas[0]):      # already 0.001 apart
            return
        between = run(middle, inner)
        connect(inner, between)
        connect(between, outer)

    midpoint = run(0.5, None)
    for end in ends:
        current = midpoint
        while current.lambdas[0] != end:
            trial = run(next_lambda(
                current, end, temperature, aim, max_step), current)
            connect(current, trial)
            current = trial

    return [windows[l1] for l1 in sorted(windows)]


def pairs(windows, temperature, target=0.5):
    """(a, b, Overlap) for every pair of neighbouring windows, in order of lambda_1."""
    return [(a, b, overlap(a, b, temperature, target)) for a, b in zip(windows, windows[1:])]
