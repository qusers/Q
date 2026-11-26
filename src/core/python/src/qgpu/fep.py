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
    # Foreign work of every frame, see work(): lambda_2 of the windows it was computed for
    # (sorted, this window among them) and U(that window) - U(this window), kcal/mol.
    foreign_lambdas: tuple = ()
    foreign_work: np.ndarray | None = None  # (n_frames, n_targets)


def _soft_coordinate(lambda2):
    """A coordinate along lambda_2 in which soft-core energies change evenly.

    The soft-core radius of a state goes as (1 - lambda)^(1/6): it changes fastest in lambda
    where a state is nearly fully on, i.e. at both ends. Equal steps of this coordinate are
    equal steps of the radii there, and nearly equal steps of lambda in the middle.
    """
    lambda2 = np.clip(lambda2, 0.0, 1.0)
    return lambda2 ** (1.0 / 6.0) - (1.0 - lambda2) ** (1.0 / 6.0)


def _from_soft_coordinate(x):
    lo, hi = 0.0, 1.0       # _soft_coordinate rises from -1 to 1
    for _ in range(60):
        middle = 0.5 * (lo + hi)
        lo, hi = (middle, hi) if _soft_coordinate(middle) < x else (lo, middle)
    return 0.5 * (lo + hi)


def foreign_targets(lambda2, reach, per_side=3):
    """lambda_2 of the windows a window at `lambda2` computes foreign work for.

    The window itself and per_side on each side, out to `reach` and not past 0 or 1; sorted.
    They are evenly spaced in _soft_coordinate, so closer together towards the ends.
    work() is exact at these and interpolates between them.
    """
    targets = {float(lambda2)}
    x0 = _soft_coordinate(lambda2)
    for end in (max(0.0, lambda2 - reach), min(1.0, lambda2 + reach)):
        x_end = _soft_coordinate(end)
        for j in range(1, per_side + 1):
            t = end if j == per_side else _from_soft_coordinate(x0 + (x_end - x0) * j / per_side)
            if abs(t - lambda2) > 1e-9:
                targets.add(round(t, 9))
    return tuple(sorted(targets))


def append_fhe(chunk_energy_file, energy_file, frames):
    """Append the foreign-work file of one run() to that of the whole window.

    Both are named after their energy file (<energy file>.fhe). `frames` is the number of
    frames already in the window's file; 0 starts it with the header. Returns the new number.
    The chunk is removed.
    """
    chunk = Path(str(chunk_energy_file) + ".fhe")
    if not chunk.exists():      # the run had no foreign targets
        return frames
    lines = chunk.read_text().splitlines()
    chunk.unlink()
    with open(str(energy_file) + ".fhe", "w" if frames == 0 else "a") as joined:
        for line in lines:
            if line.startswith("#"):
                if frames == 0:
                    joined.write(line + "\n")
            elif line.strip():
                frames += 1
                # the first column is the frame number, 10 wide, counted per run()
                joined.write(f"{frames:10d}{line[10:]}\n")
    return frames


def run_window(engine, lambdas, energy_file, block=1000, energy=10,
               min_blocks=10, max_blocks=50, error_limit=1.0, reach=0.1):
    """Run one lambda window from the engine's current state, block by block, until it is stable.

    Stable: at least min_blocks blocks, and stability() passed after each of the last two.
    The window ends at max_blocks in any case. Writes energy_file, the .en file of the
    whole window, and its foreign-work file (energy_file + ".fhe") for the windows of
    foreign_targets(lambda_2, reach). Afterwards the engine has energy output and foreign
    work switched off.
    """
    if min_blocks < 2 or max_blocks < min_blocks:
        raise InputError("need 2 <= min_blocks <= max_blocks")
    energy_file = Path(energy_file)
    if 0.0 in lambdas:
        raise InputError("run_window cannot use a lambda of exactly 0: "
                         "the engine reports no energy for that state")

    targets = foreign_targets(lambdas[1], reach)
    chunks, parts, works, passed = [], [], [], 0
    while True:
        # every run() overwrites its energy file: one chunk per block, joined below
        chunks.append(energy_file.with_name(
            f"{energy_file.name}.{len(chunks) + 1:03d}"))
        engine.update_config(steps=block, lambdas=tuple(lambdas), energy=energy,
                             energy_file=str(chunks[-1]), foreign_lambdas=targets)
        result = engine.run()
        parts.append(result.q_energy)
        works.append(result.foreign_work)
        q_energy = np.concatenate(parts)
        check = stability(q_energy[:, 0] - q_energy[:, 1], error_limit)
        passed = passed + 1 if check.stable else 0
        stable = len(chunks) >= min_blocks and passed >= 2
        if stable or len(chunks) == max_blocks:
            break

    state = engine.get_state()
    engine.update_config(energy=0, energy_file="", foreign_lambdas=())

    # .en files have no header and no frame at step 0: appending the chunks gives one file
    frames = 0
    with open(energy_file, "wb") as joined:
        for chunk in chunks:
            joined.write(chunk.read_bytes())
            frames = append_fhe(chunk, energy_file, frames)
            chunk.unlink()

    return Window(tuple(lambdas), q_energy, state, len(chunks) * block, stable, check,
                  targets, np.concatenate(works))


K_B = 0.0019872041  # Boltzmann constant, kcal/(mol K)


@dataclass(frozen=True)
class Overlap:
    fits: bool
    forward: float  # ESS fraction of the frames of a, weighted for b (0..1)
    backward: float  # ESS fraction of the frames of b, weighted for a
    mismatch: float  # kcal/mol, see overlap()


def work(window, lambdas_to):
    """U(lambdas_to) - U(window.lambdas) for every frame of the window, kcal/mol.

    The part that is linear in lambda is  q_energy @ (lambdas_to - window.lambdas), exact for
    any lambdas_to. With soft-core the Q energies themselves change with lambda; that part is
    what window.foreign_work holds beyond the linear one. It is exact at the windows in
    window.foreign_lambdas and interpolated between them; lambdas_to outside their
    range is an error. A window without foreign work gets the linear part only, which is
    the whole answer only when there is no soft-core.
    """
    q_energy = np.asarray(window.q_energy, dtype=float)
    source = np.asarray(window.lambdas, dtype=float)
    linear = q_energy @ (np.asarray(lambdas_to, dtype=float) - source)

    nodes = np.asarray(getattr(window, "foreign_lambdas", ()), dtype=float)
    if nodes.size == 0:
        return linear
    foreign = np.asarray(window.foreign_work, dtype=float)
    if foreign.shape != (len(q_energy), nodes.size):
        raise InputError("window.foreign_work must have one row per frame and one column per foreign lambda")

    t = float(lambdas_to[1])
    if not nodes[0] - 1e-9 <= t <= nodes[-1] + 1e-9:
        raise InputError(f"lambda_2 = {t:g} is outside the foreign windows of this window "
                         f"({nodes[0]:g} to {nodes[-1]:g}); run it with a larger reach")
    nearest = np.argsort(np.abs(nodes - t), kind="stable")[:3]
    if abs(nodes[nearest[0]] - t) <= 1e-9:
        return foreign[:, nearest[0]]

    # The soft-core part (foreign work minus its linear part) at the three windows nearest
    # to lambdas_to, interpolated with the parabola through them, in the coordinate in which
    # it changes evenly (see _soft_coordinate).
    x, x_nodes = _soft_coordinate(t), _soft_coordinate(nodes)
    result = linear.copy()
    for k in nearest:
        weight = np.prod([(x - x_nodes[j]) / (x_nodes[k] - x_nodes[j]) for j in nearest if j != k])
        result += weight * (foreign[:, k] - q_energy @ (np.array([1.0 - nodes[k], nodes[k]]) - source))
    return result


def _reweight(window, lambdas_to, temperature):
    """ESS fraction and mean gap of the frames of `window`, weighted for lambdas_to."""
    q_energy = np.asarray(window.q_energy, dtype=float)
    # energy change of each frame, in kT
    du = work(window, lambdas_to) / (K_B * temperature)
    # shifted: the largest weight is 1
    w = np.exp(-(du - du.min()))
    ess = w.sum() ** 2 / (len(w) * (w ** 2).sum())
    gap = q_energy[:, 0] - q_energy[:, 1]
    return float(ess), float((w * gap).sum() / w.sum())


def overlap(a, b, temperature, target=0.5):
    """Do two windows overlap enough to connect them? a and b are Window objects, or anything
    with .lambdas and .q_energy (and, with soft-core, .foreign_lambdas and .foreign_work: see work()).

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
    forward, gap_at_b = _reweight(a, b.lambdas, temperature)
    backward, gap_at_a = _reweight(b, a.lambdas, temperature)
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
        ess, _ = _reweight(window, _pair(l1 + direction * k / 1000), temperature)
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
        # reach = max_step: every window this one is compared with is within its foreign windows
        windows[lambdas[0]] = run_window(engine, lambdas, _energy_file(out, lambdas),
                                         error_limit=dg_limit / step, reach=max_step, **window_settings)
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
