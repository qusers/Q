# ABOUTME: Reduces qdyn .corr logs of the charge-correction electrostatic observable into a
# ABOUTME: per-window ensemble curve and a lambda-integrated (TI) per-leg correction.

"""Charge-correction observable reduction for the analyze pipeline.

When a production run uses ``--correction-logging``, each window writes a
``.corr`` file whose rows are ``istep  lambda(1..nstates)  U_obs(1..nstates)``.
``U_obs(state)`` is the geometry-based electrostatic observable summed with that
FEP state's Q-atom charges. Because it is linear in the Q-charge,
``U_obs(state2) - U_obs(state1)`` is the derivative of the lambda-coupled
observable, so a leg's artifact free energy is the thermodynamic-integration
integral of its ensemble average over the sampled windows:

    dG_corr(leg) = integral_0^1 <U_obs(2) - U_obs(1)>_xi  d xi   (xi = lambda of state 2)

The empirical scale and the protein-minus-water combination are applied
post-hoc; this module only produces the unscaled per-leg integral and the
per-window curve. Two-state FEP only.
"""

from pathlib import Path

import numpy as np


def read_corr_file(path):
    """Parse a ``.corr`` file into (lambdas, u_obs).

    Returns a (nstates,) lambda array from the first row (constant within a
    window) and an (nframes, nstates) observable array. Raises ValueError on an
    empty file.
    """
    path = Path(path)
    if not path.read_text().strip():
        raise ValueError(f"empty .corr file: {path}")
    data = np.atleast_2d(np.loadtxt(path))
    ncols = data.shape[1]
    nstates = (ncols - 1) // 2
    if nstates < 1 or ncols != 1 + 2 * nstates:
        raise ValueError(f"unexpected column count {ncols} in {path}")
    return data[0, 1 : 1 + nstates], data[:, 1 + nstates :]


def reduce_windows(corr_files):
    """Pool replicates per window and reduce to a sorted per-window curve.

    Files sharing a basename (e.g. ``md_0500_0500.corr`` across replicate dirs)
    are pooled: all their frames are averaged together. Returns a list of dicts
    ``{"lambda2", "du", "n_frames"}`` sorted by ``lambda2`` (the state-2 lambda),
    where ``du`` is the pooled-frame mean of ``U_obs(2) - U_obs(1)``. Two-state
    FEP only.
    """
    groups = {}
    for f in corr_files:
        lambdas, u_obs = read_corr_file(f)
        if lambdas.shape[0] != 2:
            raise ValueError(f"charge correction expects 2-state FEP, got {lambdas.shape[0]} in {f}")
        key = Path(f).name
        if key not in groups:
            groups[key] = (float(lambdas[1]), [])
        groups[key][1].append(u_obs)
    windows = []
    for lambda2, frame_blocks in groups.values():
        frames = np.vstack(frame_blocks)
        avg = frames.mean(axis=0)
        windows.append({"lambda2": lambda2, "du": float(avg[1] - avg[0]), "n_frames": int(frames.shape[0])})
    windows.sort(key=lambda w: w["lambda2"])
    return windows


def integrate_windows(windows):
    """Trapezoidal TI integral of ``du`` over ``lambda2`` (kcal/mol, unscaled)."""
    xis = np.array([w["lambda2"] for w in windows], float)
    dus = np.array([w["du"] for w in windows], float)
    return float(np.trapz(dus, xis))


def collect_leg_corr(run_dirs):
    """Reduce every ``md_*.corr`` under a leg's replicate run directories.

    Args:
        run_dirs: iterable of replicate run-directory paths for one leg.

    Returns:
        (windows, integral): the per-window curve and its lambda-integral, or
        ``(None, None)`` if no ``.corr`` files are present (logging was off).
    """
    files = []
    for d in run_dirs:
        files.extend(sorted(Path(d).glob("md_*.corr")))
    if not files:
        return None, None
    windows = reduce_windows(files)
    return windows, integrate_windows(windows)
