# ABOUTME: Reduces qdyn .corr logs of the charge-correction electrostatic observable into a
# ABOUTME: per-window ensemble curve and a lambda-integrated (TI) per-leg correction.

"""Charge-correction observable reduction for the analyze pipeline.

When a production run uses ``--correction-logging``, each window writes a
``.corr`` file whose rows are
``istep  lambda(1..nstates)  U_obs(1..nstates)  phi_cog``.

Two observables, both with the eps(r)=r screened kernel:

* ``U_obs(state2) - U_obs(state1)`` is the full per-atom charge-distribution
  interaction change (real partial charges over the actual Q-atoms);
* ``phi_cog`` is the screened protein potential at the ligand centroid, so
  ``dq * <phi_cog>`` is the net-charge monopole term, the on-the-fly analog of
  the validated static screened-Coulomb model.

A leg's artifact free energy is the thermodynamic-integration integral of the
ensemble average over the sampled windows (xi = lambda of state 2):

    dG_corr(leg)   = integral_0^1 <U_obs(2) - U_obs(1)>_xi  d xi
    dG_phi(leg)    = integral_0^1 <phi_cog>_xi  d xi   (monopole: multiply by dq)

The empirical scale, the choice of model, and the protein-minus-water
combination are applied post-hoc; this module only produces the unscaled per-leg
integrals and the per-window curve. Two-state FEP only.
"""

from pathlib import Path

import numpy as np


def read_corr_file(path):
    """Parse a ``.corr`` file into (lambdas, u_obs, phi).

    Row layout: ``step  lambda(1..n)  U_obs(1..n)  phi_cog`` (2 + 2n columns).
    Returns a (nstates,) lambda array from the first row, an (nframes, nstates)
    observable array, and an (nframes,) phi_cog array. Raises ValueError on an
    empty file.
    """
    path = Path(path)
    if not path.read_text().strip():
        raise ValueError(f"empty .corr file: {path}")
    data = np.atleast_2d(np.loadtxt(path))
    ncols = data.shape[1]
    nstates = (ncols - 2) // 2
    if nstates < 1 or ncols != 2 + 2 * nstates:
        raise ValueError(f"unexpected column count {ncols} in {path}")
    lambdas = data[0, 1 : 1 + nstates]
    u_obs = data[:, 1 + nstates : 1 + 2 * nstates]
    phi = data[:, -1]
    return lambdas, u_obs, phi


def reduce_windows(corr_files):
    """Pool replicates per window and reduce to a sorted per-window curve.

    Files sharing a basename (e.g. ``md_0500_0500.corr`` across replicate dirs)
    are pooled: all their frames are averaged together. Returns a list of dicts
    ``{"lambda2", "du", "phi", "n_frames"}`` sorted by ``lambda2`` (the state-2
    lambda), where ``du`` is the pooled-frame mean of ``U_obs(2) - U_obs(1)`` and
    ``phi`` the pooled-frame mean of ``phi_cog``. Two-state FEP only.
    """
    groups = {}
    for f in corr_files:
        lambdas, u_obs, phi = read_corr_file(f)
        if lambdas.shape[0] != 2:
            raise ValueError(f"charge correction expects 2-state FEP, got {lambdas.shape[0]} in {f}")
        key = Path(f).name
        if key not in groups:
            groups[key] = (float(lambdas[1]), [], [])
        groups[key][1].append(u_obs)
        groups[key][2].append(phi)
    windows = []
    for lambda2, u_blocks, phi_blocks in groups.values():
        u = np.vstack(u_blocks)
        p = np.concatenate(phi_blocks)
        avg = u.mean(axis=0)
        windows.append(
            {
                "lambda2": lambda2,
                "du": float(avg[1] - avg[0]),
                "phi": float(p.mean()),
                "n_frames": int(u.shape[0]),
            }
        )
    windows.sort(key=lambda w: w["lambda2"])
    return windows


def integrate_windows(windows, key="du"):
    """Trapezoidal TI integral of ``windows[key]`` over ``lambda2`` (unscaled)."""
    xis = np.array([w["lambda2"] for w in windows], float)
    vals = np.array([w[key] for w in windows], float)
    return float(np.trapz(vals, xis))


def collect_leg_corr(run_dirs):
    """Reduce every ``md_*.corr`` under a leg's replicate run directories.

    Args:
        run_dirs: iterable of replicate run-directory paths for one leg.

    Returns:
        (windows, u_integral, phi_integral): the per-window curve and the
        lambda-integrals of the full-distribution observable and of phi_cog, or
        ``(None, None, None)`` if no ``.corr`` files are present.
    """
    files = []
    for d in run_dirs:
        files.extend(sorted(Path(d).glob("md_*.corr")))
    if not files:
        return None, None, None
    windows = reduce_windows(files)
    return windows, integrate_windows(windows, "du"), integrate_windows(windows, "phi")
