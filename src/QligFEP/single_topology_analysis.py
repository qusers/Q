"""Analyse single-topology single-Hamiltonian FEP runs.

Single-Hamiltonian potentials are non-linear in lambda, so the standard qfep linear
reconstruction is invalid; instead every production window writes the true U(lambda')
at all coupling values (the .en.fl file). Only the q-atom nonbonded energy depends on
lambda, so a window's reduced potential at state l is qx(lambda_c[l])/kT and the
lambda-independent part of U cancels -- exactly what BAR/MBAR consume.

Pipeline:
  leg_dg          dG of one leg replicate (MBAR over its window .en.fl files)
  edge_ddg        per-replicate ddG = dG_protein - dG_water, with mean and SEM
  combine_network maximum-likelihood node free energies from the edge ddGs,
                  gauge-fixed to the experimental mean, vs experiment
"""
import glob
import json
import math
from pathlib import Path

import numpy as np

#: Boltzmann constant in kcal/mol (matches qfep's 0.592 at 298 K).
KB = 0.0019872041


def _load_fl(path, discard):
    """Return (lambda_c schedule, energy array [n_frames x n_lambda]) from a .en.fl."""
    sched, rows = None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith("#"):
            sched = [round(float(t), 4) for t in line.split()[2:]]
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        vals = [float(x) for x in parts[1:]]  # drop the step column
        if sched is not None and len(vals) == len(sched):
            rows.append(vals)
    return sched, np.array(rows[discard:], dtype=float)


def _windows(fl_dir, discard):
    """Load every window .en.fl in fl_dir, sorted by its coupling value."""
    files = sorted(glob.glob(str(Path(fl_dir) / "*.en.fl")))
    wins = []
    for f in files:
        lc = int(Path(f).name.split("_")[1].split(".")[0]) / 1000.0
        sched, energies = _load_fl(f, discard)
        if sched is None or energies.size == 0:
            continue
        wins.append((round(lc, 4), sched, energies))
    wins.sort(key=lambda w: w[0])
    return wins


def leg_dg_mbar(fl_dir, discard=40, temperature=298.0, max_iter=10000, tol=1e-7):
    """MBAR free energy (kcal/mol) of one leg replicate. Robust to the poor
    adjacent-window overlap of single-Hamiltonian since it uses all states jointly."""
    kt = KB * temperature
    wins = _windows(fl_dir, discard)
    if len(wins) < 2:
        return float("nan")
    sched = wins[0][1]
    col = {lc: j for j, lc in enumerate(sched)}
    K = len(sched)

    u_blocks, Nk = [], np.zeros(K)
    for lc, _, energies in wins:
        Nk[col[lc]] += len(energies)
        u_blocks.append(energies / kt)  # reduced potential at every state
    u = np.vstack(u_blocks)  # [N x K]
    with np.errstate(divide="ignore"):
        lnN = np.where(Nk > 0, np.log(Nk), -np.inf)

    f = np.zeros(K)
    for _ in range(max_iter):
        logD = _logsumexp(lnN[None, :] + f[None, :] - u, axis=1)  # [N]
        fnew = -_logsumexp(-u - logD[:, None], axis=0)  # [K]
        fnew -= fnew[0]
        if np.max(np.abs(fnew - f)) < tol:
            f = fnew
            break
        f = fnew
    return float((f[-1] - f[0]) * kt)


def leg_dg_bar(fl_dir, discard=40, temperature=298.0):
    """BAR free energy (kcal/mol) summed over adjacent windows, for comparison."""
    kt = KB * temperature
    wins = _windows(fl_dir, discard)
    if len(wins) < 2:
        return float("nan")
    col = {lc: j for j, lc in enumerate(wins[0][1])}
    total = 0.0
    for k in range(len(wins) - 1):
        lc_a, _, rows_a = wins[k]
        lc_b, _, rows_b = wins[k + 1]
        ja, jb = col[lc_a], col[lc_b]
        wf = (rows_a[:, jb] - rows_a[:, ja]) / kt
        wr = (rows_b[:, ja] - rows_b[:, jb]) / kt
        total += _bar(wf, wr) * kt
    return total


def _logsumexp(a, axis):
    m = np.max(a, axis=axis, keepdims=True)
    m = np.where(np.isfinite(m), m, 0.0)
    return (m.squeeze(axis) + np.log(np.sum(np.exp(a - m), axis=axis)))


def _bar(wf, wr, tol=1e-10):
    """Self-consistent BAR: solve Sum f(wf - df + M) = Sum f(wr + df - M)."""
    nF, nR = len(wf), len(wr)
    if nF == 0 or nR == 0:
        return float("nan")
    M = math.log(nF / nR)

    def fermi(x):
        return np.where(x > 700, 0.0, np.where(x < -700, 1.0, 1.0 / (1.0 + np.exp(np.clip(x, -700, 700)))))

    def D(df):
        return np.sum(fermi(wf - df + M)) - np.sum(fermi(wr + df - M))

    lo, hi = -1e4, 1e4
    if D(lo) * D(hi) > 0:
        return float("nan")
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        dm = D(mid)
        if abs(dm) < tol or (hi - lo) < 1e-10:
            return mid
        if D(lo) * dm < 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def _n_windows(fl_dir):
    """Number of window .en.fl files present in a replicate directory."""
    return len(glob.glob(str(Path(fl_dir) / "*.en.fl")))


def edge_ddg(root, edge, *, discard=40, temperature=298, estimator="mbar", n_windows=21):
    """Per-replicate ddG = dG_protein - dG_water for one edge, with mean and SEM.

    Expects setupFEP's layout: {root}/{2.protein,1.water}/FEP_{edge}/FEP1/{T}/{rep}/.
    Replicates are paired by directory name across the two legs. Only replicates whose
    leg has the full window schedule (n_windows .en.fl) are used -- MBAR over a
    truncated schedule misses an endpoint and gives the wrong dG; incomplete
    replicates are recorded, not silently averaged in.
    """
    leg_dg = leg_dg_mbar if estimator == "mbar" else leg_dg_bar
    legs, incomplete = {}, {}
    for leg, sub in (("protein", "2.protein"), ("water", "1.water")):
        base = Path(root) / sub / f"FEP_{edge}" / "FEP1" / str(temperature)
        legs[leg], incomplete[leg] = {}, []
        for d in sorted(glob.glob(str(base / "*"))):
            if not Path(d).is_dir():
                continue
            if _n_windows(d) < n_windows:
                incomplete[leg].append(Path(d).name)
                continue
            legs[leg][Path(d).name] = leg_dg(d, discard=discard, temperature=temperature)
    reps = sorted(set(legs["protein"]) & set(legs["water"]))
    ddgs = []
    for r in reps:
        dp, dw = legs["protein"][r], legs["water"][r]
        if math.isfinite(dp) and math.isfinite(dw):
            ddgs.append(dp - dw)
    ddgs = np.array(ddgs)
    n = len(ddgs)
    return {
        "edge": edge,
        "n_replicates": n,
        "ddg": float(np.mean(ddgs)) if n else float("nan"),
        "sem": float(np.std(ddgs, ddof=1) / math.sqrt(n)) if n > 1 else float("nan"),
        "per_replicate": {r: legs["protein"][r] - legs["water"][r] for r in reps},
        "dg_protein": legs["protein"],
        "dg_water": legs["water"],
    }


def combine_network(edges, node_exp):
    """Maximum-likelihood node free energies from edge ddGs (weighted least squares
    on g_j - g_i = ddG_ij, weight 1/sem^2), gauge-fixed so the predicted node mean
    equals the experimental node mean. Returns predicted node dG and metrics.

    Args:
        edges: list of dicts with keys 'from', 'to', 'ddg', 'sem'.
        node_exp: {node: experimental dG}.
    """
    nodes = sorted({e["from"] for e in edges} | {e["to"] for e in edges})
    idx = {n: i for i, n in enumerate(nodes)}
    K = len(nodes)

    rows, rhs, w = [], [], []
    for e in edges:
        if not math.isfinite(e["ddg"]):
            continue
        row = np.zeros(K)
        row[idx[e["to"]]] += 1.0
        row[idx[e["from"]]] -= 1.0
        rows.append(row)
        rhs.append(e["ddg"])
        sem = e.get("sem")
        w.append(1.0 / sem**2 if sem and math.isfinite(sem) and sem > 0 else 1.0)
    # gauge fix: pin sum of node free energies (soft constraint)
    rows.append(np.ones(K))
    rhs.append(0.0)
    w.append(1.0)

    A = np.array(rows)
    b = np.array(rhs)
    W = np.sqrt(np.array(w))
    g, *_ = np.linalg.lstsq(A * W[:, None], b * W, rcond=None)

    # shift predicted to the experimental mean (FEP gives relative free energies)
    common = [n for n in nodes if n in node_exp]
    if common:
        shift = np.mean([node_exp[n] for n in common]) - np.mean([g[idx[n]] for n in common])
        g = g + shift

    pred = {n: float(g[idx[n]]) for n in nodes}
    metrics = {}
    if common:
        pe = np.array([pred[n] for n in common])
        ee = np.array([node_exp[n] for n in common])
        metrics = {
            "n_nodes": len(common),
            "node_mae": float(np.mean(np.abs(pe - ee))),
            "node_rmse": float(np.sqrt(np.mean((pe - ee) ** 2))),
            "node_r": float(np.corrcoef(pe, ee)[0, 1]) if len(common) > 1 else float("nan"),
        }
    return {"node_dg_pred": pred, "metrics": metrics}


def analyze_run(root=".", mapping_json="mapping.json", *, discard=40, temperature=298, estimator="mbar"):
    """End-to-end: per-edge ddG (with replicate stats) + network combination vs
    the experimental dg_value/ddg_value in the mapping. Writes results JSON."""
    mapping = json.loads(Path(mapping_json).read_text())
    node_exp = {n: v.get("dg_value") for n, v in mapping.get("nodes", {}).items() if v.get("dg_value") is not None}

    edge_results = []
    for e in mapping["edges"]:
        edge = f"{e['from']}_{e['to']}"
        res = edge_ddg(root, edge, discard=discard, temperature=temperature, estimator=estimator)
        res.update({"from": e["from"], "to": e["to"], "ddg_exp": e.get("ddg_value")})
        edge_results.append(res)

    network = combine_network(edge_results, node_exp)
    # per-edge predicted-vs-experiment
    edge_pred = [
        {"edge": e["edge"], "ddg_pred": e["ddg"], "sem": e["sem"], "ddg_exp": e["ddg_exp"], "n": e["n_replicates"]}
        for e in edge_results
    ]
    valid = [(e["ddg_pred"], e["ddg_exp"]) for e in edge_pred
             if e["ddg_exp"] is not None and math.isfinite(e["ddg_pred"])]
    edge_metrics = {}
    if valid:
        p, x = np.array([v[0] for v in valid]), np.array([v[1] for v in valid])
        edge_metrics = {
            "n_edges": len(valid),
            "edge_mae": float(np.mean(np.abs(p - x))),
            "edge_rmse": float(np.sqrt(np.mean((p - x) ** 2))),
        }

    out = {
        "estimator": estimator,
        "discard": discard,
        "edges": edge_pred,
        "edge_metrics": edge_metrics,
        "network": network,
    }
    Path(root, "single_topology_results.json").write_text(json.dumps(out, indent=2, default=str))
    return out
