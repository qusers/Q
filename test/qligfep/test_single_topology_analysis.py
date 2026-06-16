"""Tests for the single-topology single-Hamiltonian analysis.

The network combination is exercised on synthetic edges (deterministic, no engine).
The foreign-lambda MBAR/BAR estimators are validated against the committed prototype
.en.fl files that produced the published ddG = -0.19 (matching experiment), when those
files are available; otherwise that check is skipped.
"""
from pathlib import Path

import pytest

from QligFEP.single_topology_analysis import combine_network, leg_dg_bar, leg_dg_mbar

PROTO = Path("/Users/davidararipe/projects/Q-softcore-verify/experiments/single_h")


def test_combine_network_recovers_node_free_energies():
    """A consistent set of edge ddGs recovers the node free energies (up to the
    experimental-mean gauge) by weighted least squares."""
    true_dg = {"A": -10.0, "B": -9.0, "C": -11.0}
    edges = [
        {"from": "A", "to": "B", "ddg": 1.0, "sem": 0.1},
        {"from": "B", "to": "C", "ddg": -2.0, "sem": 0.1},
        {"from": "A", "to": "C", "ddg": -1.0, "sem": 0.1},
    ]
    out = combine_network(edges, true_dg)
    for node, val in true_dg.items():
        assert out["node_dg_pred"][node] == pytest.approx(val, abs=1e-6)
    assert out["metrics"]["node_rmse"] == pytest.approx(0.0, abs=1e-6)


def test_combine_network_least_squares_balances_inconsistent_cycle():
    """An inconsistent cycle (A->B->C->A sums to 0.3, not 0) is balanced by the
    weighted least-squares fit rather than trusting any single edge."""
    edges = [
        {"from": "A", "to": "B", "ddg": 1.0, "sem": 0.1},
        {"from": "B", "to": "C", "ddg": 1.0, "sem": 0.1},
        {"from": "C", "to": "A", "ddg": -1.7, "sem": 0.1},
    ]
    out = combine_network(edges, {"A": 0.0})
    g = out["node_dg_pred"]
    # residual spread the 0.3 kcal cycle error evenly (0.1 per edge)
    assert g["B"] - g["A"] == pytest.approx(0.9, abs=1e-6)
    assert g["C"] - g["B"] == pytest.approx(0.9, abs=1e-6)


@pytest.mark.slow
@pytest.mark.skipif(not (PROTO / "singletop").exists(), reason="prototype .en.fl not present")
def test_mbar_reproduces_prototype_ddg():
    """leg_dg_mbar over the prototype .en.fl reproduces the validated single-topology
    ddG (protein ~29.6, water ~29.8, ddG ~ -0.2 matching experiment -0.16)."""
    p = leg_dg_mbar(str(PROTO / "singletop"), discard=20)
    w = leg_dg_mbar(str(PROTO / "singletop_water"), discard=20)
    assert p == pytest.approx(29.6, abs=0.5)
    assert w == pytest.approx(29.8, abs=0.5)
    assert (p - w) == pytest.approx(-0.2, abs=0.4)


@pytest.mark.slow
@pytest.mark.skipif(not (PROTO / "singletop").exists(), reason="prototype .en.fl not present")
def test_bar_and_mbar_agree_within_replicate_noise():
    """BAR and MBAR per-leg dG agree within the single-replicate spread (~1 kcal)."""
    for leg in ("singletop", "singletop_water"):
        assert leg_dg_bar(str(PROTO / leg), discard=20) == pytest.approx(
            leg_dg_mbar(str(PROTO / leg), discard=20), abs=1.0
        )
