"""Tests for the shared FEP plotting/result-shaping helpers.

These run without the Q binaries: they build a small per-edge results table and exercise
prepare_df + create_ddG_plot, which both the equilibrium and NEQ analyzers reuse.
"""

import matplotlib

matplotlib.use("Agg")  # headless backend before pyplot is imported by the module under test

import numpy as np
import pytest

from QligFEP.analysis_plotting import bootstrap_statistic, create_ddG_plot, prepare_df


def _edges():
    # A handful of edges so the bootstrap statistics have enough data.
    return {
        "edges": [
            {"from": "a", "to": "b", "ddg_value": 1.0, "Q_ddG_avg": 0.6, "Q_ddG_sem": 0.1},
            {"from": "b", "to": "c", "ddg_value": -2.0, "Q_ddG_avg": -1.4, "Q_ddG_sem": 0.2},
            {"from": "c", "to": "d", "ddg_value": 0.5, "Q_ddG_avg": 0.9, "Q_ddG_sem": 0.15},
            {"from": "d", "to": "e", "ddg_value": 2.5, "Q_ddG_avg": 2.1, "Q_ddG_sem": 0.3},
            {"from": "e", "to": "f", "ddg_value": -1.0, "Q_ddG_avg": -0.7, "Q_ddG_sem": 0.1},
        ]
    }


def test_prepare_df_adds_fep_name_and_residual():
    df = prepare_df(_edges())
    assert {"fep_name", "ddg_value", "residual"}.issubset(df.columns)
    assert "FEP_a_b" in df["fep_name"].values
    # residual = experimental - calculated
    row = df[df["fep_name"] == "FEP_a_b"].iloc[0]
    assert row["residual"] == 1.0 - 0.6


def test_prepare_df_without_experimental_data_skips_residual():
    df = prepare_df(_edges(), experimental_data=False)
    assert "fep_name" in df.columns
    assert "residual" not in df.columns


def test_create_ddG_plot_returns_fig_and_axis():
    df = prepare_df(_edges())
    fig, ax = create_ddG_plot(df, target_name="test")
    assert fig is not None and ax is not None
    assert ax.get_xlabel() and ax.get_ylabel()


def test_create_ddG_plot_axis_bounds_include_all_points():
    # Anti-correlated data: the elementwise sum of experimental+calculated collapses
    # toward zero, so axis limits derived from that sum would clip the real points off
    # the plot. The limits must instead span the combined range of both value sets.
    edges = {
        "edges": [
            {"from": "a", "to": "b", "ddg_value": 3.0, "Q_ddG_avg": -2.8, "Q_ddG_sem": 0.1},
            {"from": "b", "to": "c", "ddg_value": -3.0, "Q_ddG_avg": 2.9, "Q_ddG_sem": 0.1},
            {"from": "c", "to": "d", "ddg_value": 2.5, "Q_ddG_avg": -2.0, "Q_ddG_sem": 0.1},
            {"from": "d", "to": "e", "ddg_value": -2.5, "Q_ddG_avg": 2.2, "Q_ddG_sem": 0.1},
            {"from": "e", "to": "f", "ddg_value": 1.0, "Q_ddG_avg": -0.5, "Q_ddG_sem": 0.1},
        ]
    }
    df = prepare_df(edges)
    fig, ax = create_ddG_plot(df, target_name="test")
    xlo, xhi = ax.get_xlim()
    ylo, yhi = ax.get_ylim()
    for exp, calc in zip(df["ddg_value"], df["Q_ddG_avg"]):
        assert xlo <= exp <= xhi, f"exp {exp} outside x-range [{xlo}, {xhi}]"
        assert ylo <= calc <= yhi, f"calc {calc} outside y-range [{ylo}, {yhi}]"


def test_create_ddG_plot_saves_into_directory(tmp_path):
    # savefig=True with a directory output_path must write <target>_ddG_plot.png there.
    df = prepare_df(_edges())
    create_ddG_plot(df, target_name="mytarget", output_path=tmp_path, savefig=True)
    assert (tmp_path / "mytarget_ddG_plot.png").exists()


def test_bootstrap_statistic_point_estimate_matches_closed_form():
    a = np.array([0.6, -1.4, 0.9, 2.1, -0.7])
    b = np.array([1.0, -2.0, 0.5, 2.5, -1.0])
    rmse = bootstrap_statistic(a, b, "RMSE")
    mue = bootstrap_statistic(a, b, "MUE")
    assert rmse["mle"] == np.sqrt(np.mean((a - b) ** 2))
    assert mue["mle"] == np.mean(np.abs(a - b))
    # the confidence interval is ordered and brackets the point estimate
    assert rmse["low"] <= rmse["mle"] <= rmse["high"]


def test_bootstrap_statistic_identity_data():
    a = np.array([0.1, 0.5, 1.2, -0.3, 2.0])
    assert bootstrap_statistic(a, a, "RMSE")["mle"] == 0.0
    assert bootstrap_statistic(a, a, "MUE")["mle"] == 0.0
    # a perfectly rank-correlated series has Kendall tau == 1
    assert bootstrap_statistic(a, a, "KTAU")["mle"] == pytest.approx(1.0)
