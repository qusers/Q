"""Tests for QligFEP quality-control analysis."""

import pytest

from QligFEP.analysis_qc import (
    analyze_fep_edge,
    analyze_fep_system,
    replicate_statistics,
)


def test_replicate_statistics():
    """Statistics should be calculated across valid replicates."""
    energies = [-4.21, -4.35, -4.18, -4.30]

    result = replicate_statistics(energies)

    assert result["n_replicates"] == 4
    assert result["n_valid"] == 4
    assert result["n_failed"] == 0
    assert result["mean"] == pytest.approx(-4.26)
    assert result["std"] > 0
    assert result["sem"] > 0
    assert result["range"] == pytest.approx(0.17)


def test_replicate_statistics_ignores_nan():
    """NaN values should count as failed and be excluded from statistics."""
    energies = [-4.21, float("nan"), -4.18, -4.30]

    result = replicate_statistics(energies)

    assert result["n_replicates"] == 4
    assert result["n_valid"] == 3
    assert result["n_failed"] == 1
    assert result["mean"] == pytest.approx((-4.21 - 4.18 - 4.30) / 3)


def test_replicate_statistics_all_failed():
    """All failed replicates should return unavailable statistics."""
    energies = [float("nan"), float("nan")]

    result = replicate_statistics(energies)

    assert result["n_replicates"] == 2
    assert result["n_valid"] == 0
    assert result["n_failed"] == 2
    assert result["mean"] is None
    assert result["std"] is None
    assert result["sem"] is None
    assert result["range"] is None


def test_analyze_fep_edge():
    """QC should extract replicate energies from FepReader data."""
    data = {
        "2.protein": {
            "FEP_lig1_lig2": {
                "FEP_result": {
                    "dGbar": {
                        "energies": [-4.21, -4.35, -4.18, -4.30],
                        "avg": -4.26,
                        "sem": 0.03,
                        "std": 0.07,
                    }
                }
            }
        }
    }

    result = analyze_fep_edge(
        data,
        system="2.protein",
        fep="FEP_lig1_lig2",
    )

    assert result["system"] == "2.protein"
    assert result["fep"] == "FEP_lig1_lig2"
    assert result["method"] == "dGbar"
    assert result["n_replicates"] == 4
    assert result["n_valid"] == 4
    assert result["n_failed"] == 0
    assert result["mean"] == pytest.approx(-4.26)
    assert result["range"] == pytest.approx(0.17)


def test_analyze_fep_system():
    """QC should analyze all FEP edges in a system."""
    data = {
        "2.protein": {
            "FEP_lig1_lig2": {
                "FEP_result": {
                    "dGbar": {
                        "energies": [-4.21, -4.35, -4.18],
                    }
                }
            },
            "FEP_lig2_lig3": {
                "FEP_result": {
                    "dGbar": {
                        "energies": [-2.10, -2.20, float("nan")],
                    }
                }
            },
        }
    }

    results = analyze_fep_system(
        data,
        system="2.protein",
    )

    assert len(results) == 2

    assert results[0]["fep"] == "FEP_lig1_lig2"
    assert results[0]["n_replicates"] == 3
    assert results[0]["n_valid"] == 3
    assert results[0]["n_failed"] == 0

    assert results[1]["fep"] == "FEP_lig2_lig3"
    assert results[1]["n_replicates"] == 3
    assert results[1]["n_valid"] == 2
    assert results[1]["n_failed"] == 1
