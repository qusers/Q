"""Tests for QligFEP quality-control analysis."""

import pytest

from QligFEP.analysis_qc import (
    analyze_fep_edge,
    analyze_fep_system,
    replicate_statistics,
    check_fep_pair_consistency,
    analyze_ddg_edge,
    summarize_fep_edge_qc,
    summarize_fep_system_qc,
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


def test_check_fep_pair_consistency():
    """Matching protein and water FEP settings should pass consistency QC."""
    data = {
        "1.water": {
            "FEP_lig1_lig2": {
                "fep_stage": "production",
                "temperature": 298,
                "lambda_sum": 1.0,
            }
        },
        "2.protein": {
            "FEP_lig1_lig2": {
                "fep_stage": "production",
                "temperature": 298,
                "lambda_sum": 1.0,
            }
        },
    }

    result = check_fep_pair_consistency(
        data,
        fep="FEP_lig1_lig2",
    )

    assert result["fep_stage_match"] is True
    assert result["temperature_match"] is True
    assert result["lambda_sum_match"] is True
    assert result["consistent"] is True


def test_check_fep_pair_consistency_detects_mismatch():
    """Protein/water setting mismatches should fail consistency QC."""
    data = {
        "1.water": {
            "FEP_lig1_lig2": {
                "fep_stage": "production",
                "temperature": 298,
                "lambda_sum": 1.0,
            }
        },
        "2.protein": {
            "FEP_lig1_lig2": {
                "fep_stage": "production",
                "temperature": 310,
                "lambda_sum": 1.0,
            }
        },
    }

    result = check_fep_pair_consistency(
        data,
        fep="FEP_lig1_lig2",
    )

    assert result["temperature_match"] is False
    assert result["consistent"] is False


def test_analyze_ddg_edge():
    """QC should extract an existing ddG result."""
    data = {
        "result": {
            "ddGbar": {
                "FEP_lig1_lig2": {
                    "ddGbar_avg": 1.25,
                    "ddGbar_sem": 0.18,
                    "ddGbar_std": 0.36,
                    "from": "lig1",
                    "to": "lig2",
                }
            }
        }
    }

    result = analyze_ddg_edge(
        data,
        fep="FEP_lig1_lig2",
    )

    assert result["fep"] == "FEP_lig1_lig2"
    assert result["method"] == "ddGbar"
    assert result["from"] == "lig1"
    assert result["to"] == "lig2"
    assert result["avg"] == pytest.approx(1.25)
    assert result["sem"] == pytest.approx(0.18)
    assert result["std"] == pytest.approx(0.36)


def test_summarize_fep_edge_qc():
    """QC summary should combine ddG, replicate, and consistency information."""
    data = {
        "1.water": {
            "FEP_lig1_lig2": {
                "fep_stage": "production",
                "temperature": 298,
                "lambda_sum": 1.0,
                "FEP_result": {
                    "dGbar": {
                        "energies": [-2.0, -2.1, -1.9],
                    }
                },
            }
        },
        "2.protein": {
            "FEP_lig1_lig2": {
                "fep_stage": "production",
                "temperature": 298,
                "lambda_sum": 1.0,
                "FEP_result": {
                    "dGbar": {
                        "energies": [-3.2, -3.4, -3.3],
                    }
                },
            }
        },
        "result": {
            "ddGbar": {
                "FEP_lig1_lig2": {
                    "ddGbar_avg": -1.3,
                    "ddGbar_sem": 0.12,
                    "ddGbar_std": 0.21,
                    "from": "lig1",
                    "to": "lig2",
                }
            }
        },
    }

    result = summarize_fep_edge_qc(
        data,
        fep="FEP_lig1_lig2",
    )

    assert result["fep"] == "FEP_lig1_lig2"
    assert result["from"] == "lig1"
    assert result["to"] == "lig2"

    assert result["ddg"] == pytest.approx(-1.3)
    assert result["ddg_sem"] == pytest.approx(0.12)

    assert result["protein_n_valid"] == 3
    assert result["protein_n_failed"] == 0

    assert result["water_n_valid"] == 3
    assert result["water_n_failed"] == 0

    assert result["consistent"] is True


def test_summarize_fep_system_qc():
    """QC summary should include every matching FEP edge."""
    data = {
        "1.water": {
            "FEP_lig1_lig2": {
                "fep_stage": "production",
                "temperature": 298,
                "lambda_sum": 1.0,
                "FEP_result": {
                    "dGbar": {
                        "energies": [-2.0, -2.1, -1.9],
                    }
                },
            },
            "FEP_lig2_lig3": {
                "fep_stage": "production",
                "temperature": 298,
                "lambda_sum": 1.0,
                "FEP_result": {
                    "dGbar": {
                        "energies": [-1.0, -1.2, -1.1],
                    }
                },
            },
        },
        "2.protein": {
            "FEP_lig1_lig2": {
                "fep_stage": "production",
                "temperature": 298,
                "lambda_sum": 1.0,
                "FEP_result": {
                    "dGbar": {
                        "energies": [-3.2, -3.4, -3.3],
                    }
                },
            },
            "FEP_lig2_lig3": {
                "fep_stage": "production",
                "temperature": 298,
                "lambda_sum": 1.0,
                "FEP_result": {
                    "dGbar": {
                        "energies": [-2.0, -2.2, -2.1],
                    }
                },
            },
        },
        "result": {
            "ddGbar": {
                "FEP_lig1_lig2": {
                    "ddGbar_avg": -1.3,
                    "ddGbar_sem": 0.12,
                    "ddGbar_std": 0.21,
                    "from": "lig1",
                    "to": "lig2",
                },
                "FEP_lig2_lig3": {
                    "ddGbar_avg": -1.0,
                    "ddGbar_sem": 0.10,
                    "ddGbar_std": 0.18,
                    "from": "lig2",
                    "to": "lig3",
                },
            }
        },
    }

    results = summarize_fep_system_qc(data)

    assert len(results) == 2
    assert results[0]["fep"] == "FEP_lig1_lig2"
    assert results[1]["fep"] == "FEP_lig2_lig3"
    assert results[0]["ddg"] == pytest.approx(-1.3)
    assert results[1]["ddg"] == pytest.approx(-1.0)
