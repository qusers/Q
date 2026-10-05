"""Tests for QligFEP quality-control analysis."""

import pytest

from QligFEP.analysis_qc import (
    analyze_ddg_edge,
    analyze_fep_edge,
    analyze_fep_system,
    check_fep_pair_consistency,
    cycle_closure_error,
    find_cycle_basis,
    replicate_statistics,
    summarize_cycle_closure_qc,
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


def test_replicate_statistics_excludes_all_nonfinite_values():
    result = replicate_statistics([1.0, float("nan"), float("inf"), -float("inf"), 3.0])

    assert result["n_replicates"] == 5
    assert result["n_valid"] == 2
    assert result["n_failed"] == 3
    assert result["mean"] == pytest.approx(2.0)
    assert result["std"] == pytest.approx(1.0)
    assert result["sem"] == pytest.approx(1.0 / 2**0.5)
    assert result["range"] == pytest.approx(2.0)


def test_replicate_statistics_all_nonfinite():
    result = replicate_statistics([float("inf"), -float("inf")])

    assert result["n_valid"] == 0
    assert result["n_failed"] == 2
    assert all(result[key] is None for key in ("mean", "std", "sem", "range"))


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
                "lambda_sum": 100,
                "input_n_lambdas": 100,
            }
        },
        "2.protein": {
            "FEP_lig1_lig2": {
                "fep_stage": "production",
                "temperature": 298,
                "lambda_sum": 100,
                "input_n_lambdas": 100,
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
    assert result["consistency_status"] == "consistent"


def test_check_fep_pair_consistency_detects_mismatch():
    """Protein/water setting mismatches should fail consistency QC."""
    data = {
        "1.water": {
            "FEP_lig1_lig2": {
                "fep_stage": "production",
                "temperature": 298,
                "lambda_sum": 100,
                "input_n_lambdas": 100,
            }
        },
        "2.protein": {
            "FEP_lig1_lig2": {
                "fep_stage": "production",
                "temperature": 310,
                "lambda_sum": 100,
                "input_n_lambdas": 100,
            }
        },
    }

    result = check_fep_pair_consistency(
        data,
        fep="FEP_lig1_lig2",
    )

    assert result["temperature_match"] is False
    assert result["consistent"] is False
    assert result["consistency_status"] == "mismatch"


@pytest.mark.parametrize("mismatch", [False, True])
def test_missing_window_metadata_is_unverified_unless_another_check_fails(mismatch):
    water = {"fep_stage": "FEP1", "temperature": "298", "lambda_sum": 100}
    protein = {**water, "temperature": "310" if mismatch else "298", "input_n_lambdas": 100}
    data = {"1.water": {"FEP_A_B": water}, "2.protein": {"FEP_A_B": protein}}

    result = check_fep_pair_consistency(data, "FEP_A_B")

    assert result["lambda_sum_match"] is None
    assert result["consistent"] is (False if mismatch else None)
    assert result["consistency_status"] == ("mismatch" if mismatch else "unverified")


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


@pytest.mark.parametrize("value", [None, float("nan"), float("inf"), -float("inf")])
def test_analyze_ddg_edge_normalizes_unavailable_statistics(value):
    data = {
        "result": {
            "ddGbar": {
                "FEP_A_B": {
                    "from": "A",
                    "to": "B",
                    "ddGbar_avg": value,
                    "ddGbar_sem": value,
                    "ddGbar_std": value,
                }
            }
        }
    }

    result = analyze_ddg_edge(data, "FEP_A_B")

    assert all(result[key] is None for key in ("avg", "sem", "std"))


def test_summarize_fep_edge_qc():
    """QC summary should combine ddG, replicate, and consistency information."""
    data = {
        "1.water": {
            "FEP_lig1_lig2": {
                "fep_stage": "production",
                "temperature": 298,
                "lambda_sum": 100,
                "input_n_lambdas": 100,
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
                "lambda_sum": 100,
                "input_n_lambdas": 100,
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
                "lambda_sum": 100,
                "input_n_lambdas": 100,
                "FEP_result": {
                    "dGbar": {
                        "energies": [-2.0, -2.1, -1.9],
                    }
                },
            },
            "FEP_lig2_lig3": {
                "fep_stage": "production",
                "temperature": 298,
                "lambda_sum": 100,
                "input_n_lambdas": 100,
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
                "lambda_sum": 100,
                "input_n_lambdas": 100,
                "FEP_result": {
                    "dGbar": {
                        "energies": [-3.2, -3.4, -3.3],
                    }
                },
            },
            "FEP_lig2_lig3": {
                "fep_stage": "production",
                "temperature": 298,
                "lambda_sum": 100,
                "input_n_lambdas": 100,
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


def test_cycle_closure_error():
    """Cycle closure should sum directed ddG values around a closed cycle."""
    data = {
        "result": {
            "ddGbar": {
                "FEP_lig1_lig2": {
                    "ddGbar_avg": 1.2,
                    "from": "lig1",
                    "to": "lig2",
                },
                "FEP_lig2_lig3": {
                    "ddGbar_avg": -0.4,
                    "from": "lig2",
                    "to": "lig3",
                },
                "FEP_lig3_lig1": {
                    "ddGbar_avg": -0.3,
                    "from": "lig3",
                    "to": "lig1",
                },
            }
        }
    }

    result = cycle_closure_error(
        data,
        cycle=["lig1", "lig2", "lig3", "lig1"],
    )

    assert result == pytest.approx(0.5)


def test_cycle_closure_error_handles_reverse_edge():
    """Cycle closure should invert ddG when traversing an edge backwards."""
    data = {
        "result": {
            "ddGbar": {
                "FEP_lig1_lig2": {
                    "ddGbar_avg": 1.2,
                    "from": "lig1",
                    "to": "lig2",
                },
                "FEP_lig2_lig3": {
                    "ddGbar_avg": -0.4,
                    "from": "lig2",
                    "to": "lig3",
                },
                "FEP_lig1_lig3": {
                    "ddGbar_avg": 0.3,
                    "from": "lig1",
                    "to": "lig3",
                },
            }
        }
    }

    result = cycle_closure_error(
        data,
        cycle=["lig1", "lig2", "lig3", "lig1"],
    )

    assert result == pytest.approx(0.5)


def test_summarize_cycle_closure_qc():
    """Cycle QC should report closure errors for all detected cycles."""
    data = {
        "result": {
            "ddGbar": {
                "FEP_lig1_lig2": {
                    "ddGbar_avg": 1.2,
                    "from": "lig1",
                    "to": "lig2",
                },
                "FEP_lig2_lig3": {
                    "ddGbar_avg": -0.4,
                    "from": "lig2",
                    "to": "lig3",
                },
                "FEP_lig3_lig1": {
                    "ddGbar_avg": -0.3,
                    "from": "lig3",
                    "to": "lig1",
                },
            }
        }
    }

    results = summarize_cycle_closure_qc(data)

    assert len(results) == 1
    assert results[0]["cycle"] == ["lig1", "lig2", "lig3", "lig1"]
    assert results[0]["n_edges"] == 3
    assert results[0]["closure_error"] == pytest.approx(0.5)
    assert results[0]["abs_closure_error"] == pytest.approx(0.5)
    assert results[0]["status"] == "ok"
    assert results[0]["reason"] == ""


@pytest.mark.parametrize("value", [None, float("nan"), float("inf"), -float("inf")])
def test_cycle_qc_retains_unavailable_cycles_and_reports_usable_cycles(value):
    edges = {
        "ab": {"from": "A", "to": "B", "ddGbar_avg": 1.0},
        "bc": {"from": "B", "to": "C", "ddGbar_avg": value},
        "ca": {"from": "C", "to": "A", "ddGbar_avg": -2.0},
        "de": {"from": "D", "to": "E", "ddGbar_avg": 1.0},
        "ef": {"from": "E", "to": "F", "ddGbar_avg": 2.0},
        "fd": {"from": "F", "to": "D", "ddGbar_avg": -2.5},
    }

    results = summarize_cycle_closure_qc({"result": {"ddGbar": edges}})

    assert len(results) == 2
    unavailable, usable = results
    assert unavailable["cycle"] == ["A", "B", "C", "A"]
    assert unavailable["status"] == "unavailable"
    assert unavailable["closure_error"] is None
    assert unavailable["abs_closure_error"] is None
    assert unavailable["reason"] == "Missing FEP edge: B -> C"
    assert usable["status"] == "ok"
    assert usable["closure_error"] == pytest.approx(0.5)


def test_cycle_qc_for_network_without_cycles():
    data = {"result": {"ddGbar": {"ab": {"from": "A", "to": "B", "ddGbar_avg": 1.0}}}}

    assert summarize_cycle_closure_qc(data) == []


def test_find_cycle_basis():
    """Cycle basis should contain only independent cycles."""
    data = {
        "result": {
            "ddGbar": {
                "ab": {"ddGbar_avg": 1.0, "from": "A", "to": "B"},
                "bc": {"ddGbar_avg": 1.0, "from": "B", "to": "C"},
                "ca": {"ddGbar_avg": -2.0, "from": "C", "to": "A"},
                "cd": {"ddGbar_avg": 1.0, "from": "C", "to": "D"},
                "da": {"ddGbar_avg": -1.0, "from": "D", "to": "A"},
            }
        }
    }

    cycles = find_cycle_basis(data)

    assert len(cycles) == 2
    assert all(cycle[0] == cycle[-1] for cycle in cycles)


def test_cycle_closure_error_rejects_invalid_ddg():
    """Cycle closure should fail clearly when a required ddG is unavailable."""
    data = {
        "result": {
            "ddGbar": {
                "FEP_lig1_lig2": {
                    "ddGbar_avg": 1.2,
                    "from": "lig1",
                    "to": "lig2",
                },
                "FEP_lig2_lig3": {
                    "ddGbar_avg": None,
                    "from": "lig2",
                    "to": "lig3",
                },
                "FEP_lig3_lig1": {
                    "ddGbar_avg": -0.3,
                    "from": "lig3",
                    "to": "lig1",
                },
            }
        }
    }

    with pytest.raises(ValueError, match="Missing FEP edge: lig2 -> lig3"):
        cycle_closure_error(
            data,
            cycle=["lig1", "lig2", "lig3", "lig1"],
        )
