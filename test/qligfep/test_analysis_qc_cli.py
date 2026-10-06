import json
from argparse import Namespace
from pathlib import Path

import pandas as pd
import pytest

from QligFEP import analyze_FEP


def test_qc_csv_output(monkeypatch, tmp_path):
    """QC CLI output should create edge and cycle CSV files."""

    class FakeFepReader:
        def __init__(self, *args, **kwargs):
            self.data = {
                "1.water": {
                    "FEP_A_B": {
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
                    "FEP_B_C": {
                        "fep_stage": "production",
                        "temperature": 298,
                        "lambda_sum": 100,
                        "input_n_lambdas": 100,
                        "FEP_result": {
                            "dGbar": {
                                "energies": [-1.0, -1.1, -0.9],
                            }
                        },
                    },
                    "FEP_C_A": {
                        "fep_stage": "production",
                        "temperature": 298,
                        "lambda_sum": 100,
                        "input_n_lambdas": 100,
                        "FEP_result": {
                            "dGbar": {
                                "energies": [-0.5, -0.6, -0.4],
                            }
                        },
                    },
                },
                "2.protein": {
                    "FEP_A_B": {
                        "fep_stage": "production",
                        "temperature": 298,
                        "lambda_sum": 100,
                        "input_n_lambdas": 100,
                        "FEP_result": {
                            "dGbar": {
                                "energies": [-3.0, -3.1, -2.9],
                            }
                        },
                    },
                    "FEP_B_C": {
                        "fep_stage": "production",
                        "temperature": 298,
                        "lambda_sum": 100,
                        "input_n_lambdas": 100,
                        "FEP_result": {
                            "dGbar": {
                                "energies": [-1.5, -1.6, -1.4],
                            }
                        },
                    },
                    "FEP_C_A": {
                        "fep_stage": "production",
                        "temperature": 298,
                        "lambda_sum": 100,
                        "input_n_lambdas": 100,
                        "FEP_result": {
                            "dGbar": {
                                "energies": [-0.8, -0.9, -0.7],
                            }
                        },
                    },
                },
                "result": {
                    "ddGbar": {
                        "FEP_A_B": {
                            "ddGbar_avg": -1.0,
                            "ddGbar_sem": 0.1,
                            "ddGbar_std": 0.2,
                            "from": "A",
                            "to": "B",
                        },
                        "FEP_B_C": {
                            "ddGbar_avg": -0.5,
                            "ddGbar_sem": 0.1,
                            "ddGbar_std": 0.2,
                            "from": "B",
                            "to": "C",
                        },
                        "FEP_C_A": {
                            "ddGbar_avg": 0.7,
                            "ddGbar_sem": 0.1,
                            "ddGbar_std": 0.2,
                            "from": "C",
                            "to": "A",
                        },
                    }
                },
            }

            self.ignored_edges = []
            self.verbose_qEnergies = []
            self.verbose_dgBar = []
            self.run_data = []

        def read_perturbations(self, *args, **kwargs):
            pass

        def load_new_system(self, *args, **kwargs):
            pass

        def calculate_ddG(self):
            pass

        def save_json_data(self):
            pass

        def populate_mapping_dictionary(self, *args, **kwargs):
            output_file = kwargs["output_file"]

            Path(output_file).write_text(
                json.dumps(
                    {
                        "edges": [],
                    }
                )
            )

    monkeypatch.setattr(analyze_FEP, "FepReader", FakeFepReader)
    monkeypatch.setattr(
        analyze_FEP,
        "prepare_df",
        lambda *args, **kwargs: pd.DataFrame(),
    )
    monkeypatch.chdir(tmp_path)

    args = Namespace(
        log="info",
        water_dir="1.water",
        protein_dir="2.protein",
        target="qc_test",
        json_file="mapping.json",
        n_lambdas=None,
        allow_missing=False,
        no_run_data=True,
        method="ddGbar",
        qc=True,
        experimental_key=None,
        save_verbose=False,
    )

    analyze_FEP.main(args)

    edge_file = tmp_path / "qc_test_fep_qc.csv"
    cycle_file = tmp_path / "qc_test_cycle_qc.csv"

    assert edge_file.exists()
    assert cycle_file.exists()

    edge_df = pd.read_csv(edge_file)
    cycle_df = pd.read_csv(cycle_file)

    assert len(edge_df) == 3
    assert "max_leg_std" in edge_df.columns
    assert "max_leg_range" in edge_df.columns

    assert len(cycle_df) == 1
    assert "abs_closure_error" in cycle_df.columns
    assert cycle_df.iloc[0]["n_edges"] == 3


def test_no_qc_csv_output_without_flag(monkeypatch, tmp_path):
    """QC CSV files should not be created when --qc is not requested."""

    class FakeFepReader:
        def __init__(self, *args, **kwargs):
            self.data = {
                "1.water": {},
                "2.protein": {},
                "result": {"ddGbar": {}},
            }
            self.ignored_edges = []
            self.verbose_qEnergies = []
            self.verbose_dgBar = []
            self.run_data = []

        def read_perturbations(self, *args, **kwargs):
            pass

        def load_new_system(self, *args, **kwargs):
            pass

        def calculate_ddG(self):
            pass

        def save_json_data(self):
            pass

        def populate_mapping_dictionary(self, *args, **kwargs):
            output_file = kwargs["output_file"]
            Path(output_file).write_text(json.dumps({"edges": [{"ddg": 0.0}]}))

    monkeypatch.setattr(analyze_FEP, "FepReader", FakeFepReader)
    monkeypatch.setattr(
        analyze_FEP,
        "prepare_df",
        lambda *args, **kwargs: pd.DataFrame(),
    )
    monkeypatch.chdir(tmp_path)

    args = Namespace(
        log="info",
        water_dir="1.water",
        protein_dir="2.protein",
        target="no_qc_test",
        json_file="mapping.json",
        n_lambdas=None,
        allow_missing=False,
        no_run_data=True,
        method="ddGbar",
        qc=False,
        experimental_key=None,
        save_verbose=False,
    )

    analyze_FEP.main(args)

    assert not (tmp_path / "no_qc_test_fep_qc.csv").exists()
    assert not (tmp_path / "no_qc_test_cycle_qc.csv").exists()


def _write_edge(system, fep, n_inputs):
    root = system / fep
    (root / "inputfiles").mkdir(parents=True)
    for index in range(n_inputs):
        (root / "inputfiles" / f"md_{index:04d}.inp").touch()
    for replicate in (1, 2):
        directory = root / "FEP1" / "298" / str(replicate)
        directory.mkdir(parents=True)
        (directory / "qfep.out").write_text(str(replicate / 10))


@pytest.fixture
def calculation(monkeypatch, tmp_path):
    """Exercise the real reader and CLI, replacing only the Qfep text parsers."""
    edges = [{"from": source, "to": target} for source, target in [("A", "B"), ("B", "C"), ("C", "A")]]
    (tmp_path / "mapping.json").write_text(json.dumps({"edges": edges}))
    for system in ("1.water", "2.protein"):
        for edge in edges:
            _write_edge(tmp_path / system, f"FEP_{edge['from']}_{edge['to']}", n_inputs=3)

    def read_qfep(path):
        text = path.read_text()
        if text == "ERROR":
            raise OSError("Failed Qfep calculation")
        value = float(text) + (0.5 if "2.protein" in path.parts else 0)
        return [value] * 5

    monkeypatch.setattr(analyze_FEP, "read_qfep", read_qfep)
    monkeypatch.setattr(analyze_FEP, "read_qfep_verbose", lambda path: (None, None))
    monkeypatch.chdir(tmp_path)
    return Namespace(
        log="error",
        water_dir="1.water",
        protein_dir="2.protein",
        target="integration",
        json_file="mapping.json",
        n_lambdas=None,
        allow_missing=False,
        no_run_data=True,
        method="ddGbar",
        qc=True,
        experimental_key=None,
        save_verbose=False,
    )


@pytest.mark.parametrize("result_text", ["ERROR", "nan", "inf", "-inf"])
def test_qc_cli_preserves_primary_outputs_with_unavailable_cycle(calculation, tmp_path, result_text):
    for file in (tmp_path / "1.water" / "FEP_B_C").glob("FEP1/298/*/qfep.out"):
        file.write_text(result_text)

    analyze_FEP.main(calculation)

    assert (tmp_path / "integration_FEP_results.json").exists()
    mapping = json.loads((tmp_path / "mapping_ddG.json").read_text())
    assert len(mapping["edges"]) == 3
    edges = pd.read_csv(tmp_path / "integration_fep_qc.csv").set_index("fep")
    assert edges.loc["FEP_B_C", "water_n_failed"] == 2
    assert edges.loc["FEP_B_C", "water_n_valid"] == 0
    assert pd.isna(edges.loc["FEP_B_C", "ddg"])
    cycles = pd.read_csv(tmp_path / "integration_cycle_qc.csv")
    assert len(cycles) == 1
    assert cycles.loc[0, "status"] == "unavailable"
    assert pd.isna(cycles.loc[0, "closure_error"])
    assert cycles.loc[0, "reason"] == "Missing FEP edge: B -> C"


@pytest.mark.parametrize("override", [None, 100])
def test_window_counts_are_read_per_edge_and_system(calculation, tmp_path, override):
    (tmp_path / "2.protein" / "FEP_B_C" / "inputfiles" / "md_extra.inp").touch()
    calculation.n_lambdas = override

    analyze_FEP.main(calculation)

    data = json.loads((tmp_path / "integration_FEP_results.json").read_text())
    assert data["1.water"]["FEP_B_C"]["input_n_lambdas"] == 2
    assert data["2.protein"]["FEP_B_C"]["input_n_lambdas"] == 3
    assert data["2.protein"]["FEP_C_A"]["input_n_lambdas"] == 2
    assert data["2.protein"]["FEP_B_C"]["lambda_sum"] == (3 if override is None else override)
    qc = pd.read_csv(tmp_path / "integration_fep_qc.csv").set_index("fep")
    assert qc.loc["FEP_B_C", "consistency_status"] == "mismatch"
    assert not qc.loc["FEP_B_C", "lambda_sum_match"]
    assert qc.loc["FEP_C_A", "consistency_status"] == "consistent"


def test_archived_results_have_unverified_input_count_check(calculation, tmp_path):
    for file in tmp_path.glob("*.*/FEP_*/inputfiles/*.inp"):
        file.unlink()
    calculation.n_lambdas = 100

    analyze_FEP.main(calculation)

    qc = pd.read_csv(tmp_path / "integration_fep_qc.csv")
    assert (qc["consistency_status"] == "unverified").all()
    assert qc["lambda_sum_match"].isna().all()
    assert qc["consistent"].isna().all()


def test_primary_results_are_saved_before_optional_report_write(calculation, monkeypatch, tmp_path):
    def fail_write(*args, **kwargs):
        raise PermissionError("Cannot write QC report")

    monkeypatch.setattr(pd.DataFrame, "to_csv", fail_write)

    with pytest.raises(PermissionError, match="Cannot write QC report"):
        analyze_FEP.main(calculation)

    assert (tmp_path / "integration_FEP_results.json").exists()
    assert (tmp_path / "mapping_ddG.json").exists()


def test_acyclic_network_writes_empty_cycle_csv_with_headers(calculation, tmp_path):
    import shutil

    for system in ("1.water", "2.protein"):
        shutil.rmtree(tmp_path / system / "FEP_C_A")
    mapping = json.loads((tmp_path / "mapping.json").read_text())
    mapping["edges"] = mapping["edges"][:2]
    (tmp_path / "mapping.json").write_text(json.dumps(mapping))

    analyze_FEP.main(calculation)

    cycles = pd.read_csv(tmp_path / "integration_cycle_qc.csv")
    assert cycles.empty
    assert list(cycles.columns) == [
        "cycle",
        "n_edges",
        "closure_error",
        "abs_closure_error",
        "status",
        "reason",
    ]
