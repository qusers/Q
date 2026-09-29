import json
from argparse import Namespace
from pathlib import Path

import pandas as pd

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
                        "lambda_sum": 1.0,
                        "FEP_result": {
                            "dGbar": {
                                "energies": [-2.0, -2.1, -1.9],
                            }
                        },
                    },
                    "FEP_B_C": {
                        "fep_stage": "production",
                        "temperature": 298,
                        "lambda_sum": 1.0,
                        "FEP_result": {
                            "dGbar": {
                                "energies": [-1.0, -1.1, -0.9],
                            }
                        },
                    },
                    "FEP_C_A": {
                        "fep_stage": "production",
                        "temperature": 298,
                        "lambda_sum": 1.0,
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
                        "lambda_sum": 1.0,
                        "FEP_result": {
                            "dGbar": {
                                "energies": [-3.0, -3.1, -2.9],
                            }
                        },
                    },
                    "FEP_B_C": {
                        "fep_stage": "production",
                        "temperature": 298,
                        "lambda_sum": 1.0,
                        "FEP_result": {
                            "dGbar": {
                                "energies": [-1.5, -1.6, -1.4],
                            }
                        },
                    },
                    "FEP_C_A": {
                        "fep_stage": "production",
                        "temperature": 298,
                        "lambda_sum": 1.0,
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
