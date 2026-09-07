"""Small parser tests plus an opt-in integration against pinned HPC data."""
import importlib.util
import json
import os
from pathlib import Path
import struct

import pytest


ROOT = Path(__file__).resolve().parents[2]
CASE = ROOT / "experiments/historical-cmet-reproduction"
SPEC = importlib.util.spec_from_file_location("historical_cmet", CASE / "reproduce.py")
REPRO = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(REPRO)


def record(payload):
    marker = struct.pack("<i", len(payload))
    return marker + payload + marker


def test_fixed_width_bar_reader_handles_adjacent_negative_column(tmp_path):
    path = tmp_path / "qfep.out"
    path.write_text("# Part 6: BAR Bennet:\n# lambda dG total\n"+
                    f"   {0.999:8.6f}{-1234.123:9.3f}{-1234.123:9.3f}\n")
    assert REPRO.bar_rows(path) == [(0.999, -1234.123, -1234.123)]


@pytest.mark.parametrize("content", [b"x", struct.pack("<i", -1), record(b"abc")[:-1]])
def test_damaged_binary_records_fail(tmp_path, content):
    path = tmp_path / "bad.re"
    path.write_bytes(content)
    with pytest.raises(ValueError):
        list(REPRO.records(path))


def test_restart_offsets_decode_single_precision(tmp_path):
    path = tmp_path / "eq5.re"
    path.write_bytes(record(b"coords")+record(b"velocities")+record(struct.pack("<i3f", 3, .1, -.2, .3)))
    assert REPRO.restart_offsets(path) == pytest.approx([.1, -.2, .3])


def test_archived_log_charge_and_radius_convention(tmp_path):
    path = tmp_path / "md.log"
    path.write_text("Eff. solvent radius = 23.300\n SUM 1.000 -0.004\n"
                    "Total charge of non-Q atoms = 3.00\n SUM -1 -2 -3\n")
    assert REPRO.boundary_parameters(path) == {"radius": 23.3, "q_from": 1., "q_to": -.004, "q_env": 3.}


@pytest.mark.parametrize("stored_lambda", [.999, .998])
def test_binary_mapping_verification(tmp_path, stored_lambda):
    path = tmp_path / "window.en"
    path.write_bytes(record(struct.pack("<i15d", 1, stored_lambda, *([0.]*14)))+
                     record(struct.pack("<i15d", 2, 1-stored_lambda, *([0.]*14)))+record(b""))
    manifest = {"retained": [{"source_energy": str(path), "lambdas": [".999", ".001"]}]}
    if stored_lambda == .999:
        assert REPRO.check_energy_mappings(manifest)["all_stored_lambda_pairs_match_md_inputs"]
    else:
        with pytest.raises(ValueError, match="Stored mapping"):
            REPRO.check_energy_mappings(manifest)


def test_historical_cmet_native_reproduction(tmp_path):
    if os.environ.get("Q_RUN_HISTORICAL_CMET") != "1":
        pytest.skip("Set Q_RUN_HISTORICAL_CMET=1 after fetching and building the pinned reference")
    result = REPRO.run(tmp_path / "analysis", CASE / "engine-build/qfep")
    frozen = json.loads((CASE / "result.json").read_text())
    assert result["gate"] == frozen["gate"]
    assert result["download_verification"] == frozen["download_verification"]
    assert len(result["reductions"]) == 4
    assert result["forward_reverse_offsets_identical"] is False
    for new, reference in zip(result["reductions"], frozen["reductions"]):
        assert new["label"] == reference["label"]
        assert new["bar_kcal_mol"] == pytest.approx(reference["bar_kcal_mol"], abs=.002)
        assert new["born_interval_shift_kcal_mol"] == pytest.approx(reference["born_interval_shift_kcal_mol"], abs=1e-10)
        assert new["maximum_cumulative_row_error_kcal_mol"] <= .002
        assert new["retained_pairs"] == 98
