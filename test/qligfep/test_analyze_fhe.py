from pathlib import Path
import json

import numpy as np
import pytest

from QligFEP.analyze_fhe import _bar_delta_f, analyze_fhe, read_fhe, write_results


def _write_fhe(
    path: Path,
    source_l1: float,
    source_l2: float,
    rows: list[tuple[float, float]],
    *,
    work_format: bool = False,
) -> None:
    lines = [
        "# Q_FOREIGN_WORKS_V2" if work_format else "# Q_FOREIGN_HAMILTONIANS_V1",
        "# units kcal/mol",
        "# method gapsys",
        "# target_state2_weights 0.0000000 1.0000000",
        "# columns frame source_lambda1 source_lambda2 H(target_1) H(target_2)",
    ]
    lines.extend(
        f"{frame:10d} {source_l1:10.7f} {source_l2:10.7f} {h0:24.16e} {h1:24.16e}"
        for frame, (h0, h1) in enumerate(rows, start=1)
    )
    path.write_text("\n".join(lines) + "\n")


def test_read_fhe_and_adjacent_bar(tmp_path):
    state0 = tmp_path / "md_1000_0000.en.fhe"
    state1 = tmp_path / "md_0000_1000.en.fhe"
    _write_fhe(state0, 1.0, 0.0, [(0.0, 2.0), (0.0, 2.0)])
    _write_fhe(state1, 0.0, 1.0, [(-2.0, 0.0), (-2.0, 0.0)])

    parsed = read_fhe(state0)
    assert parsed.source_lambda2 == 0.0
    assert parsed.energy_at(1.0).tolist() == [2.0, 2.0]
    assert parsed.work_to(1.0).tolist() == [2.0, 2.0]

    result = analyze_fhe([state0, state1], temperature=298.0)
    assert result["total_delta_g_kcal"] == pytest.approx(2.0)
    assert result["poor_windows"] == 0
    assert result["adjacent"][0]["overlap_gap_kT"] == pytest.approx(0.0)
    assert result["adjacent"][0]["fermi_ess_forward"] == pytest.approx(2.0)
    assert result["adjacent"][0]["fermi_ess_reverse"] == pytest.approx(2.0)
    assert result["adjacent"][0]["bar_residual"] == pytest.approx(0.0)


def test_work_format_and_legacy_totals_give_identical_bar(tmp_path):
    legacy0 = tmp_path / "legacy0.fhe"
    legacy1 = tmp_path / "legacy1.fhe"
    work0 = tmp_path / "work0.fhe"
    work1 = tmp_path / "work1.fhe"
    _write_fhe(legacy0, 1.0, 0.0, [(10.0, 12.0), (20.0, 22.0)])
    _write_fhe(legacy1, 0.0, 1.0, [(8.0, 10.0), (18.0, 20.0)])
    _write_fhe(work0, 1.0, 0.0, [(0.0, 2.0), (0.0, 2.0)], work_format=True)
    _write_fhe(work1, 0.0, 1.0, [(-2.0, 0.0), (-2.0, 0.0)], work_format=True)

    parsed = read_fhe(work0)
    assert parsed.format == "Q_FOREIGN_WORKS_V2"
    assert parsed.work_to(1.0).tolist() == [2.0, 2.0]
    with pytest.raises(ValueError, match="stores work"):
        parsed.energy_at(1.0)

    legacy = analyze_fhe([legacy0, legacy1], temperature=298.0)
    work = analyze_fhe([work0, work1], temperature=298.0)
    assert work["format"] == "Q_FOREIGN_WORKS_V2"
    assert work["total_delta_g_kcal"] == pytest.approx(legacy["total_delta_g_kcal"])
    assert work["adjacent"][0]["delta_g_kcal"] == pytest.approx(
        legacy["adjacent"][0]["delta_g_kcal"]
    )
    assert work["adjacent"][0]["mean_forward_work_kcal"] == pytest.approx(
        legacy["adjacent"][0]["mean_forward_work_kcal"]
    )


def test_frame_selection_applies_per_window(tmp_path):
    state0 = tmp_path / "state0.fhe"
    state1 = tmp_path / "state1.fhe"
    rows0 = [(0.0, 2.0)] * 6
    rows1 = [(-2.0, 0.0)] * 6
    _write_fhe(state0, 1.0, 0.0, rows0, work_format=True)
    _write_fhe(state1, 0.0, 1.0, rows1, work_format=True)

    result = analyze_fhe([state0, state1], skip=1, stride=2)
    assert result["frame_selection"] == {"skip": 1, "stride": 2}
    assert result["adjacent"][0]["frames"] == 3
    assert result["total_delta_g_kcal"] == pytest.approx(2.0)

    with pytest.raises(ValueError, match="skip must be non-negative"):
        analyze_fhe([state0, state1], skip=-1)
    with pytest.raises(ValueError, match="stride must be positive"):
        analyze_fhe([state0, state1], stride=0)
    with pytest.raises(ValueError, match="no frames remain"):
        analyze_fhe([state0, state1], skip=6)


def test_read_fhe_rejects_nonfinite(tmp_path):
    path = tmp_path / "bad.fhe"
    _write_fhe(path, 1.0, 0.0, [(0.0, float("inf"))])
    with pytest.raises(ValueError, match="non-finite"):
        read_fhe(path)


def test_unequal_counts_and_large_work_are_stable(tmp_path):
    state0, state1 = tmp_path / "a.fhe", tmp_path / "b.fhe"
    _write_fhe(state0, 1, 0, [(0, 2)] * 7, work_format=True)
    _write_fhe(state1, 0, 1, [(-2, 0)] * 3, work_format=True)
    result = analyze_fhe([state0, state1])
    assert result["total_delta_g_kcal"] == pytest.approx(2)
    assert result["adjacent"][0]["frames_reverse"] == 3
    assert result["adjacent"][0]["bar_residual"] == pytest.approx(0)
    # Both acceptances underflow in ordinary floating-point arithmetic. The
    # root must still reflect the individual works, not a spurious 0 == 0.
    with np.errstate(over="raise", invalid="raise"):
        assert _bar_delta_f(np.array([1000., 2000.]), np.array([900., 3000.]), 1) == pytest.approx(50)


def test_bar_count_correction_against_exact_discrete_ensembles():
    # Partition functions: ZA=4, ZB=6. Integer Boltzmann weights let us
    # enumerate both canonical ensembles exactly, at unequal sample counts.
    forward = np.repeat(-np.log(np.array([4., 2.]) / [1., 3.]), [1, 3])
    reverse = np.repeat(np.log(np.array([4., 2.]) / [1., 3.]), [4, 2])
    assert _bar_delta_f(forward, reverse, 1) == pytest.approx(-np.log(6/4))


@pytest.mark.parametrize("old,new,message", [
    ("0.0000000 1.0000000", "nan 1.0000000", "target lambdas"),
    ("0.0000000 1.0000000", "0.0000000 0.0000000", "duplicate target"),
    ("0.0000000 1.0000000", "0.5000000 1.0000000", "missing source"),
])
def test_rejects_malformed_targets(tmp_path, old, new, message):
    path = tmp_path / "bad.fhe"
    _write_fhe(path, 1, 0, [(0, 2)], work_format=True)
    path.write_text(path.read_text().replace(old, new))
    with pytest.raises(ValueError, match=message):
        read_fhe(path)


def test_rejects_invalid_sources_and_nonzero_self_work(tmp_path):
    path = tmp_path / "bad.fhe"
    _write_fhe(path, .5, 0, [(0, 2)], work_format=True)
    with pytest.raises(ValueError, match="sum to one"):
        read_fhe(path)
    _write_fhe(path, 1, 0, [(1, 2)], work_format=True)
    with pytest.raises(ValueError, match="source-to-source"):
        read_fhe(path)


def test_rejects_duplicate_paths_and_nonfinite_temperature(tmp_path):
    path = tmp_path / "a.fhe"
    _write_fhe(path, 1, 0, [(0, 2)], work_format=True)
    with pytest.raises(ValueError, match="duplicate sidecar"):
        analyze_fhe([path, path])
    with pytest.raises(ValueError, match="temperature must be positive"):
        analyze_fhe([path], temperature=float("nan"))


@pytest.mark.parametrize("change", ["method", "parameters", "missing"])
def test_rejects_incompatible_softcore_hamiltonians(tmp_path, change):
    a,b = tmp_path / "a.fhe", tmp_path / "b.fhe"
    _write_fhe(a, 1, 0, [(0, 2)], work_format=True)
    _write_fhe(b, 0, 1, [(-2, 0)], work_format=True)
    if change == "method":
        b.write_text(b.read_text().replace("# method gapsys", "# method beutler"))
    else:
        a.write_text(a.read_text() + "# softcore_parameters 0.5 100\n")
        if change == "parameters":
            b.write_text(b.read_text() + "# softcore_parameters 0.5 42\n")
    with pytest.raises(ValueError, match="different softcore methods/parameters"):
        analyze_fhe([a,b])


def test_equilibrium_reader_uses_foreign_bar_without_qfep(tmp_path, monkeypatch):
    from QligFEP.analyze_FEP import FepReader

    directory = tmp_path / "298" / "1"
    directory.mkdir(parents=True)
    state0, state1 = directory / "a.fhe", directory / "b.fhe"
    _write_fhe(state0, 1, 0, [(0, 2)] * 2, work_format=True)
    _write_fhe(state1, 0, 1, [(-2, 0)] * 2, work_format=True)
    result = analyze_fhe([state0, state1])
    write_results(result, directory / "foreign_bar")
    reader = FepReader.__new__(FepReader)
    reader.methods_list = ["dG", "dGf", "dGr", "dGos", "dGbar"]
    monkeypatch.setattr(reader, "run_qfep", lambda _: pytest.fail("legacy QFEP invoked"))
    values, failed, _, _ = reader.read_single_replicate(directory / "qfep.out")
    assert not failed
    assert values[-1] == pytest.approx(2)
    assert np.isnan(values[:4]).all()
    # Failed runs must not reuse an old successful result.
    (directory / "FAILED.log").write_text("MD failed")
    values, failed, _, _ = reader.read_single_replicate(directory / "qfep.out")
    assert failed == [1]
    assert np.isnan(values).all()


@pytest.mark.parametrize("method", ["gapsys", "beutler", "gapsys_v2", "beutler_coul"])
def test_nonlinear_method_cannot_fall_back_to_linear_qfep(tmp_path, method):
    from QligFEP.analyze_FEP import FepReader

    (tmp_path / "FEP1.fep").write_text(f"[FEP]\nsoftcore_method {method}\n")
    reader = FepReader.__new__(FepReader)
    with pytest.raises(ValueError, match="requires foreign-work BAR"):
        reader.run_qfep(tmp_path / "qfep.inp")


def test_result_collection_discovers_foreign_only_replicates(tmp_path, monkeypatch):
    from QligFEP.analyze_FEP import FepReader

    monkeypatch.chdir(tmp_path)
    directory = tmp_path / "water/FEP_a_b/FEP1/298/1"
    directory.mkdir(parents=True)
    state0, state1 = directory / "a.fhe", directory / "b.fhe"
    _write_fhe(state0, 1, 0, [(0, 2)] * 2, work_format=True)
    _write_fhe(state1, 0, 1, [(-2, 0)] * 2, work_format=True)
    write_results(analyze_fhe([state0, state1]), directory / "foreign_bar")
    mapping = tmp_path / "mapping.json"
    mapping.write_text(json.dumps({"edges": [{"from": "a", "to": "b"}]}))
    reader = FepReader("water", "test", str(mapping), n_lambdas=2)
    monkeypatch.setattr(reader, "run_qfep", lambda _: pytest.fail("legacy QFEP invoked"))
    reader.read_perturbations(add_run_data=False)
    result = reader.data["water"]["FEP_a_b"]
    assert result["FEP_result"]["dGbar"]["avg"] == pytest.approx(2)
    assert result["CrashedReplicates"] == []
    assert not (directory / "qfep.out").exists()
