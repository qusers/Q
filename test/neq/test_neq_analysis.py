"""Tests for the NEQ analysis extras: switch-completion diagnostics, run diagnostics from
slurm logs, experimental-data matching, correlation metrics, and the end-to-end main(). All
run without Q binaries.
"""

import matplotlib

matplotlib.use("Agg")  # headless backend before analyze_neq imports pyplot via analysis_plotting

import argparse
import json

import numpy as np
import pytest

from QligFEP.analyze_neq import (
    KB_KCAL,
    analyze_edge,
    bar_delta_f,
    beta_from_units,
    collect_works_with_counts,
    correlation_metrics,
    load_experimental_ddG,
    main,
    parse_arguments,
    parse_run_diagnostics,
    read_final_work,
    works_by_replicate,
)

COMPLETE = "qdyn version 6.0.1 terminated normally.\n"


def test_work_units_default_is_kcal_with_physical_beta(monkeypatch):
    # the analyzer defaults to the physically consistent kcal convention (beta = 1/kT),
    # not kT (beta = 1)
    monkeypatch.setattr("sys.argv", ["qligfep_neq_analyze"])
    assert parse_arguments().work_units == "kcal"
    assert beta_from_units("kcal", 300.0) == pytest.approx(1.0 / (KB_KCAL * 300.0))
    assert beta_from_units("kT", 300.0) == 1.0


def test_leg_dir_flags_match_equilibrium_analyzer(monkeypatch):
    # leg directories use -p/-w (same as qligfep_analyze); single-edge mode uses -pe/-we
    monkeypatch.setattr("sys.argv", ["qligfep_neq_analyze", "-p", "prot", "-w", "wat"])
    a = parse_arguments()
    assert (a.protein_dir, a.water_dir) == ("prot", "wat")
    assert a.protein_edge is None and a.water_edge is None

    monkeypatch.setattr("sys.argv", ["qligfep_neq_analyze", "-pe", "e/prot", "-we", "e/wat"])
    a = parse_arguments()
    assert (a.protein_edge, a.water_edge) == ("e/prot", "e/wat")
    assert (a.protein_dir, a.water_dir) == ("2.protein", "1.water")


def _switch_log(path, work, complete=True):
    body = f"At step 20000, work accumulated was   {work:.16E} and dU and dlambda were 0.1 0.1\n"
    path.write_text(body + (COMPLETE if complete else "ABNORMAL TERMINATION of qdyn\n"))


def test_collect_works_with_counts_tracks_attempted_and_completed(tmp_path):
    _switch_log(tmp_path / "neq_1_0.log", 2.0)  # forward ok
    _switch_log(tmp_path / "neq_1_1.log", 2.5)  # forward ok
    _switch_log(tmp_path / "neq_1_2.log", 9.9, complete=False)  # forward failed/incomplete
    _switch_log(tmp_path / "neq_0_0.log", -1.0)  # reverse ok
    forward, reverse, counts = collect_works_with_counts(str(tmp_path))
    assert sorted(forward) == pytest.approx([2.0, 2.5])
    assert sorted(reverse) == pytest.approx([-1.0])
    assert counts["forward_attempted"] == 3
    assert counts["forward_completed"] == 2
    assert counts["reverse_attempted"] == 1
    assert counts["reverse_completed"] == 1
    assert counts["failed"] == 1


def _crooks_consistent_works(true_df, variance, n_forward, n_reverse):
    """Deterministic samples of two Crooks-consistent Gaussian work distributions.

    In reduced units (beta=1) the forward work ~ N(a, s^2) and reverse work ~ N(s^2 - a, s^2)
    with a = true_df + s^2 / 2 satisfy P_F(x) / P_R(-x) = exp(x - true_df), so the count-corrected
    BAR estimate recovers ``true_df`` for any forward/reverse sample sizes. Samples sit on
    inverse-CDF midpoints so the empirical distributions match the Gaussians with negligible noise.
    """
    from scipy.stats import norm

    a = true_df + variance / 2.0

    def quantile_samples(mean, n):
        q = (np.arange(n) + 0.5) / n
        return mean + np.sqrt(variance) * norm.ppf(q)

    return quantile_samples(a, n_forward), quantile_samples(variance - a, n_reverse)


def test_bar_equal_counts_recovers_true_df():
    # Crooks-consistent Gaussians, equal counts -> the ln(N_F/N_R) term is zero.
    forward, reverse = _crooks_consistent_works(true_df=2.0, variance=2.0, n_forward=6000, n_reverse=6000)
    assert bar_delta_f(forward, reverse, beta=1.0) == pytest.approx(2.0, abs=0.02)


def _slurm(path, seed, replicate, runtime="1h:2m:3s", failure=None):
    lines = [
        f"Parameters T=300, replicate={replicate}, seed={seed}, neq_reps=5",
        "Running NEQ ...",
    ]
    if failure:
        lines.append(failure)
    lines += [
        "#    EXPRESS LOG for jobid: 12345",
        f"#    Runtime: {runtime}",
        f"#    Random seed: {seed}",
        f"#    Replicate Number: {replicate}",
    ]
    path.write_text("\n".join(lines) + "\n")


def test_parse_run_diagnostics_reads_runtime_seed_status(tmp_path):
    _slurm(tmp_path / "slurm.run1.node.1.out", seed=42, replicate=1)
    _slurm(tmp_path / "slurm.run2.node.2.out", seed=43, replicate=2, failure="DUE TO TIME LIMIT")
    diags = parse_run_diagnostics(str(tmp_path))
    by_rep = {d["replicate"]: d for d in diags}
    assert by_rep["1"]["seed"] == "42"
    assert by_rep["1"]["runtime"] == "1h:2m:3s"
    assert by_rep["1"]["status"] == "SUCCESS"
    assert by_rep["2"]["status"] == "TIMEOUT"


def test_load_experimental_ddG_maps_by_fep_name(tmp_path):
    mapping = {
        "edges": [
            {"from": "ejm_31", "to": "ejm_42", "ddg_value": 1.23},
            {"from": "ejm_42", "to": "ejm_43", "ddg_value": -0.5},
        ]
    }
    mfile = tmp_path / "mapping.json"
    mfile.write_text(json.dumps(mapping))
    exp = load_experimental_ddG(str(mfile), "ddg_value")
    assert exp["FEP_ejm_31_ejm_42"] == 1.23
    assert exp["FEP_ejm_42_ejm_43"] == -0.5


def test_correlation_metrics_perfect_correlation():
    m = correlation_metrics([1.0, 2.0, 3.0, 4.0], [1.0, 2.0, 3.0, 4.0])
    assert m["n"] == 4
    assert m["r2"] == pytest.approx(1.0)
    assert m["pearson"] == pytest.approx(1.0)
    assert m["rmse"] == pytest.approx(0.0)
    assert m["mae"] == pytest.approx(0.0)


def test_correlation_metrics_known_offset():
    # predicted is exp shifted by +1: perfect correlation, RMSE/MAE = 1
    m = correlation_metrics([2.0, 3.0, 4.0], [1.0, 2.0, 3.0])
    assert m["pearson"] == pytest.approx(1.0)
    assert m["rmse"] == pytest.approx(1.0)
    assert m["mae"] == pytest.approx(1.0)


# --- end-to-end main() over a synthetic setupFEP tree -------------------------------------

import pandas as pd  # noqa: E402


def _make_leg(leg_dir, replicates):
    """replicates: list of (fwd_works, rev_works, n_failed_fwd), one per replicate dir 298/<i>/."""
    for rep_idx, (fwd_works, rev_works, n_failed_fwd) in enumerate(replicates, start=1):
        rundir = leg_dir / "298" / str(rep_idx)
        rundir.mkdir(parents=True)
        for i, w in enumerate(fwd_works):
            _switch_log(rundir / f"neq_1_{i}.log", w)
        for j in range(n_failed_fwd):
            _switch_log(rundir / f"neq_1_{len(fwd_works) + j}.log", 9.9, complete=False)
        for i, w in enumerate(rev_works):
            _switch_log(rundir / f"neq_0_{i}.log", w)
        _slurm(leg_dir / f"slurm.run{rep_idx}.node.{rep_idx}.out", seed=6 + rep_idx, replicate=rep_idx)


def _args(**over):
    base = dict(
        protein_dir="2.protein",
        water_dir="1.water",
        protein_edge=None,
        water_edge=None,
        temperature=298.0,
        work_units="kT",
        output="neq_results.csv",
        json_file=None,
        experimental_key=None,
        target="neq",
        no_run_data=False,
        log="warning",
    )
    base.update(over)
    return argparse.Namespace(**base)


def test_main_end_to_end_with_experiment_and_diagnostics(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    # three edges with distinct protein signals so the predicted ddG varies (correlation defined)
    signals = {"FEP_a_b": 2.0, "FEP_c_d": -1.0, "FEP_e_f": 0.5}
    for name, s in signals.items():
        fwd = [s - 0.2, s, s + 0.2]
        rev = [-(s - 0.2), -s, -(s + 0.2)]
        water_fwd, water_rev = [0.1, 0.0, -0.1], [-0.1, 0.0, 0.1]
        # two replicates (so a SEM is defined); a failed forward switch in rep 1 of the first edge
        failed = 1 if name == "FEP_a_b" else 0
        _make_leg(tmp_path / "2.protein" / name, [(fwd, rev, failed), (fwd, rev, 0)])
        _make_leg(tmp_path / "1.water" / name, [(water_fwd, water_rev, 0), (water_fwd, water_rev, 0)])

    mapping = {
        "edges": [
            {"from": "a", "to": "b", "ddg_value": 1.8},
            {"from": "c", "to": "d", "ddg_value": -0.9},
            {"from": "e", "to": "f", "ddg_value": 0.6},
        ]
    }
    mapping_file = tmp_path / "mapping.json"
    mapping_file.write_text(json.dumps(mapping))

    df = main(_args(json_file=str(mapping_file), experimental_key="ddg_value", target="tyk2"))

    # results frame carries experiment + failed-switch diagnostics
    assert set(df["edge"]) == set(signals)
    assert "ddG_exp" in df.columns
    a_b = df[df["edge"] == "FEP_a_b"].iloc[0]
    assert a_b["ddG_exp"] == 1.8
    assert a_b["n_failed_protein"] >= 1
    assert a_b["n_replicates"] == 2
    assert "ddG_sem_kcal" in df.columns

    # outputs written to disk
    results = pd.read_csv(tmp_path / "neq_results.csv")
    assert "ddG_exp" in results.columns and len(results) == 3
    run_data = pd.read_csv(tmp_path / "neq_results_run_data.csv")
    assert "status" in run_data.columns and len(run_data) >= 1
    assert (tmp_path / "tyk2_neq_ddG_plot.png").exists()


def test_main_injects_ddG_into_mapping_json(tmp_path, monkeypatch):
    # matching analyze_FEP: --json-file (no experimental key needed) writes a <name>_ddG.json
    # copy with Q_ddG_avg/sem/std injected per edge, leaving the original mapping untouched.
    monkeypatch.chdir(tmp_path)
    # one fwd + one rev switch per replicate -> BAR dF = (fwd - rev)/2. protein dF = [1, 3],
    # water dF = 0 -> ddG = [1, 3]: mean 2, std sqrt(2), sem 1.
    protein = [([2.0], [0.0], 0), ([6.0], [0.0], 0)]
    water = [([0.0], [0.0], 0), ([0.0], [0.0], 0)]
    _make_leg(tmp_path / "2.protein" / "FEP_a_b", protein)
    _make_leg(tmp_path / "1.water" / "FEP_a_b", water)

    mapping = {"edges": [{"from": "a", "to": "b", "ddg_value": 1.8}]}
    mapping_file = tmp_path / "mapping.json"
    mapping_file.write_text(json.dumps(mapping))

    main(_args(work_units="kcal", json_file=str(mapping_file), no_run_data=True))

    ddG_file = tmp_path / "mapping_ddG.json"
    assert ddG_file.exists()
    edge = json.loads(ddG_file.read_text())["edges"][0]
    assert edge["Q_ddG_avg"] == pytest.approx(2.0)
    assert edge["Q_ddG_std"] == pytest.approx(np.sqrt(2))
    assert edge["Q_ddG_sem"] == pytest.approx(1.0)
    assert edge["ddg_value"] == 1.8  # existing edge fields are preserved

    # the original mapping file is not mutated
    assert "Q_ddG_avg" not in json.loads(mapping_file.read_text())["edges"][0]


def test_populate_mapping_json_writes_null_for_nan(tmp_path):
    # an edge whose ddG could not be estimated (NaN) is serialized as JSON null, not NaN.
    from QligFEP.analyze_neq import populate_mapping_json

    df = pd.DataFrame(
        [{"edge": "FEP_a_b", "ddG_kcal": np.nan, "ddG_sem_kcal": np.nan, "ddG_std_kcal": np.nan}]
    )
    mapping_file = tmp_path / "mapping.json"
    mapping_file.write_text(json.dumps({"edges": [{"from": "a", "to": "b"}]}))
    out_file = tmp_path / "mapping_ddG.json"
    populate_mapping_json(df, str(mapping_file), str(out_file))
    edge = json.loads(out_file.read_text())["edges"][0]
    assert edge["Q_ddG_avg"] is None
    assert edge["Q_ddG_sem"] is None
    assert edge["Q_ddG_std"] is None


def test_main_without_experiment_skips_plot(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    _make_leg(tmp_path / "2.protein" / "FEP_a_b", [([2.0, 2.2, 1.8], [-2.0, -2.2, -1.8], 0)])
    _make_leg(tmp_path / "1.water" / "FEP_a_b", [([0.1, 0.0, -0.1], [-0.1, 0.0, 0.1], 0)])
    df = main(_args(no_run_data=True))
    assert "ddG_exp" not in df.columns
    assert (tmp_path / "neq_results.csv").exists()
    assert not (tmp_path / "neq_neq_ddG_plot.png").exists()


def test_works_by_replicate_groups_by_dir(tmp_path):
    # switches are grouped by their replicate directory (298/1, 298/2, ...)
    leg = tmp_path / "2.protein" / "FEP_a_b"
    _make_leg(leg, [([1.0, 1.2], [-1.0], 0), ([2.0], [-2.0, -2.2], 0)])
    by_rep = works_by_replicate(str(leg))
    assert set(by_rep) == {"1", "2"}
    fwd1, rev1 = by_rep["1"]
    assert sorted(fwd1) == pytest.approx([1.0, 1.2])
    assert rev1 == pytest.approx([-1.0])
    fwd2, rev2 = by_rep["2"]
    assert fwd2 == pytest.approx([2.0])
    assert sorted(rev2) == pytest.approx([-2.2, -2.0])


class _CountingHandle:
    """File-handle proxy that tallies bytes served through read()/readlines()/iteration."""

    def __init__(self, fh, counter):
        self._fh = fh
        self._counter = counter

    def read(self, *a, **k):
        data = self._fh.read(*a, **k)
        self._counter["n"] += len(data)
        return data

    def readlines(self, *a, **k):
        lines = self._fh.readlines(*a, **k)
        self._counter["n"] += sum(len(line) for line in lines)
        return lines

    def __iter__(self):
        for line in self._fh:
            self._counter["n"] += len(line)
            yield line

    def __enter__(self):
        self._fh.__enter__()
        return self

    def __exit__(self, *a):
        return self._fh.__exit__(*a)

    def __getattr__(self, name):
        return getattr(self._fh, name)


def test_read_final_work_reads_only_the_tail_of_large_logs(tmp_path, monkeypatch):
    # Real switch logs are ~10 MB; the final work line and the "terminated normally" footer are
    # both at the very end. read_final_work must read only the tail, not stream the whole file --
    # otherwise a leg's ~100 logs cost gigabytes on Lustre.
    log = tmp_path / "neq_1_0.log"
    filler = ("At step 1, work accumulated was   0.1 and dU and dlambda were 0.1 0.1\n") * 30000
    log.write_text(
        filler + "At step 20000, work accumulated was   3.5 and dU and dlambda were 0.1 0.1\n" + COMPLETE
    )
    size = log.stat().st_size
    assert size > 2_000_000  # a genuinely large log

    counter = {"n": 0}
    real_open = open

    def counting_open(path, *a, **k):
        fh = real_open(path, *a, **k)
        if str(path) == str(log):
            return _CountingHandle(fh, counter)
        return fh

    monkeypatch.setattr("builtins.open", counting_open)
    assert read_final_work(str(log)) == pytest.approx(3.5)
    assert counter["n"] < size / 2  # tail only, not the whole multi-MB file


def test_read_final_work_incomplete_switch_returns_none(tmp_path):
    # a switch that never printed "terminated normally" (crash / disk-full) is unfinished -> None
    log = tmp_path / "neq_1_0.log"
    log.write_text(
        ("At step 1, work accumulated was   0.1 and dU and dlambda were 0.1 0.1\n") * 30000
        + "ABNORMAL TERMINATION of qdyn\n"
    )
    assert read_final_work(str(log)) is None


def test_each_switch_log_is_read_once_per_edge(tmp_path, monkeypatch):
    # I/O contract: analyzing an edge reads every switch log exactly once. The per-replicate
    # grouping and the attempted/completed switch counts both come from that single pass, so a
    # log must not be re-opened to gather counts and then re-opened to group by replicate.
    protein = [([1.0, 1.2], [-1.0, -1.1], 1), ([2.0], [-2.0], 0)]  # rep 1 has a failed forward
    water = [([0.1], [-0.1], 0), ([0.0], [0.0], 0)]
    p_dir = tmp_path / "2.protein" / "FEP_a_b"
    w_dir = tmp_path / "1.water" / "FEP_a_b"
    _make_leg(p_dir, protein)
    _make_leg(w_dir, water)

    import QligFEP.analyze_neq as mod

    n_logs = len(list(p_dir.rglob("neq_*.log"))) + len(list(w_dir.rglob("neq_*.log")))
    calls = {"n": 0}
    real = mod.read_final_work

    def counting(path):
        calls["n"] += 1
        return real(path)

    monkeypatch.setattr(mod, "read_final_work", counting)
    beta = mod.beta_from_units("kcal", 300.0)
    row = analyze_edge("FEP_a_b", p_dir, w_dir, beta, "kcal", 300.0)

    assert calls["n"] == n_logs  # each log opened exactly once, not twice
    assert row["n_failed_protein"] == 1  # counts still gathered in the single pass
    assert row["n_replicates"] == 2  # per-replicate grouping still intact


def test_per_replicate_ddG_is_mean_and_sem_over_replicates(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    # For one forward + one reverse switch, BAR dF = (fwd - rev) / 2 exactly. Build 3 replicates
    # whose protein dF = g and water dF = 0, so ddG(rep) = g. g = [1, 2, 3] -> mean 2, std 1,
    # sem 1/sqrt(3). The manuscript estimator averages per-replicate ddG and reports the SEM.
    g = [1.0, 2.0, 3.0]
    protein = [([2 * gr], [0.0], 0) for gr in g]
    water = [([0.0], [0.0], 0) for _ in g]
    _make_leg(tmp_path / "2.protein" / "FEP_a_b", protein)
    _make_leg(tmp_path / "1.water" / "FEP_a_b", water)
    df = main(_args(work_units="kcal", no_run_data=True))
    row = df[df["edge"] == "FEP_a_b"].iloc[0]
    assert row["n_replicates"] == 3
    assert row["ddG_kcal"] == pytest.approx(2.0)
    assert row["ddG_std_kcal"] == pytest.approx(1.0)
    assert row["ddG_sem_kcal"] == pytest.approx(1.0 / np.sqrt(3))
