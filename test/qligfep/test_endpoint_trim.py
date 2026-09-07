"""Explicitly truncated BAR must not masquerade as a full-endpoint result."""
from decimal import Decimal
from pathlib import Path
import math
import os
import struct
import subprocess

import pytest

from QligFEP.endpoint_trim import prepare, summarize


def campaign(tmp_path, labels=("1", ".9999", ".5", ".0001", "0")):
    md = tmp_path / "md"
    md.mkdir()
    energies = []
    for i, label in enumerate(labels):
        name = f"arbitrary_name_{i}.en"
        (md / f"window_{i}.inp").write_text(
            f"[files]\nenergy {name}\n[lambdas]\n{label} {1-Decimal(label)} ! mapping\n"
        )
        (tmp_path / name).write_bytes(b"unchanged energy fixture")
        energies.append(name)
    inp = tmp_path / "qfep.inp"
    inp.write_text("\n".join([str(len(labels)), "2 0", ".596 0 0", "5", "1", "0", "0", "1 -1", *energies])+"\n")
    return inp, md, tmp_path / "trimmed"


def write_summary(out, rows=None):
    rows = rows or [".999900 0.000 0.000", ".500000 1.000 1.000", ".000100 2.000 3.000"]
    (out / "qfep.out").write_text("# Part 6: BAR Bennet:\n# lambda(1) dG sum(dG)\n"+"\n".join(rows)+"\n")


def test_default_retains_nearest_interiors_and_preserves_originals(tmp_path):
    inp, md, out = campaign(tmp_path)
    before = {p: p.read_bytes() for p in tmp_path.rglob("*") if p.is_file()}
    result = prepare(inp, md, out)
    assert [w["lambdas"][0] for w in result["retained"]] == ["0.9999", "0.5", "0.0001"]
    assert [w["excluded_reason"] for w in result["excluded"]] == ["exact_endpoint"] * 2
    assert result["state_constant_gap_multiplier"] == "0.9998"
    assert result["full_endpoint_free_energy"] is False
    assert all(p.read_bytes() == data for p, data in before.items())
    assert (out / "window_0000.en").is_symlink()
    assert (out / "qfep.inp").read_text().splitlines()[1:8] == inp.read_text().splitlines()[1:8]
    assert (out / "qfep.inp").read_text().splitlines()[0] == "3"
    write_summary(out)
    summary = summarize(out)
    assert summary["bar_kcal_mol"] == 3
    assert summary["retained_pairs"] == 2
    assert summary["full_endpoint_free_energy"] is False


@pytest.mark.parametrize("labels,multiplier", [
    (("0", ".001", ".999", "1"), "-0.998"),
    (("1", ".99", ".003", "0"), "0.987"),
])
def test_direction_and_asymmetric_constant_shift(tmp_path, labels, multiplier):
    inp, md, out = campaign(tmp_path, labels)
    assert prepare(inp, md, out)["state_constant_gap_multiplier"] == multiplier


def test_explicit_shared_interval(tmp_path):
    inp, md, out = campaign(tmp_path, ("1", ".9999", ".999", ".5", ".001", ".0001", "0"))
    result = prepare(inp, md, out, ".001", ".999")
    assert len(result["retained"]) == 3
    assert result["state_constant_gap_multiplier"] == "0.998"
    assert len(result["excluded"]) == 4


def test_excluded_files_need_not_be_available(tmp_path):
    inp, md, out = campaign(tmp_path)
    (tmp_path / "arbitrary_name_0.en").unlink()
    (tmp_path / "arbitrary_name_4.en").unlink()
    assert len(prepare(inp, md, out)["retained"]) == 3


@pytest.mark.parametrize("lower,upper", [(".001", None), ("0", "1"), (".7", ".3"), (".001", ".999"), ("NaN", ".999")])
def test_invalid_or_unsampled_interval_fails_without_creating_output(tmp_path, lower, upper):
    inp, md, out = campaign(tmp_path)
    with pytest.raises(ValueError):
        prepare(inp, md, out, lower, upper)
    assert not out.exists()


@pytest.mark.parametrize("labels", [("1", ".5", "0"), ("1", ".1", ".8", "0"), ("1", ".5", ".5", "0")])
def test_invalid_ladder(tmp_path, labels):
    inp, md, out = campaign(tmp_path, labels)
    with pytest.raises(ValueError):
        prepare(inp, md, out)
    assert not out.exists()


@pytest.mark.parametrize("mapping", [".1 .1", "-.1 1.1", ".3333333 .6666667", ".5 .5 0"])
def test_invalid_md_mapping(tmp_path, mapping):
    inp, md, out = campaign(tmp_path)
    (md / "window_1.inp").write_text(f"[files]\nenergy arbitrary_name_1.en\n[lambdas]\n{mapping}\n")
    with pytest.raises(ValueError):
        prepare(inp, md, out)


@pytest.mark.parametrize("problem", ["missing_md", "duplicate_md", "missing_energy", "empty_energy", "gas", "coupling", "empty_input"])
def test_preparation_failures(tmp_path, problem):
    inp, md, out = campaign(tmp_path)
    if problem == "missing_md":
        (md / "window_2.inp").unlink()
    elif problem == "duplicate_md":
        (md / "copy.inp").write_bytes((md / "window_2.inp").read_bytes())
    elif problem == "missing_energy":
        (tmp_path / "arbitrary_name_2.en").unlink()
    elif problem == "empty_energy":
        (tmp_path / "arbitrary_name_2.en").write_bytes(b"")
    elif problem == "gas":
        inp.write_text(inp.read_text().replace(".596 0 0", ".596 0 1"))
    elif problem == "coupling":
        inp.write_text(inp.read_text().replace("2 0", "2 1", 1))
    else:
        inp.write_text("")
    with pytest.raises(ValueError):
        prepare(inp, md, out)
    assert not out.exists()


def test_no_overwrite(tmp_path):
    inp, md, out = campaign(tmp_path)
    out.mkdir()
    with pytest.raises(FileExistsError):
        prepare(inp, md, out)


def test_separate_energy_directory_preserves_original_input(tmp_path):
    inp, md, out = campaign(tmp_path)
    original = inp.read_bytes()
    nested_input = md / "qfep.inp"
    nested_input.write_bytes(original)
    with pytest.raises(ValueError, match="Missing or empty retained"):
        prepare(nested_input, md, out)
    result = prepare(nested_input, md, out, energy_dir=tmp_path)
    assert nested_input.read_bytes() == original
    assert result["energy_base"] == str(tmp_path.resolve())
    assert (out / "window_0000.en").resolve() == tmp_path / "arbitrary_name_1.en"


@pytest.mark.parametrize("problem", ["truncated", "wrong_lambda", "inconsistent", "nonfinite", "overflow", "nonzero_origin", "duplicate", "input_changed", "energy_changed"])
def test_summary_rejects_misleading_results(tmp_path, problem):
    inp, md, out = campaign(tmp_path)
    prepare(inp, md, out)
    rows = [".999900 0 0", ".500000 1 1", ".000100 2 3"]
    if problem == "truncated":
        rows.pop()
    elif problem == "wrong_lambda":
        rows[-1] = ".000000 2 3"
    elif problem == "inconsistent":
        rows[-1] = ".000100 2 4"
    elif problem == "nonfinite":
        rows[-1] = ".000100 NaN NaN"
    elif problem == "overflow":
        rows[-1] = ".000100 ********* *********"
    elif problem == "nonzero_origin":
        rows[0] = ".999900 1 1"
    elif problem == "duplicate":
        rows += ["# Part 6: BAR Bennet:", *rows]
    elif problem == "input_changed":
        (out / "qfep.inp").write_text("modified")
    else:
        (out / "window_0000.en").write_bytes(b"changed")
    write_summary(out, rows)
    with pytest.raises(ValueError):
        summarize(out)


def test_real_qfep_discrete_model(tmp_path):
    """Optional native-binary integration against an analytic partition function.

    Two configurations have U1=0, U2=c +/- d. Each mapping's empirical counts
    approximate its exact Boltzmann distribution; no fitting to QFEP output.
    Set QFEP_ENDPOINT_TEST_BINARY to opt in (uses Q's native record format).
    """
    binary = os.environ.get("QFEP_ENDPOINT_TEST_BINARY")
    if not binary:
        pytest.skip("Set QFEP_ENDPOINT_TEST_BINARY for native integration")
    inp, md, out = campaign(tmp_path)
    kt, c, d, count = .596, 1.2, .2, 10000

    def record(payload):
        marker = struct.pack("=i", len(payload))
        return marker+payload+marker

    for i, lambda1 in enumerate((1., .9999, .5, .0001, 0.)):
        if i in (0, 4):
            # Deliberately invalid exact-endpoint files must never be opened.
            continue
        lambda2 = 1-lambda1
        nplus = round(count/(1+math.exp(2*lambda2*d/kt)))
        data = bytearray()
        for frame in range(count):
            u2 = c + (d if frame < nplus else -d)
            for state, weight, potential in [(1, lambda1, 0.), (2, lambda2, u2)]:
                data.extend(record(struct.pack("=i15d", state, weight, potential, *([0.]*13))))
            data.extend(record(b""))
        (tmp_path / f"arbitrary_name_{i}.en").write_bytes(data)
    prepare(inp, md, out)
    with (out / "qfep.inp").open() as source:
        result = subprocess.run([str(Path(binary).resolve())], stdin=source, cwd=out,
                                capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stderr
    (out / "qfep.out").write_text(result.stdout)
    summary = summarize(out)
    free_energy = lambda weight: c*weight-kt*math.log(math.cosh(weight*d/kt))
    expected = free_energy(.9999)-free_energy(.0001)
    assert summary["bar_kcal_mol"] == pytest.approx(expected, abs=.002)
