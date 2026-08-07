"""Tests for qdyn execution, rolling checkpoints, and resume behaviour."""

import json
import os
from pathlib import Path

import pytest

from QligFEP.production_runner import ProductionRunner, StageExecutionError


def _write_input(
    root: Path,
    name: str,
    final: str,
    restart: str | None,
    energy: str | None = None,
) -> None:
    restart_line = f"restart {restart}\n" if restart else ""
    energy_line = f"energy {energy}\n" if energy else ""
    seed_line = "random_seed SEED_VAR\ninitial_temperature 1\n" if restart is None else ""
    (root / f"{name}.inp").write_text(
        f"""[MD]
steps 10
temperature T_VAR
{seed_line}[intervals]
energy {10 if energy else 0}
trajectory 100
[files]
topology dualtop.top
trajectory {name}.dcd
{restart_line}{energy_line}final {final}
fep FEP_VAR
[lambdas]
0.5 0.5
"""
    )


def _write_inputs(root: Path) -> None:
    previous = None
    for number in range(1, 6):
        final = f"eq{number}.re"
        _write_input(root, f"eq{number}", final, previous)
        previous = final
    _write_input(root, "md_1000_0000", "md_1000_0000.re", "eq5.re", "md_1000_0000.en")
    _write_input(
        root,
        "md_0500_0500",
        "md_0500_0500.re",
        "md_1000_0000.re",
        "md_0500_0500.en",
    )
    _write_input(
        root,
        "md_0000_1000",
        "md_0000_1000.re",
        "md_0500_0500.re",
        "md_0000_1000.en",
    )
    (root / "dualtop.top").write_text("topology")
    (root / "FEP1.fep").write_text("fep")
    (root / "qfep.inp").write_text(
        "\n".join(
            [
                "3",
                "2 0",
                "0.592 0",
                "10",
                "10",
                "10",
                "0",
                "0",
                "1 0",
                "md_1000_0000.en",
                "md_0500_0500.en",
                "md_0000_1000.en",
            ]
        )
        + "\n"
    )


def _write_fake_qdyn(path: Path) -> None:
    path.write_text(
        """#!/usr/bin/env python3
import os
import struct
import sys
from pathlib import Path

input_path = Path(sys.argv[1])
text = input_path.read_text()
stage = text.splitlines()[0].split(':', 1)[1].strip()
final = None
energy = None
lambdas = []
section = ''
for raw in text.splitlines():
    line = raw.strip()
    if line.startswith('[') and line.endswith(']'):
        section = line[1:-1].lower()
        continue
    fields = line.split(maxsplit=1)
    if section == 'files' and fields:
        if fields[0].lower() == 'final':
            final = fields[1]
        elif fields[0].lower() == 'energy':
            energy = fields[1]
    elif section == 'lambdas' and fields:
        lambdas = [float(value) for value in line.split()]

failure_stage = os.environ.get('FAKE_FAIL_STAGE')
failure_marker = Path('.fake-failure-used')
if stage == failure_stage and not failure_marker.exists():
    failure_marker.write_text(stage)
    print('ABNORMAL TERMINATION of qdyn')
    sys.exit(0)

if final is None:
    print('missing final file', file=sys.stderr)
    sys.exit(2)
Path(final).write_text('checkpoint for ' + stage)
if energy:
    with Path(energy).open('wb') as stream:
        for state, lambda_value in enumerate(lambdas, 1):
            payload = struct.pack('<i15d', state, lambda_value, *[float(state + i) for i in range(14)])
            marker = struct.pack('<i', len(payload))
            stream.write(marker + payload + marker)
        empty = struct.pack('<i', 0)
        stream.write(empty + empty)
print('qdyn terminated normally.')
"""
    )
    path.chmod(0o755)


def _write_fake_qfep(path: Path) -> None:
    path.write_text(
        """#!/usr/bin/env python3
import os
import sys

if os.environ.get('FAKE_QFEP_FAIL'):
    print('qfep terminated abnormally')
    sys.exit(2)
print('# Part 1: Free energy perturbation summary:')
print('fake free-energy result')
"""
    )
    path.chmod(0o755)


def _runner(input_dir: Path, work_dir: Path, qdyn: Path) -> ProductionRunner:
    return ProductionRunner(
        input_dir=input_dir,
        work_dir=work_dir,
        temperature=298,
        seed=1234,
        fep_file="FEP1.fep",
        qdyn=str(qdyn),
    )


def test_executes_with_one_input_log_and_bounded_restart_files(tmp_path):
    input_dir = tmp_path / "inputs"
    work_dir = tmp_path / "replicate"
    input_dir.mkdir()
    _write_inputs(input_dir)
    qdyn = tmp_path / "fake_qdyn.py"
    _write_fake_qdyn(qdyn)

    state = _runner(input_dir, work_dir, qdyn).run_dynamics()

    assert state["status"] == "dynamics_complete"
    assert len(state["completed_stages"]) == 8
    assert (work_dir / "md_0000_1000.re").read_text() == "checkpoint for md_0000_1000"
    assert not list(work_dir.glob("checkpoint.*.re"))
    assert not list(work_dir.glob("checkpoint.*.next.re"))
    assert [path.name for path in work_dir.glob("*.inp")] == ["current.inp"]
    log = (work_dir / "run.log").read_text()
    assert log.count("===== START ") == 8
    assert log.count("normal_termination=True") == 8
    persisted = json.loads((work_dir / "run-state.json").read_text())
    assert persisted == state


def test_resumes_after_missing_normal_termination_marker(tmp_path, monkeypatch):
    input_dir = tmp_path / "inputs"
    work_dir = tmp_path / "replicate"
    input_dir.mkdir()
    _write_inputs(input_dir)
    qdyn = tmp_path / "fake_qdyn.py"
    _write_fake_qdyn(qdyn)
    runner = _runner(input_dir, work_dir, qdyn)
    monkeypatch.setenv("FAKE_FAIL_STAGE", "md_0500_0500")

    with pytest.raises(StageExecutionError, match="normal termination marker=False"):
        runner.run_dynamics()

    failed = json.loads((work_dir / "run-state.json").read_text())
    assert failed["status"] == "failed"
    assert failed["completed_stages"][-1] == "md_1000_0000"
    assert "md_0500_0500" not in failed["completed_stages"]
    completed_before_resume = len(failed["completed_stages"])

    monkeypatch.delenv("FAKE_FAIL_STAGE")
    resumed = runner.run_dynamics()

    assert resumed["status"] == "dynamics_complete"
    log = (work_dir / "run.log").read_text()
    assert log.count("===== START ") == len(resumed["completed_stages"]) + 1
    for stage in failed["completed_stages"]:
        assert log.count(f"===== START {stage} =====") == 1
    assert log.count("===== START md_0500_0500 =====") == 2
    assert len(resumed["completed_stages"]) > completed_before_resume


def test_rejects_resume_with_different_scientific_parameters(tmp_path):
    input_dir = tmp_path / "inputs"
    work_dir = tmp_path / "replicate"
    input_dir.mkdir()
    _write_inputs(input_dir)
    qdyn = tmp_path / "fake_qdyn.py"
    _write_fake_qdyn(qdyn)
    _runner(input_dir, work_dir, qdyn).run_dynamics()

    changed = ProductionRunner(
        input_dir=input_dir,
        work_dir=work_dir,
        temperature=310,
        seed=1234,
        fep_file="FEP1.fep",
        qdyn=str(qdyn),
    )
    with pytest.raises(ValueError, match="parameters differ"):
        changed.run_dynamics()


def test_full_run_converts_then_removes_native_energies(tmp_path):
    input_dir = tmp_path / "inputs"
    work_dir = tmp_path / "replicate"
    input_dir.mkdir()
    _write_inputs(input_dir)
    qdyn = tmp_path / "fake_qdyn.py"
    qfep = tmp_path / "fake_qfep.py"
    _write_fake_qdyn(qdyn)
    _write_fake_qfep(qfep)

    state = _runner(input_dir, work_dir, qdyn).run(str(qfep))

    assert state["status"] == "complete"
    assert state["qfep_complete"] is True
    assert state["binary_energies_removed"] is True
    assert state["energy_conversion"]["frames"] == 3
    assert state["energy_conversion"]["rows"] == 6
    assert (work_dir / "qfep.out").read_text().startswith("# Part 1:")
    energies_csv = (work_dir / "energies.csv").read_text()
    assert "source_file,frame,record_type" in energies_csv
    assert "md_1000_0000.en" in energies_csv
    assert "md_0000_1000.en" in energies_csv
    assert not list(work_dir.glob("*.en"))
    assert not (work_dir / "current.inp").exists()


def test_qfep_failure_preserves_native_energies(tmp_path, monkeypatch):
    input_dir = tmp_path / "inputs"
    work_dir = tmp_path / "replicate"
    input_dir.mkdir()
    _write_inputs(input_dir)
    qdyn = tmp_path / "fake_qdyn.py"
    qfep = tmp_path / "fake_qfep.py"
    _write_fake_qdyn(qdyn)
    _write_fake_qfep(qfep)
    monkeypatch.setenv("FAKE_QFEP_FAIL", "1")

    with pytest.raises(StageExecutionError, match="qfep failed"):
        _runner(input_dir, work_dir, qdyn).run(str(qfep))

    assert len(list(work_dir.glob("*.en"))) == 3
    assert not (work_dir / "energies.csv").exists()
    state = json.loads((work_dir / "run-state.json").read_text())
    assert state["status"] == "failed"
    assert state["qfep_complete"] is False
    assert state["binary_energies_removed"] is False


def test_keep_en_retains_native_energies_after_conversion(tmp_path):
    input_dir = tmp_path / "inputs"
    work_dir = tmp_path / "replicate"
    input_dir.mkdir()
    _write_inputs(input_dir)
    qdyn = tmp_path / "fake_qdyn.py"
    qfep = tmp_path / "fake_qfep.py"
    _write_fake_qdyn(qdyn)
    _write_fake_qfep(qfep)

    state = _runner(input_dir, work_dir, qdyn).run(str(qfep), keep_binary_energies=True)

    assert state["status"] == "complete"
    assert state["binary_energies_removed"] is False
    assert len(list(work_dir.glob("*.en"))) == 3
    assert (work_dir / "energies.csv").is_file()


def test_fake_executables_are_runnable(tmp_path):
    qdyn = tmp_path / "fake_qdyn.py"
    qfep = tmp_path / "fake_qfep.py"
    _write_fake_qdyn(qdyn)
    _write_fake_qfep(qfep)
    assert os.access(qdyn, os.X_OK)
    assert os.access(qfep, os.X_OK)
