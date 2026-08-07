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


def _write_fake_qdyn(path: Path) -> None:
    path.write_text(
        """#!/usr/bin/env python3
import os
import sys
from pathlib import Path

input_path = Path(sys.argv[1])
text = input_path.read_text()
stage = text.splitlines()[0].split(':', 1)[1].strip()
final = None
section = ''
for raw in text.splitlines():
    line = raw.strip()
    if line.startswith('[') and line.endswith(']'):
        section = line[1:-1].lower()
        continue
    fields = line.split(maxsplit=1)
    if section == 'files' and fields and fields[0].lower() == 'final':
        final = fields[1]

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
print('qdyn terminated normally.')
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


def test_fake_qdyn_is_executable(tmp_path):
    qdyn = tmp_path / "fake_qdyn.py"
    _write_fake_qdyn(qdyn)
    assert os.access(qdyn, os.X_OK)
