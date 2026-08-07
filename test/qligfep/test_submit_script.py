"""Tests for the POSIX edge-level Slurm submission adapter."""

import os
import subprocess
from pathlib import Path

import QligFEP


def test_submit_script_records_job_and_prevents_duplicate_submission(tmp_path):
    edge_dir = tmp_path / "FEP_lig1_lig2"
    input_dir = edge_dir / "inputfiles"
    input_dir.mkdir(parents=True)
    (input_dir / "runHABROK.sh").write_text("#!/bin/bash\n")

    template = Path(QligFEP.__file__).parent / "INPUTS" / "FEP_submit.sh"
    submit_script = edge_dir / "FEP_submit.sh"
    submit_script.write_text(template.read_text().replace("RUNFILE", "runHABROK.sh"))

    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    fake_sbatch = fake_bin / "sbatch"
    fake_sbatch.write_text(
        "#!/bin/sh\n"
        'printf "%s\\n" "$*" >> "$FAKE_SBATCH_CALLS"\n'
        "printf '123456;cluster\\n'\n"
    )
    fake_sbatch.chmod(0o755)

    calls = tmp_path / "sbatch-calls.txt"
    env = {
        **os.environ,
        "PATH": f"{fake_bin}:{os.environ['PATH']}",
        "FAKE_SBATCH_CALLS": str(calls),
    }
    first = subprocess.run(
        ["sh", str(submit_script)],
        cwd=tmp_path,
        env=env,
        check=True,
        capture_output=True,
        text=True,
    )
    second = subprocess.run(
        ["sh", str(submit_script)],
        cwd=tmp_path,
        env=env,
        check=True,
        capture_output=True,
        text=True,
    )
    selected_retry = subprocess.run(
        ["sh", str(submit_script), "2", "5"],
        cwd=tmp_path,
        env=env,
        check=True,
        capture_output=True,
        text=True,
    )

    assert "Submitted Slurm job 123456;cluster" in first.stdout
    assert "Already submitted as job 123456;cluster" in second.stdout
    assert "Submitted Slurm job 123456;cluster" in selected_retry.stdout
    assert (edge_dir / "submission-jobid.txt").read_text() == "123456;cluster\n"
    assert calls.read_text().splitlines() == [
        f"--parsable {input_dir / 'runHABROK.sh'}",
        f"--parsable --array=2,5 {input_dir / 'runHABROK.sh'}",
    ]
    history = (edge_dir / "submission-history.tsv").read_text().splitlines()
    assert [line.split("\t")[1:] for line in history] == [
        ["all", "123456;cluster"],
        ["2,5", "123456;cluster"],
    ]
