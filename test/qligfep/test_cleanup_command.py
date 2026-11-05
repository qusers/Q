"""Cleanup must not fail when a requested pattern matches nothing (e.g. -clean dcd
with trajectories disabled); run scripts use `set -e`."""

import subprocess

from QligFEP.qligfep import cleanup_command


def test_cleanup_tolerates_missing_patterns_and_removes_matches(tmp_path):
    for name in ("md.re", "eq1.inp", "md.en", "md.log"):
        (tmp_path / name).write_text("x")
    script = "set -eo pipefail\n" + cleanup_command(["dcd", "re", "inp"]) + "echo done\n"
    result = subprocess.run(["bash", "-c", script], cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert sorted(p.name for p in tmp_path.iterdir()) == ["md.en", "md.log"]
