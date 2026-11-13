"""Tests for the provenance headers stamped on generated Q .prm/.lib files."""

import getpass
import re
from pathlib import Path

import pytest

from QligFEP import provenance
from QligFEP.IO import get_force_field_paths, read_prm
from QligFEP.provenance import (
    BEGIN_MARKER,
    copy_with_provenance,
    is_comment,
    provenance_header,
    provenance_record,
)
from QligFEP.qresfep import QresFEP


def _labels(header: str) -> list[str]:
    return [m.group(1) for m in re.finditer(r"^\* (\S[^:]*?)\s*:", header, flags=re.M)]


def test_header_is_entirely_q_comments():
    header = provenance_header("tests.generator", details={"Note": "x " * 200})
    lines = header.splitlines()
    assert lines[0] == BEGIN_MARKER
    for line in lines:
        assert line.startswith("*"), line
        assert len(line) <= 100, line  # Q's parser reads at most 400 characters per line


def test_header_records_version_generator_date_and_details(monkeypatch):
    monkeypatch.setenv("SOURCE_DATE_EPOCH", "0")
    header = provenance_header("tests.generator", details={"Force field": "FF-X"})
    assert f"by QligFEP {provenance.__version__}" in header
    assert "* Generator : tests.generator" in header
    assert "* Date      : 1970-01-01, 00:00:00 UTC" in header
    assert _labels(header) == ["Generated", "Generator", "Date", "Force field"]


def test_header_has_no_sensitive_information(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    header = provenance_header("tests.generator")
    assert "end QligFEP provenance" not in header
    assert str(tmp_path) not in header
    assert str(Path.home()) not in header
    assert getpass.getuser() not in header.replace("tests.generator", "")
    assert "sha256" not in header.lower()


@pytest.mark.parametrize(
    "line, expected",
    [("* hi", True), ("  ! hi", True), ("# hi", True), ("{LIG}  ! c", False), ("", False), ("C1 C2", False)],
)
def test_is_comment(line, expected):
    assert is_comment(line) is expected


def test_copy_keeps_original_contents_after_header(tmp_path):
    lib_path, _ = get_force_field_paths("AMBER14sb")
    target = copy_with_provenance(lib_path, tmp_path / "copy.lib", "tests.generator")
    text = target.read_text()
    original = Path(lib_path).read_text()
    assert text.startswith(BEGIN_MARKER)
    assert text.endswith(original)
    assert "* Copied from: AMBER14sb.lib" in text
    assert str(Path(lib_path).parent) not in text  # only the file name, never the path


def test_python_readers_ignore_the_header(tmp_path):
    lib_path, prm_path = get_force_field_paths("AMBER14sb")
    stamped_prm = copy_with_provenance(prm_path, tmp_path / "stamped.prm", "tests.generator")
    assert read_prm([str(stamped_prm)]) == read_prm([prm_path])

    stamped_lib = copy_with_provenance(lib_path, tmp_path / "stamped.lib", "tests.generator")
    wanted = {"ALA", "LEU"}
    assert QresFEP._read_library_entries(stamped_lib, wanted) == QresFEP._read_library_entries(
        Path(lib_path), wanted
    )


def test_git_description_is_safe_outside_a_checkout(monkeypatch):
    provenance._git_description.cache_clear()

    def boom(*args, **kwargs):
        raise OSError("no git")

    monkeypatch.setattr(provenance.subprocess, "run", boom)
    try:
        assert provenance._git_description() is None
        assert "QligFEP" in provenance_header("tests.generator")
    finally:
        provenance._git_description.cache_clear()


def test_record_matches_the_header_and_is_json_serialisable():
    import json

    provenance._git_description.cache_clear()
    record = provenance_record()
    assert json.loads(json.dumps(record)) == record
    assert record["QligFEP_version"] == provenance.__version__
    git = provenance._git_description()
    if git is None:
        assert record["git_commit"] is None and record["uncommitted_changes"] is None
    else:
        assert re.fullmatch(r"[0-9a-f]{12}", record["git_commit"])
        assert record["uncommitted_changes"] is ("uncommitted" in git)
        assert f"git {record['git_commit']}" in provenance_header("tests.generator")


def test_record_is_safe_outside_a_checkout(monkeypatch):
    provenance._git_description.cache_clear()
    monkeypatch.setattr(provenance, "_git_description", lambda: None)
    assert provenance_record() == {
        "QligFEP_version": provenance.__version__,
        "git_commit": None,
        "uncommitted_changes": None,
    }
