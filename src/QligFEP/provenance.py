"""Provenance headers for the Q parameter (.prm) and library (.lib) files QligFEP writes.

Every .prm/.lib file produced by QligFEP starts with a short comment block that
records what generated it: the QligFEP version (and git commit), the function that
wrote the file, the date, and a few non-sensitive details such as the force field.
No user names, paths or checksums are recorded. The block uses ``*`` comment lines, which Q's
Fortran parser (``prmfile.f90``) and QligFEP's own readers already skip, so it has
no effect on how the file is read.

``provenance_record`` returns the version and git commit as a dict, for the JSON run
configuration files (``fep_config.json``, ``resfep_config.json``).

Set ``SOURCE_DATE_EPOCH`` to pin the recorded timestamp (e.g. for reproducible
output in tests).
"""

from __future__ import annotations

import datetime as _dt
import os
import subprocess
import textwrap
from collections.abc import Mapping
from functools import lru_cache
from pathlib import Path
from typing import Optional, Union

from . import __version__

PathLike = Union[str, "os.PathLike[str]"]

BEGIN_MARKER = "*---- QligFEP provenance " + "-" * 55
_UNCOMMITTED = "(uncommitted changes)"  # the source tree had modified tracked files
_WIDTH = 100  # Q reads 400-character lines; stay far below that.


@lru_cache(maxsize=1)
def _git_description() -> Optional[str]:
    """Return ``<short-hash>[ (uncommitted changes)]`` for the QligFEP source tree, if it is a git checkout."""
    package_dir = Path(__file__).resolve().parent
    try:
        commit = subprocess.run(
            ["git", "-C", str(package_dir), "rev-parse", "--short=12", "HEAD"],
            capture_output=True,
            text=True,
            timeout=5,
            check=True,
        ).stdout.strip()
        status = subprocess.run(
            ["git", "-C", str(package_dir), "status", "--porcelain", "--untracked-files=no"],
            capture_output=True,
            text=True,
            timeout=5,
            check=True,
        ).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return None
    if not commit:
        return None
    return f"{commit} {_UNCOMMITTED}" if status else commit


def _timestamp() -> str:
    epoch = os.environ.get("SOURCE_DATE_EPOCH")
    if epoch is not None:
        try:
            moment = _dt.datetime.fromtimestamp(int(epoch), tz=_dt.timezone.utc)
        except ValueError:
            moment = _dt.datetime.now(tz=_dt.timezone.utc)
    else:
        moment = _dt.datetime.now(tz=_dt.timezone.utc)
    return moment.strftime("%Y-%m-%d, %H:%M:%S UTC")


def _field(label: str, value: str) -> list[str]:
    """Format one ``label: value`` entry, wrapping long values onto continuation lines."""
    prefix = f"* {label:<10}: "
    value = " ".join(str(value).split())  # no embedded newlines or tabs
    lines = textwrap.wrap(
        value,
        width=_WIDTH,
        initial_indent=prefix,
        subsequent_indent="*" + " " * (len(prefix) - 1),
        break_long_words=True,
        break_on_hyphens=False,
    )
    return lines or [prefix.rstrip()]


def provenance_header(generator: str, details: Optional[Mapping[str, object]] = None) -> str:
    """Return the provenance comment block for a generated Q .prm/.lib file.

    Args:
        generator: the function that wrote the file, e.g. ``"QligFEP.openff2Q.OpenFF2Q.write_lib_Q"``.
        details: extra ``label -> value`` entries (force field, charge model, ...). Keep these
            free of user names and file system paths.

    Returns:
        The header text, ending with a newline.
    """
    git = _git_description()
    version = f"QligFEP {__version__}" + (f", git {git}" if git else "")
    lines = [BEGIN_MARKER]
    lines += _field("Generated", f"by {version}")
    lines += _field("Generator", generator)
    lines += _field("Date", _timestamp())
    for label, value in (details or {}).items():
        lines += _field(str(label), str(value))
    return "\n".join(lines) + "\n"


def provenance_record() -> dict[str, object]:
    """Return the same provenance as a JSON-serialisable dict, for run configuration files.

    ``git_commit`` and ``uncommitted_changes`` are ``None`` when QligFEP is not running from a git
    checkout (e.g. a regular, non-editable install); ``QligFEP_version`` is always set.
    """
    git = _git_description()
    return {
        "QligFEP_version": __version__,
        "git_commit": git.split()[0] if git else None,
        "uncommitted_changes": git.endswith(_UNCOMMITTED) if git else None,
    }


def is_comment(line: str) -> bool:
    """True for lines Q treats as whole-line comments (first non-blank char is ``!``, ``#`` or ``*``)."""
    stripped = line.lstrip()
    return bool(stripped) and stripped[0] in "!#*"


def copy_with_provenance(
    source: PathLike,
    destination: PathLike,
    generator: str,
    details: Optional[Mapping[str, object]] = None,
) -> Path:
    """Copy a .prm/.lib file, prepending a provenance header.

    The original contents (including any header the source already carries) follow
    unchanged, so the file's lineage is preserved.
    """
    source, destination = Path(source), Path(destination)
    contents = source.read_text()
    header = provenance_header(generator, details={"Copied from": source.name, **(details or {})})
    destination.write_text(header + contents)
    return destination
