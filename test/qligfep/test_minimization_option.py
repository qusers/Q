"""Tests for opt-in FIRE minimization in the FEP setup CLIs."""

import sys

import pytest

from QligFEP.CLI.parser_base import parse_arguments
from QligFEP.CLI.setupFEP import create_call


@pytest.mark.parametrize(("extra_args", "expected"), [([], False), (["--minimize"], True)])
def test_qligfep_minimize_flag(monkeypatch, extra_args, expected):
    """The qligfep CLI leaves minimization off unless explicitly requested."""
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "qligfep",
            "-l1",
            "lig1",
            "-l2",
            "lig2",
            "-FF",
            "AMBER14sb",
            "-s",
            "water",
            "-c",
            "LOCAL",
            *extra_args,
        ],
    )

    assert parse_arguments("QligFEP").minimize is expected


def test_setupfep_forwards_minimize_flag():
    """setupFEP includes the opt-in flag in each generated qligfep call."""
    kwargs = {
        "lig1": "lig1",
        "lig2": "lig2",
        "FF": "AMBER14sb",
        "system": "water",
        "cluster": "LOCAL",
        "replicates": "1",
        "sampling": "sigmoidal",
        "sphereradius": "25",
        "start": "0.5",
        "windows": "10",
        "temperature": "298",
        "timestep": "2fs",
        "rest": "heavyatom_p",
        "dr_force": 0.5,
        "log": "info",
    }

    assert "--minimize" not in create_call(**kwargs, minimize=False)
    assert "--minimize" in create_call(**kwargs, minimize=True)
