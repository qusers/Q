"""Tests for opt-in FIRE minimization in the FEP setup CLIs."""

import sys

import pytest

from QligFEP.CLI.parser_base import parse_arguments
from QligFEP.CLI.setupFEP import create_call
from QligFEP.templates.equilibration import get_equilibration_configs


@pytest.mark.parametrize('timestep', ['1fs', '2fs'])
@pytest.mark.parametrize('minimize', [False, True])
def test_minimization_and_boundary_overrides_survive_integration(timestep, minimize):
    configs = get_equilibration_configs(
        timestep, 14, minimize=minimize, perstate_polarization=True,
        polarization_adaptation=False, perstate_born=True)
    assert len(configs) == 5
    assert configs[0].params.minimize is minimize
    assert configs[0].params.constrain_hydrogens is minimize
    for config in configs:
        assert config.params.perstate_polarization is True
        assert config.params.polarization_adaptation is False
        assert config.params.perstate_born is True
        assert config.params.shell_radius == 14


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
