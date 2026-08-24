"""CLI regression tests for residue-FEP protocol controls."""

from argparse import Namespace

import pytest

from QligFEP.CLI import qresfep_cli, setup_resfep_cli
from QligFEP.resfep_protocols import (
    MANUSCRIPT_PREPARATION_SETTINGS,
    MANUSCRIPT_SEEDS,
    MANUSCRIPT_SETTINGS,
)


def test_qresfep_accepts_an_explicit_seed_per_replicate(monkeypatch):
    monkeypatch.setattr(
        "sys.argv",
        [
            "qresfep",
            "-m",
            "LYS27GLY",
            "-mc",
            "A",
            "-S",
            "protein",
            "-c",
            "SNELLIUS",
            "-r",
            "10",
            "--seeds",
            *(str(seed) for seed in MANUSCRIPT_SEEDS),
        ],
    )

    args = qresfep_cli.parse_arguments()

    assert args.replicates == 10
    assert args.seeds == list(MANUSCRIPT_SEEDS)
    assert args.random_state is None


@pytest.mark.parametrize("value", ["on", "off"])
def test_qresfep_accepts_the_fortran_separate_scaling_values(monkeypatch, value):
    monkeypatch.setattr(
        "sys.argv",
        [
            "qresfep",
            "-m",
            "LYS27GLY",
            "-mc",
            "A",
            "-S",
            "protein",
            "-c",
            "SNELLIUS",
            "--separate-scaling",
            value,
        ],
    )

    assert qresfep_cli.parse_arguments().separate_scaling == value


def test_qresfep_manuscript_settings_force_the_protocol(monkeypatch):
    monkeypatch.setattr(
        "sys.argv",
        [
            "qresfep",
            "-m",
            "LYS27GLY",
            "-mc",
            "A",
            "-S",
            "protein",
            "-c",
            "SNELLIUS",
            "-w",
            "12",
            "-r",
            "3",
            "-eqs",
            "100",
            "-ps",
            "200",
            "--separate-scaling",
            "on",
            "--seeds",
            "1",
            "2",
            "3",
            "--manuscript-settings",
        ],
    )

    args = qresfep_cli.parse_arguments()

    for name, expected in MANUSCRIPT_SETTINGS.items():
        if isinstance(expected, tuple):
            expected = list(expected)
        assert getattr(args, name) == expected


def test_setup_manuscript_settings_include_preparation(monkeypatch):
    monkeypatch.setattr(
        "sys.argv",
        [
            "setup_resFEP",
            "-i",
            "1A2P.pdb",
            "-M",
            "mutations.txt",
            "-mc",
            "A",
            "-c",
            "HABROK",
            "-r",
            "18",
            "--legs",
            "protein",
            "--manuscript-settings",
        ],
    )

    args = setup_resfep_cli.parse_arguments()

    for name, expected in MANUSCRIPT_PREPARATION_SETTINGS.items():
        if isinstance(expected, tuple):
            expected = list(expected)
        assert getattr(args, name) == expected
    assert args.seeds == list(MANUSCRIPT_SEEDS)


def test_series_passes_flexible_protocol_options_to_every_qresfep_setup():
    args = Namespace(
        tripeptide_flanks="A",
        windows=50,
        sampling="exponential",
        start="1",
        timestep="2fs",
        temperature="298",
        replicates=10,
        log="info",
        eq5_steps=2_500_000,
        production_steps=10_000,
        to_clean=["inp", "re", "top", "dcd"],
        write_trajectories=False,
        separate_scaling="off",
        manuscript_settings=False,
        random_state=None,
        seeds=list(MANUSCRIPT_SEEDS),
    )

    options = setup_resfep_cli._qresfep_options(args)

    seed_index = options.index("--seeds") + 1
    assert options[seed_index : seed_index + len(MANUSCRIPT_SEEDS)] == [
        str(seed) for seed in MANUSCRIPT_SEEDS
    ]
    assert "--no-trajectories" in options
    assert "--manuscript-settings" not in options
    assert options[options.index("--separate-scaling") + 1] == "off"
    assert options[options.index("-eqs") + 1] == "2500000"
    assert options[options.index("-ps") + 1] == "10000"


def test_series_passes_the_manuscript_preset_without_repeating_its_values():
    args = Namespace(
        log="info",
        manuscript_settings=True,
        to_clean=["dcd"],
        write_trajectories=False,
    )

    options = setup_resfep_cli._qresfep_options(args)

    assert options == [
        "-log",
        "info",
        "--manuscript-settings",
        "-clean",
        "dcd",
        "--no-trajectories",
    ]
