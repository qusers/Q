"""CLI regression tests for reproducible residue-FEP seed vectors."""

from argparse import Namespace

from QligFEP.CLI import qresfep_cli, setup_resfep_cli


MAY_SEEDS = [971, 2856, 8253, 1405, 5342, 7486, 5374, 1123, 8397, 8884]


def test_qresfep_accepts_an_explicit_seed_per_replicate(monkeypatch):
    monkeypatch.setattr(
        "sys.argv",
        [
            "qresfep",
            "-m", "LYS27GLY",
            "-mc", "A",
            "-S", "protein",
            "-c", "SNELLIUS",
            "-r", "10",
            "--seeds", *(str(seed) for seed in MAY_SEEDS),
        ],
    )

    args = qresfep_cli.parse_arguments()

    assert args.replicates == 10
    assert args.seeds == MAY_SEEDS
    assert args.random_state is None


def test_series_passes_the_explicit_seed_vector_to_every_qresfep_setup():
    args = Namespace(
        tripeptide_flanks="A",
        windows=50,
        sampling="exponential",
        start="1",
        timestep="2fs",
        temperature="298",
        replicates=10,
        shell_restraint=0.0,
        log="info",
        eq5_steps=2_500_000,
        to_clean=["inp", "re", "top", "dcd"],
        write_trajectories=False,
        separate_scaling=False,
        random_state=None,
        seeds=MAY_SEEDS,
    )

    options = setup_resfep_cli._qresfep_options(args)

    assert options[options.index("--seeds") + 1 :] == [str(seed) for seed in MAY_SEEDS]
    assert "--no-trajectories" in options
    assert "--coupled-thermostat" in options
    assert options[options.index("-eqs") + 1] == "2500000"


