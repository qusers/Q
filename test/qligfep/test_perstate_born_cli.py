"""CLI plumbing for the per-state Born correction: setupFEP.create_call must render
the --perstate-born-correction / --born-coefficient flags only when requested."""

from QligFEP.CLI.setupFEP import create_call

BASE = dict(
    lig1="a", lig2="b", FF="AMBER14sb", system="water", cluster="SNELLIUS",
    replicates=10, sampling="sigmoidal", sphereradius=25, start=0.5, windows=10,
    temperature=298, timestep="2fs", rest="hybridization_p", dr_force=0.5, log="info",
)


def test_no_born_flag_by_default():
    assert "--perstate-born-correction" not in create_call(**BASE, perstate_born=False)


def test_perstate_born_flag_rendered_when_enabled():
    assert "--perstate-born-correction" in create_call(**BASE, perstate_born=True)


def test_no_coefficient_flag_without_override():
    call = create_call(**BASE, perstate_born=True, born_coefficient=None)
    assert "--born-coefficient" not in call


def test_coefficient_override_rendered():
    call = create_call(**BASE, perstate_born=True, born_coefficient=7.26)
    assert "--born-coefficient 7.26" in call
