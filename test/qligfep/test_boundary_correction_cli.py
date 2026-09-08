"""CLI plumbing for endpoint-resolved spherical-boundary corrections."""

from QligFEP.CLI.setupFEP import create_call

BASE = dict(
    lig1="a",
    lig2="b",
    FF="AMBER14sb",
    system="water",
    cluster="SNELLIUS",
    replicates=10,
    sampling="sigmoidal",
    sphereradius=25,
    start=0.5,
    windows=10,
    temperature=298,
    timestep="2fs",
    rest="hybridization_p",
    dr_force=0.5,
    log="info",
)


def test_boundary_flags_absent_by_default():
    call = create_call(**BASE)
    assert "--perstate-polarization" not in call
    assert "--perstate-born-correction" not in call
    assert "--born-dielectric" not in call
    assert "--born-coefficient" not in call


def test_boundary_flags_render_when_requested():
    call = create_call(
        **BASE,
        perstate_polarization=True,
        perstate_born=True,
        born_dielectric=60.0,
        born_coefficient=7.1,
    )
    assert "--perstate-polarization" in call
    assert "--perstate-born-correction" in call
    assert "--born-dielectric 60.0" in call
    assert "--born-coefficient 7.1" in call
