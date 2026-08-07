"""Tests for platform-independent production-run planning and input rendering."""

from pathlib import Path

import pytest

from QligFEP.production_runner import (
    RunnerConfigurationError,
    build_run_plan,
    launcher_command,
    render_runtime_input,
)


def _write_input(
    root: Path,
    name: str,
    final: str,
    restart: str | None,
    energy: str | None = None,
    lambdas: tuple[float, float] = (0.5, 0.5),
) -> None:
    restart_line = f"restart                   {restart}\n" if restart else ""
    energy_line = f"energy                    {energy}\n" if energy else ""
    seed_lines = "random_seed               SEED_VAR\ninitial_temperature       1\n" if restart is None else ""
    (root / f"{name}.inp").write_text(
        f"""[MD]
steps                     10
temperature               T_VAR
{seed_lines}[intervals]
energy                    {10 if energy else 0}
trajectory                100

[files]
topology                  dualtop.top
trajectory                {name}.dcd
{restart_line}{energy_line}final                     {final}
fep                       FEP_VAR

[trajectory_atoms]
not excluded

[lambdas]
{lambdas[0]} {lambdas[1]}
"""
    )


def _write_equilibration(root: Path) -> None:
    previous = None
    for number in range(1, 6):
        final = f"eq{number}.re"
        _write_input(root, f"eq{number}", final, previous)
        previous = final


def test_plans_sequential_endpoint_run_in_one_production_lane(tmp_path):
    _write_equilibration(tmp_path)
    _write_input(tmp_path, "md_1000_0000", "md_1000_0000.re", "eq5.re", "md_1000_0000.en", (1, 0))
    _write_input(
        tmp_path,
        "md_0500_0500",
        "md_0500_0500.re",
        "md_1000_0000.re",
        "md_0500_0500.en",
    )
    _write_input(
        tmp_path,
        "md_0000_1000",
        "md_0000_1000.re",
        "md_0500_0500.re",
        "md_0000_1000.en",
        (0, 1),
    )

    plan = build_run_plan(tmp_path)

    assert [item.stage.name for item in plan.stages] == [
        "eq1",
        "eq2",
        "eq3",
        "eq4",
        "eq5",
        "md_1000_0000",
        "md_0500_0500",
        "md_0000_1000",
    ]
    assert [item.lane for item in plan.stages[:5]] == ["equilibration"] * 5
    assert [item.lane for item in plan.stages[5:]] == ["production"] * 3
    assert plan.terminal_restarts == ("md_0000_1000.re",)


def test_plans_centre_start_branches_depth_first_with_bounded_lanes(tmp_path):
    _write_equilibration(tmp_path)
    _write_input(tmp_path, "md_0500_0500", "md_0500_0500.re", "eq5.re", "md_0500_0500.en")
    _write_input(
        tmp_path,
        "md_0750_0250",
        "md_0750_0250.re",
        "md_0500_0500.re",
        "md_0750_0250.en",
        (0.75, 0.25),
    )
    _write_input(
        tmp_path,
        "md_1000_0000",
        "md_1000_0000.re",
        "md_0750_0250.re",
        "md_1000_0000.en",
        (1, 0),
    )
    _write_input(
        tmp_path,
        "md_0250_0750",
        "md_0250_0750.re",
        "md_0500_0500.re",
        "md_0250_0750.en",
        (0.25, 0.75),
    )
    _write_input(
        tmp_path,
        "md_0000_1000",
        "md_0000_1000.re",
        "md_0250_0750.re",
        "md_0000_1000.en",
        (0, 1),
    )

    plan = build_run_plan(tmp_path)
    production = plan.stages[5:]

    assert [item.stage.name for item in production] == [
        "md_0500_0500",
        "md_0250_0750",
        "md_0000_1000",
        "md_0750_0250",
        "md_1000_0000",
    ]
    assert [item.lane for item in production] == [
        "production",
        "branch-1",
        "branch-1",
        "branch-2",
        "branch-2",
    ]
    assert set(plan.terminal_restarts) == {"md_0000_1000.re", "md_1000_0000.re"}
    assert {item.checkpoint for item in production} == {
        "checkpoint.production.re",
        "checkpoint.branch-1.re",
        "checkpoint.branch-2.re",
    }


def test_runtime_input_uses_shared_inputs_and_disables_trajectory(tmp_path):
    _write_equilibration(tmp_path)
    _write_input(tmp_path, "md_1000_0000", "md_1000_0000.re", "eq5.re", "md_1000_0000.en", (1, 0))
    (tmp_path / "dualtop.top").write_text("topology")
    (tmp_path / "FEP1.fep").write_text("fep")
    stage = build_run_plan(tmp_path).stages[-1]

    rendered = render_runtime_input(
        stage,
        tmp_path,
        temperature=298,
        seed=12345,
        fep_file="FEP1.fep",
        restart_file="checkpoint.equilibration.re",
    )

    assert "T_VAR" not in rendered
    assert "SEED_VAR" not in rendered
    assert "FEP_VAR" not in rendered
    assert f"topology                  {(tmp_path / 'dualtop.top').resolve()}" in rendered
    assert f"fep                       {(tmp_path / 'FEP1.fep').resolve()}" in rendered
    assert "restart                   checkpoint.equilibration.re" in rendered
    assert "final                     checkpoint.production.next.re" in rendered
    assert "trajectory                0" in rendered
    assert "md_1000_0000.dcd" not in rendered
    assert "energy                    md_1000_0000.en" in rendered


def test_runtime_input_can_retain_trajectories(tmp_path):
    _write_equilibration(tmp_path)
    _write_input(tmp_path, "md_1000_0000", "md_1000_0000.re", "eq5.re", "md_1000_0000.en", (1, 0))
    stage = build_run_plan(tmp_path).stages[-1]

    rendered = render_runtime_input(
        stage,
        tmp_path,
        temperature=298,
        seed=1,
        fep_file="FEP1.fep",
        restart_file="checkpoint.equilibration.re",
        trajectory=True,
    )

    assert "trajectory                100" in rendered
    assert "trajectory                md_1000_0000.dcd" in rendered


def test_rejects_missing_restart_dependency(tmp_path):
    _write_equilibration(tmp_path)
    _write_input(tmp_path, "md_1000_0000", "md_1000_0000.re", "missing.re", "md.en", (1, 0))

    with pytest.raises(RunnerConfigurationError, match="not produced by any stage"):
        build_run_plan(tmp_path)


def test_launcher_is_tokenized_without_shell_interpolation():
    assert launcher_command("mpirun -n 8 --map-by core", "/opt/Q/qdynp", "current.inp") == [
        "mpirun",
        "-n",
        "8",
        "--map-by",
        "core",
        "/opt/Q/qdynp",
        "current.inp",
    ]
