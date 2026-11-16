"""Tests for the YAML configuration layer: profiles, protocols and stage overrides."""

from __future__ import annotations

import argparse
import textwrap
from pathlib import Path

import pytest
import yaml

from QligFEP.CLI import config_cli
from QligFEP.CLI.config_args import resolve_run_settings
from QligFEP.config import ConfigError, StageOverride, load_config
from QligFEP.config.detect import render_user_config, requested_modules, write_build_info
from QligFEP.config.schema import parse_duration_fs
from QligFEP.qligfep import QligFEP
from QligFEP.templates import get_production_config, render_md_input
from QligFEP.templates.overrides import apply_stage_override


@pytest.fixture(autouse=True)
def isolated_config(tmp_path, monkeypatch):
    """Keep the developer's own user file and environment out of every test."""
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "xdg"))
    monkeypatch.delenv("QLIGFEP_CONFIG", raising=False)
    monkeypatch.delenv("QLIGFEP_Q_BIN", raising=False)
    monkeypatch.chdir(tmp_path)


def write_yaml(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(textwrap.dedent(text))
    return path


def user_file(tmp_path: Path) -> Path:
    return tmp_path / "xdg" / "qligfep" / "config.yaml"


# ---------------------------------------------------------------------------------
# Profiles
# ---------------------------------------------------------------------------------


class TestPackagedProfiles:
    def test_historical_clusters_are_still_available(self):
        names = load_config().profile_names()
        for name in ["CSB", "SNELLIUS", "TETRA", "DARDEL", "HABROK", "LOCAL", "LOCALP"]:
            assert name in names

    def test_site_quirks_are_data_not_code(self):
        config = load_config()
        assert config.get_profile("SNELLIUS").slurm.partition == "rome"
        assert config.get_profile("DARDEL").slurm.partition == "shared"
        # CSB does not accept --mem-per-cpu
        assert "--mem-per-cpu" not in config.get_profile("CSB").sbatch_header("equilibrium", "j", 3)

    def test_no_packaged_profile_carries_a_project_account(self):
        config = load_config()
        for name in config.profile_names():
            assert config.get_profile(name).slurm.account is None, name

    def test_local_profiles(self):
        config = load_config()
        local, localp = config.get_profile("LOCAL"), config.get_profile("LOCALP")
        assert local.is_local and not local.use_mpi and local.resources.ntasks == 1
        assert localp.is_local and localp.use_mpi and localp.resources.ntasks >= 1

    def test_no_profile_and_no_default_is_explained(self):
        with pytest.raises(ConfigError, match="No execution profile selected"):
            load_config().get_profile(None)

    def test_unknown_profile_suggests_the_closest_name(self):
        with pytest.raises(ConfigError, match="Did you mean 'SNELLIUS'"):
            load_config().get_profile("SNELIUS")


class TestProfileLayering:
    def test_user_file_changes_only_the_keys_it_states(self, tmp_path):
        write_yaml(
            user_file(tmp_path),
            """
            default_profile: TETRA
            profiles:
              TETRA: {slurm: {account: my-project}}
            """,
        )
        profile = load_config().get_profile()
        assert profile.name == "TETRA"
        assert profile.slurm.account == "my-project"
        assert profile.resources.ntasks == 8  # untouched packaged value

    def test_project_file_overrides_user_file(self, tmp_path):
        write_yaml(user_file(tmp_path), "profiles: {HABROK: {resources: {ntasks: 12}}}")
        write_yaml(tmp_path / "qligfep.yaml", "profiles: {HABROK: {resources: {ntasks: 4}}}")
        assert load_config().get_profile("HABROK").resources.ntasks == 4

    def test_extends_and_per_system_overrides(self, tmp_path):
        write_yaml(
            tmp_path / "qligfep.yaml",
            """
            profiles:
              mine:
                extends: HABROK
                resources: {ntasks: 16}
                systems:
                  water: {resources: {ntasks: 4}}
            """,
        )
        config = load_config()
        assert config.get_profile("mine").modules == ("OpenMPI/4.1.4-GCC-11.3.0",)
        assert config.get_profile("mine", system="protein").resources.ntasks == 16
        assert config.get_profile("mine", system="water").resources.ntasks == 4

    def test_circular_extends_is_reported(self, tmp_path):
        write_yaml(tmp_path / "qligfep.yaml", "profiles: {a: {extends: b}, b: {extends: a}}")
        with pytest.raises(ConfigError, match="circular"):
            load_config().get_profile("a")

    def test_typo_in_a_key_fails_loudly_and_names_the_file(self, tmp_path):
        path = write_yaml(tmp_path / "qligfep.yaml", "profiles: {x: {resources: {ntaks: 4}}}")
        with pytest.raises(ConfigError, match=r"Did you mean 'ntasks'") as error:
            load_config()
        assert str(path) in str(error.value)

    def test_auto_cores_need_local_execution(self, tmp_path):
        write_yaml(tmp_path / "qligfep.yaml", "profiles: {x: {resources: {ntasks: auto}}}")
        with pytest.raises(ConfigError, match="only valid for local execution"):
            load_config().get_profile("x")

    def test_user_file_must_not_define_protocols(self, tmp_path):
        write_yaml(user_file(tmp_path), "protocols: {fast: {replicates: 2}}")
        with pytest.raises(ConfigError, match="not read from the user file"):
            load_config()


class TestEngineLocation:
    def test_package_binaries_by_default(self):
        profile = load_config().get_profile("HABROK")
        assert profile.q_bin_source == "package"
        assert profile.q_bin.parts[-3:] == ("q6", "bin", "q6")

    def test_profile_env_and_cli_precedence(self, tmp_path, monkeypatch):
        write_yaml(tmp_path / "qligfep.yaml", f"profiles: {{dev: {{extends: HABROK, q_bin: {tmp_path}/a}}}}")
        config = load_config()
        assert config.get_profile("dev").q_bin == tmp_path / "a"
        monkeypatch.setenv("QLIGFEP_Q_BIN", str(tmp_path / "b"))
        assert config.get_profile("dev").q_bin == (tmp_path / "b").resolve()
        assert config.get_profile("dev", q_bin=tmp_path / "c").q_bin == (tmp_path / "c").resolve()

    def test_modules_from_build_follow_the_binaries(self, tmp_path):
        q_bin = tmp_path / "clone" / "bin"
        q_bin.mkdir(parents=True)
        write_build_info(q_bin, modules=["foss/2023a"])
        write_yaml(
            tmp_path / "qligfep.yaml",
            f"profiles: {{dev: {{extends: HABROK, q_bin: {q_bin}, modules: from-build}}}}",
        )
        profile = load_config().get_profile("dev")
        assert profile.modules == ("foss/2023a",)
        assert profile.environment_block() == "module purge\nmodule load foss/2023a\n"

    def test_only_requested_modules_are_recorded_with_lmod(self, monkeypatch):
        import base64

        table = """_ModuleTable_ = { mT = {
        GCCcore = { fn = "x", fullName = "GCCcore/11.3.0", loadOrder = 2, propT = {}, stackDepth = 2, status = "active", },
        OpenMPI = { fn = "x", fullName = "OpenMPI/4.1.4-GCC-11.3.0", loadOrder = 3, propT = {}, stackDepth = 0, },
        StdEnv = { fn = "x", fullName = "StdEnv", loadOrder = 1, propT = {}, stackDepth = 0, },
        }, }"""
        encoded = base64.b64encode(table.encode()).decode()
        half = len(encoded) // 2
        monkeypatch.setenv("_ModuleTable_Sz_", "2")
        monkeypatch.setenv("_ModuleTable001_", encoded[:half])
        monkeypatch.setenv("_ModuleTable002_", encoded[half:])
        monkeypatch.setenv("LOADEDMODULES", "StdEnv:GCCcore/11.3.0:OpenMPI/4.1.4-GCC-11.3.0")
        assert requested_modules() == ["StdEnv", "OpenMPI/4.1.4-GCC-11.3.0"]

    def test_all_loaded_modules_without_lmod(self, monkeypatch):
        monkeypatch.delenv("_ModuleTable_Sz_", raising=False)
        monkeypatch.setenv("LOADEDMODULES", "a/1:b/2")
        assert requested_modules() == ["a/1", "b/2"]

    def test_from_build_without_a_manifest_says_what_to_do(self, tmp_path):
        write_yaml(
            tmp_path / "qligfep.yaml",
            f"profiles: {{dev: {{extends: HABROK, q_bin: {tmp_path}, modules: from-build}}}}",
        )
        with pytest.raises(ConfigError, match="qligfep_config build-info"):
            load_config().get_profile("dev")


class TestSbatchHeader:
    def test_directives(self, tmp_path):
        write_yaml(
            tmp_path / "qligfep.yaml",
            """
            profiles:
              mine:
                resources: {ntasks: 12, time: "0-01:00:00", mem_per_cpu: 900}
                slurm: {account: acc, partition: short, extra_directives: ["--qos=fast", "#SBATCH --constraint=zen3"]}
            """,
        )
        header = load_config().get_profile("mine").sbatch_header("equilibrium", "w_a_b", 10)
        for line in [
            "#SBATCH --ntasks-per-node=12",
            "#SBATCH --mem-per-cpu=900",
            "#SBATCH -A acc",
            "#SBATCH -p short",
            "#SBATCH --export=ALL",
            "#SBATCH --time=0-01:00:00",
            "#SBATCH -J w_a_b",
            "#SBATCH --array=1-10",
            "#SBATCH --qos=fast",
            "#SBATCH --constraint=zen3",
        ]:
            assert line in header.splitlines()

    def test_memory_default_depends_on_the_method(self):
        profile = load_config().get_profile("HABROK")
        assert "--mem-per-cpu=512" in profile.sbatch_header("equilibrium", "j", 1)
        assert "--mem-per-cpu=1700" in profile.sbatch_header("nonequilibrium", "j", 1)

    def test_export_can_be_left_to_the_site(self, tmp_path):
        write_yaml(tmp_path / "qligfep.yaml", "profiles: {mine: {slurm: {export: null}}}")
        assert "--export" not in load_config().get_profile("mine").sbatch_header("equilibrium", "j", 1)

    def test_local_neq_cores_stay_within_the_profile(self, tmp_path):
        write_yaml(tmp_path / "qligfep.yaml", "profiles: {box: {extends: LOCALP, resources: {ntasks: 4}}}")
        assert load_config().get_profile("box").neq_cores(20) == 4

    def test_neq_cores_policy(self, tmp_path):
        write_yaml(tmp_path / "qligfep.yaml", "profiles: {waves: {extends: HABROK, resources: {neq_cores: 4}}}")
        config = load_config()
        auto = config.get_profile("HABROK")  # ntasks 8
        assert auto.neq_cores(6) == 8  # never below ntasks
        assert auto.neq_cores(20) == 20  # grows to hold every switch
        assert config.get_profile("waves").neq_cores(20) == 4  # explicit is respected


# ---------------------------------------------------------------------------------
# Protocols
# ---------------------------------------------------------------------------------


class TestProtocols:
    def test_default_protocol_matches_the_historical_defaults(self):
        protocol = load_config().get_protocol()
        assert (protocol.method, protocol.timestep, protocol.temperature) == ("equilibrium", "2fs", "298")
        assert (protocol.replicates, protocol.sphere_radius) == (10, 25)
        assert (protocol.windows, protocol.sampling, protocol.start) == (100, "sigmoidal", "0.5")
        assert (protocol.switch_pairs, protocol.schedule, protocol.steepness) == (5, "sigmoidal", 8.0)
        assert protocol.stages == {}

    def test_project_protocol_inherits_from_default(self, tmp_path):
        write_yaml(
            tmp_path / "qligfep.yaml",
            """
            protocols:
              fast-ligand-simulation:
                replicates: 5
                temperature: [298, 310]
                equilibrium: {windows: 50}
                stages:
                  eq5: {time: 100ps, inp: {MD: {bath_coupling: 5}}}
            """,
        )
        protocol = load_config().get_protocol("fast-ligand-simulation")
        assert protocol.replicates == 5 and protocol.windows == 50
        assert protocol.sampling == "sigmoidal"  # inherited
        assert protocol.temperature == "298,310"
        assert protocol.stages["eq5"].time_fs == 100_000
        assert protocol.stages["eq5"].inp == {"MD": {"bath_coupling": 5}}

    def test_options_of_the_other_method_are_an_error(self, tmp_path):
        write_yaml(
            tmp_path / "qligfep.yaml",
            "protocols: {p: {method: nonequilibrium, equilibrium: {windows: 50}}}",
        )
        with pytest.raises(ConfigError, match="would be ignored"):
            load_config().get_protocol("p")

    def test_stage_of_the_other_method_is_an_error(self, tmp_path):
        write_yaml(tmp_path / "qligfep.yaml", "protocols: {p: {stages: {switch: {steps: 10}}}}")
        with pytest.raises(ConfigError, match="not a stage of the 'equilibrium' method"):
            load_config().get_protocol("p")

    def test_unknown_stage_and_double_length_are_errors(self, tmp_path):
        write_yaml(tmp_path / "qligfep.yaml", "protocols: {p: {stages: {eq6: {steps: 10}}}}")
        with pytest.raises(ConfigError, match="unknown key 'eq6'"):
            load_config()
        write_yaml(tmp_path / "qligfep.yaml", "protocols: {p: {stages: {eq5: {steps: 10, time: 1ps}}}}")
        with pytest.raises(ConfigError, match="either 'time' or 'steps'"):
            load_config().get_protocol("p")

    @pytest.mark.parametrize("text, fs", [("500fs", 500), ("10ps", 10_000), ("2.5 ns", 2_500_000)])
    def test_durations(self, text, fs):
        assert parse_duration_fs(text, "x") == fs

    def test_duration_needs_a_unit(self):
        with pytest.raises(ConfigError, match="fs, ps or ns"):
            parse_duration_fs("100", "x")


class TestStageOverrides:
    def params(self):
        return get_production_config("2fs", shell_radius=25).params

    def test_time_is_converted_with_the_stage_stepsize(self):
        params = apply_stage_override(self.params(), StageOverride(time_fs=20_000), "production")
        assert params.steps == 10_000

    def test_time_that_is_not_a_whole_number_of_steps_is_refused(self):
        with pytest.raises(ConfigError, match="not a whole number"):
            apply_stage_override(self.params(), StageOverride(time_fs=3), "production")

    def test_modelled_keywords_replace_the_value_in_place(self):
        override = StageOverride(inp={"cut-offs": {"solute_solvent": 12}, "MD": {"shake_hydrogens": False}})
        params = apply_stage_override(self.params(), override, "production")
        assert params.cutoff_solute_solvent == 12 and params.constrain_hydrogens is False
        text = render_md_input(params, "0.5", "0.5", "md.dcd", "md.re")
        assert text.count("solute_solvent") == 1
        assert "solute_solvent            12" in text

    def test_other_keywords_are_passed_through_into_their_section(self):
        override = StageOverride(inp={"intervals": {"temperature": 10}, "MD": {"excluded_thermostat": True}})
        params = apply_stage_override(self.params(), override, "production")
        text = render_md_input(params, "0.5", "0.5", "md.dcd", "md.re")
        intervals = text.split("[intervals]")[1].split("\n[")[0]
        assert "temperature               10" in intervals
        md = text.split("[MD]")[1].split("\n[")[0]
        assert "excluded_thermostat       on" in md

    def test_generated_sections_cannot_be_overridden(self):
        with pytest.raises(ConfigError, match="generated by QligFEP"):
            apply_stage_override(self.params(), StageOverride(inp={"lambdas": {"x": 1}}), "production")
        with pytest.raises(ConfigError, match="set per replicate"):
            apply_stage_override(self.params(), StageOverride(inp={"MD": {"random_seed": 1}}), "eq1")

    def test_no_override_leaves_the_input_untouched(self):
        params = self.params()
        assert apply_stage_override(params, None, "production") is params


# ---------------------------------------------------------------------------------
# Command line precedence
# ---------------------------------------------------------------------------------


def cli_args(**given) -> argparse.Namespace:
    names = [
        "cluster", "protocol", "config", "q_bin", "sphereradius", "start", "temperature", "replicates",
        "sampling", "windows", "timestep", "neq", "neq_reps", "neq_steps", "neq_eq_steps",
        "neq_relax_steps", "neq_L", "neq_schedule",
    ]  # fmt: skip
    return argparse.Namespace(**{**dict.fromkeys(names), **given})


class TestCommandLinePrecedence:
    def test_protocol_fills_what_the_command_line_left_unset(self, tmp_path):
        write_yaml(
            tmp_path / "qligfep.yaml",
            """
            default_profile: HABROK
            default_protocol: fast
            protocols:
              fast: {replicates: 3, equilibrium: {windows: 20}}
            """,
        )
        args = cli_args(windows="40")
        profile, protocol = resolve_run_settings(args)
        assert (profile.name, protocol.name) == ("HABROK", "fast")
        assert args.replicates == "3"  # from the protocol
        assert args.windows == "40"  # the flag wins
        assert args.sphereradius == "25" and args.neq is False  # packaged defaults

    def test_protocol_selects_the_nonequilibrium_method(self, tmp_path):
        write_yaml(
            tmp_path / "qligfep.yaml",
            """
            protocols:
              neq:
                method: nonequilibrium
                nonequilibrium: {switch_pairs: 8}
                stages: {switch: {time: 40ps}}
            """,
        )
        args = cli_args(cluster="HABROK", protocol="neq")
        resolve_run_settings(args)
        assert args.neq is True and args.neq_reps == 8
        assert args.neq_steps == 20_000  # 40 ps at 2 fs, reported in the run record
        assert args.neq_eq_steps == 1000 and args.explicit_neq_steps == []

    def test_explicit_length_flag_beats_the_protocol_stage(self, tmp_path):
        write_yaml(
            tmp_path / "qligfep.yaml",
            "protocols: {neq: {method: nonequilibrium, stages: {switch: {time: 40ps, inp: {MD: {lrf: off}}}}}}",
        )
        args = cli_args(cluster="HABROK", protocol="neq", neq_steps=123)
        resolve_run_settings(args)
        assert args.neq_steps == 123 and args.explicit_neq_steps == ["neq_steps"]
        switch = args.stage_overrides["switch"]
        assert switch.time_fs is None and switch.inp == {"MD": {"lrf": False}}

    def test_account_and_partition_flags_override_the_profile(self):
        args = cli_args(cluster="DARDEL", account="naiss-1", partition=None)
        profile, _ = resolve_run_settings(args)
        header = profile.sbatch_header("equilibrium", "j", 1)
        assert "#SBATCH -A naiss-1" in header and "#SBATCH -p shared" in header
        profile, _ = resolve_run_settings(cli_args(cluster="DARDEL", partition="main"))
        header = profile.sbatch_header("equilibrium", "j", 1)
        assert "#SBATCH -p main" in header and "#SBATCH -A" not in header

    def test_configuration_errors_exit_with_a_message(self):
        with pytest.raises(SystemExit, match="Configuration error: Unknown profile"):
            resolve_run_settings(cli_args(cluster="NOPE"))


# ---------------------------------------------------------------------------------
# Generated job scripts
# ---------------------------------------------------------------------------------


def make_run(tmp_path: Path, cluster: str, **attributes) -> QligFEP:
    """A bare QligFEP holding only what the script writers read."""
    run = QligFEP.__new__(QligFEP)
    defaults = dict(
        cluster=cluster, system="water", lig1="a", lig2="b", replicates="3", temperature="298",
        seeds=[1, 2, 3], start="0.5", to_clean=None, production=False, minimize=False, neq=False, neq_reps=3,
    )  # fmt: skip
    for key, value in {**defaults, **attributes}.items():
        setattr(run, key, value)
    for name in ["eq1", "eq2", "md_0500_0500", "md_0750_0250", "md_0250_0750"]:
        (tmp_path / f"{name}.inp").write_text("")
    return run


WINDOW_FILES = [[], ["md_0750_0250.inp"], ["md_0250_0750.inp"]]


class TestRunScripts:
    def test_slurm_script(self, tmp_path):
        run = make_run(tmp_path, "SNELLIUS")
        run.write_runfile(str(tmp_path), WINDOW_FILES)
        script = (tmp_path / "runSNELLIUS.sh").read_text()
        assert "#SBATCH --array=1-3" in script and "#SBATCH -p rome" in script
        assert "module load OpenMPI/4.1.5-GCC-12.3.0" in script
        assert "time mpirun -n $ncores --map-by core --bind-to core $qdyn eq1.inp > eq1.log" in script
        assert script.index("md_0500_0500.inp") < script.index("md_0750_0250.inp") < script.index("md_0250_0750.inp")
        # no template placeholder may survive
        for placeholder in ["SBATCH_HEADER", "ENGINE_PATHS", "ENGINE_VARS", "ENVIRONMENT_SETUP", "DEFAULT_NTASKS"]:
            assert placeholder not in script

    def test_script_body_does_not_require_slurm(self, tmp_path):
        run = make_run(tmp_path, "SNELLIUS")
        run.write_runfile(str(tmp_path), WINDOW_FILES)
        script = (tmp_path / "runSNELLIUS.sh").read_text()
        assert "task_id=${QLIGFEP_TASK_ID:-${SLURM_ARRAY_TASK_ID:-}}" in script
        assert "ncores=${QLIGFEP_NTASKS:-${SLURM_NTASKS:-16}}" in script
        assert "$SLURM_ARRAY_TASK_ID" not in script and "$SLURM_NTASKS" not in script

    def test_serial_local_script_and_launcher(self, tmp_path):
        run = make_run(tmp_path, "LOCAL")
        run.write_runfile(str(tmp_path), WINDOW_FILES)
        run.write_submitfile(str(tmp_path))
        script = (tmp_path / "runLOCAL.sh").read_text()
        assert "\n#SBATCH" not in script and "time mpirun" not in script and "module load" not in script
        assert "time $qdyn eq1.inp > eq1.log" in script
        submit = (tmp_path / "FEP_submit.sh").read_text()
        assert "sbatch" not in submit
        assert "runfile=inputfiles/runLOCAL.sh" in submit and "total_tasks=3" in submit
        assert 'QLIGFEP_TASK_ID=$task bash "$runfile"' in submit

    def test_slurm_launcher_submits_an_array(self, tmp_path):
        run = make_run(tmp_path, "HABROK")
        run.write_submitfile(str(tmp_path))
        submit = (tmp_path / "FEP_submit.sh").read_text()
        assert "sbatch $submitfile" in submit and "inputfiles/runHABROK.sh" in submit

    def test_one_task_per_temperature_and_replicate(self, tmp_path):
        run = make_run(tmp_path, "HABROK", temperature="298,310")
        run.write_runfile(str(tmp_path), WINDOW_FILES)
        script = (tmp_path / "runHABROK.sh").read_text()
        assert "temperatures=(298 310)" in script and "#SBATCH --array=1-6" in script

    def test_local_mpi_launch_has_no_core_binding(self, tmp_path):
        run = make_run(tmp_path, "LOCALP")
        run.write_runfile(str(tmp_path), WINDOW_FILES)
        script = (tmp_path / "runLOCALP.sh").read_text()
        assert "time mpirun -n $ncores $qdyn eq1.inp" in script and "--bind-to" not in script

    def test_custom_mpi_command(self, tmp_path):
        write_yaml(tmp_path / "qligfep.yaml", 'profiles: {srun: {extends: HABROK, mpi_command: "srun -n $ncores"}}')
        run = make_run(tmp_path, "srun")
        run.write_runfile(str(tmp_path), WINDOW_FILES)
        assert "time srun -n $ncores $qdyn eq1.inp" in (tmp_path / "runsrun.sh").read_text()

    def test_neq_runs_locally_with_an_mpi_profile(self, tmp_path):
        write_yaml(tmp_path / "qligfep.yaml", "profiles: {box: {extends: LOCALP, resources: {ntasks: 4, neq_cores: 2}}}")
        run = make_run(tmp_path, "box", neq=True)
        run.write_neq_runfile(str(tmp_path), [])
        run.write_submitfile(str(tmp_path))
        script = (tmp_path / "runbox.sh").read_text()
        assert 'mpi_map=""' in script and "ncores=${QLIGFEP_NTASKS:-${SLURM_NTASKS:-2}}" in script
        assert "export QLIGFEP_NTASKS=${QLIGFEP_NTASKS:-2}" in (tmp_path / "FEP_submit.sh").read_text()

    def test_neq_is_refused_for_a_serial_profile(self):
        with pytest.raises(ValueError, match="engine is the serial"):
            QligFEP("a", "b", "AMBER14sb", "water", cluster="LOCAL", neq=True)


# ---------------------------------------------------------------------------------
# qligfep_config
# ---------------------------------------------------------------------------------


class TestConfigCommand:
    def test_init_writes_a_commented_file_that_loads(self, tmp_path, capsys):
        assert config_cli.main(["init", "--execution", "slurm-array", "--name", "site", "--account", "acc"]) == 0
        path = user_file(tmp_path)
        text = path.read_text()
        assert "# slurm-array: each replicate is a task of a Slurm job array" in text
        profile = load_config().get_profile()
        assert profile.name == "site" and profile.slurm.account == "acc" and not profile.is_local

    def test_init_does_not_guess_the_slurm_account(self, tmp_path, capsys, monkeypatch):
        detect = config_cli.detect_environment
        monkeypatch.setattr(config_cli, "detect_environment", lambda q_bin: {**detect(q_bin), "accounts": ["users"]})
        assert config_cli.main(["init", "--execution", "slurm-array"]) == 0
        profile = load_config().get_profile()
        assert profile.slurm.account is None
        assert "#SBATCH -A" not in profile.sbatch_header("equilibrium", "j", 1)
        assert "default account" in capsys.readouterr().err

    def test_init_never_overwrites(self, tmp_path, capsys):
        path = write_yaml(user_file(tmp_path), "# mine\ndefault_profile: LOCAL\n")
        assert config_cli.main(["init"]) == 1
        assert path.read_text() == "# mine\ndefault_profile: LOCAL\n"
        assert "left untouched" in capsys.readouterr().err

    def test_generated_local_file_is_valid_yaml(self):
        env = {
            "hostname": "h", "has_slurm": False, "cpu_count": 8, "mpirun": None, "accounts": [],
            "partitions": [], "binaries": {"qdynp": True},
        }  # fmt: skip
        data = yaml.safe_load(render_user_config(env, "local", "local", "auto", "0-24:00:00", [], None, None))
        # no mpirun on the machine: fall back to the serial engine on one core
        assert data["profiles"]["local"]["engine"] == "qdyn"
        assert data["profiles"]["local"]["resources"]["ntasks"] == 1

    def test_show_lists_and_resolves(self, tmp_path, capsys):
        write_yaml(tmp_path / "qligfep.yaml", "protocols: {fast: {description: quick look, replicates: 2}}")
        assert config_cli.main(["show"]) == 0
        listing = capsys.readouterr().out
        assert "SNELLIUS" in listing and "fast" in listing and "quick look" in listing
        assert config_cli.main(["show", "--protocol", "fast"]) == 0
        assert yaml.safe_load(capsys.readouterr().out)["protocol"]["replicates"] == 2

    def test_check_reports_an_invalid_file(self, tmp_path, capsys):
        write_yaml(tmp_path / "qligfep.yaml", "profiles: {x: {engine: qdynx}}")
        assert config_cli.main(["check"]) == 1
        assert "qdynx" in capsys.readouterr().out
