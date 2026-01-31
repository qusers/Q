"""Unit tests for the FEP input file template system."""

import pytest
from conftest import GOLDEN_FEP_SETUP_DIR, parse_inp_file

from QligFEP.templates import (
    MDParameters,
    get_equilibration_configs,
    get_production_config,
    render_md_input,
)
from QligFEP.templates.sections import (
    format_distance_restraints,
    format_sequence_restraint,
    format_water_restraint,
)

GOLDEN_DIR = GOLDEN_FEP_SETUP_DIR / "tyk2_ejm_31_ejm_45"


class TestMDParameters:
    """Tests for MDParameters dataclass."""

    def test_default_values(self):
        """Verify default parameter values."""
        params = MDParameters(
            steps=5000,
            stepsize=2.0,
            temperature=298,
            bath_coupling=10.0,
        )

        assert params.shake_solvent is True
        assert params.shake_hydrogens is True
        assert params.shake_solute is False
        assert params.lrf is True
        assert params.cutoff_q_atom == 99
        assert params.shell_force == 10.0
        assert params.polarisation is True
        assert params.topology == "dualtop.top"
        assert params.fep_file == "FEP_VAR"

    def test_temperature_placeholder(self):
        """Verify temperature can hold T_VAR placeholder string."""
        params = MDParameters(
            steps=5000,
            stepsize=2.0,
            temperature="T_VAR",
            bath_coupling=10.0,
        )

        assert params.temperature == "T_VAR"


class TestRenderMdInput:
    """Tests for render_md_input function."""

    def test_render_basic_file(self):
        """Verify basic MD input file rendering."""
        params = MDParameters(
            steps=5000,
            stepsize=2.0,
            temperature=298,
            bath_coupling=10.0,
            shell_radius=25,
        )

        content = render_md_input(
            params=params,
            lambda1="0.500",
            lambda2="0.500",
            trajectory_file="test.dcd",
            final_file="test.re",
        )

        assert "[MD]" in content
        assert "steps                     5000" in content
        assert "stepsize                  2.0" in content
        assert "temperature               298" in content
        assert "[cut-offs]" in content
        assert "[sphere]" in content
        assert "shell_radius              25" in content
        assert "[lambdas]" in content
        assert "0.500 0.500" in content
        assert "[distance_restraints]" in content

    def test_render_eq1_with_random_seed(self):
        """Verify eq1 includes random_seed and initial_temperature."""
        params = MDParameters(
            steps=10000,
            stepsize=0.1,
            temperature=1,
            bath_coupling=0.2,
            shell_radius=25,
        )

        content = render_md_input(
            params=params,
            lambda1="0.500",
            lambda2="0.500",
            trajectory_file="eq1.dcd",
            final_file="eq1.re",
            is_eq1=True,
        )

        assert "random_seed               SEED_VAR" in content
        assert "initial_temperature       1" in content

    def test_render_without_random_seed(self):
        """Verify non-eq1 files don't have random_seed."""
        params = MDParameters(
            steps=5000,
            stepsize=2.0,
            temperature=50,
            bath_coupling=2.0,
            shell_radius=25,
        )

        content = render_md_input(
            params=params,
            lambda1="0.500",
            lambda2="0.500",
            trajectory_file="eq2.dcd",
            final_file="eq2.re",
            restart_file="eq1.re",
        )

        assert "random_seed" not in content
        assert "initial_temperature" not in content
        assert "restart                   eq1.re" in content

    def test_render_with_energy_interval(self):
        """Verify production MD includes energy section."""
        params = MDParameters(
            steps=5000,
            stepsize=2.0,
            temperature="T_VAR",
            bath_coupling=10.0,
            shell_radius=25,
            interval_energy=10,
        )

        content = render_md_input(
            params=params,
            lambda1="0.500",
            lambda2="0.500",
            trajectory_file="md_0500_0500.dcd",
            final_file="md_0500_0500.re",
            restart_file="eq5.re",
            energy_file="md_0500_0500.en",
        )

        # Check intervals section has energy
        assert "energy                    10" in content
        # Check files section has energy file
        assert "energy                    md_0500_0500.en" in content

    def test_render_with_restraints(self):
        """Verify distance and sequence restraints are included."""
        params = MDParameters(
            steps=5000,
            stepsize=2.0,
            temperature=298,
            bath_coupling=10.0,
            shell_radius=25,
        )

        distance_restraints = "4681 4713 0.0 0.1 0.5 0\n4682 4714 0.0 0.1 0.5 0"
        sequence_restraints = "4681 4751 10.0 0  0"

        content = render_md_input(
            params=params,
            lambda1="0.500",
            lambda2="0.500",
            trajectory_file="test.dcd",
            final_file="test.re",
            distance_restraints=distance_restraints,
            sequence_restraints=sequence_restraints,
        )

        assert "4681 4713 0.0 0.1 0.5 0" in content
        assert "4682 4714 0.0 0.1 0.5 0" in content
        assert "4681 4751 10.0 0  0" in content


class TestEquilibrationConfigs:
    """Tests for equilibration configuration generation."""

    def test_configs_count(self):
        """Verify 5 equilibration configs are generated."""
        configs = get_equilibration_configs("2fs", shell_radius=25)
        assert len(configs) == 5
        assert [c.name for c in configs] == ["eq1", "eq2", "eq3", "eq4", "eq5"]

    def test_eq1_config(self):
        """Verify eq1 has fixed timestep and temperature."""
        configs = get_equilibration_configs("2fs", shell_radius=25)
        eq1 = configs[0]

        assert eq1.params.steps == 10000
        assert eq1.params.stepsize == 0.1
        assert eq1.params.temperature == 1
        assert eq1.params.bath_coupling == 0.2
        assert eq1.params.shake_hydrogens is True
        assert eq1.sequence_restraint_force == 10.0
        assert eq1.distance_restraint_force == 1.5

    def test_eq5_uses_water_restraint(self):
        """Verify eq5 is configured for WATER_RESTRAINT."""
        configs = get_equilibration_configs("2fs", shell_radius=25)
        eq5 = configs[4]

        assert eq5.params.temperature == "T_VAR"
        assert eq5.use_water_restraint is True

    def test_timestep_1fs_vs_2fs(self):
        """Verify timestep affects step counts and settings."""
        configs_2fs = get_equilibration_configs("2fs", shell_radius=25)
        configs_1fs = get_equilibration_configs("1fs", shell_radius=25)

        # eq2 should have different step counts
        eq2_2fs = configs_2fs[1]
        eq2_1fs = configs_1fs[1]

        assert eq2_2fs.params.steps == 5000
        assert eq2_1fs.params.steps == 10000

        assert eq2_2fs.params.stepsize == 2.0
        assert eq2_1fs.params.stepsize == 1.0

        # shake_hydrogens differs
        assert eq2_2fs.params.shake_hydrogens is True
        assert eq2_1fs.params.shake_hydrogens is False

    def test_temperature_progression(self):
        """Verify temperatures increase through equilibration."""
        configs = get_equilibration_configs("2fs", shell_radius=25)

        temps = [c.params.temperature for c in configs]
        assert temps[0] == 1  # eq1
        assert temps[1] == 50  # eq2
        assert temps[2] == 150  # eq3
        assert temps[3] == 275  # eq4
        assert temps[4] == "T_VAR"  # eq5

    def test_restraint_force_progression(self):
        """Verify sequence restraint forces decrease through eq1-eq4."""
        configs = get_equilibration_configs("2fs", shell_radius=25)

        forces = [c.sequence_restraint_force for c in configs[:4]]
        assert forces == [10.0, 10.0, 5.0, 2.5]


class TestProductionConfig:
    """Tests for production MD configuration."""

    def test_production_config_2fs(self):
        """Verify production config for 2fs timestep."""
        config = get_production_config("2fs", shell_radius=25)

        assert config.params.steps == 5000
        assert config.params.stepsize == 2.0
        assert config.params.temperature == "T_VAR"
        assert config.params.interval_energy == 10
        assert config.distance_restraint_force == 0.5

    def test_production_config_1fs(self):
        """Verify production config for 1fs timestep."""
        config = get_production_config("1fs", shell_radius=25)

        assert config.params.steps == 10000
        assert config.params.stepsize == 1.0

    def test_custom_dr_force(self):
        """Verify custom distance restraint force."""
        config = get_production_config("2fs", shell_radius=25, distance_restraint_force=1.0)
        assert config.distance_restraint_force == 1.0

    def test_extended_sampling_steps(self):
        """Verify custom step count for extended sampling."""
        config = get_production_config("2fs", shell_radius=25, steps=10000)
        assert config.params.steps == 10000

    def test_custom_intervals(self):
        """Verify custom interval settings."""
        config = get_production_config(
            "2fs",
            shell_radius=25,
            interval_output=50,
            interval_energy=5,
            interval_trajectory=200,
        )
        assert config.params.interval_output == 50
        assert config.params.interval_energy == 5
        assert config.params.interval_trajectory == 200

    def test_custom_bath_coupling(self):
        """Verify custom bath coupling."""
        config = get_production_config("2fs", shell_radius=25, bath_coupling=5.0)
        assert config.params.bath_coupling == 5.0


class TestSectionFormatters:
    """Tests for restraint section formatters."""

    def test_format_distance_restraints(self):
        """Verify distance restraints formatting."""
        pairs = [(4681, 4713), (4682, 4714)]
        result = format_distance_restraints(pairs, force=0.5)

        assert "4681 4713 0.0 0.1 0.5 0" in result
        assert "4682 4714 0.0 0.1 0.5 0" in result

    def test_format_sequence_restraint(self):
        """Verify sequence restraint formatting."""
        result = format_sequence_restraint(4681, 4751, force=10.0)
        # Format: atom_start (6 chars left), space, atom_end (6 chars left), space, force (4 chars right)
        assert result == "4681   4751   10.0 0  0"

    def test_format_water_restraint(self):
        """Verify water restraint formatting."""
        result = format_water_restraint(4681, 4751, force=1.0)
        # Format: atom_start (7 chars left-aligned), atom_end (7 chars), force, 0, 1
        assert "4681" in result
        assert "4751" in result
        assert "1.0 0 1" in result


class TestGoldenFileComparison:
    """Tests comparing rendered output against golden files."""

    @pytest.fixture
    def golden_eq1_parsed(self):
        """Load and parse golden eq1.inp."""
        if not GOLDEN_DIR.exists():
            pytest.skip("Golden files not available")
        content = (GOLDEN_DIR / "eq1.inp").read_text()
        return parse_inp_file(content)

    @pytest.fixture
    def golden_md_parsed(self):
        """Load and parse golden md_0500_0500.inp."""
        if not GOLDEN_DIR.exists():
            pytest.skip("Golden files not available")
        content = (GOLDEN_DIR / "md_0500_0500.inp").read_text()
        return parse_inp_file(content)

    def test_eq1_md_section_matches_golden(self, golden_eq1_parsed):
        """Verify eq1 [MD] section matches golden file."""
        configs = get_equilibration_configs("2fs", shell_radius=25)
        eq1_config = configs[0]

        content = render_md_input(
            params=eq1_config.params,
            lambda1="0.500",
            lambda2="0.500",
            trajectory_file="eq1.dcd",
            final_file="eq1.re",
            is_eq1=True,
        )
        rendered_parsed = parse_inp_file(content)

        golden_md = golden_eq1_parsed["MD"]
        rendered_md = rendered_parsed["MD"]

        # Compare key values (ignoring placeholders that vary)
        assert rendered_md["steps"] == golden_md["steps"]
        assert rendered_md["stepsize"] == golden_md["stepsize"]
        assert rendered_md["temperature"] == golden_md["temperature"]
        assert rendered_md["bath_coupling"] == golden_md["bath_coupling"]
        assert rendered_md["random_seed"] == golden_md["random_seed"]

    def test_md_section_matches_golden(self, golden_md_parsed):
        """Verify production MD [MD] section matches golden file."""
        config = get_production_config("2fs", shell_radius=25)

        content = render_md_input(
            params=config.params,
            lambda1="0.500",
            lambda2="0.500",
            trajectory_file="md_0500_0500.dcd",
            final_file="md_0500_0500.re",
            restart_file="eq5.re",
            energy_file="md_0500_0500.en",
        )
        rendered_parsed = parse_inp_file(content)

        golden_md = golden_md_parsed["MD"]
        rendered_md = rendered_parsed["MD"]

        assert rendered_md["steps"] == golden_md["steps"]
        assert rendered_md["stepsize"] == golden_md["stepsize"]
        assert rendered_md["temperature"] == golden_md["temperature"]
        assert rendered_md["shake_hydrogens"] == golden_md["shake_hydrogens"]


class TestQfepTemplate:
    """Tests for qfep.inp template functions."""

    def test_format_energy_files_basic(self):
        """Verify energy filenames are generated from lambda values."""
        from QligFEP.templates.qfep import format_energy_files

        lambdas = ["1.000", "0.500", "0.000"]
        result = format_energy_files(lambdas)

        assert len(result) == 3
        assert result[0] == "md_1000_0000.en"
        assert result[1] == "md_0500_0500.en"
        assert result[2] == "md_0000_1000.en"

    def test_format_energy_files_single(self):
        """Verify single lambda pair."""
        from QligFEP.templates.qfep import format_energy_files

        lambdas = ["0.500"]
        result = format_energy_files(lambdas)

        assert len(result) == 1
        assert result[0] == "md_0500_0500.en"

    def test_render_qfep_input_structure(self):
        """Verify qfep.inp has correct structure."""
        from QligFEP.templates.qfep import render_qfep_input

        content = render_qfep_input(
            total_lambdas=5,
            temperature=298.0,
            windows=5,
            energy_files=["md_1000_0000.en", "md_0750_0250.en", "md_0500_0500.en"],
        )

        lines = content.strip().split("\n")
        # Line 0: total_lambdas
        assert lines[0] == "5"
        # Line 1: 2  0
        assert lines[1] == "2  0"
        # Line 2: kT  windows
        assert "0.592" in lines[2]  # kT at 298K
        assert "5" in lines[2]
        # Lines 3,4,5: windows repeated
        assert lines[3] == "5"
        assert lines[4] == "5"
        assert lines[5] == "5"
        # Lines 6,7: 0, 0
        assert lines[6] == "0"
        assert lines[7] == "0"
        # Line 8: 1 0
        assert lines[8] == "1 0"
        # Energy files follow
        assert "md_1000_0000.en" in content
        assert "md_0750_0250.en" in content
        assert "md_0500_0500.en" in content

    def test_render_qfep_kT_calculation(self):
        """Verify kT is calculated correctly for different temperatures."""
        from QligFEP.templates.qfep import render_qfep_input

        # At 300K, kT = 0.00198720 * 300 = 0.5962
        content = render_qfep_input(
            total_lambdas=3,
            temperature=300.0,
            windows=3,
            energy_files=["test.en"],
        )

        assert "0.596" in content  # kT at 300K

        # At 310K, kT = 0.00198720 * 310 = 0.6160
        content = render_qfep_input(
            total_lambdas=3,
            temperature=310.0,
            windows=3,
            energy_files=["test.en"],
        )

        assert "0.616" in content  # kT at 310K


class TestQprepFEPTemplate:
    """Tests for qprep.inp FEP template functions."""

    def test_render_water_system(self):
        """Verify qprep.inp for water system."""
        from QligFEP.templates.qprep import QprepFEPParameters, render_qprep_fep_input

        params = QprepFEPParameters(
            ff_lib="qamber14.lib",
            lig1_lib="lig1.lib",
            lig2_lib="lig2_renumber.lib",
            ligand_prm="AMBER14sb_lig1_lig2_merged.prm",
            ligand_pdb="lig1_lig2.pdb",
            center="10.0 20.0 30.0",
            sphere_radius="25",
            solute_density="0.05794",
            solvent="1 HOH",
        )

        content = render_qprep_fep_input(params)

        assert "rl qamber14.lib" in content
        assert "rl lig1.lib" in content
        assert "rl lig2_renumber.lib" in content
        assert "rprm AMBER14sb_lig1_lig2_merged.prm" in content
        assert "rp lig1_lig2.pdb" in content
        assert "boundary 1 10.0 20.0 30.0 25" in content
        assert "solvate 10.0 20.0 30.0 25 1 HOH" in content
        assert "set solute_density 0.05794" in content
        assert "maketop MKC_p" in content
        assert "writetop dualtop.top" in content

    def test_render_protein_system(self):
        """Verify qprep.inp for protein system."""
        from QligFEP.templates.qprep import QprepFEPParameters, render_qprep_fep_input

        params = QprepFEPParameters(
            ff_lib="qamber14.lib",
            lig1_lib="lig1.lib",
            lig2_lib="lig2_renumber.lib",
            ligand_prm="AMBER14sb_lig1_lig2_merged.prm",
            ligand_pdb="lig1_lig2.pdb",
            center="-5.123 10.456 0.789",
            sphere_radius="25",
            solute_density="0.06123",
            solvent="4 water.pdb",
        )

        content = render_qprep_fep_input(params)

        assert "solvate -5.123 10.456 0.789 25 4 water.pdb" in content
        assert "set solute_density 0.06123" in content

    def test_render_vacuum_system(self):
        """Verify qprep.inp for vacuum system (no solvation)."""
        from QligFEP.templates.qprep import QprepFEPParameters, render_qprep_fep_input

        params = QprepFEPParameters(
            ff_lib="qamber14.lib",
            lig1_lib="lig1.lib",
            lig2_lib="lig2_renumber.lib",
            ligand_prm="merged.prm",
            ligand_pdb="combined.pdb",
            center="0.0 0.0 0.0",
            sphere_radius="25",
            solute_density="0.05794",
            solvent="",
            solvate=False,
        )

        content = render_qprep_fep_input(params)

        assert "!solvate" in content or "solvate" not in content.replace("!solvate", "")

    def test_render_with_cysbonds(self):
        """Verify qprep.inp includes cysbond lines."""
        from QligFEP.templates.qprep import QprepFEPParameters, render_qprep_fep_input

        params = QprepFEPParameters(
            ff_lib="qamber14.lib",
            lig1_lib="lig1.lib",
            lig2_lib="lig2_renumber.lib",
            ligand_prm="merged.prm",
            ligand_pdb="combined.pdb",
            center="0.0 0.0 0.0",
            sphere_radius="25",
            solute_density="0.05794",
            solvent="1 HOH",
            cysbonds="addbond 22:SG 45:SG y\n",
        )

        content = render_qprep_fep_input(params)

        assert "addbond 22:SG 45:SG y" in content
        # cysbonds should appear before maketop
        assert content.index("addbond 22:SG 45:SG y") < content.index("maketop")


class TestQprepProteinTemplate:
    """Tests for qprep.inp protein template functions."""

    def test_render_basic(self):
        """Verify basic protein qprep.inp rendering."""
        from QligFEP.templates.qprep import QprepProteinParameters, render_qprep_protein_input

        params = QprepProteinParameters(
            ff_lib_path="/path/to/qamber14.lib",
            ff_prm_path="/path/to/qamber14.prm",
            pdb_file_path="protein.pdb",
            cog="10.0 20.0 30.0",
            sphere_radius="25",
            solvent_pack="2.3",
        )

        content = render_qprep_protein_input(params)

        assert "rl /path/to/qamber14.lib" in content
        assert "rprm /path/to/qamber14.prm" in content
        assert "rp protein.pdb" in content
        assert "boundary 1 10.0 20.0 30.0 25" in content
        assert "solvate 10.0 20.0 30.0 25 1 HOH" in content
        assert "set solvent_pack 2.3" in content
        assert "maketop MKC_p" in content
        assert "writetop dualtop.top" in content
        assert "q" in content

    def test_render_with_cysbonds(self):
        """Verify protein qprep.inp includes cysbond lines."""
        from QligFEP.templates.qprep import QprepProteinParameters, render_qprep_protein_input

        params = QprepProteinParameters(
            ff_lib_path="qamber14.lib",
            ff_prm_path="qamber14.prm",
            pdb_file_path="protein.pdb",
            cog="0.0 0.0 0.0",
            sphere_radius="25",
            solvent_pack="2.3",
            cysbonds="addbond 10:SG 25:SG y\naddbond 50:SG 80:SG y\n",
        )

        content = render_qprep_protein_input(params)

        assert "addbond 10:SG 25:SG y" in content
        assert "addbond 50:SG 80:SG y" in content
        # cysbonds should appear before maketop
        assert content.index("addbond 10:SG 25:SG y") < content.index("maketop")
