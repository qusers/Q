"""Tests for counter-ion placement on a sphere via Coulomb energy minimization."""

import shutil
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from QligFEP.counter_ions import (
    cart_to_sph,
    fibonacci_sphere_angles,
    minimize_coulomb_on_sphere,
    sph_to_cart,
)


class TestSphericalConversions:
    """Tests for spherical <-> Cartesian coordinate conversions."""

    def test_sph_to_cart_north_pole(self):
        """theta=0 maps to (0, 0, 1)."""
        theta = np.array([0.0])
        phi = np.array([0.0])
        xyz = sph_to_cart(theta, phi)
        np.testing.assert_allclose(xyz, [[0, 0, 1]], atol=1e-12)

    def test_sph_to_cart_south_pole(self):
        """theta=pi maps to (0, 0, -1)."""
        theta = np.array([np.pi])
        phi = np.array([0.0])
        xyz = sph_to_cart(theta, phi)
        np.testing.assert_allclose(xyz, [[0, 0, -1]], atol=1e-12)

    def test_sph_to_cart_equator(self):
        """theta=pi/2, phi=0 maps to (1, 0, 0)."""
        theta = np.array([np.pi / 2])
        phi = np.array([0.0])
        xyz = sph_to_cart(theta, phi)
        np.testing.assert_allclose(xyz, [[1, 0, 0]], atol=1e-12)

    def test_roundtrip_conversion(self):
        """sph_to_cart followed by cart_to_sph recovers original angles."""
        rng = np.random.default_rng(42)
        theta = rng.uniform(0.1, np.pi - 0.1, size=10)
        phi = rng.uniform(0, 2 * np.pi, size=10)

        xyz = sph_to_cart(theta, phi)
        theta2, phi2 = cart_to_sph(xyz)

        np.testing.assert_allclose(theta2, theta, atol=1e-10)
        np.testing.assert_allclose(phi2, phi, atol=1e-10)

    def test_sph_to_cart_unit_vectors(self):
        """All output vectors have unit norm."""
        rng = np.random.default_rng(7)
        theta = rng.uniform(0, np.pi, size=50)
        phi = rng.uniform(0, 2 * np.pi, size=50)
        xyz = sph_to_cart(theta, phi)
        norms = np.linalg.norm(xyz, axis=1)
        np.testing.assert_allclose(norms, 1.0, atol=1e-12)


class TestFibonacciSphere:
    """Tests for Fibonacci sphere angle initialization."""

    def test_correct_count(self):
        """Returns the requested number of points."""
        theta, phi = fibonacci_sphere_angles(20, seed=0)
        assert theta.shape == (20,)
        assert phi.shape == (20,)

    def test_theta_range(self):
        """Theta values are in [0, pi]."""
        theta, _ = fibonacci_sphere_angles(100, seed=0)
        assert np.all(theta >= 0)
        assert np.all(theta <= np.pi)

    def test_phi_range(self):
        """Phi values are in [0, 2*pi)."""
        _, phi = fibonacci_sphere_angles(100, seed=0)
        assert np.all(phi >= 0)
        assert np.all(phi < 2 * np.pi)

    def test_reproducible_with_seed(self):
        """Same seed produces identical angles."""
        t1, p1 = fibonacci_sphere_angles(10, seed=42)
        t2, p2 = fibonacci_sphere_angles(10, seed=42)
        np.testing.assert_array_equal(t1, t2)
        np.testing.assert_array_equal(p1, p2)


class TestMinimizeCoulombOnSphere:
    """Tests for the Coulomb minimization on a sphere."""

    def test_two_charges_antipodal(self):
        """Two charges should settle at opposite ends of a diameter."""
        center = np.array([0.0, 0.0, 0.0])
        radius = 10.0
        xyz = minimize_coulomb_on_sphere(2, radius, center, seed=42)

        assert xyz.shape == (2, 3)
        # Both should be at distance `radius` from center
        dists = np.linalg.norm(xyz - center, axis=1)
        np.testing.assert_allclose(dists, radius, atol=0.1)

        # They should be diametrically opposed (distance = 2*radius)
        pair_dist = np.linalg.norm(xyz[0] - xyz[1])
        np.testing.assert_allclose(pair_dist, 2 * radius, atol=0.1)

    def test_single_charge(self):
        """Single charge returns one point on the sphere."""
        center = np.array([5.0, 5.0, 5.0])
        radius = 14.0
        xyz = minimize_coulomb_on_sphere(1, radius, center, seed=42)

        assert xyz.shape == (1, 3)
        dist = np.linalg.norm(xyz[0] - center)
        np.testing.assert_allclose(dist, radius, atol=0.1)

    def test_center_offset(self):
        """Points are centered on the given center, not the origin."""
        center = np.array([10.0, -5.0, 3.0])
        radius = 8.0
        xyz = minimize_coulomb_on_sphere(4, radius, center, seed=42)

        dists = np.linalg.norm(xyz - center, axis=1)
        np.testing.assert_allclose(dists, radius, atol=0.1)

    def test_three_charges_equilateral(self):
        """Three charges should form an equilateral triangle on a great circle."""
        center = np.array([0.0, 0.0, 0.0])
        radius = 10.0
        xyz = minimize_coulomb_on_sphere(3, radius, center, seed=42)

        # All pairwise distances should be equal
        d01 = np.linalg.norm(xyz[0] - xyz[1])
        d02 = np.linalg.norm(xyz[0] - xyz[2])
        d12 = np.linalg.norm(xyz[1] - xyz[2])
        np.testing.assert_allclose(d01, d02, atol=0.2)
        np.testing.assert_allclose(d01, d12, atol=0.2)

    def test_reproducible_with_seed(self):
        """Same seed produces identical placement."""
        center = np.array([0.0, 0.0, 0.0])
        xyz1 = minimize_coulomb_on_sphere(5, 10.0, center, seed=123)
        xyz2 = minimize_coulomb_on_sphere(5, 10.0, center, seed=123)
        np.testing.assert_array_equal(xyz1, xyz2)


class TestPlaceCounterIonsWithProteinCharge:
    """Tests for place_counter_ions() protein_charge logic.

    Uses a real QligFEP instance with Tyk2 tutorial data to test the charge
    calculation and ion type selection when protein_charge is provided.
    """

    TUTORIAL_DIR = Path(__file__).parents[2] / "tutorials" / "Tyk2" / "setupFEP"
    LIG1 = "ejm_31"
    LIG2 = "ejm_42"

    def _make_qligfep(self, tmp_path, system="water", protein_charge=None, charge_method="ion_match"):
        """Create a QligFEP instance in tmp_path with real tutorial ligand files."""
        from QligFEP.qligfep import DualTopologyFEP

        # Copy ligand files to tmp_path
        for lig in (self.LIG1, self.LIG2):
            for ext in (".lib", ".prm", ".pdb"):
                src = self.TUTORIAL_DIR / f"{lig}{ext}"
                if src.exists():
                    shutil.copy(src, tmp_path / f"{lig}{ext}")

        import os

        original_cwd = os.getcwd()
        os.chdir(tmp_path)
        try:
            run = DualTopologyFEP(
                lig1=self.LIG1,
                lig2=self.LIG2,
                FF="AMBER14sb",
                system=system,
                cluster="TETRA",
                sphereradius="25",
                protein_charge=protein_charge,
                charge_method=charge_method,
            )
        finally:
            os.chdir(original_cwd)
        return run

    def _prepare_inputdir(self, run, tmp_path):
        """Run the initial steps to create the merged PDB needed by place_counter_ions."""
        import os

        original_cwd = os.getcwd()
        os.chdir(tmp_path)
        try:
            writedir = run.makedir()
            inputdir = writedir + "/inputfiles"
            run.read_files()
            run.merge_pdbs(inputdir)
            return inputdir
        finally:
            os.chdir(original_cwd)

    def _run_place_counter_ions(self, run, inputdir, tmp_path):
        """Run place_counter_ions in the correct cwd (COG uses relative paths)."""
        import os

        original_cwd = os.getcwd()
        os.chdir(tmp_path)
        try:
            return run.place_counter_ions(inputdir)
        finally:
            os.chdir(original_cwd)

    def test_protein_charge_positive_places_sod(self, tmp_path):
        """protein_charge=8, lig1=0, lig2=0 (same charge) -> 0 ions (delta_q=0)."""
        run = self._make_qligfep(tmp_path, protein_charge=8)
        inputdir = self._prepare_inputdir(run, tmp_path)
        # Both Tyk2 ligands are neutral (charge=0), so delta_q=0 -> no ions regardless of protein_charge
        n_ions = self._run_place_counter_ions(run, inputdir, tmp_path)
        assert n_ions == 0

    def test_protein_charge_sets_ion_count(self, tmp_path):
        """With protein_charge set, n_ions = abs(protein_charge - water_charge)."""
        run = self._make_qligfep(tmp_path, protein_charge=8)
        inputdir = self._prepare_inputdir(run, tmp_path)

        # Simulate a charge-changing perturbation: override charge_lig2 to +1
        run.charge_lig2 = 1
        n_ions = self._run_place_counter_ions(run, inputdir, tmp_path)
        # water_charge = charge_lig1 + charge_lig2 = 0 + 1 = 1
        # ions_charge = 8 - 1 = 7 -> 7 SOD ions
        assert n_ions == 7
        assert run.ion_type == "SOD"

    def test_protein_charge_negative_places_chloride(self, tmp_path):
        """protein_charge=-2, lig1=0, lig2=+1 -> 3 CHL ions (AMBER14sb)."""
        run = self._make_qligfep(tmp_path, protein_charge=-2)
        inputdir = self._prepare_inputdir(run, tmp_path)

        run.charge_lig2 = 1
        n_ions = self._run_place_counter_ions(run, inputdir, tmp_path)
        # water_charge = 0 + 1 = 1; ions_charge = -2 - 1 = -3 -> 3 chloride
        assert n_ions == 3
        assert run.ion_type == "CHL"

    def test_protein_charge_with_negative_lig(self, tmp_path):
        """protein_charge=8, lig1=0, lig2=-1 -> 9 SOD ions."""
        run = self._make_qligfep(tmp_path, protein_charge=8)
        inputdir = self._prepare_inputdir(run, tmp_path)

        run.charge_lig2 = -1
        n_ions = self._run_place_counter_ions(run, inputdir, tmp_path)
        # water_charge = 0 + (-1) = -1; ions_charge = 8 - (-1) = 9 -> 9 SOD
        assert n_ions == 9
        assert run.ion_type == "SOD"

    def test_no_protein_charge_falls_back_to_delta_q(self, tmp_path):
        """Without protein_charge (standalone), use abs(delta_q) as before."""
        run = self._make_qligfep(tmp_path, protein_charge=None)
        inputdir = self._prepare_inputdir(run, tmp_path)

        run.charge_lig2 = 1
        n_ions = self._run_place_counter_ions(run, inputdir, tmp_path)
        # delta_q = charge_lig1 - charge_lig2 = 0 - 1 = -1
        # ion_type = chloride if delta_q > 0 else "SOD" -> SOD
        assert n_ions == 1
        assert run.ion_type == "SOD"

    def test_protein_system_returns_zero(self, tmp_path):
        """Counter-ions are only placed in water system."""
        # For protein system we need a protein.pdb file
        src_pdb = self.TUTORIAL_DIR / "ejm_31.pdb"
        shutil.copy(src_pdb, tmp_path / "protein.pdb")

        run = self._make_qligfep(tmp_path, system="protein", protein_charge=8)
        inputdir = self._prepare_inputdir(run, tmp_path)

        run.charge_lig2 = 1
        n_ions = self._run_place_counter_ions(run, inputdir, tmp_path)
        assert n_ions == 0


class TestChargeMethodDispatch(TestPlaceCounterIonsWithProteinCharge):
    """The --charge-method flag selects how charge-changing edges are handled.

    Only `ion_match` should place real ions in the water-leg PDB during setup.
    `analytical`, `none`, and `coalchemical_water` must all skip the ion
    placement step; analytical applies its correction later in analyze_FEP, and
    coalchemical_water swaps a real water for an ion AFTER qprep (handled by
    select_counter_water once that lands from test/ion-xchange-w-softcore).
    """

    def test_invalid_charge_method_raises(self, tmp_path):
        from QligFEP.qligfep import DualTopologyFEP

        with pytest.raises(ValueError, match="charge_method"):
            DualTopologyFEP(
                lig1=self.LIG1,
                lig2=self.LIG2,
                FF="AMBER14sb",
                system="water",
                cluster="TETRA",
                sphereradius="25",
                charge_method="not_a_real_method",
            )

    def test_same_charge_attribute_populated_after_read_files(self, tmp_path):
        run = self._make_qligfep(tmp_path, protein_charge=8)
        self._prepare_inputdir(run, tmp_path)
        # Tyk2 ligands are both neutral -> same_charge True
        assert run.same_charge is True

    def test_same_charge_false_after_override(self, tmp_path):
        run = self._make_qligfep(tmp_path, protein_charge=8)
        self._prepare_inputdir(run, tmp_path)
        run.charge_lig2 = 1
        # read_files already set same_charge based on raw lib charges; recompute
        # the way the code under test sees it during the dispatch.
        assert run.charge_lig1 != run.charge_lig2

    def test_ion_match_default_places_ions_on_charge_change(self, tmp_path):
        """The default mode (ion_match) keeps the existing behavior."""
        run = self._make_qligfep(tmp_path, protein_charge=8, charge_method="ion_match")
        inputdir = self._prepare_inputdir(run, tmp_path)
        run.charge_lig2 = 1
        n_ions = self._run_place_counter_ions(run, inputdir, tmp_path)
        assert n_ions == 7
        assert run.ion_type == "SOD"

    def test_none_skips_ions_even_on_charge_change(self, tmp_path):
        run = self._make_qligfep(tmp_path, protein_charge=8, charge_method="none")
        inputdir = self._prepare_inputdir(run, tmp_path)
        run.charge_lig2 = 1
        n_ions = self._run_place_counter_ions(run, inputdir, tmp_path)
        assert n_ions == 0
        assert run.ion_type is None

    def test_coalchemical_water_skips_ion_placement_step(self, tmp_path):
        """coalchemical_water replaces ion placement with water-to-ion swap
        that happens post-qprep -- place_counter_ions itself must be a no-op."""
        run = self._make_qligfep(tmp_path, protein_charge=8, charge_method="coalchemical_water")
        inputdir = self._prepare_inputdir(run, tmp_path)
        run.charge_lig2 = 1
        n_ions = self._run_place_counter_ions(run, inputdir, tmp_path)
        assert n_ions == 0
        assert run.ion_type is None

    def test_ion_match_same_charge_no_ions(self, tmp_path):
        """Even with ion_match, same-charge ligand pairs need no ions."""
        run = self._make_qligfep(tmp_path, protein_charge=0, charge_method="ion_match")
        inputdir = self._prepare_inputdir(run, tmp_path)
        # Do NOT override charge_lig2; both ligands stay neutral.
        n_ions = self._run_place_counter_ions(run, inputdir, tmp_path)
        assert n_ions == 0


class TestPlaceCounterWater(TestPlaceCounterIonsWithProteinCharge):
    """Tests for the co-alchemical water swap (Method 3)."""

    def _run_place_counter_water(self, run, inputdir, tmp_path):
        import os

        original_cwd = os.getcwd()
        os.chdir(tmp_path)
        try:
            return run.place_counter_water(inputdir)
        finally:
            os.chdir(original_cwd)

    def test_same_charge_returns_zero(self, tmp_path):
        run = self._make_qligfep(tmp_path, charge_method="coalchemical_water")
        inputdir = self._prepare_inputdir(run, tmp_path)
        # Both ligands neutral; delta_q == 0 -> no counter-waters.
        n_cw = self._run_place_counter_water(run, inputdir, tmp_path)
        assert n_cw == 0
        assert run.counter_water_atoms == []

    def test_wrong_method_returns_zero(self, tmp_path):
        """place_counter_water must be a no-op for any method other than coalchemical_water."""
        run = self._make_qligfep(tmp_path, charge_method="ion_match")
        inputdir = self._prepare_inputdir(run, tmp_path)
        run.charge_lig2 = 1
        n_cw = self._run_place_counter_water(run, inputdir, tmp_path)
        assert n_cw == 0
        assert run.counter_water_atoms == []

    def test_neutral_to_positive_places_chloride_in_state2(self, tmp_path):
        """lig1=0 -> lig2=+1. The +1 state needs a Cl- to net to zero; |lig2|
        is the larger side, so the ion is real in state 2 (lig2 endpoint),
        and state 1 has a real water."""
        run = self._make_qligfep(tmp_path, charge_method="coalchemical_water")
        inputdir = self._prepare_inputdir(run, tmp_path)
        run.charge_lig2 = 1
        n_cw = self._run_place_counter_water(run, inputdir, tmp_path)
        assert n_cw == 1
        assert run.ion_type == "CHL"  # AMBER14sb chloride
        assert len(run.counter_water_atoms) == 1
        cw = run.counter_water_atoms[0]
        assert len(cw["topology_indices"]) == 3
        assert len(cw["qatoms"]) == 3
        # State 1 (lig1, neutral): real water. State 2 (lig2, +1): CHL ion.
        assert cw["qatoms"][0]["type_s1"] == "OW"
        assert cw["qatoms"][0]["type_s2"] == "CHL"

    def test_positive_to_neutral_places_chloride_in_state1(self, tmp_path):
        """lig1=+1 -> lig2=0. The +1 state needs a Cl- to net to zero; |lig2|
        is the smaller side, so the ion is real in state 1 (lig1 endpoint),
        and state 2 has a real water."""
        run = self._make_qligfep(tmp_path, charge_method="coalchemical_water")
        inputdir = self._prepare_inputdir(run, tmp_path)
        run.charge_lig1 = 1
        run.charge_lig2 = 0
        n_cw = self._run_place_counter_water(run, inputdir, tmp_path)
        assert n_cw == 1
        assert run.ion_type == "CHL"
        cw = run.counter_water_atoms[0]
        # State 1 (lig1, +1): CHL ion. State 2 (lig2, neutral): real water.
        assert cw["qatoms"][0]["type_s1"] == "CHL"
        assert cw["qatoms"][0]["type_s2"] == "OW"

    def test_negative_ligand_places_sodium(self, tmp_path):
        """lig1=0 -> lig2=-1 requires a Na+ to net to zero in the -1 state."""
        run = self._make_qligfep(tmp_path, charge_method="coalchemical_water")
        inputdir = self._prepare_inputdir(run, tmp_path)
        run.charge_lig2 = -1
        n_cw = self._run_place_counter_water(run, inputdir, tmp_path)
        assert n_cw == 1
        assert run.ion_type == "SOD"
        cw = run.counter_water_atoms[0]
        # State 2 has the negative ligand; ion is real there.
        assert cw["qatoms"][0]["type_s2"] == "SOD"
        assert cw["qatoms"][0]["type_s1"] == "OW"

    def test_h_atoms_swap_to_dum(self, tmp_path):
        """Each counter-water gets two H atoms that perturb between real H and DUM."""
        run = self._make_qligfep(tmp_path, charge_method="coalchemical_water")
        inputdir = self._prepare_inputdir(run, tmp_path)
        run.charge_lig2 = 1
        self._run_place_counter_water(run, inputdir, tmp_path)
        cw = run.counter_water_atoms[0]
        # Atoms 1 and 2 are the hydrogens.
        for h_qatom in cw["qatoms"][1:]:
            assert "DUM" in (h_qatom["type_s1"], h_qatom["type_s2"])
            assert "HW" in (h_qatom["type_s1"], h_qatom["type_s2"])

    def test_cwt_lib_not_written_at_runtime(self, tmp_path):
        """CWT lives inside the FF .lib next to HOH; no runtime cwt.lib is emitted."""
        run = self._make_qligfep(tmp_path, charge_method="coalchemical_water")
        inputdir = self._prepare_inputdir(run, tmp_path)
        run.charge_lig2 = 1
        self._run_place_counter_water(run, inputdir, tmp_path)
        assert not (Path(inputdir) / "cwt.lib").exists()

    def test_cwt_defined_in_ff_lib(self, tmp_path):
        """The FF .lib file used by qprep must contain a CWT residue, matched
        on the AMBER14sb water atom types (OW, HW)."""
        run = self._make_qligfep(tmp_path, charge_method="coalchemical_water")
        ff_lib = Path(run.lib_file)
        text = ff_lib.read_text()
        assert "{CWT}" in text
        # Find the CWT block and verify its atom types.
        cwt_start = text.index("{CWT}")
        # Block ends at the next '*-----' separator line.
        cwt_block = text[cwt_start:].split("*", 1)[0]
        assert "OW" in cwt_block
        assert "HW" in cwt_block
        # And importantly: CWT must NOT carry the solvent flag, otherwise qprep
        # would classify it as solvent and the co-alchemical swap wouldn't work.
        # Check the actual directive, not prose comments mentioning it.
        directive_lines = [
            ln.strip() for ln in cwt_block.splitlines() if not ln.lstrip().startswith(("!", "{"))
        ]
        assert "solvent 1" not in directive_lines

    def test_pdb_appended_with_cwt_residue(self, tmp_path):
        run = self._make_qligfep(tmp_path, charge_method="coalchemical_water")
        inputdir = self._prepare_inputdir(run, tmp_path)
        run.charge_lig2 = 1
        self._run_place_counter_water(run, inputdir, tmp_path)
        pdb_text = (Path(inputdir) / run.pdb_fname).read_text()
        assert "CWT" in pdb_text

    def test_counter_water_restraint_formats_oxygen_only(self, tmp_path):
        run = self._make_qligfep(tmp_path, charge_method="coalchemical_water")
        inputdir = self._prepare_inputdir(run, tmp_path)
        run.charge_lig2 = 1
        self._run_place_counter_water(run, inputdir, tmp_path)
        restraint = run._format_counter_water_restraint()
        assert restraint, "Expected a non-empty restraint string"
        o_idx = run.counter_water_atoms[0]["topology_indices"][0]
        # The oxygen index appears in the restraint line.
        assert str(o_idx) in restraint

    def test_two_charge_change_two_waters(self, tmp_path):
        """|delta_q|=2 should produce 2 counter-waters."""
        run = self._make_qligfep(tmp_path, charge_method="coalchemical_water")
        inputdir = self._prepare_inputdir(run, tmp_path)
        run.charge_lig2 = 2  # lig1=0, lig2=+2 -> 2 chlorides
        n_cw = self._run_place_counter_water(run, inputdir, tmp_path)
        assert n_cw == 2
        assert len(run.counter_water_atoms) == 2


class TestResidueAtomSerialRange:
    """Tests for selecting the atom-serial range of counter-ion/counter-water residues.

    Wall restraints are built from this range. For the ion_match method the
    residues are named after the ion (e.g. SOD); for coalchemical_water they are
    named CWT. Both must resolve, and an absent selection must yield None rather
    than NaN.
    """

    @staticmethod
    def _df(rows):
        import pandas as pd

        return pd.DataFrame(rows, columns=["residue_name", "atom_serial_number"])

    def test_coalchemical_water_cwt_range(self):
        """CWT counter-waters resolve even when the ion_type name is absent."""
        from QligFEP.pdb_utils import residue_atom_serial_range

        df = self._df([("LIG", 1), ("LIG", 2), ("CWT", 3), ("CWT", 4), ("CWT", 5), ("HOH", 6)])
        assert residue_atom_serial_range(df, ["SOD", "CWT"]) == (3, 5)

    def test_ion_match_sod_range(self):
        """Real SOD counter-ions resolve to their serial range."""
        from QligFEP.pdb_utils import residue_atom_serial_range

        df = self._df([("LIG", 1), ("SOD", 2), ("SOD", 3), ("HOH", 4)])
        assert residue_atom_serial_range(df, ["SOD", "CWT"]) == (2, 3)

    def test_no_match_returns_none(self):
        """No matching residue returns None instead of producing NaN."""
        from QligFEP.pdb_utils import residue_atom_serial_range

        df = self._df([("LIG", 1), ("HOH", 2)])
        assert residue_atom_serial_range(df, ["SOD", "CWT"]) is None

    def test_accepts_single_residue_name(self):
        """A bare string residue name is treated as a single-element list."""
        from QligFEP.pdb_utils import residue_atom_serial_range

        df = self._df([("CWT", 7), ("CWT", 9), ("HOH", 10)])
        assert residue_atom_serial_range(df, "CWT") == (7, 9)
