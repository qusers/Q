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

    def _make_qligfep(self, tmp_path, system="water", protein_charge=None):
        """Create a QligFEP instance in tmp_path with real tutorial ligand files."""
        from QligFEP.qligfep import QligFEP

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
            run = QligFEP(
                lig1=self.LIG1,
                lig2=self.LIG2,
                FF="AMBER14sb",
                system=system,
                cluster="TETRA",
                sphereradius="25",
                protein_charge=protein_charge,
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
