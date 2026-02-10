"""Tests for center of geometry calculation across SDF formats (V2000 and V3000)."""

from pathlib import Path

import pytest

from QligFEP.CLI.cog_cli import MolecularCOG

RESOURCES = Path(__file__).parent / "resources"


@pytest.fixture
def tyk2_v2000_sdf():
    return RESOURCES / "tyk2_ligands.sdf"


@pytest.fixture
def v3000_sdf():
    return Path("/Users/davidararipe/projects/Q/qparams-hydrogen-handling/ligands.sdf")


def _parse_cog_string(cog_str):
    """Parse '[x y z]' string into list of floats."""
    return [float(v) for v in cog_str.strip("[]").split()]


class TestSdfCog:
    def test_v2000_sdf_returns_valid_cog(self, tyk2_v2000_sdf):
        cog = MolecularCOG(tyk2_v2000_sdf)
        result = cog()
        coords = _parse_cog_string(result)
        assert len(coords) == 3
        assert all(isinstance(c, float) for c in coords)

    def test_v2000_sdf_heavy_atom_cog(self, tyk2_v2000_sdf):
        """COG should be computed from heavy atoms only."""
        cog = MolecularCOG(tyk2_v2000_sdf)
        result = cog()
        coords = _parse_cog_string(result)
        assert coords == pytest.approx([-4.589, 25.490, -30.631], abs=0.01)

    def test_v3000_sdf_returns_valid_cog(self, v3000_sdf):
        cog = MolecularCOG(v3000_sdf)
        result = cog()
        coords = _parse_cog_string(result)
        assert len(coords) == 3
        assert all(isinstance(c, float) for c in coords)

    def test_v3000_sdf_heavy_atom_cog(self, v3000_sdf):
        """COG should be computed from heavy atoms only, even for V3000 format."""
        cog = MolecularCOG(v3000_sdf)
        result = cog()
        coords = _parse_cog_string(result)
        assert coords == pytest.approx([-1.128, 0.150, 0.004], abs=0.01)
