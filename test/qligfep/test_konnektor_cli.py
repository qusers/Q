"""Tests for the konnektor-based perturbation network CLI."""

import json
import shutil
import tempfile
from pathlib import Path

import pytest
from rdkit import Chem

from QligFEP.CLI.konnektor_cli import KonnektorWrap

RESOURCES_DIR = Path(__file__).parent / "resources"


@pytest.fixture
def tyk2_sdf() -> Path:
    return RESOURCES_DIR / "tyk2_ligands.sdf"


@pytest.fixture
def cmet_sdf() -> Path:
    return RESOURCES_DIR / "cmet_ligands.sdf"


@pytest.fixture
def tyk2_output_dir():
    tmpdir = Path(tempfile.mkdtemp())
    yield tmpdir
    shutil.rmtree(tmpdir, ignore_errors=True)


@pytest.fixture
def tyk2_mst_result(tyk2_sdf, tyk2_output_dir):
    """Run KonnektorWrap with MST on tyk2 and return the result dict."""
    kw = KonnektorWrap(inp=str(tyk2_sdf), out=str(tyk2_output_dir), network="mst")
    return kw.run()


@pytest.fixture
def tyk2_star_result(tyk2_sdf):
    """Run KonnektorWrap with star network on tyk2."""
    tmpdir = Path(tempfile.mkdtemp())
    kw = KonnektorWrap(inp=str(tyk2_sdf), out=str(tmpdir), network="star")
    result = kw.run()
    shutil.rmtree(tmpdir, ignore_errors=True)
    return result


class TestOutputStructure:
    def test_output_json_has_correct_structure(self, tyk2_mst_result):
        """Output should contain 'nodes' and 'edges' top-level keys."""
        assert "nodes" in tyk2_mst_result
        assert "edges" in tyk2_mst_result

    def test_edges_have_required_fields(self, tyk2_mst_result):
        """Each edge should have weight, source, target, from, to, same_charge."""
        required_fields = {"weight", "source", "target", "from", "to", "same_charge"}
        for edge in tyk2_mst_result["edges"]:
            assert required_fields.issubset(edge.keys()), f"Missing fields in edge: {required_fields - edge.keys()}"

    def test_nodes_have_required_fields(self, tyk2_mst_result):
        """Each node should have smiles and formal_charge."""
        for name, node in tyk2_mst_result["nodes"].items():
            assert "smiles" in node, f"Node '{name}' missing 'smiles'"
            assert "formal_charge" in node, f"Node '{name}' missing 'formal_charge'"

    def test_weights_are_valid_scores(self, tyk2_mst_result):
        """Edge weights should be floats in [0, 1]."""
        for edge in tyk2_mst_result["edges"]:
            assert isinstance(edge["weight"], float)
            assert 0.0 <= edge["weight"] <= 1.0

    def test_edge_names_reference_valid_nodes(self, tyk2_mst_result):
        """Edge 'from' and 'to' should reference existing node names."""
        node_names = set(tyk2_mst_result["nodes"].keys())
        for edge in tyk2_mst_result["edges"]:
            assert edge["from"] in node_names, f"Edge 'from' ({edge['from']}) not in nodes"
            assert edge["to"] in node_names, f"Edge 'to' ({edge['to']}) not in nodes"


class TestNetworkTopologies:
    def test_mst_has_n_minus_1_edges(self, tyk2_mst_result):
        """MST should have exactly N-1 edges for N nodes."""
        n_nodes = len(tyk2_mst_result["nodes"])
        n_edges = len(tyk2_mst_result["edges"])
        assert n_edges == n_nodes - 1

    def test_star_has_n_minus_1_edges(self, tyk2_star_result):
        """Star network should have N-1 edges (all connected to central node)."""
        n_nodes = len(tyk2_star_result["nodes"])
        n_edges = len(tyk2_star_result["edges"])
        assert n_edges == n_nodes - 1

    def test_all_nodes_present(self, tyk2_sdf, tyk2_mst_result):
        """All ligands from input SDF should appear as nodes."""
        supplier = Chem.SDMolSupplier(str(tyk2_sdf))
        expected_names = {m.GetProp("_Name") for m in supplier if m is not None}
        actual_names = set(tyk2_mst_result["nodes"].keys())
        assert actual_names == expected_names


class TestChargeHandling:
    def test_charged_perturbations_flagged(self, cmet_sdf):
        """Edges between differently charged ligands should have same_charge=False."""
        tmpdir = Path(tempfile.mkdtemp())
        kw = KonnektorWrap(inp=str(cmet_sdf), out=str(tmpdir), network="mst")
        result = kw.run()
        shutil.rmtree(tmpdir, ignore_errors=True)

        # cmet has both neutral and charged ligands, so some edges should have same_charge=False
        same_charges = [e["same_charge"] for e in result["edges"]]
        assert not all(same_charges), "Expected some edges with same_charge=False in cmet dataset"


class TestSetupFEPCompatibility:
    def test_output_compatible_with_setupFEP(self, tyk2_sdf, tyk2_output_dir):
        """Output JSON should be parseable by setupFEP's ligpairs_from_json."""
        from QligFEP.CLI.setupFEP import ligpairs_from_json

        kw = KonnektorWrap(inp=str(tyk2_sdf), out=str(tyk2_output_dir), network="mst")
        kw.run()

        json_path = tyk2_output_dir / "mapping.json"
        assert json_path.exists(), "mapping.json should be written to output directory"

        pairs = ligpairs_from_json(str(json_path))
        assert len(pairs) > 0
        for lig1, lig2 in pairs:
            assert isinstance(lig1, str)
            assert isinstance(lig2, str)


class TestJsonOutput:
    def test_json_file_written(self, tyk2_sdf, tyk2_output_dir):
        """Run should write mapping.json to the output directory."""
        kw = KonnektorWrap(inp=str(tyk2_sdf), out=str(tyk2_output_dir), network="mst")
        kw.run()
        assert (tyk2_output_dir / "mapping.json").exists()

    def test_json_roundtrips(self, tyk2_sdf, tyk2_output_dir):
        """Written JSON should be valid and match the returned dict."""
        kw = KonnektorWrap(inp=str(tyk2_sdf), out=str(tyk2_output_dir), network="mst")
        result = kw.run()

        with open(tyk2_output_dir / "mapping.json") as f:
            loaded = json.load(f)

        assert loaded["nodes"] == result["nodes"]
        assert len(loaded["edges"]) == len(result["edges"])


class TestUserDefinedProperties:
    def test_user_defined_properties_in_nodes(self, tyk2_sdf, tyk2_output_dir):
        """SDF properties (like r_exp_dg) should flow through to nodes."""
        kw = KonnektorWrap(inp=str(tyk2_sdf), out=str(tyk2_output_dir), network="mst")
        result = kw.run()

        # tyk2 SDF has user-defined properties; check they're present
        has_extra_props = False
        for name, node in result["nodes"].items():
            extra_keys = set(node.keys()) - {"smiles", "formal_charge"}
            if extra_keys:
                has_extra_props = True
                break
        assert has_extra_props, "Expected user-defined SDF properties in nodes"

    def test_delta_properties_on_edges(self, tyk2_sdf, tyk2_output_dir):
        """Numerical user-defined properties should produce delta_* fields on edges."""
        kw = KonnektorWrap(inp=str(tyk2_sdf), out=str(tyk2_output_dir), network="mst")
        result = kw.run()

        # Check if any edge has delta_ prefixed keys
        has_deltas = any(any(k.startswith("delta_") for k in e.keys()) for e in result["edges"])
        assert has_deltas, "Expected delta_ properties on edges from numerical SDF properties"
