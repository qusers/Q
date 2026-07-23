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

    def test_star_with_central_ligand(self, tyk2_sdf):
        """Star with central_ligand should connect all edges to that ligand."""
        tmpdir = Path(tempfile.mkdtemp())
        central = "ejm_31"
        kw = KonnektorWrap(inp=str(tyk2_sdf), out=str(tmpdir), network="star", central_ligand=central)
        result = kw.run()
        shutil.rmtree(tmpdir, ignore_errors=True)

        # Every edge should have ejm_31 as either 'from' or 'to'
        for edge in result["edges"]:
            assert central in (edge["from"], edge["to"]), (
                f"Edge {edge['from']} -> {edge['to']} does not involve central ligand '{central}'"
            )
        # Should still have N-1 edges
        assert len(result["edges"]) == len(result["nodes"]) - 1

    def test_star_with_invalid_central_ligand_raises(self, tyk2_sdf):
        """Passing a non-existent central ligand name should raise ValueError."""
        tmpdir = Path(tempfile.mkdtemp())
        kw = KonnektorWrap(inp=str(tyk2_sdf), out=str(tmpdir), network="star", central_ligand="nonexistent")
        with pytest.raises(ValueError, match="not found"):
            kw.run()
        shutil.rmtree(tmpdir, ignore_errors=True)

    def test_redundant_mst_has_more_edges_than_mst(self, tyk2_sdf):
        """Redundant MST should have more edges than a plain MST."""
        tmpdir = Path(tempfile.mkdtemp())
        kw = KonnektorWrap(inp=str(tyk2_sdf), out=str(tmpdir), network="rmst")
        result = kw.run()
        shutil.rmtree(tmpdir, ignore_errors=True)

        n_nodes = len(result["nodes"])
        n_edges = len(result["edges"])
        # Must have more edges than MST (N-1) due to redundancy
        assert n_edges > n_nodes - 1

    def test_redundant_mst_custom_redundancy(self, tyk2_sdf):
        """Redundant MST with n_redundancy=3 should have more edges than n_redundancy=2."""
        tmpdir2 = Path(tempfile.mkdtemp())
        tmpdir3 = Path(tempfile.mkdtemp())

        kw2 = KonnektorWrap(inp=str(tyk2_sdf), out=str(tmpdir2), network="rmst", n_redundancy=2)
        result2 = kw2.run()
        kw3 = KonnektorWrap(inp=str(tyk2_sdf), out=str(tmpdir3), network="rmst", n_redundancy=3)
        result3 = kw3.run()

        shutil.rmtree(tmpdir2, ignore_errors=True)
        shutil.rmtree(tmpdir3, ignore_errors=True)

        assert len(result3["edges"]) >= len(result2["edges"])

    def test_nedges_produces_connected_network(self, tyk2_sdf):
        """N-edges network should produce a valid connected graph."""
        tmpdir = Path(tempfile.mkdtemp())
        kw = KonnektorWrap(inp=str(tyk2_sdf), out=str(tmpdir), network="nedges", connectivity=3)
        result = kw.run()
        shutil.rmtree(tmpdir, ignore_errors=True)

        n_nodes = len(result["nodes"])
        n_edges = len(result["edges"])
        # At minimum a connected graph has N-1 edges
        assert n_edges >= n_nodes - 1

    def test_cyclic_has_more_edges_than_mst(self, tyk2_sdf):
        """Cyclic network should have more edges than MST due to cycle closure."""
        tmpdir = Path(tempfile.mkdtemp())
        kw = KonnektorWrap(inp=str(tyk2_sdf), out=str(tmpdir), network="cyclic")
        result = kw.run()
        shutil.rmtree(tmpdir, ignore_errors=True)

        n_nodes = len(result["nodes"])
        n_edges = len(result["edges"])
        assert n_edges > n_nodes - 1


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

    def test_separate_charges_no_cross_charge_edges(self, cmet_sdf):
        """With separate_charges=True, no edges should cross charge boundaries."""
        tmpdir = Path(tempfile.mkdtemp())
        kw = KonnektorWrap(inp=str(cmet_sdf), out=str(tmpdir), network="mst", separate_charges=True)
        result = kw.run()
        shutil.rmtree(tmpdir, ignore_errors=True)

        for edge in result["edges"]:
            assert edge["same_charge"], (
                f"Edge {edge['from']} -> {edge['to']} crosses charge boundary with separate_charges=True"
            )

    def test_separate_charges_all_nodes_present(self, cmet_sdf):
        """With separate_charges=True, all ligands should still appear in nodes."""
        tmpdir = Path(tempfile.mkdtemp())
        kw = KonnektorWrap(inp=str(cmet_sdf), out=str(tmpdir), network="mst", separate_charges=True)
        result = kw.run()
        shutil.rmtree(tmpdir, ignore_errors=True)

        supplier = Chem.SDMolSupplier(str(cmet_sdf))
        expected_names = {m.GetProp("_Name") for m in supplier if m is not None}
        actual_names = set(result["nodes"].keys())
        assert actual_names == expected_names


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
        for lig1, lig2, same_charge in pairs:
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

    def test_no_exp_key_no_deltas(self, tyk2_sdf, tyk2_output_dir):
        """Without exp_key, nodes keep original props and edges have no deltas."""
        kw = KonnektorWrap(inp=str(tyk2_sdf), out=str(tyk2_output_dir), network="mst")
        result = kw.run()

        # Nodes should still have original property names
        has_original = any("r_exp_dg" in node for node in result["nodes"].values())
        assert has_original, "Expected original 'r_exp_dg' property on nodes"

        # Edges should have no delta_ or ddg_value keys
        for edge in result["edges"]:
            delta_keys = [k for k in edge if k.startswith("delta_") or k == "ddg_value"]
            assert delta_keys == [], f"Edge has unexpected keys: {delta_keys}"


class TestExpKey:
    @pytest.fixture
    def tyk2_result_with_exp_key(self, tyk2_sdf):
        """Run KonnektorWrap with exp_key on tyk2."""
        tmpdir = Path(tempfile.mkdtemp())
        kw = KonnektorWrap(inp=str(tyk2_sdf), out=str(tmpdir), network="mst", exp_key="r_exp_dg")
        result = kw.run()
        yield result
        shutil.rmtree(tmpdir, ignore_errors=True)

    def test_exp_key_renames_to_dg_value_on_nodes(self, tyk2_result_with_exp_key):
        """When exp_key is set, nodes should have dg_value, not the original key."""
        for name, node in tyk2_result_with_exp_key["nodes"].items():
            assert "dg_value" in node, f"Node '{name}' missing 'dg_value'"
            assert "r_exp_dg" not in node, f"Node '{name}' still has original 'r_exp_dg'"

    def test_exp_key_computes_ddg_value_on_edges(self, tyk2_result_with_exp_key):
        """Edges should have ddg_value = dg_value(to) - dg_value(from)."""
        nodes = tyk2_result_with_exp_key["nodes"]
        for edge in tyk2_result_with_exp_key["edges"]:
            assert "ddg_value" in edge, f"Edge {edge['from']}->{edge['to']} missing 'ddg_value'"
            expected = nodes[edge["to"]]["dg_value"] - nodes[edge["from"]]["dg_value"]
            assert edge["ddg_value"] == pytest.approx(expected), (
                f"Edge {edge['from']}->{edge['to']}: ddg_value={edge['ddg_value']}, "
                f"expected {expected}"
            )

    def test_exp_key_no_delta_prefix_keys(self, tyk2_result_with_exp_key):
        """No edge should have any key starting with delta_."""
        for edge in tyk2_result_with_exp_key["edges"]:
            delta_keys = [k for k in edge if k.startswith("delta_")]
            assert delta_keys == [], f"Edge {edge['from']}->{edge['to']} has delta_ keys: {delta_keys}"
