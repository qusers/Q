"""Tests for the lomap-based perturbation network CLI."""

import shutil
import tempfile
from pathlib import Path

import pytest

from QligFEP.CLI.lomap_wrap_cli import LomapWrap

RESOURCES_DIR = Path(__file__).parent / "resources"


@pytest.fixture
def tyk2_sdf() -> Path:
    return RESOURCES_DIR / "tyk2_ligands.sdf"


@pytest.fixture
def tyk2_lomap_result(tyk2_sdf):
    """Run LomapWrap on tyk2 and return the result dict."""
    outdir = Path(tempfile.mkdtemp()) / "lomap_out"
    lw = LomapWrap(inp=str(tyk2_sdf), out=str(outdir))
    result = lw.run_lomap()
    yield result
    shutil.rmtree(outdir, ignore_errors=True)


class TestExpKey:
    @pytest.fixture
    def tyk2_result_with_exp_key(self, tyk2_sdf):
        """Run LomapWrap with exp_key on tyk2."""
        outdir = Path(tempfile.mkdtemp()) / "lomap_exp"
        lw = LomapWrap(inp=str(tyk2_sdf), out=str(outdir), exp_key="r_exp_dg")
        result = lw.run_lomap()
        yield result
        shutil.rmtree(outdir, ignore_errors=True)

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

    def test_no_exp_key_no_deltas(self, tyk2_lomap_result):
        """Without exp_key, nodes keep original props and edges have no deltas."""
        # Nodes should still have original property names
        has_original = any("r_exp_dg" in node for node in tyk2_lomap_result["nodes"].values())
        assert has_original, "Expected original 'r_exp_dg' property on nodes"

        # Edges should have no delta_ or ddg_value keys
        for edge in tyk2_lomap_result["edges"]:
            delta_keys = [k for k in edge if k.startswith("delta_") or k == "ddg_value"]
            assert delta_keys == [], f"Edge has unexpected keys: {delta_keys}"
