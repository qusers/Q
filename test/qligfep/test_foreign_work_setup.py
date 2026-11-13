"""Setup/runner coverage without external tutorial preparation artifacts."""

import subprocess

import pytest
import pandas as pd

from QligFEP.qligfep import QligFEP, nonlinear_softcore_activation, zero_lj_atom_types


@pytest.mark.parametrize("start", ["1", "0.5"])
@pytest.mark.parametrize("cluster", ["LOCAL", "SNELLIUS"])
@pytest.mark.parametrize("method", ["gapsys", "beutler"])
def test_generated_workflow_preserves_modernization(tmp_path, start, cluster, method):
    run = QligFEP("a", "b", "AMBER14sb", "water", cluster=cluster, start=start,
                  softcore_method=method, minimize=True)
    lambdas = ["1.000", "0.750", "0.500", "0.250", "0.000"]
    files = run.write_md_files(lambdas, str(tmp_path), 2, 2, [])
    run.write_runfile(str(tmp_path), files)
    eq1 = (tmp_path / "eq1.inp").read_text()
    assert "minimize                  on" in eq1
    assert "[foreign_lambdas]" not in eq1
    for path in tmp_path.glob("md*.inp"):
        content = path.read_text()
        assert "constraint_algorithm      lincs settle" in content
        assert "foreign_hamiltonian_mode  correction" in content
        source = content.split("[lambdas]\n")[1].splitlines()[0].split()[1]
        targets = content.split("[foreign_lambdas]\n")[1].split("[")[0].split()
        assert source in targets
        index = list(reversed(lambdas)).index(source)
        assert targets == list(reversed(lambdas))[max(0, index-1):index+2]
    script_path = tmp_path / f"run{cluster}.sh"
    script = script_path.read_text()
    assert 'qligfep_fhe_analyze . -T "$temperature"' in script
    assert "*.fhe" not in script
    assert subprocess.run(["bash", "-n", str(script_path)], capture_output=True).returncode == 0


def test_activation_is_derived_from_endpoint_types():
    assert nonlinear_softcore_activation([(1, "C", "DUM"), (2, "DUM", "H"), (3, "C", "N"), (4, "DUM", "DUM")]) == {
        1: (1, 0), 2: (0, 1), 3: (0, 0), 4: (0, 0),
    }


def test_activation_covers_atoms_whose_lj_wall_vanishes():
    vdw = ["C  1.9080  0.0860  0  0  1.9080  0.0430  12.01", "HO  0.0000  0.0000  0  0  0.0000  0.0000  1.008",
           "N  1.8240  0.1700  0  0  1.8240  0.0850  14.01"]
    zero = zero_lj_atom_types(vdw)
    assert zero == {"DUM", "HO"}
    assert nonlinear_softcore_activation(
        [(1, "C", "HO"), (2, "HO", "N"), (3, "HO", "DUM"), (4, "C", "N"), (5, "HO", "HO")], zero
    ) == {1: (1, 0), 2: (0, 1), 3: (1, 0), 4: (0, 0), 5: (0, 0)}


def test_neq_is_rejected_at_setup():
    with pytest.raises(ValueError, match="equilibrium FEP only"):
        QligFEP("a", "b", "AMBER14sb", "water", neq=True, softcore_method="gapsys")


@pytest.mark.parametrize("method", ["standard", "beutler", "gapsys"])
def test_fep_writer_and_md_writer_share_method(tmp_path, monkeypatch, method):
    monkeypatch.setattr("QligFEP.qligfep.read_pdb_to_dataframe",
                        lambda _: pd.DataFrame({"residue_name": ["LIG"]}))
    run = QligFEP("a", "b", "AMBER14sb", "water", cluster="LOCAL")
    changes = [(1, "C", "DUM"), (2, "H", "DUM"), (3, "DUM", "C"), (4, "DUM", "H")]
    charges = [(1, -.1, 0), (2, .1, 0), (3, 0, -.1), (4, 0, .1)]
    run.write_FEP_file(charges, changes, ["C 1.9 0.1 0 0 1.9 0.05 12.011",
                                        "H 0 0 0 0 0 0 1.008"], str(tmp_path), 2, 2,
                       softcore_method=method)
    content = (tmp_path / "FEP1.fep").read_text()
    flags = [line.split() for line in content.split("[softcore]")[1].split("[")[0].splitlines() if line.strip()]
    if method in {"gapsys", "beutler"}:
        assert [row[1:] for row in flags] == [["1", "0"], ["1", "0"], ["0", "1"], ["0", "1"]]
        assert "softcore_use_max_potential off" in content
    else:
        assert [row[1:] for row in flags] == [["0", "20"], ["0", "20"], ["20", "0"], ["20", "0"]]
        assert "softcore_use_max_potential on" in content
    assert ("softcore_method" in content) == (method != "standard")
    run.write_md_files(["1.000", "0.500", "0.000"], str(tmp_path), 2, 2, [])
    assert ("[foreign_lambdas]" in (tmp_path / "md_1000_0000.inp").read_text()) == (method in {"gapsys", "beutler"})
