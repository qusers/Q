"""Public softcore names and fail-closed setup, without external tutorial data.

Generated FEP flags and production scripts are tested in test_foreign_work_setup.
"""

import pytest

from QligFEP.qligfep import QligFEP


@pytest.mark.parametrize("method", ["standard", "gapsys", "beutler"])
def test_supported_method_names(method):
    run = QligFEP("a", "b", "AMBER14sb", "water", softcore_method=method)
    assert run.softcore_method == method


@pytest.mark.parametrize("method", ["gapsys_v2", "beutler_coul", "invalid_method"])
def test_removed_names_are_not_aliases(method):
    with pytest.raises(ValueError, match="softcore"):
        QligFEP("a", "b", "AMBER14sb", "water", softcore_method=method)


@pytest.mark.parametrize("method", ["gapsys", "beutler"])
def test_nonlinear_neq_is_rejected(method):
    with pytest.raises(ValueError, match="equilibrium FEP only"):
        QligFEP("a", "b", "AMBER14sb", "water", neq=True, softcore_method=method)


@pytest.mark.parametrize("method", ["standard", "gapsys", "beutler", "gapsys_v2", "beutler_coul"])
def test_cli_method_choices(monkeypatch, method):
    from QligFEP.CLI.parser_base import parse_arguments

    monkeypatch.setattr("sys.argv", ["qligfep", "-l1", "a", "-l2", "b", "-FF", "AMBER14sb",
                                     "-s", "water", "-c", "LOCAL", "-sc", method])
    if method in {"standard", "gapsys", "beutler"}:
        assert parse_arguments("QligFEP").softcore_method == method
    else:
        with pytest.raises(SystemExit) as error:
            parse_arguments("QligFEP")
        assert error.value.code == 2
