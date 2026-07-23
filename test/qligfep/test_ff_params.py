"""Force-field parameter file invariants for AMBER14sb.

The [atom_types] block of a Q AMBER .prm stores six values per type:

    R*1   R*2   epsilon1   R*3   epsilon2&3   mass

Column 5 (epsilon2&3) is the 1-4 Lennard-Jones well depth. AMBER scales 1-4
vdW by 1/SCNB = 1/2 with the same R*, so column 5 must equal column 3 / 2.
Hand-added blocks (phospho-AA, POPC Lipid21, ions) once stored R*/2 here by
mistake; this test guards against that class of regression.
"""

from decimal import Decimal
from pathlib import Path

import pytest

PRM = Path(__file__).resolve().parents[2] / "src" / "QligFEP" / "FF" / "AMBER14sb.prm"


def _atom_type_rows():
    rows = []
    in_block = False
    for line in PRM.read_text().splitlines():
        s = line.strip()
        if s == "[atom_types]":
            in_block = True
            continue
        if in_block and s.startswith("[") and s != "[atom_types]":
            break
        if not in_block or not s or s.startswith(("*", "!")):
            continue
        parts = s.split()
        if len(parts) < 7:
            continue
        name = parts[0]
        try:
            cols = [Decimal(p) for p in parts[1:7]]
        except Exception:
            continue
        rows.append((name, cols))
    return rows


def test_prm_exists_and_has_atom_types():
    rows = _atom_type_rows()
    assert len(rows) > 50, f"expected many atom types, parsed {len(rows)}"


@pytest.mark.parametrize("name,cols", _atom_type_rows(), ids=lambda v: v if isinstance(v, str) else "")
def test_one_four_epsilon_is_half_of_epsilon(name, cols):
    eps, eps14 = cols[2], cols[4]
    expected = eps / 2
    # tolerance absorbs rounding (e.g. 0.015191 vs 0.0151915); the historical bug
    # stored R*/2, which is off by orders of magnitude and always caught.
    tol = abs(expected) * Decimal("0.01") + Decimal("1e-6")
    assert abs(eps14 - expected) <= tol, (
        f"{name}: column 5 (1-4 epsilon) = {eps14}, expected eps/2 = {expected} "
        f"(eps={eps}, R*/2={cols[0] / 2})"
    )
