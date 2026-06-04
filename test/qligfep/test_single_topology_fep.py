"""Golden-file tests for the single-topology hybrid builder (SingleTopologyFEP).

Reproduces the validated ejm_31 -> ejm_47 single-topology build (the one that gave
ddG = -0.19, matching experiment) from committed input ligands, comparing the
generated hybrid .lib and FEP file against committed golden references. Uses the
water-leg layout (FEP atom offset 0, ligand is the only solute), so no qprep /
engine run is needed -- this exercises the pure file-generation logic.
"""
import os
import shutil
from pathlib import Path

import pytest

RES = Path(__file__).parent / "resources" / "single_topology"
INPUTS = RES / "inputs"
GOLDEN = RES / "golden"
LIG_A, LIG_B = "ejm_31", "ejm_47"


# ---------------------------------------------------------------- parsers ----
def parse_lib_atoms(path):
    """[(name, type, charge)] from a Q .lib [atoms] block."""
    atoms, blk = [], None
    for ln in Path(path).read_text().splitlines():
        s = ln.split()
        if not s:
            continue
        if s[0].startswith("["):
            blk = s[0]
            continue
        if blk == "[atoms]" and s[0].isdigit():
            atoms.append((s[1], s[2], round(float(s[3]), 4)))
    return atoms


def parse_lib_pairs(path, section):
    """List of token-tuples for a [bonds]/[impropers] section."""
    out, blk = [], None
    for ln in Path(path).read_text().splitlines():
        s = ln.split()
        if not s:
            continue
        if s[0].startswith("["):
            blk = s[0]
            continue
        if blk == section:
            out.append(tuple(s))
    return out


def parse_fep_blocks(path):
    """{'[block]': [[tokens], ...]} for a Q FEP file."""
    blocks, cur = {}, None
    for ln in Path(path).read_text().splitlines():
        st = ln.strip()
        if st.startswith("[") and st.endswith("]"):
            cur = st
            blocks[cur] = []
            continue
        if cur and st and not st.startswith("!"):
            blocks[cur].append(st.split())
    return blocks


# ---------------------------------------------------------------- fixture ----
@pytest.fixture(scope="module")
def built(tmp_path_factory):
    """Build the ejm_31 -> ejm_47 hybrid once and expose the generated files."""
    from QligFEP.single_topology_fep import SingleTopologyFEP

    work = tmp_path_factory.mktemp("singletop_build")
    for lig in (LIG_A, LIG_B):
        for ext in ("lib", "prm", "pdb", "sdf"):
            shutil.copy(INPUTS / f"{lig}.{ext}", work / f"{lig}.{ext}")

    cwd = os.getcwd()
    os.chdir(work)
    try:
        run = SingleTopologyFEP(
            lig1=LIG_A, lig2=LIG_B, FF="AMBER14sb",
            system="water", cluster="TETRA", sphereradius="15",
        )
        writedir = run.makedir()
        inputdir = writedir + "/inputfiles"
        a = run.read_files()
        run.change_lib(a[0][1], inputdir)
        fep_vdw = run.change_prm(a[0][1], inputdir)
        run.write_FEP_file(
            a[1][0], a[1][1], fep_vdw, inputdir, a[2][0], a[2][1],
            softcore_method="gapsys", single_hamiltonian=True,
        )
    finally:
        os.chdir(cwd)

    return {
        "run": run,
        "lib": Path(inputdir) / "singletop.lib",
        "fep": Path(inputdir) / "FEP1.fep",
    }


# ------------------------------------------------------------------ tests ----
def test_hybrid_partition(built):
    """41 hybrid atoms: 30 shared (morph) + 2 vanish + 9 appear; endpoints 32/39."""
    roles = [h["role"] for h in built["run"].hyb]
    assert len(roles) == 41
    assert roles.count("shared") == 30
    assert roles.count("vanish") == 2
    assert roles.count("appear") == 9
    n_real_a = sum(r in ("shared", "vanish") for r in roles)
    n_real_b = sum(r in ("shared", "appear") for r in roles)
    assert (n_real_a, n_real_b) == (32, 39)


def test_endpoint_net_charges_neutral(built):
    hyb = built["run"].hyb
    q_a = sum(h["sA"][1] for h in hyb)
    q_b = sum(h["sB"][1] for h in hyb)
    assert abs(q_a) < 0.02
    assert abs(q_b) < 0.02


def test_lib_atoms_match_golden(built):
    assert parse_lib_atoms(built["lib"]) == parse_lib_atoms(GOLDEN / "singletop.lib")


def test_lib_bonds_match_golden(built):
    gen = parse_lib_pairs(built["lib"], "[bonds]")
    gold = parse_lib_pairs(GOLDEN / "singletop.lib", "[bonds]")
    assert sorted(gen) == sorted(gold)


def test_lib_impropers_match_golden(built):
    gen = parse_lib_pairs(built["lib"], "[impropers]")
    gold = parse_lib_pairs(GOLDEN / "singletop.lib", "[impropers]")
    assert sorted(gen) == sorted(gold)


def test_fep_change_charges_match_golden(built):
    def norm(block):
        return [(int(r[0]), round(float(r[1]), 4), round(float(r[2]), 4)) for r in block]

    gen = parse_fep_blocks(built["fep"])["[change_charges]"]
    gold = parse_fep_blocks(GOLDEN / "singletop.fep")["[change_charges]"]
    assert norm(gen) == norm(gold)


def test_fep_change_atoms_match_golden(built):
    gen = parse_fep_blocks(built["fep"])["[change_atoms]"]
    gold = parse_fep_blocks(GOLDEN / "singletop.fep")["[change_atoms]"]
    assert gen == gold


def test_fep_softcore_match_golden(built):
    gen = parse_fep_blocks(built["fep"])["[softcore]"]
    gold = parse_fep_blocks(GOLDEN / "singletop.fep")["[softcore]"]
    assert gen == gold


def test_fep_atom_types_match_golden(built):
    gen = {r[0]: r[1:] for r in parse_fep_blocks(built["fep"])["[atom_types]"]}
    gold = {r[0]: r[1:] for r in parse_fep_blocks(GOLDEN / "singletop.fep")["[atom_types]"]}
    assert set(gen) == set(gold)
    for t, vals in gold.items():
        assert gen[t] == vals, t


def test_fep_header(built):
    text = built["fep"].read_text()
    assert "states 2" in text
    assert "softcore_use_max_potential on" in text
    assert "softcore_method gapsys" in text
    assert "single_hamiltonian on" in text


def test_fep_atoms_offset_zero_for_water(built):
    """Water leg: ligand is the only solute, so topology index == FEP index."""
    atoms = parse_fep_blocks(built["fep"])["[atoms]"]
    assert [(int(r[0]), int(r[1])) for r in atoms] == [(k, k) for k in range(1, 42)]
