"""Amino-acid side-chain composition, used to build hybrid residues for QresFEP.

A dual-topology mutation carries both side chains at once: the wild-type atoms
keep their names and the mutant atoms are lower-cased, so a single residue holds
``CB``/``cb``, ``CG``/``cg`` and so on. Deciding which atoms belong to the side
chain, and which wild-type atom each mutant atom should be restrained onto,
needs an explicit table -- the force field library alone does not say where the
backbone ends.
"""

from __future__ import annotations

from .functions import euclidian_overlap

#: Atoms shared by both topologies and therefore never duplicated or softcored.
#: Glycine mutations are the exception: ``HA`` becomes part of the perturbed set,
#: so :func:`backbone_atoms` drops it.
BACKBONE_ATOMS = ("C", "O", "CA", "N", "H", "HA")

#: Every side-chain atom name the supported force fields use. Membership in this
#: set is how atoms of the hybrid residue are recognised in the Q topology.
SIDE_CHAIN_ATOM_NAMES = (
    "N", "H", "C", "O",
    "CA",
    "HA", "HA2", "HA3",
    "CB",
    "HB", "HB1", "HB2", "HB3",
    "CG", "CG1", "CG2",
    "OG", "OG1",
    "SG",
    "HG", "HG1", "HG2", "HG3", "HG11", "HG12", "HG13", "HG21", "HG22", "HG23",
    "CD", "CD1", "CD2",
    "ND1", "ND2",
    "OD1", "OD2",
    "SD",
    "HD", "HD1", "HD2", "HD3", "HD11", "HD12", "HD13", "HD21", "HD22", "HD23",
    "CE", "CE1", "CE2", "CE3",
    "OE1", "OE2",
    "NE", "NE1", "NE2",
    "HE", "HE1", "HE2", "HE3", "HE21", "HE22",
    "CZ", "CZ2", "CZ3",
    "NZ",
    "HZ", "HZ1", "HZ2", "HZ3",
    "CH2",
    "OH",
    "NH1", "NH2",
    "HH", "HH2", "HH11", "HH12", "HH21", "HH22",
)

#: Side-chain atoms of each residue, in library order. Protonation variants are
#: listed separately because they differ by one atom and QresFEP has to restrain
#: the right one.
SIDE_CHAINS: dict[str, list[str]] = {
    "ALA": ["CB", "HB1", "HB2", "HB3"],
    "ARG": ["CB", "HB2", "HB3", "CG", "HG2", "HG3", "CD", "HD2", "HD3", "NE", "HE",
            "CZ", "NH1", "HH11", "HH12", "NH2", "HH21", "HH22"],
    "ARN": ["CB", "HB2", "HB3", "CG", "HG2", "HG3", "CD", "HD2", "HD3", "NE", "HE",
            "CZ", "NH1", "HH11", "NH2", "HH21", "HH22"],
    "ASH": ["CB", "HB2", "HB3", "CG", "OD1", "OD2", "HD1"],
    "ASN": ["CB", "HB2", "HB3", "CG", "OD1", "ND2", "HD21", "HD22"],
    "ASP": ["CB", "HB2", "HB3", "CG", "OD1", "OD2"],
    "CYS": ["CB", "HB2", "HB3", "SG", "HG"],
    "GLH": ["CB", "HB2", "HB3", "CG", "HG2", "HG3", "CD", "OE1", "OE2", "HE1"],
    "GLN": ["CB", "HB2", "HB3", "CG", "HG2", "HG3", "CD", "OE1", "NE2", "HE21", "HE22"],
    "GLU": ["CB", "HB2", "HB3", "CG", "HG2", "HG3", "CD", "OE1", "OE2"],
    "GLY": [],
    "HID": ["CB", "HB2", "HB3", "CG", "ND1", "HD1", "CD2", "HD2", "CE1", "HE1", "NE2"],
    "HIE": ["CB", "HB2", "HB3", "CG", "ND1", "CD2", "HD2", "CE1", "HE1", "NE2", "HE2"],
    "HIP": ["CB", "HB2", "HB3", "CG", "ND1", "HD1", "CD2", "HD2", "CE1", "HE1", "NE2", "HE2"],
    "ILE": ["CB", "HB", "CG1", "HG12", "HG13", "CG2", "HG21", "HG22", "HG23",
            "CD1", "HD11", "HD12", "HD13"],
    "LEU": ["CB", "HB2", "HB3", "CG", "HG", "CD1", "HD11", "HD12", "HD13",
            "CD2", "HD21", "HD22", "HD23"],
    "LYN": ["CB", "HB2", "HB3", "CG", "HG2", "HG3", "CD", "HD2", "HD3", "CE", "HE2", "HE3",
            "NZ", "HZ1", "HZ2"],
    "LYS": ["CB", "HB2", "HB3", "CG", "HG2", "HG3", "CD", "HD2", "HD3", "CE", "HE2", "HE3",
            "NZ", "HZ1", "HZ2", "HZ3"],
    "MET": ["CB", "HB2", "HB3", "CG", "HG2", "HG3", "SD", "CE", "HE1", "HE2", "HE3"],
    "PHE": ["CB", "HB2", "HB3", "CG", "CD1", "HD1", "CD2", "HD2", "CE1", "HE1",
            "CE2", "HE2", "CZ", "HZ"],
    "SER": ["CB", "HB2", "HB3", "OG", "HG"],
    "THR": ["CB", "HB", "CG2", "HG21", "HG22", "HG23", "OG1", "HG1"],
    "TRP": ["CB", "HB2", "HB3", "CG", "CD1", "HD1", "CD2", "NE1", "HE1", "CE2", "CE3",
            "HE3", "CZ2", "HZ2", "CZ3", "HZ3", "CH2", "HH2"],
    "TYR": ["CB", "HB2", "HB3", "CG", "CD1", "HD1", "CD2", "HD2", "CE1", "HE1",
            "CE2", "HE2", "CZ", "OH", "HH"],
    "VAL": ["CB", "HB", "CG1", "HG11", "HG12", "HG13", "CG2", "HG21", "HG22", "HG23"],
}

#: The first side-chain torsion of each residue, keyed by the atoms whose
#: connection to the *other* topology's CB has to be zeroed. Residues are grouped
#: by which chi1 atoms they carry.
_CHI1_ATOMS: dict[str, list[str]] = {
    **{res: ["HB2", "HB3", "CG"] for res in
       ("ASP", "ASH", "GLU", "GLH", "HID", "HIE", "HIP", "ARG", "ARN",
        "LYS", "LYN", "PHE", "LEU", "MET", "TRP", "TYR", "ASN", "GLN")},
    "ILE": ["HB", "CG2", "CG1"],
    "VAL": ["HB", "CG2", "CG1"],
    "THR": ["HB", "CG2", "OG1"],
    "CYS": ["HB2", "HB3", "SG"],
    "SER": ["HB2", "HB3", "OG"],
    "ALA": ["HB2", "HB3", "HB1"],
}


def normalize_pdb_atom_name(line: str) -> str:
    """Rewrite a PDB line's atom name from PyMOL's convention to the Q libraries'.

    PyMOL writes branched hydrogens with the branch index in front (``1HB``,
    ``2HG1``) while the force field libraries put it last (``HB1``, ``HG12``).
    Glycine's two alpha hydrogens are a separate case: PyMOL leaves the first
    unnumbered.

    Args:
        line: A PDB ``ATOM``/``HETATM`` line.

    Returns:
        The line with its atom-name field rewritten; unchanged if it already
        follows the library convention.
    """
    if not line.startswith(("ATOM", "HETATM")) or len(line) < 20:
        return line

    name = line[12:16].strip()
    residue = line[17:20].strip()

    if len(name) > 1 and name[0] in "123" and name[1] == "H":
        name = name[1:] + name[0]
    elif residue == "GLY" and name == "HA":
        name = "HA2"
    else:
        return line

    # Four-character names fill the field; shorter ones are indented by one, which
    # is what separates an atom name from an element symbol in the PDB format.
    field = name if len(name) == 4 else f" {name:<3}"
    return line[:12] + field + line[16:]


def backbone_atoms(from_gly: bool, to_gly: bool) -> tuple[str, ...]:
    """Return the atoms shared by both topologies for this mutation.

    ``HA`` is shared for every mutation except those involving glycine, where the
    alpha carbon changes its hydrogen count and ``HA`` has to be perturbed.

    Args:
        from_gly: True when the wild-type residue is glycine.
        to_gly: True when the mutant residue is glycine.
    """
    if from_gly or to_gly:
        return tuple(atom for atom in BACKBONE_ATOMS if atom != "HA")
    return BACKBONE_ATOMS


def chi1_atoms(residue: str) -> list[str]:
    """Return the chi1 side-chain atoms of a residue, or an empty list for glycine."""
    return list(_CHI1_ATOMS.get(residue, []))


def side_chain_pair(wild_type: str, mutant: str) -> tuple[list[str], list[str]]:
    """Return the side-chain atom names of both topologies of a mutation.

    Mutant names are lower-cased, matching how they are written into the hybrid
    residue library and PDB.

    Args:
        wild_type: Three-letter code of the wild-type residue.
        mutant: Three-letter code of the mutant residue.

    Raises:
        KeyError: If either residue has no side-chain entry.
    """
    return list(SIDE_CHAINS[wild_type]), [atom.lower() for atom in SIDE_CHAINS[mutant]]


def match_heavy_atoms(
    wild_type_atoms: list[str],
    mutant_atoms: list[str],
    hybrid_residue: list[list],
    tolerance: float = 0.5,
) -> list[tuple[str, str]]:
    """Pair wild-type and mutant heavy atoms that occupy the same position.

    Both side chains sit in the pocket simultaneously, so the pairs found here
    become distance restraints that hold the mutant topology on top of the
    wild-type one. Only atoms whose names agree on the first character (``CB`` with
    ``cb``, ``CG`` with ``cg``) are considered, which keeps chemically unrelated
    atoms that happen to overlap from being tied together.

    Args:
        wild_type_atoms: Wild-type side-chain atom names.
        mutant_atoms: Mutant side-chain atom names, lower-cased.
        hybrid_residue: Parsed PDB atom records of the hybrid residue.
        tolerance: Maximum separation in Angstrom for a pair to count.

    Returns:
        ``(wild_type_atom, mutant_atom)`` name pairs.
    """
    heavy_wt = [atom for atom in wild_type_atoms if not atom.startswith("H")]
    heavy_mut = [atom for atom in mutant_atoms if not atom.startswith("h")]

    positions_wt, positions_mut = {}, {}
    for record in hybrid_residue:
        name = record[2]
        coordinates = [float(record[8]), float(record[9]), float(record[10])]
        if name in heavy_wt:
            positions_wt[name] = coordinates
        if name in heavy_mut:
            positions_mut[name] = coordinates

    pairs = []
    for wt_name, wt_xyz in positions_wt.items():
        for mut_name, mut_xyz in positions_mut.items():
            if wt_name[0] != mut_name[0].upper():
                continue
            if euclidian_overlap(wt_xyz, mut_xyz, tolerance):
                pairs.append((wt_name, mut_name))
    return pairs
