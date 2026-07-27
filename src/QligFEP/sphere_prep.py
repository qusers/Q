"""Record of what ``qprep_prot`` did to a protein, serialised as ``prep.json``.

Setups that run downstream of the sphere preparation need facts that the output
PDB files do not carry on their own: where the sphere sits, what net charge it
encloses, which cysteines are bridged, and how the residue numbers of the input
PDB relate to the sequential numbering Q assigns in the topology.

``QresFEP`` depends on all four. A mutation is requested against the *input* PDB
numbering (``LEU39ALA`` on chain A), but every atom index it later writes into a
FEP file is in Q numbering, and the reference (tripeptide) leg has to be charge
matched against the protein sphere.

Q numbers residues positionally: ``qprep`` discards the residue numbers in the
PDB it reads and counts residues from 1 in file order. :func:`build_residue_map`
therefore derives the mapping from the order of the prepared PDB rather than
from its residue-number column, which may carry gaps after out-of-sphere
fragments are dropped.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path

import pandas as pd

from .logger import logger

PREP_JSON = "prep.json"

#: Residue names that never take part in the protein residue numbering Q assigns
#: to the solute -- waters are appended by ``solvate`` after the solute is read.
_SOLVENT_RESIDUES = frozenset({"HOH", "SOL", "WAT"})


@dataclass(frozen=True)
class ResidueMapping:
    """One residue's identity in both numbering schemes.

    Attributes:
        q_number: Residue number in the Q topology (1-based, positional).
        pdb_number: Residue number in the PDB handed to ``qprep_prot``.
        chain: Chain identifier in that PDB (empty string when unset).
        name: Residue name as prepared (may differ from the input, e.g. a
            neutralised ``GLU`` prepared as ``GLH``).
    """

    q_number: int
    pdb_number: int
    chain: str
    name: str


@dataclass
class SpherePrep:
    """Everything ``qprep_prot`` learned about the sphere it built.

    Attributes:
        input_pdb: Name of the PDB the user passed in.
        prepared_pdb: Name of the PDB actually handed to ``qprep``.
        force_field: Force field used for the preparation.
        center: Sphere centre as ``[x, y, z]``.
        radius: Sphere radius in Angstrom.
        total_charge: Net charge enclosed by the sphere, from ``qprep.out``.
        disulfides: Bridged cysteine pairs, as Q residue numbers.
        residues: Per-residue mapping between PDB and Q numbering.
    """

    input_pdb: str
    prepared_pdb: str
    force_field: str
    center: list[float]
    radius: float
    total_charge: int
    disulfides: list[tuple[int, int]] = field(default_factory=list)
    residues: list[ResidueMapping] = field(default_factory=list)

    def q_number(self, pdb_number: int, chain: str | None = None) -> int:
        """Return the Q residue number for a residue of the input PDB.

        Args:
            pdb_number: Residue number in the input PDB.
            chain: Chain identifier. Optional when the residue number is unique
                across chains, required otherwise.

        Returns:
            The residue's number in the Q topology.

        Raises:
            KeyError: If no residue matches, or if several do and no chain was
                given to disambiguate.
        """
        hits = [r for r in self.residues if r.pdb_number == pdb_number]
        if chain is not None:
            hits = [r for r in hits if r.chain == chain]
        if not hits:
            where = f"residue {pdb_number}" + (f" of chain {chain!r}" if chain else "")
            raise KeyError(f"{where} is not in the prepared sphere ({self.prepared_pdb})")
        if len(hits) > 1:
            chains = sorted({r.chain for r in hits})
            raise KeyError(
                f"residue {pdb_number} occurs in chains {chains}; pass a chain to disambiguate"
            )
        return hits[0].q_number

    def residue_name(self, pdb_number: int, chain: str | None = None) -> str:
        """Return the prepared residue name for a residue of the input PDB."""
        q = self.q_number(pdb_number, chain)
        return next(r.name for r in self.residues if r.q_number == q)

    def write(self, directory: Path) -> Path:
        """Write this record to ``prep.json`` in ``directory`` and return its path."""
        path = Path(directory) / PREP_JSON
        payload = asdict(self)
        payload["disulfides"] = [list(pair) for pair in self.disulfides]
        path.write_text(json.dumps(payload, indent=2) + "\n")
        return path

    @classmethod
    def read(cls, directory_or_file: Path) -> SpherePrep:
        """Load a record from a ``prep.json`` file or the directory holding one.

        Raises:
            FileNotFoundError: If no ``prep.json`` is found.
        """
        path = Path(directory_or_file)
        if path.is_dir():
            path = path / PREP_JSON
        if not path.exists():
            raise FileNotFoundError(
                f"{path} not found. Run `qprep_prot` first -- it writes {PREP_JSON} "
                "alongside water.pdb."
            )
        payload = json.loads(path.read_text())
        residues = [ResidueMapping(**entry) for entry in payload.pop("residues", [])]
        disulfides = [tuple(pair) for pair in payload.pop("disulfides", [])]
        return cls(residues=residues, disulfides=disulfides, **payload)


def _solute_residues_in_order(prepared_pdb: Path) -> list[tuple[int, str, str]]:
    """Return ``(residue_number, chain, name)`` per solute residue, in file order.

    Waters are skipped: ``qprep`` adds the sphere's solvent itself, after the solute.
    """
    from .pdb_utils import read_pdb_to_dataframe

    frame = read_pdb_to_dataframe(prepared_pdb)
    solute = frame[~frame["residue_name"].isin(_SOLVENT_RESIDUES)]
    keys = ["residue_seq_number", "residue_name", "chain_id", "insertion_code"]
    # drop_duplicates keeps file order, which is what Q counts.
    residues = solute[keys].drop_duplicates(keep="first")
    return [
        (int(row.residue_seq_number), str(row.chain_id).strip(), str(row.residue_name).strip())
        for row in residues.itertuples(index=False)
    ]


def build_residue_map(
    prepared_pdb: Path,
    original_numbering: dict[int, tuple[int, str, str]] | None = None,
) -> list[ResidueMapping]:
    """Map the residues of a prepared PDB onto the numbering Q will assign.

    Q numbers solute residues from 1 in the order they appear, so the mapping is
    positional.

    Args:
        prepared_pdb: The PDB handed to ``qprep``.
        original_numbering: The mapping returned by
            :func:`~QligFEP.pdb_utils.reindex_pdb_residues`, which is the only
            record of the residue numbers the user's PDB carried. Without it the
            prepared file's own numbers are reported instead.

    Returns:
        One :class:`ResidueMapping` per solute residue, in Q order.
    """
    mapping = []
    for q, (prepared_number, chain, name) in enumerate(
        _solute_residues_in_order(prepared_pdb), start=1
    ):
        pdb_number = prepared_number
        if original_numbering and prepared_number in original_numbering:
            pdb_number, original_chain, _icode = original_numbering[prepared_number]
            chain = original_chain or chain
        mapping.append(
            ResidueMapping(q_number=q, pdb_number=pdb_number, chain=chain, name=name)
        )
    return mapping


def verify_residue_map(residues: list[ResidueMapping], topology_pdb: Path) -> bool:
    """Check a residue map against the topology PDB ``qprep`` wrote.

    ``top_p.pdb`` carries Q's own numbering, so it is the authority on whether
    the positional mapping held. A mismatch is logged rather than raised: the
    rest of the preparation is still usable, and the caller may not need the map.

    Args:
        residues: Mapping produced by :func:`build_residue_map`.
        topology_pdb: Path to ``top_p.pdb``.

    Returns:
        True when every solute residue matches by number and name.
    """
    from .pdb_utils import read_pdb_to_dataframe

    if not Path(topology_pdb).exists():
        logger.debug(f"{topology_pdb} not found; skipping residue map verification")
        return False

    frame = read_pdb_to_dataframe(topology_pdb)
    solute = frame[~frame["residue_name"].isin(_SOLVENT_RESIDUES)]
    topology = (
        solute[["residue_seq_number", "residue_name"]]
        .drop_duplicates(keep="first")
        .itertuples(index=False)
    )
    expected = {r.q_number: r.name for r in residues}

    mismatches = []
    seen = 0
    for row in topology:
        seen += 1
        q = int(row.residue_seq_number)
        name = str(row.residue_name).strip()
        if expected.get(q) != name:
            mismatches.append(f"Q{q}: topology says {name}, map says {expected.get(q)}")

    if seen != len(residues):
        mismatches.append(f"topology has {seen} solute residues, map has {len(residues)}")

    if mismatches:
        logger.warning(
            "Residue map disagrees with the topology Q wrote -- residue numbering "
            "downstream may be wrong:\n  " + "\n  ".join(mismatches[:10])
        )
        return False
    logger.debug(f"Residue map verified against {topology_pdb} ({seen} residues)")
    return True


def parse_disulfide_pairs(
    cysbond_lines: str, prepared_to_q: dict[int, int]
) -> list[tuple[int, int]]:
    """Convert ``addbond`` lines into Q residue-number pairs.

    ``cysbonds_for_qprep`` reports residue numbers as written in the prepared PDB.
    Those coincide with Q numbering only when nothing was dropped from that file,
    so they are translated here.

    Args:
        cysbond_lines: The ``[!]addbond N:SG M:SG y`` block, one bond per line.
        prepared_to_q: Prepared-PDB residue number to Q residue number.

    Returns:
        Bridged pairs as Q residue numbers, in the order they were reported.
    """
    pairs = []
    for line in cysbond_lines.splitlines():
        tokens = line.replace("!", "").split()
        if not tokens or tokens[0] != "addbond":
            continue
        try:
            numbers = [int(token.split(":")[0]) for token in tokens[1:3]]
        except ValueError:
            logger.debug(f"Skipping unparsable cysbond line: {line!r}")
            continue
        translated = []
        for number in numbers:
            if number not in prepared_to_q:
                logger.warning(
                    f"Disulfide partner {number} is not a solute residue of the "
                    "prepared PDB; recording it unchanged"
                )
                translated.append(number)
            else:
                translated.append(prepared_to_q[number])
        pairs.append((translated[0], translated[1]))
    return pairs


def collect(
    input_pdb: Path,
    prepared_pdb: Path,
    force_field: str,
    center: list[float],
    radius: float,
    total_charge: int,
    cysbond_lines: str,
    topology_pdb: Path | None = None,
    original_numbering: dict[int, tuple[int, str, str]] | None = None,
) -> SpherePrep:
    """Assemble a :class:`SpherePrep` from a finished preparation run.

    Args:
        input_pdb: The PDB the user passed to ``qprep_prot``.
        prepared_pdb: The PDB actually handed to ``qprep``.
        force_field: Force field name used.
        center: Sphere centre as ``[x, y, z]``.
        radius: Sphere radius in Angstrom.
        total_charge: Net charge in the sphere, parsed from ``qprep.out``.
        cysbond_lines: The ``addbond`` block produced during preparation.
        topology_pdb: ``top_p.pdb``, used to verify the residue map when present.
        original_numbering: Residue renumbering recorded during preparation, so
            the map can report the numbers the user's PDB carried.
    """
    residues = build_residue_map(Path(prepared_pdb), original_numbering)
    if topology_pdb is not None:
        verify_residue_map(residues, Path(topology_pdb))

    prepared_to_q = {
        prepared_number: q
        for q, (prepared_number, _chain, _name) in enumerate(
            _solute_residues_in_order(Path(prepared_pdb)), start=1
        )
    }

    return SpherePrep(
        input_pdb=Path(input_pdb).name,
        prepared_pdb=Path(prepared_pdb).name,
        force_field=force_field,
        center=[float(c) for c in center],
        radius=float(radius),
        total_charge=int(total_charge),
        disulfides=parse_disulfide_pairs(cysbond_lines, prepared_to_q),
        residues=residues,
    )


def residue_map_frame(prep: SpherePrep) -> pd.DataFrame:
    """Return the residue mapping as a DataFrame, for inspection in notebooks."""
    return pd.DataFrame(asdict(entry) for entry in prep.residues)
