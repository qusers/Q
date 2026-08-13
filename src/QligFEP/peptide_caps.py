"""Build the capped reference peptide of a mutation FEP, without PyMOL.

The reference leg reproduces a mutation with no protein around it: the mutated
residue is cut out, given flanking residues, and capped with ACE and NMA. What
that capping produces is an *idealised* peptide -- the flanks and caps are built
in a fixed extended conformation, and the anchoring residue's own carbonyl oxygen
and amide hydrogen are moved to match it. Only the residue itself carries
coordinates that mean anything, and only its side chain is being perturbed.

Because the built part is rigid and always the same, it does not need a molecular
editor to construct: it is one fragment, placed by superimposing it on the
backbone of the residue it attaches to. The templates in
``INPUTS/tripeptide_templates`` hold that fragment for each choice of flanks,
recorded from the PyMOL mutagenesis wizard that used to do this at run time. The
placement reproduces it to about 0.2 A, well inside what the 2.5 ns of
equilibration that follows will move it.

Flanks are ``A`` (alanine), ``G`` (glycine), ``X`` (none) or ``Z`` (the residue's
native neighbours). ``Z`` brings its own flanks from the protein, so it caps them
with the same fragment ``X`` uses.
"""

from __future__ import annotations

import copy
from functools import lru_cache
from pathlib import Path

import numpy as np

from .settings.settings import CONFIGS

#: Backbone atoms the fragment is placed on. Three points fix the frame.
FRAME_ATOMS = ("N", "CA", "C")

#: Residue name marking the template's central residue -- the one whose backbone
#: the fragment is placed against, and whose O and H the capping moves.
CENTRAL_RESIDUE = "XXX"

#: Which template each choice of flanks uses. `Z` keeps the native neighbours, so
#: it needs the caps alone, which is what `X` holds.
_TEMPLATE_FOR_FLANKS = {"A": "A", "G": "G", "X": "X", "Z": "X"}


class PeptideBuildError(Exception):
    """Raised when the reference peptide cannot be built."""


@lru_cache(maxsize=None)
def _load_template(name: str) -> tuple[tuple, ...]:
    """Return a template as ``(residue_number, residue_name, atom_name, x, y, z)``."""
    path = Path(CONFIGS["INPUT_DIR"]) / "tripeptide_templates" / f"flanks_{name}.pdb"
    if not path.exists():
        raise PeptideBuildError(f"Missing reference peptide template: {path}")
    atoms = []
    for line in path.read_text().splitlines():
        if line.startswith("ATOM"):
            atoms.append(
                (
                    int(line[22:26]),
                    line[17:20].strip(),
                    line[12:16].strip(),
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                )
            )
    return tuple(atoms)


def _superimpose(mobile: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return the rotation and translation taking `mobile` onto `target`.

    Ordinary Kabsch superposition. With three points it is exact up to the
    difference in the two triangles, which here is the difference between the
    template residue's backbone geometry and the real one.
    """
    mobile_center, target_center = mobile.mean(axis=0), target.mean(axis=0)
    correlation = (mobile - mobile_center).T @ (target - target_center)
    u, _, vt = np.linalg.svd(correlation)
    # Guard against the reflection SVD is free to return: a mirrored peptide
    # would have the wrong chirality everywhere.
    handedness = np.sign(np.linalg.det(vt.T @ u.T))
    rotation = vt.T @ np.diag([1.0, 1.0, handedness]) @ u.T
    return rotation, target_center - rotation @ mobile_center


def _backbone_frame(residue: list[list], where: str) -> np.ndarray:
    """Return the N/CA/C coordinates of one residue.

    Raises:
        PeptideBuildError: If the residue is missing a backbone atom, which
            leaves nothing to attach the cap to.
    """
    coordinates = {atom[2]: atom[8:11] for atom in residue}
    missing = [name for name in FRAME_ATOMS if name not in coordinates]
    if missing:
        raise PeptideBuildError(
            f"The {where} residue of the reference peptide has no {', '.join(missing)} atom, "
            "so there is no backbone to cap."
        )
    return np.array([coordinates[name] for name in FRAME_ATOMS], dtype=float)


def _template_frame(template: tuple[tuple, ...]) -> np.ndarray:
    coordinates = {
        atom[2]: atom[3:6] for atom in template if atom[1] == CENTRAL_RESIDUE
    }
    return np.array([coordinates[name] for name in FRAME_ATOMS], dtype=float)


def _record(atom_name: str, residue_name: str, residue_number: int, xyz) -> list:
    """Return a PDB atom record in the form :func:`~QligFEP.pdb_utils.pdb_parse_out` wants."""
    return [
        "ATOM  ", 0, atom_name, " ", residue_name, "", residue_number, "",
        float(xyz[0]), float(xyz[1]), float(xyz[2]), 1.00, 0.00, "  ", "  ",
    ]


def build_reference_peptide(residues: list[list[list]], flanks: str) -> list[list]:
    """Cap a stretch of residues into the reference peptide of a mutation FEP.

    Args:
        residues: The residues to cap, N-terminus first, each a list of atom
            records. One residue for ``A``/``G``/``X`` flanks, three for ``Z``.
        flanks: ``A``, ``G``, ``X`` or ``Z``.

    Returns:
        Atom records of the whole capped peptide, in order, with residue numbers
        assigned and serial numbers left for the caller to fill.

    Raises:
        PeptideBuildError: If the flanks are not one of the four, or a residue
            has no backbone to attach to.
    """
    if flanks not in _TEMPLATE_FOR_FLANKS:
        raise PeptideBuildError(
            f"Unknown flanking residue choice {flanks!r}; expected one of "
            f"{', '.join(sorted(_TEMPLATE_FOR_FLANKS))}"
        )
    if not residues:
        raise PeptideBuildError("No residues to cap")

    template = _load_template(_TEMPLATE_FOR_FLANKS[flanks])
    template_frame = _template_frame(template)
    central = next(atom[0] for atom in template if atom[1] == CENTRAL_RESIDUE)

    # The fragment is placed twice: once on the first residue, which carries the
    # N-terminal cap, and once on the last, which carries the C-terminal one.
    # For a single residue the two placements are the same.
    n_rotation, n_shift = _superimpose(template_frame, _backbone_frame(residues[0], "first"))
    c_rotation, c_shift = _superimpose(template_frame, _backbone_frame(residues[-1], "last"))

    first_number = central
    last_number = central + len(residues) - 1

    leading, trailing, moved = [], [], {}
    for number, residue_name, atom_name, x, y, z in template:
        position = np.array([x, y, z])
        if residue_name == CENTRAL_RESIDUE:
            # Capping moves the anchoring residue's carbonyl oxygen and amide
            # hydrogen; its N, CA and C are the frame and stay where they are.
            if atom_name == "H":
                moved[("first", "H")] = n_rotation @ position + n_shift
            elif atom_name == "O":
                moved[("last", "O")] = c_rotation @ position + c_shift
            continue
        if number < central:
            leading.append(
                _record(atom_name, residue_name, number, n_rotation @ position + n_shift)
            )
        else:
            trailing.append(
                _record(
                    atom_name,
                    residue_name,
                    last_number + (number - central),
                    c_rotation @ position + c_shift,
                )
            )

    peptide = list(leading)
    for offset, residue in enumerate(residues):
        number = first_number + offset
        is_first, is_last = offset == 0, offset == len(residues) - 1
        # Proline has no amide hydrogen; nothing should be put on its nitrogen.
        keep_h = not (is_first and residue and residue[0][4].upper().endswith("PRO"))
        for atom in residue:
            atom = copy.deepcopy(atom)
            atom[5], atom[6] = "", number
            if is_first and atom[2] == "H" and keep_h:
                atom[8], atom[9], atom[10] = moved[("first", "H")]
            elif is_last and atom[2] == "O":
                atom[8], atom[9], atom[10] = moved[("last", "O")]
            peptide.append(atom)

        if is_first and keep_h and not any(atom[2] == "H" for atom in residue):
            peptide.append(_record("H", residue[0][4], number, moved[("first", "H")]))

    return peptide + trailing
