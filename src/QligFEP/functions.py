import numpy as np

from .pdb_utils import pdb_parse_in


def sigmoid(t, k):
    return k * t / (k - t + 1.0)


def COG(pdbfile, include="ATOM,HETATM"):
    """
    Calculates center of geometry of a protein and/or ligand structure.
    Returns:
        center (list): List of float coordinates [x,y,z] that represent the
        center of geometry (precision 3).
    """

    center = [None, None, None]
    include = tuple(include.split(","))

    with open(pdbfile) as pdb:

        # extract coordinates [ [x1,y1,z1], [x2,y2,z2], ... ]
        coordinates = []
        for line in pdb:
            if line.startswith(include):
                coordinates.append(
                    [
                        float(line[30:38]),  # x_coord
                        float(line[38:46]),  # y_coord
                        float(line[46:54]),  # z_coord
                    ]
                )

        # calculate center of geometry
        center = [
            sum([coordinates[i][j] / (len(coordinates)) for i in range(len(coordinates))]) for j in range(3)
        ]
        center = [round(center[i], 3) for i in range(3)]
    return center


def euclidian_overlap(coord1, coord2, distance):
    """
    Calculates whether two points in space overlap within a certain distance
    Returns:
        Boolean
    """
    return (coord1[0] - coord2[0]) ** 2 + (coord1[1] - coord2[1]) ** 2 + (
        coord1[2] - coord2[2]
    ) ** 2 < distance**2


def overlapping_pairs(pdbfile, reslist, include=("ATOM", "HETATM")):
    """
    Calculates whether input pdb has overlaying atoms, based on provided residue names
    Returns:
        dictionary of overlaying atoms based on atom number
    """
    coordinates = []
    overlapping_atoms = []
    atomlist = []
    # Parse the input pdbfile
    with open(pdbfile) as infile:
        for line in infile:
            if line.startswith(include):
                line_parse = pdb_parse_in(line)
                if line_parse[4] in reslist:
                    coordinates.append(
                        [line_parse[1], line_parse[8], line_parse[9], line_parse[10], line_parse[13]]
                    )
    for at1 in coordinates:
        for at2 in coordinates:
            condit1 = at1[0] != at2[0]  # don't restrain atoms that aren't the same
            condit2 = ((at1[1] - at2[1]) ** 2 + (at1[2] - at2[2]) ** 2 + (at1[3] - at2[3]) ** 2) < 0.8  # dist
            condit3 = at1[4] == at2[4] and at1[4].strip() != "H"  # don't restraint hydrogens
            if all([condit1, condit2, condit3]):
                overlapping_atoms.append([at1[0], at2[0]])

    total = len(overlapping_atoms)
    for i in range(0, (int(total / 2))):
        atomlist.append(overlapping_atoms[i])

    return atomlist


def resfep_lambda_ladder(windows: int, sampling: str) -> list[str]:
    """Return the lambda values of one QresFEP stage, from 1.000 down to 0.000.

    Distinct from :func:`QligFEP.IO.get_lambdas`, which the ligand protocol uses:
    that one is built on a rational sigmoid and returns ``windows + 1`` values,
    whereas a mutation stage needs exactly ``windows``, spaced as published.

    The ``exponential`` and ``reverse_exponential`` ladders are mirror images.
    Pairing them across the two stages puts the closely-spaced windows at the end
    of each stage, which is where that stage's own topology is being switched and
    the free energy changes fastest.

    Args:
        windows: Number of lambda values to produce.
        sampling: ``linear``, ``sigmoidal``, ``exponential`` or ``reverse_exponential``.

    Returns:
        Lambda values as strings with three decimals, descending from 1.000.

    Raises:
        ValueError: If `sampling` is not a known scheme, or `windows` is below 2.
    """
    if windows < 2:
        raise ValueError(f"A FEP stage needs at least 2 windows, got {windows}")

    fraction = np.linspace(0.0, 1.0, windows)

    if sampling == "linear":
        values = fraction
    elif sampling == "sigmoidal":
        curve = 1 / (1 + np.exp(np.linspace(-6, 6, windows)))
        values = 1 - (curve - curve.min()) / (curve.max() - curve.min())
    elif sampling in ("exponential", "reverse_exponential"):
        k = -3.0 if sampling == "exponential" else 3.0
        values = (np.exp(k * fraction) - 1) / (np.exp(k) - 1)
    else:
        raise ValueError(f"Unknown lambda sampling scheme: {sampling!r}")

    # Rounding error at the endpoints would otherwise render as "-0.000", which
    # ends up in file names and in Q's [lambdas] section. np.clip preserves the
    # sign bit of negative zero, so normalize values that round to zero explicitly.
    clipped = np.clip(values[::-1], 0.0, 1.0)
    return ["0.000" if abs(value) < 0.0005 else f"{value:.3f}" for value in clipped]


#: Reciprocal molecular volumes used to set qprep's `solute_density`, in A^-3.
PROTEIN_VOLUME = 0.05794
LIPID_VOLUME = 0.03431  # from octane


def sphere_solute_density(pdb_file, center, radius) -> float:
    """Return the solute density qprep should pack solvent against.

    Lipids are looser than protein, so a membrane system inside the sphere needs a
    density between the two, weighted by how many heavy atoms of each are enclosed.
    Waters and hydrogens are ignored -- qprep's figure counts heavy solute atoms.

    Args:
        pdb_file: Path to the PDB holding the solute.
        center: Sphere center as ``[x, y, z]``.
        radius: Sphere radius in Angstrom.

    Returns:
        Density in A^-3; ``PROTEIN_VOLUME`` when the sphere holds no lipid.
    """
    center = np.array([float(c) for c in center], dtype=float)
    radius = float(radius)

    n_solute = 0
    n_lipid = 0
    with open(pdb_file) as infile:
        for line in infile:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            residue_name = line[17:20].strip()
            atom_name = line[12:16].strip()
            if residue_name == "HOH" or not atom_name or atom_name[0] == "H":
                continue
            try:
                position = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            except ValueError:
                continue
            if np.linalg.norm(position - center) <= radius:
                n_solute += 1
                if residue_name == "POP" and atom_name[0] == "C":
                    n_lipid += 1

    if not n_solute:
        return PROTEIN_VOLUME

    n_protein = n_solute - n_lipid
    return (n_protein * PROTEIN_VOLUME + n_lipid * LIPID_VOLUME) / n_solute


def kT(T):
    k = 0.0019872041  # kcal/(mol.K)
    kT = k * T
    kT = f"{kT:.3f}"
    return kT


def avg_sem(array):
    FEP_sum = array.sum(axis=0)
    dG = np.nanmean(FEP_sum)
    cnt = len(FEP_sum)
    sem = np.nanstd(FEP_sum, ddof=1) / np.sqrt(cnt)
    return [dG, sem]
