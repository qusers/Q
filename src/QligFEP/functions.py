import numpy as np
from openff.units import unit

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


def convert_value(value, original_type, final_type, temperature=298.0, out_unit=None):
    """Adapted from: https://github.com/openforcefield/protein-ligand-benchmark/blob/main/plbenchmark/utils.py
    Converts an experimental value into another derived quantity with specified unit. Please note that
    ic50 and ki values should be in nanomolar.

    Args:
        value: float, numerical value
        original_type: string, code for the original observable. Can be `dg`, `ki`, `ic50`, `pic50`
        final_type: string, code for the desired derived quantity. Can be `dg`, `ki`, `ic50`, `pic50`
        temperature: float, temperature in kelvin. Defaults to 298.0 (default temperature in QligFEP)
        out_unit: unit of type :py:class:`pint`, output unit of final_type, needs to fit to the requested final_type. Defaults to None.

    Raises:
        NotImplementedError: whenever the conversion is not possible.

    Returns:
        float with desired unit
    """
    # define default units
    if out_unit is None:
        if final_type == "dg":
            out_unit = unit("kilocalories / mole")
        elif final_type == "ki" or final_type == "ic50":
            out_unit = unit("nanomolar")
        elif final_type == "pic50":
            out_unit = unit("")

    if original_type == "dg":
        if final_type == "dg":
            return (value.to(out_unit)).magnitude
        elif final_type == "ki" or final_type == "ic50":
            result = np.exp(-value / (unit.molar_gas_constant * temperature * unit.kelvin)) * unit.molar
            return (result.to(out_unit)).magnitude
        elif final_type == "pic50":
            result = value / (unit.molar_gas_constant * temperature * unit.kelvin) / np.log(10)
            return (result.to(out_unit)).magnitude
        else:
            raise NotImplementedError(
                f"Conversion to observable {final_type} not possible. "
                f"Observable must be any of: dg, ki, ic50 or pic50."
            )
    elif original_type == "ki":
        if final_type == "dg":
            if value < 1e-15 * unit("molar"):
                return (0.0 * out_unit).magnitude
            else:
                result = unit.molar_gas_constant * temperature * unit.kelvin * np.log(value / unit.molar)
                return (result.to(out_unit).round(2)).magnitude
        elif final_type == "ki" or final_type == "ic50":
            return (value.to(out_unit)).magnitude
        elif final_type == "pic50":
            if value < 1e-15 * unit("molar"):
                return (-1e15 * out_unit).magnitude
            else:
                result = -np.log(value / unit.molar) / np.log(10)
                return (result).magnitude
        else:
            raise NotImplementedError(
                f"Conversion to observable {final_type} not possible. "
                f"Observable must be any of: dg, ki, ic50 or pic50."
            )
    elif original_type == "ic50":
        if final_type == "dg":
            if value < 1e-15 * unit("molar"):
                return (0.0 * out_unit).magnitude
            else:
                result = (
                    unit.molar_gas_constant
                    * temperature
                    * unit.kelvin
                    * np.log(value.to("molar") / unit.molar)
                )
                return (result.to(out_unit).round(2)).magnitude
        elif final_type == "ki" or final_type == "ic50":
            return (value.to(out_unit)).magnitude
        elif final_type == "pic50":
            if value.to("molar") < 1e-15 * unit("molar"):
                return (-1e15 * out_unit).magnitude
            else:
                result = -np.log(value / unit.molar) / np.log(10)
                return (result).magnitude
        else:
            raise NotImplementedError(
                f"Conversion to observable {final_type} not possible. "
                f"Observable must be any of: dg, ki, ic50 or pic50."
            )
    elif original_type == "pic50":
        if final_type == "dg":
            result = -1 * unit.molar_gas_constant * temperature * unit.kelvin * value * np.log(10)
            return (result.to(out_unit).round(2)).magnitude
        elif final_type == "ki" or final_type == "ic50":
            result = 10 ** (-value) * unit("molar")
            return (result.to(out_unit)).magnitude
        elif final_type == "pic50":
            return (value.to(out_unit)).magnitude
        else:
            raise NotImplementedError(
                f"Conversion to observable {final_type} not possible. "
                f"Observable must be any of: dg, ki, ic50 or pic50."
            )
