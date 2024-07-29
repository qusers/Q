"""Module containing the command line interface for the center of geometry calculation."""

import argparse
import re
from pathlib import Path

from loguru import logger


class MolecularCOG:
    """Class implemented to calculate the center of geometry of a molecule. It supports
    PDB and SDF file formats. The center of geometry is calculated as the average of all
    the coordinates of the atoms in the molecule.

    In case of SDF files, the class will calculate the center of geometry for each ligand
    in the file and return the average of all the centers.
    """

    def __init__(self, filepath):
        self.filepath = Path(filepath).absolute()
        self.filetype = self.filepath.suffix

    def COG(self, include="ATOM,HETATM"):
        if self.filetype == ".pdb":
            return self._cog_pdb(include)
        elif self.filetype == ".sdf":
            return self._cog_sdf()
        else:
            raise ValueError(f"Unsupported file type: {self.filetype}")

    def _cog_pdb(self, include):
        """
        Calculates center of geometry for PDB files.
        """
        center = [0.0, 0.0, 0.0]
        include = tuple(include.split(","))
        count = 0
        with self.filepath.open() as file:
            for line in file:
                if line.startswith(include):
                    count += 1
                    center[0] += float(line[30:38])
                    center[1] += float(line[38:46])
                    center[2] += float(line[46:54])
        if count > 0:
            center = [x / count for x in center]
        return f"[{round(center[0], 3)} {round(center[1], 3)} {round(center[2], 3)}]"

    def _cog_sdf(self):
        coordinate_regex = re.compile(r"^\s*(-?\d+\.\d{4})\s*(-?\d+\.\d{4})\s*(-?\d+\.\d{4})\s+\S")
        centers = []
        coordinates = []
        with self.filepath.open() as file:
            for line in file:
                if coordinate_regex.match(line):
                    coords = [float(coordinate_regex.match(line).group(i)) for i in range(1, 4)]
                    coordinates.append(coords)
                elif line.startswith("$$$$") and coordinates:
                    centers.append(self._calculate_center(coordinates))
                    coordinates = []

        if coordinates:  # For the last set of coordinates if the file does not end with '$$$$'
            centers.append(self._calculate_center(coordinates))

        if len(centers) > 1:
            logger.warning("Calculating for all ligands in the file.")
            for i, center in enumerate(centers):
                logger.debug(f"Ligand {i+1} center: {center}")

        overall_center = [sum(x) / len(centers) for x in zip(*centers)]
        return f"[{round(overall_center[0], 3)} {round(overall_center[1], 3)} {round(overall_center[2], 3)}]"

    def _calculate_center(self, coordinates):
        count = len(coordinates)
        center = [sum(x[i] for x in coordinates) / count for i in range(3)]
        return [round(x, 3) for x in center]

    def __call__(self, include="ATOM,HETATM"):
        """When the object is called, it will calculate the center of geometry of the molecule.

        Args:
            include: if working with PDB files, which atom types to include. Defaults to 'ATOM,HETATM'.

        Returns:
            center of geometry as a string (e.g. '[x y z]')
        """
        return self.COG(include)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Calculate the center of geometry for molecular structures in PDB or SDF format."
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="input",
        required=True,
        help="Path to the molecular structure file (PDB or SDF format).",
        type=str,
    )
    parser.add_argument(
        "-incl",
        "--include",
        dest="include",
        default="ATOM,HETATM",
        help="For PDB files, specify the atom types to include (e.g., 'ATOM,HETATM').",
        type=str,
    )
    return parser.parse_args()


def main(args: argparse.Namespace) -> None:
    cog = MolecularCOG(args.input)
    center = cog(args.include)
    logger.info(f"Center of geometry: {center}")


def main_exe():
    args = parse_arguments()
    main(args)


if __name__ == "__main__":
    main_exe()
