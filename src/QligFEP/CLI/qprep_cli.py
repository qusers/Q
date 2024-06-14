"""Module containing the command line interface for the qprep fortran program."""

import argparse
from pathlib import Path
from typing import Optional

import numpy as np

from ..IO import run_command
from ..logger import logger, setup_logger
from ..pdb_utils import pdb_HOH_nn, read_pdb_to_dataframe, write_dataframe_to_pdb
from ..settings.settings import CONFIGS
from .utils import handle_cysbonds

# NOTE: cysbonds will have \n after each bond -> `maketop MKC_p` is in a different line
qprep_inp_content = """rl {ff_lib_path}
rprm {ff_prm_path}
! TO DO Change if protein system is used
rp {pdb_file_path}
! set solute_density 0.05794
! NOTE, this is now large for water system, change for protein system
set solvent_pack {solvent_pack}
boundary 1 {cog} {sphereradius}
solvate {cog} {sphereradius} 1 HOH
{cysbond}maketop MKC_p
writetop dualtop.top
wp top_p.pdb y
rt dualtop.top
mask none
mask not excluded
wp complexnotexcluded.pdb y
q
"""


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Runs qprep to generate the water.pdb, the dualtop.top, and the complexnotexcluded.pdb files for "
            "an input protein.pdb file. In case of failure, please inspect the qprep.log file. "
            "Depending on the forcefield you're using, you might need some additional preparation steps "
            "(e.g.: running `pdb2amber` if you're using that forcefield)."
        )
    )
    parser.add_argument(
        "-i",
        "--input_pdb_file",
        dest="input_pdb_file",
        required=True,
        help="input protein to run qprep",
    )
    parser.add_argument(
        "-FF",
        "--forcefield",
        dest="FF",
        default="AMBER14sb",
        choices=["OPLS2015", "AMBER14sb", "CHARMM36"],
        help="Forcefield to be used. Defaults to Amber14sb.",
    )
    parser.add_argument(
        "-cog",
        "--center_of_geometry",
        dest="cog",
        help=(
            "Center of geometry for the protein. The format is 'x y z', where all numbers "
            "contain 3 decimal cases. This center of geometry can be obtained using the `qcog ` "
            "command, but if you include ligands using the `-lig` option, this COG will be "
            "automatically calculated. If you want to calculate the COG manually, you can use "
            "this option. Defaults to None."
        ),
        required=False,
        default=None,
        nargs=3,
        type=str,
    )
    parser.add_argument(
        "-r",
        "--sphereradius",
        dest="sphereradius",
        required=False,
        default=25,
        help="Size of the simulation sphere. If float, only one decimal case will be used. Defaults to 25.",
        type=float,
    )
    parser.add_argument(
        "-b",
        "--cysbond",
        dest="cysbond",
        default="auto",
        help=(
            "Add cystein bonds. Input should be formatted with the atom numbers"
            "(participating in the Cys bond) connected by `_` and with different bonds "
            "separated by `,` as in: `atom1_atom2,atom3_atom4`. Defaults to `auto`, where "
            "cystein bonds will be automatically detected within distance of 1.8 to 2.2 A."
        ),
        type=str,
    )
    parser.add_argument(
        "-sp",
        "--solvent_pack",
        dest="solvent_pack",
        default=3.0,
        help=(
            "Parameter to qprep.inp `set solvent_pack`. According to Q's manual, this value "
            "corresponds to the minimum distance between solute and solvent heavy atoms when "
            "adding solvent (e.g.: HOH) and defaults to 2.4. In QligFEP we use a value of 2.3 "
            "on the water leg. For the protein leg, however, this value often creates water "
            "molecules that are too close to the structure, causing calculations to crash. "
            "close to the ligands, causing calculations to crash. Defaults to 3.0."
        ),
        type=float,
    )
    parser.add_argument(
        "-log",
        "--log-level",
        dest="log_level",
        default="info",
        choices=["info", "debug", "warning", "error", "critical"],
        help=(
            "Set the logging level. Defaults to info. "
            "Choose between: info, debug, warning, error, critical."
        ),
        type=str,
    )
    return parser.parse_args()


def main(args: Optional[argparse.Namespace] = None, **kwargs) -> None:
    """Either runs the qprep program with the given arguments via **kwargs or parses
    the arguments and runs the program.

    Args:
        args: argparse Namespace containing the arguments for the qprep program. Defaults to None
    """
    cwd = Path.cwd()
    if args.log_level != "INFO":
        setup_logger(args.log_level.upper())
    qprep_path = CONFIGS["QPREP"]
    logger.debug(f"Running qprep from path: {qprep_path}")
    sphereradius = f"{args.sphereradius:.1f}"
    formatted_solvent_pack = f"{args.solvent_pack:.1f}"

    cog = " ".join(args.cog)
    logger.debug(f"COG is {cog}")

    pdb_file = str(cwd / args.input_pdb_file)

    ff_lib_path = str(Path(CONFIGS["FF_DIR"]) / f"{args.FF}.lib")
    ff_prm_path = str(Path(CONFIGS["FF_DIR"]) / f"{args.FF}.prm")

    qprep_inp_path = cwd / "qprep.inp"

    cysbonds = handle_cysbonds(args.cysbond, pdb_file, comment_out=True)

    if args is not None:
        param_dict = {
            "pdb_file_path": pdb_file,
            "cog": cog,
            "ff_lib_path": ff_lib_path,
            "ff_prm_path": ff_prm_path,
            "sphereradius": sphereradius,
            "cysbond": cysbonds,
            "solvent_pack": formatted_solvent_pack,
        }
    else:
        param_dict = {**kwargs}

    if qprep_inp_path.exists():
        logger.warning("qprep.inp already exists!! Overwriting...")
    with qprep_inp_path.open("w") as qprep_inp_f:
        qprep_inp_f.write(qprep_inp_content.format(**param_dict))

    options = " < qprep.inp > qprep.out"
    logger.debug(f"Running command {qprep_path} {options}")
    run_command(qprep_path, options, string=True)
    logger.info("qprep run finished. Check the output `qprep.out` for more information.")

    waterfile = Path(cwd / "water.pdb")
    # Write water file and deal with possible errors
    if not Path("complexnotexcluded.pdb").exists():
        logger.error(
            "`complexnotexcluded.pdb` file not found. This is as sign qprep didn't "
            "run correctly. Check the outoput in your console and try again..."
        )
        logger.info(
            "If your console contains something like `libgfortran.so.5: cannot "
            "open shared object file: No such file or directory`, you might need to load "
            "some module in your HPC system that you used to compile Q."
        )
        raise FileNotFoundError("complexnotexcluded.pdb file not found. Something went wrong")
    with open("complexnotexcluded.pdb") as f:
        lines = f.readlines()
    with open("water.pdb", "w") as f:
        for line in lines:
            if line.startswith("ATOM") and line[17:20] == "HOH":
                f.write(line)
        logger.info("water.pdb file created.")

    # Now that the water file is created, we remove water molecules outside the sphere radius
    cog = [float(i) for i in cog.split()]
    water_df = read_pdb_to_dataframe(waterfile)
    oxygen_subset = water_df.query('atom_name == "O"')
    euclidean_distances = oxygen_subset[["x", "y", "z"]].sub(cog).pow(2).sum(1).apply(np.sqrt)
    outside = np.where(euclidean_distances > args.sphereradius)[0]
    outside_HOH_residues = oxygen_subset.iloc[outside].residue_seq_number.unique()  # noqa: F841
    if outside.shape[0] > 0:
        logger.warning(f"Found {outside.shape[0]} water molecules outside the sphere radius.")
        logger.warning("Removing these water molecules from the water.pdb file.")
        todrop_idxs = water_df.query("residue_seq_number in @outside_HOH_residues").index
        water_df.drop(index=todrop_idxs, inplace=True)
        new_distances = (
            water_df.query('atom_name == "O"')[["x", "y", "z"]].sub(cog).pow(2).sum(1).apply(np.sqrt)
        )
        logger.debug(f"Final highest distance is {new_distances.max():.2f} A")
        write_dataframe_to_pdb(water_df, waterfile)
    else:
        logger.info("All water molecules are inside the sphere radius.")
        logger.debug(f"Final highest distance to COG is {euclidean_distances.max():.2f} A")


def main_exe():
    args = parse_arguments()
    main(args)


if __name__ == "__main__":
    main_exe()
