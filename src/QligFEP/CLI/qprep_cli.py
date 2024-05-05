"""Module containing the command line interface for the qprep fortran program."""

import argparse
from io import StringIO
from pathlib import Path
from typing import Optional

from ..chemIO import MoleculeIO
from ..IO import run_command
from ..logger import logger, setup_logger
from ..pdb_utils import (
    disulfide_search,
    nest_pdb,
    read_pdb_to_dataframe,
    write_dataframe_to_pdb,
)
from ..settings.settings import CONFIGS

# NOTE: cysbonds will have \n after each bond -> `maketop MKC_p` is in a different line
qprep_inp_content = """rl {ff_lib_path}
rprm {ff_prm_path}
! TO DO Change if protein system is used
rp {pdb_file_path}
! set solute_density 0.05794
! NOTE, this is now large for water system, change for protein system
set solvent_pack {solvent_pack}
boundary 1 {cog} 25.0
solvate {cog} {sphereradius} {qprep_type}
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
        "-pdb",
        "--pdb_file",
        dest="pdb_file",
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
        help="Center of geometry for the protein. The format is 'x y z', where all numbers contain 3 decimal cases.",
        required=True,
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
        default=2.3,
        help=(
            "Parameter to qprep.inp `set solvent_pack`. According to Q's manual, this value "
            "corresponds to the minimum distance between solute and solvent heavy atoms when "
            "adding solvent (e.g.: HOH) and defaults to 2.4. In QligFEP we use a value of 2.3, "
            "but while generating a water.pdb that is suitable for input in `qligfep` and "
            "`setupFEP` (with `-t protein`), this can create water molecules that are too "
            "close to the ligands, causing calculations to crash. Please increase this value "
            "when you call this program with the `-t protein` argument. Defaults to 2.3."
        ),
        type=float,
    )
    parser.add_argument(
        "-il",
        "--include_ligands",
        dest="include_ligands",
        default=None,
        help=(
            "Optional argument to include ligands while preparing the water sphere with your `-pdb` file. "
            "This is useful if you want to include ligands to `mask` your water sphere, leaving a blank space "
            "to be filled in when you call `qligfep` or `setupFEP`. The input for this argument can be either "
            "pdb or sdf files. If you include a sdf file that contains multiple ligands, the program will "
            "include the coordinates for all of them."
        ),
        type=float,
    )
    parser.add_argument(
        "-t",
        "--qprep_type",
        dest="qprep_type",
        default="water",
        choices=["protein", "water"],
        help=(
            "Type of system to be solvated. If water, the cysbond argument will be ignored. "
            "But qprep will require a water.pdb file in the same directory. Defaults to water."
        ),
        type=str,
    )
    parser.add_argument(
        "-log",
        "--log-level",
        dest="log_level",
        default="INFO",
        choices=["INFO", "DEBUG", "WARNING", "ERROR", "CRITICAL"],
        help=(
            "Set the logging level. Defaults to INFO. "
            "Choose between: INFO, DEBUG, WARNING, ERROR, CRITICAL."
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

    # write temporary pdb_file depending on args.include_ligands parameter
    if args.include_ligands is not None:
        lig_path = Path(args.include_ligands)
        assert lig_path.exists(), f"File {lig_path} does not exist."
        if lig_path.suffix == ".pdb":
            raise ValueError(".pdb file format not yet implemented for ligands. Please use .sdf format.")
        molio = MoleculeIO(args.include_ligands)
        input_pdb_df = read_pdb_to_dataframe(args.pdb_file)
        offset = input_pdb_df["atom_serial_number"].astype(int).max()
        output = StringIO()  # just to prevent writing to disk
        ligands_df = molio.write_to_single_pdb(output)
        del output
        ligands_df = ligands_df.assign(
            atom_serial_number=lambda x: x["atom_serial_number"].astype(int) + offset,
        )
        # make an output suffix that will be always unique, including time information:
        write_dataframe_to_pdb(ligands_df, cwd / f"{args.pdb_file}_withLigands.pdb")
        pdb_file = str(cwd / f"{args.pdb_file}_withLigands.pdb")
    else:
        pdb_file = str(cwd / args.pdb_file)

    ff_lib_path = str(Path(CONFIGS["FF_DIR"]) / f"{args.FF}.lib")
    ff_prm_path = str(Path(CONFIGS["FF_DIR"]) / f"{args.FF}.prm")

    if args.qprep_type == "protein":
        annotation_type = "4 water.pdb"
    if args.qprep_type == "water":
        annotation_type = "1 HOH"

    qprep_filename = f"qprep_{args.qprep_type}"
    qprep_inp_path = cwd / f"{qprep_filename}.inp"

    cysbonds = "" if args.qprep_type == "water" else args.cysbond
    if cysbonds == "auto":
        with open(pdb_file) as f:
            pdb_lines = f.readlines()
            npdb = nest_pdb(pdb_lines)
            npdb, cysbonds = disulfide_search(npdb)
            del npdb
        cysbonds = "".join([f"!addbond {atomN[0]} {atomN[1]} y\n" for atomN in cysbonds])
    elif cysbonds != "":
        cysbonds = cysbonds.split(",")
        cysbonds = "".join([f"!addbond {b.split('_')[0]} {b.split('_')[1]} y\n" for b in cysbonds])
    elif cysbonds == "":  # TODO: do this in a smarter way...
        pass
    else:
        raise ValueError("Invalid cysbond input. Please check the input format.")

    # format the cysbonds = addbond at1 at2 y
    if args is not None:
        param_dict = {
            "pdb_file_path": pdb_file,
            "cog": cog,
            "ff_lib_path": ff_lib_path,
            "ff_prm_path": ff_prm_path,
            "sphereradius": sphereradius,
            "cysbond": cysbonds,
            "qprep_type": annotation_type,
            "solvent_pack": formatted_solvent_pack,
        }
    else:
        param_dict = {**kwargs}

    # write qprep.inp with the formatted qprep_inp_content using Path
    if qprep_inp_path.exists():
        logger.warning("qprep.inp already exists!! Skipping qprep.inp file generation...")
    else:
        with qprep_inp_path.open("w") as qprep_inp_f:
            qprep_inp_f.write(qprep_inp_content.format(**param_dict))

    options = f" < {qprep_filename}.inp > {qprep_filename}.out"
    logger.debug(f"Running command {qprep_path} {options}")
    run_command(qprep_path, options, string=True)
    logger.info(f"qprep run finished. Check the output `{qprep_filename}.out` for more information.")
    # if the -t argument is water, we should read the pdb files and write a water.pdb file with only water molecules
    if args.qprep_type == "water":
        waterfile = Path(cwd / "water.pdb")
        if waterfile.exists():
            logger.warning("water.pdb already exists!! Skipping water.pdb file generation...")
        else:
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


def main_exe():
    args = parse_arguments()
    main(args)


if __name__ == "__main__":
    main_exe()
