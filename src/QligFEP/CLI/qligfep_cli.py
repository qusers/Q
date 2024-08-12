"""Module containing the QligFEP command line interface."""

import argparse
from typing import Optional

from ..logger import logger
from ..qligfep import QligFEP
from .parser_base import parse_arguments


def main(args: Optional[argparse.Namespace] = None, **kwargs) -> None:
    """Main function for qligfep_cli.py. Takes arguments from argparse and passes them
    to QligFEP class. If no arguments are given, the function will use the keyword arguments
    that are passed to it.

    Args:
        args: argparse.Namespace object containing all the arguments.
        kwargs: keyword arguments that will be passed to QligFEP class.
    """
    if args is not None:
        param_dict = {
            "lig1": args.lig1,
            "lig2": args.lig2,
            "FF": args.FF,
            "system": args.system,
            "cluster": args.cluster,
            "sphereradius": args.sphereradius,
            "cysbond": args.cysbond,  # TODO: add automatic detection of cysbond as in qprep_cli.py
            "start": args.start,
            "temperature": args.temperature,
            "replicates": args.replicates,
            "sampling": args.sampling,
            "timestep": args.timestep,
            "softcore": args.softcore,
            "to_clean": args.to_clean,
            "random_state": args.random_state,
        }
    else:
        param_dict = {}
    param_dict.update(kwargs)
    run = QligFEP(**param_dict)
    run.set_timestep()

    writedir = run.makedir()
    inputdir = writedir + "/inputfiles"
    a = run.read_files()
    changes_for_libfiles = a[0][1]
    changes_for_prmfiles = a[0][1]
    change_charges = a[1][0]
    change_vdw = a[1][1]
    lig_size1, lig_size2 = a[2][0], a[2][1]

    # Write the merged files
    logger.debug("Writing changes on files")
    run.change_lib(changes_for_libfiles, inputdir)
    logger.debug("Changing parameters")
    FEP_vdw = run.change_prm(changes_for_prmfiles, inputdir)
    logger.debug("Writing PDB files")
    run.merge_pdbs(inputdir)

    run.write_water_pdb(inputdir)

    logger.debug("Getting the lambdas")
    lambdas = run.get_lambdas(args.windows, args.sampling)
    logger.debug("Writing atom mapping for distance restraints")

    run.avoid_water_protein_clashes(
        inputdir, header=f"{run.sphereradius}.0 SPHERE", save_removed=(args.log in ["trace", "debug"])
    )

    logger.debug("Writing the QPREP files & running qprep")
    run.write_qprep(inputdir)
    run.qprep(inputdir)
    logger.debug("Writing FEP files")
    run.write_FEP_file(change_charges, change_vdw, FEP_vdw, inputdir, lig_size1, lig_size2)
    overlapping_atoms = run.set_restraints(writedir, args.restraint_method, strict_check=True)

    # Handling the correct offset here
    logger.debug("Writing the MD files")
    if args.start == "0.5":
        file_list = run.write_MD_05(lambdas, inputdir, lig_size1, lig_size2, overlapping_atoms)
        run.write_runfile(inputdir, file_list)

    if args.start == "1":
        file_list = run.write_MD_1(lambdas, inputdir, lig_size1, lig_size2, overlapping_atoms)
        run.write_runfile(inputdir, file_list)
    logger.debug(f"Generated files: {file_list}")
    logger.debug("Writing the submit files")
    run.write_submitfile(writedir)
    logger.debug("Writing the QFEP files")
    run.write_qfep(args.windows, lambdas)


def main_exe():
    args = parse_arguments(program="QligFEP")
    main(args)


if __name__ == "__main__":
    main_exe()
