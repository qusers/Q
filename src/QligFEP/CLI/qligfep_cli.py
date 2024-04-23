"""Module containing the QligFEP command line interface."""

import argparse
from QligFEP.qligfep import QligFEP
from typing import Optional
from ..logger import logger

def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog='QligFEP',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '       == Generate FEP files for dual topology ligand FEP == ')

    parser.add_argument('-l1', '--lig_1',
                        dest = "lig1",
                        required = True,
                        help = "name of ligand 1",
                        type=str)

    parser.add_argument('-l2', '--lig_2',
                        dest = "lig2",
                        required = True,
                        help = "name of ligand 2",
                        type=str)

    parser.add_argument('-FF', '--forcefield',
                        dest = "FF",
                        required = True,
                        choices = ['OPLS2005', 'OPLS2015', 'AMBER14sb', 'CHARMM36', 'CHARMM22', 'CHARMM_TEST'],
                        help = "Forcefield to be used.")

    parser.add_argument('-s', '--system',
                        dest = "system",
                        required = True,
                        choices = ['water', 'protein', 'vacuum'],
                        help = "what type of system we are setting up")

    parser.add_argument('-c', '--cluster',
                        dest = "cluster",
                        required = True,
                        help = "cluster you want to submit to, cluster specific parameters added to settings."
                       )

    parser.add_argument('-r', '--sphereradius',
                        dest = "sphereradius",
                        required = False,
                        default = '25',
                        help = "Size of the simulation sphere. Defaults to 25."
                       )

    parser.add_argument('-b', '--cysbond',
                        dest = "cysbond",
                        default = None,
                        help = (
                            "Add cystein bonds. Input should be formatted with the atom numbers"
                            "(participating in the Cys bond) connected by `_` and with different bonds "
                            "separated by `,` as in: `atom1_atom2,atom3_atom4`"
                            )
                        )

    parser.add_argument('-l', '--start',
                        dest = "start",
                        default = '0.5',
                        choices = ['1', '0.5'],
                        help = "Starting FEP in the middle or endpoint. Defaults to 0.5."
                       )

    parser.add_argument('-T', '--temperature',
                        dest = "temperature",
                        default = '298',
                        help = "Temperature(s), mutliple tempereratures given as 'T1,T2,...,TN'. Defaults to 298K"
                       )

    parser.add_argument('-R', '--replicates',
                        dest = "replicates",
                        default = '10',
                        help = "How many repeats should be run. Defaults to 10."
                       )

    parser.add_argument('-S', '--sampling',
                        dest = "sampling",
                        default = 'sigmoidal',
                        choices = ['linear', 'sigmoidal', 'exponential', 'reverse_exponential'],
                        help = "Lambda spacing type to be used"
                       )

    parser.add_argument('-w', '--windows',
                        dest = "windows",
                        help = "Total number of windows that will be run. Defaults to 100.",
                        type=str,
                       )
    parser.add_argument('-sc', '--softcore',
                        dest = "softcore",
                        default = False,
                        action="store_true",
                        help = "Turn on if you want to use softcore"
                       )
    parser.add_argument('-ts', '--timestep',
                        dest = "timestep",
                        choices = ['1fs','2fs'],
                        default = "2fs",
                        help = "Simulation timestep, default 2fs"
                       )
    parser.add_argument('-clean', '--files-to-clean',
                        dest="to_clean",
                        nargs="+",
                        default=None,
                        help=(
                            "Files to clean after the simulation. The arguments are given as a list of strings "
                            "and the cleaning is done by adding the command `rm -rf *{arg1} *{arg2}` to the job submission. "
                            "Usage example: `-clean dcd` will remove all dcd files after the simulation. If left as None, won't clean any files."
                            )
                        )
    return parser.parse_args()


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
            "lig1" : args.lig1,
            "lig2" : args.lig2,
            "FF": args.FF,
            "system": args.system,
            "cluster": args.cluster,
            "sphereradius": args.sphereradius,
            "cysbond": args.cysbond,
            "start": args.start,
            "temperature": args.temperature,
            "replicates": args.replicates,
            "sampling": args.sampling,
            "timestep": args.timestep,
            "softcore": args.softcore,
            "to_clean": args.to_clean,
        }
    else: 
        param_dict = {}
    param_dict.update(kwargs)
    run = QligFEP(**param_dict)
    run.set_timestep()

    writedir = run.makedir()
    inputdir = writedir + '/inputfiles'
    a = run.read_files()
    changes_for_libfiles = a[0][1]
    changes_for_prmfiles = a[0][1]
    change_charges       = a[1][0]
    change_vdw           = a[1][1]
    changes_for_pdbfiles = a[0][0]
    lig_size1, lig_size2 = a[2][0], a[2][1]

    # Write the merged files
    logger.debug("Writing changes on files")
    run.change_lib(changes_for_libfiles, inputdir)
    logger.debug("Changing parameters")
    FEP_vdw = run.change_prm(changes_for_prmfiles, inputdir)
    logger.debug("Writing FEP files")
    run.write_FEP_file(change_charges, change_vdw, FEP_vdw, inputdir, lig_size1, lig_size2)
    logger.debug('Writing PDB files')
    run.merge_pdbs(inputdir)
    if args.system == 'protein':
        run.write_water_pdb(inputdir)
    logger.debug('Getting the lambdas')
    lambdas = run.get_lambdas(args.windows, args.sampling)
    logger.debug('Run the overlapping atoms')
    overlapping_atoms = run.overlapping_atoms(writedir)
    
    # Handling the correct offset here
    logger.debug('Writing the MD files')
    if args.start == '0.5':
        file_list = run.write_MD_05(lambdas, inputdir, lig_size1, lig_size2, overlapping_atoms)
        run.write_runfile(inputdir, file_list)    
        
    if args.start == '1':
        file_list = run.write_MD_1(lambdas, inputdir, lig_size1, lig_size2, overlapping_atoms)
        run.write_runfile(inputdir, file_list)    
    logger.debug(f'Generated files: {file_list}')
    logger.debug('Writing the submit files')
    run.write_submitfile(writedir)
    logger.debug('Writing the QFEP files')
    run.write_qfep(args.windows, lambdas)
    logger.debug('Writing the QPREP files')
    run.write_qprep(inputdir)
    logger.debug('Running QFEP')
    run.qprep(inputdir)

def main_exe():
    args = parse_arguments()
    main(args)

if __name__ == '__main__':
    main_exe()