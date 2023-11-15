import argparse
from QligFEP.qligfep import QligFEP

def parse_arguments():
    parser = argparse.ArgumentParser(
        prog='QligFEP',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '       == Generate FEP files for dual topology ligand FEP == ')

    parser.add_argument('-l1', '--lig_1',
                        dest = "lig1",
                        required = True,
                        help = "name of ligand 1")

    parser.add_argument('-l2', '--lig_2',
                        dest = "lig2",
                        required = True,
                        help = "name of ligand 2")

    parser.add_argument('-FF', '--forcefield',
                        dest = "FF",
                        required = True,
                        choices = ['OPLS2005', 'OPLS2015', 'AMBER14sb', 'CHARMM36', 'CHARMM22', 'CHARMM_TEST'],
                        help = "Forcefield to be used")

    parser.add_argument('-s', '--system',
                        dest = "system",
                        required = True,
                        choices = ['water', 'protein', 'vacuum'],
                        help = "what type of system we are setting up")

    parser.add_argument('-c', '--cluster',
                        dest = "cluster",
                        required = True,
                        help = "cluster you want to submit to, cluster specific parameters added to settings"
                       )

    parser.add_argument('-r', '--sphereradius',
                        dest = "sphereradius",
                        required = False,
                        default = '15',
                        help = "size of the simulation sphere"
                       )

    parser.add_argument('-b', '--cysbond',
                        dest = "cysbond",
                        default = None,
                        help = "Temporary function to add cysbonds at1:at2,at3:at4 etc."
                       )

    parser.add_argument('-l', '--start',
                        dest = "start",
                        default = '0.5',
                        choices = ['1', '0.5'],
                        help = "Starting FEP in the middle or endpoint"
                       )

    parser.add_argument('-T', '--temperature',
                        dest = "temperature",
                        default = '298',
                        help = "Temperature(s), mutliple tempereratures given as 'T1,T2,...,TN'"
                       )

    parser.add_argument('-R', '--replicates',
                        dest = "replicates",
                        default = '10',
                        help = "How many repeats should be run"
                       )

    parser.add_argument('-S', '--sampling',
                        dest = "sampling",
                        default = 'linear',
                        choices = ['linear', 'sigmoidal', 'exponential', 'reverse_exponential'],
                        help = "Lambda spacing type to be used"
                       )

    parser.add_argument('-w', '--windows',
                        dest = "windows",
                        default = '50',
                        help = "Total number of windows that will be run"
                       )
    return parser.parse_args()


def main():
    args = parse_arguments()
    run = QligFEP(lig1 = args.lig1,
              lig2 = args.lig2,
              FF= args.FF,
              system = args.system,
              cluster = args.cluster,
              sphereradius = args.sphereradius,
              cysbond = args.cysbond,
              start = args.start,
              temperature = args.temperature,
              replicates = args.replicates,
              sampling =args.sampling,
              windows= args.windows
             )

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
    run.change_lib(changes_for_libfiles, inputdir)
    FEP_vdw = run.change_prm(changes_for_prmfiles, inputdir)
    run.write_FEP_file(change_charges, change_vdw, FEP_vdw, inputdir, lig_size1, lig_size2)
    run.merge_pdbs(inputdir)
    if args.system == 'protein':
        run.write_water_pdb(inputdir)
    lambdas = run.get_lambdas(args.windows, args.sampling)
    overlapping_atoms = run.overlapping_atoms(writedir)
    
    # Handling the correct offset here
    if args.start == '0.5':
        file_list = run.write_MD_05(lambdas, inputdir, lig_size1, lig_size2)
        run.write_runfile(inputdir, file_list)    
        
    if args.start == '1':
        file_list = run.write_MD_1(lambdas, inputdir, lig_size1, lig_size2, overlapping_atoms)
        run.write_runfile(inputdir, file_list)    
    
    run.write_submitfile(writedir)
    run.write_qfep(inputdir, args.windows, lambdas)
    run.write_qprep(inputdir)
    run.qprep(inputdir)