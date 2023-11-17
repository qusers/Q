import argparse
from QligFEP.qmapfep import Init


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog='QmapFEP',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='FEP map generator based on selected distance metrics.')
    parser.add_argument('-i', '--insdf',
                        dest="isdf",
                        required=True,
                        help=".sdf file name")
    parser.add_argument('-m', '--metric',
                        dest="metric",
                        default='MFP',
                        choices=['MFP', 'Tanimoto', 'MCS', 'SMILES'],
                        required=False,
                        help="Distance metric for ligand pairwairse comparison")
    parser.add_argument('-o', '--otxt',
                        dest="o",
                        required=False,
                        default=None,
                        help=(
                            "Name for output .json file. Default will be the lowercase "
                            "name of the input sdf."
                        ))
    parser.add_argument('-wd', '--workdir',
                        dest="wd",
                        required=False,
                        default=None,
                        help=(
                            "Name of the directory to store the output files. Default "
                            "will create a directory with the same name as the input "
                            "sdf file. "
                        ))
    
    return parser.parse_args()
    
def main():
    args = parse_arguments()
    if args.wd is None:
        args.wd = args.isdf.split('.')[0]
    if args.o is None:
        args.o = args.isdf.split('.')[0].lower()
    data = {
        'isdf': args.isdf,
        'metric': args.metric,
        'o': args.o,
        'wd': args.wd
    }
    Init(data)
