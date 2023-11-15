import argparse
from QligFEP.qmapfep import Init

class Startup(object):
    """
    Create dual topology FEP files based on two ligands
    """
    def __init__(self, data, *args, **kwargs):
        START = Init(data)
        
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
                        required=True,
                        help="Name for output file")

    parser.add_argument('-wd', '--workdir',
                        dest="wd",
                        default='workdir',
                        help="Name for the working directory")                        
    
    return parser.parse_args()
    
def main(args):
    parser = argparse.ArgumentParser(
        prog='QmapFEP',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='FEP map generator based on selected distance metrics.')
    Startup(vars(args))
    

if __name__ == "__main__":
    args = parse_arguments()
    main(args)
