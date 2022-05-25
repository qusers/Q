import argparse
import io
import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src/')))

import QmapFEP

class Startup(object):
    """
    Create dual topology FEP files based on two ligands
    """
    def __init__(self, data, *args, **kwargs):
        START = QmapFEP.Init(data)

if __name__ == "__main__":
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
    
    args = parser.parse_args()

    Startup(vars(args))
