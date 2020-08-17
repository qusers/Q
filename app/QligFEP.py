import argparse
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src/')))

import QligFEP

class Startup(object):
    """
    Create dual topology FEP files based on two ligands
    """
    def __init__(self, oldFEP, wd, cluster, overwrite,*args, **kwargs):
        
        data = {'oldFEP'    : oldFEP,
                'wd'        : wd,
                'cluster'   : cluster,
                'overwrite' : overwrite
               }
        
        START = QligFEP.Init(data)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='QligFEP',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '       == Generate FEP files for dual topology ligand FEP == ')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')    
    parser.add_argument('-of', '--oldFEP',
                        dest = "oldFEP",
                        required = False,
                        help = "If specified, read old QligFEP input file format")
        
    parser.add_argument('-wd', '--workdir',
                        dest = "wd",
                        required = True,
                        help = "Workdir folder")
            
    parser.add_argument('-c', '--cluster',
                        dest = "cluster",
                        required = False,
                        help = "Workdir folder")
                
    parser.add_argument('-ow', '--overwrite',
                        dest = "overwrite",
                        required = False,
                        default = False,
                        action = 'store_true',
                        help = "Set flag to overwrite old folder completely")
    
    args = parser.parse_args()
    Startup(oldFEP = args.oldFEP,
            wd   = args.wd,
            cluster   = args.cluster,
            overwrite   = args.overwrite,
            )
