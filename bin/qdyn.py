import argparse
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src/')))

import qdyn

class Startup(object):
    def __init__(self,top,fep,md,re,wd,verbose):
        data = {'top'       :   top,
                'fep'       :   fep,
                'md'        :   md,
                're'        :   re,
                'wd'        :   wd,
                'verbose'   :   verbose
               }
        START = qdyn.Init(data)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='Qdyn',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '       == Qdyn == ')
    
    parser.add_argument('--version', 
                        action='version', 
                        version='%(prog)s 0.1.0')
    
    parser.add_argument('-t', '--top',
                        dest = "top",
                        default = None,
                        required = True,
                        help = "Q topology file")
    
    parser.add_argument('-m', '--md',
                        dest = "md",
                        default = None,
                        required = True,                        
                        help = "MD input file, use .inp for Q file format, .json for Python \n"   \
                               "type input"
                        )
           
    parser.add_argument('-f', '--fep',
                        dest = "fep",
                        default = None,
                        required = False,                        
                        help = "FEPfile, use .fep for Q file format, .json for Python")
       
    parser.add_argument('-r', '--re',
                        dest = "re",
                        default = None,
                        required = False,                                                
                        help = "Restart file")
         
    parser.add_argument('-d', '--workdir',
                        dest = "workdir",
                        default = None,
                        required = True,                                                
                        help = "Working directory")
         
    parser.add_argument('--verbose',
                        dest = "verbose",
                        default = False,
                        required = False,                                                
                        action = 'store_true',
                        help = "Working directory")

    args = parser.parse_args()
    
    Startup(top = args.top,
            fep = args.fep,            
            md = args.md,
            re = args.re,
            wd = args.workdir,
            verbose = args.verbose,
           )