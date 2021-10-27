#!/usr/bin/env python3

import argparse
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src/')))

import qdyn

class Startup(object):
    def __init__(self,top,fep,md,re,wd,verbose,gpu,clean):
        data = {'top'       :   top,
                'fep'       :   fep,
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

    parser.add_argument('-w', '--workdir',
                        dest = "wd",
                        default = None,
                        required = False,
                        help = "Working directory")
    
    parser.add_argument('-F', '--oldFEP',
                        dest = "oldFEP",
                        default = None,
                        required = True,                        
                        help = "Directory with QligFEP v1.0 inputfiles, implemented for backward compatability."
                        )

    parser.add_argument('--overwrite',
                        dest = "overwrite",
                        default = False,
                        required = False,                                                
                        action = 'store_true',
                        help = "OVerwrite files in specified working directory")    

    args = parser.parse_args()
    
    Startup(vars(args))
