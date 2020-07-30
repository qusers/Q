import argparse
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src/')))

import qdyn

class Startup(object):
    def __init__(self,top):
        start = qdyn.Init(top)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='Q',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '       == Q == ')

    parser.add_argument('-t', '--top',
                        dest = "top",
                        default = None,
                        required = True,
                        help = "Q topology file")
    
    parser.add_argument('-i', '--inp',
                        dest = "inp",
                        default = None,
                        required = False,                        
                        help = "Python type MD input file")
        
    parser.add_argument('-d', '--qdir',
                        dest = "qdir",
                        default = None,
                        required = False,                        
                        help = "Directory with standard Q inputfiles")
            
    parser.add_argument('-wp', '--wpython',
                        dest = "wpython",
                        default = None,
                        help = "Toggle to write out python readable json files of the md")
                
    parser.add_argument('-r', '--rundir',
                        dest = "rundir",
                        default = None,
                        help = "Toggle to write out python readable json files of the md")
    
    args = parser.parse_args()
    
    Startup(top = args.top)