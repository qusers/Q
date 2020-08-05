import argparse
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src/')))

import qfep

class Startup(object):
    def __init__(self,ener,wd,inp):
        if inp == None:
            if ener == None or wd ==None :
                print(">>> FATAL: need input arguments or inputfile")
                
            else:
                data = { 'ener'  : ener,
                         'wd'   : wd,
                       }
                
        else:
            print("Reading data from input file")
            self.inp = inp
            data = self.read_input()
        
        START = qfep.Init(data)
    
    def read_input(self):
        print(self.inp)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='Qdyn',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '       == Qdyn == ')
    
    parser.add_argument('--version', 
                        action='version', 
                        version='%(prog)s 0.1.0')
    
    parser.add_argument('-e', '--ener',
                        dest = "ener",
                        default = None,
                        required = False,
                        help = "Energy file from Qdyn")
         
    parser.add_argument('-d', '--workdir',
                        dest = "wd",
                        default = None,
                        required = False,                                                
                        help = "Working directory")

          
    parser.add_argument('-i', '--input',
                        dest = "inp",
                        default = None,
                        required = False,                                                
                        help = "Input file, other arguments ignored ")

    args = parser.parse_args()
        
    Startup(ener = args.ener,
            wd   = args.wd,
            inp  = args.inp,
           )
