import argparse
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src/')))

import qcalc

class Startup(object):
    def __init__(self,top,otrj,itrj,ilib,wtraj,calc):
        data = { 'top' : top,
                 'otrj' : otrj,
                 'itrj' : itrj,
                 'ilib' : ilib,
                 'wraj' : wtraj,
                 'calc' : calc,
               }
        START = qcalc.Init(data)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='Qdyn',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '       == Qdyn == ')
    
    parser.add_argument('--version', 
                        action='version', 
                        version='%(prog)s 0.1.0')
    
            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='Qcalc,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '       == Python based MD engine == ')

    
    parser.add_argument('-t', '--top',
                        dest = "top",
                        default = None,
                        required = True,
                        help = "Q topology file")
        
    parser.add_argument('-ot', '--otrj',
                        dest = "otrj",
                        default = None,
                        help = " Output trajectory file")
            
    parser.add_argument('-it', '--itrj',
                        dest = "itrj",
                        default = None,
                        help = " Input trajectory file")
                
    parser.add_argument('-il', '--ilib',
                        dest = "ilib",
                        default = None,
                        help = "Library files for used topology")
    
    parser.add_argument("-wt", "--wtraj",
                        default = False,
                        action = 'store_true',
                        help="Write trajecotry")
    
    parser.add_argument('-c', '--calc',
                        dest = "calc",
                        required = False,
                        help = "Comma seperated list of calculations to be performed, requires a *.calc inputfile")
    
    args = parser.parse_args()
    
    
    
    Startup(top = args.top,
            otrj = args.otrj,
            itrj = args.itrj,
            ilib = args.ilib,
            wtraj = args.wtraj,
            calc = args.calc,
           )