#!/usr/bin/env python3

import argparse
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src/')))

import qcalc

class Startup(object):
    def __init__(self,top,wd,itrj,ilib,otrj,calc):
        data = { 'top'  : top,
                 'wd'   : wd,
                 'itrj' : itrj,
                 'ilib' : ilib,
                 'otrj' : otrj,
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
    
    parser.add_argument('-t', '--top',
                        dest = "top",
                        default = None,
                        required = True,
                        help = "Q .json topology file")
        
    parser.add_argument('-wd', '--workdir',
                        dest = "wd",
                        default = None,
                        help = " Output folder")
            
    parser.add_argument('-it', '--itrj',
                        dest = "itrj",
                        default = None,
                        help = " Input trajectory file")
                
    parser.add_argument('-il', '--ilib',
                        dest = "ilib",
                        default = None,
                        help = "Library files for used topology")
    
    parser.add_argument("-ot", "--otrj",
                        dest = 'otrj',
                        choices = ['.pdb'],
                        help="Write trajecotry")
    
    parser.add_argument('-c', '--calc',
                        dest = "calc",
                        required = False,
                        help = "Comma seperated list of calculations to be performed, requires a *.calc inputfile")
    
    args = parser.parse_args()
    
    
    
    Startup(top = args.top,
            wd = args.wd,
            itrj = args.itrj,
            ilib = args.ilib,
            otrj = args.otrj,
            calc = args.calc,
           )
