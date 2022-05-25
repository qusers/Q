#!/usr/bin/env python3

import argparse
import os
import sys
import json

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src/')))

import qfep

class Startup(object):
    def __init__(self,ener,wd,inp):
        self.data = {
                    'workdir'           : None,
                    'states'            : None,
                    'offdiag_elements'  : None,
                    'kT'                : None,
                    'points_to_skip'    : None,
                    'only_QQ'           : None,
                    'gap_bins'          : None,
                    'points_per_bin'    : None,
                    'alpha_state'       : None,
                    'linear_combination': [],
                    'energy_files'      : []
        
        }
        
        if inp == None:
            if ener == None or wd == None:
                print(">>> FATAL: need input arguments or inputfile")
                
            else:
                self.data['energy_files'] = ener
        
        else:
            print("Reading data from input file")
            self.inp = inp
            self.read_input()
        
        self.data['workdir']      = wd

        START = qfep.Init(self.data)
    
    def read_input(self):
        block = 0
        if self.inp.split('.')[-1] != 'json':
            with open(self.inp) as infile:
                for line in infile:
                    line = line.strip()
                    line = line.split()
                    if line[0] == 'energy_files':
                        block = 1
                        continue

                    if block == 0:
                        if len(line) == 2:
                            if not line[0] in self.data:
                                print(">>> FATAL: keyword {} in inputfile invalid".format(key))
                            self.data[line[0]] = line[1]

                    if block == 1:
                        self.data['energy_files'].append([line[0],line[1]])
                            
        else:
            with open(self.inp) as json_file:
                self.data = json.load(json_file)
                self.data['energy_files'] = self.data['energy_files']
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='Qfep',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '       == Qfep == ')
    
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
