# Standard Python libraries
import os
import itertools
from os import path
import json

# Q-GPU libraries
import IO

class Energy(object):
    """ Trajectory object class, one nparray per timeframe"""
    def __init__(self):
        self.data = {
                      'frames'      : [],
                      'Temp'        : [],
                      'q_total'     : [],
                      'energies'    : { 'Ubond':[],
                                        'Uangle':[],
                                        'Utor':[],
                                        'Uradx':[],
                                        'Upolx':[],
                                        'Ushell':[],
                                        'Ufix':[],
                                        'Upres':[],
                                        'Urestr':[],
                                        'Ucoul':[],
                                        'Uvdw':[],
                                        'Ukin':[],
                                        'Upot':[],
                                        'Utot':[],                                     
                                      }
                    }        
        
class Read_Energy(object):
    def __init__(self,ener,states):
        data = Energy()
        self.data = data.data
        self.ener = ener
        self.states = states
        
        for state in range(0, int(self.states)):
            self.data['q_total'].append([])
        
    def QDYN(self):
        with open (self.ener) as infile:
            block = 0
            for line in infile:
                line = line.strip()
                line = line.split('=')
                line[0] = line[0].strip()
                
                # read in block
                if line[0].split()[0] == '>>>':
                    line = line[0].split()
                    if line[1] == 'State':
                        state = int(line[2])
                        block = 1
                        
                    if line[1] == 'Total':
                        block = 2
                
                    continue
                
                if block == 1:
                    if line[0] == 'Q-SUM':
                        state = int(state)
                        self.data['q_total'][state].append(line[1])
                
                if block == 2:
                    if len(line) == 1:
                        self.data['frames'].append(int(line[0]))

                    else:
                        if line[0] == 'Temp':
                            self.data['Temp'].append(float(line[1]))

                        else:
                            self.data['energies'][line[0]].append(float(line[1]))

        return(self.data)
    
    def JSON(self,outfile):
        with open(outfile) as json_file:
            self.data = json.load(json_file)
        return(self.data)
    
class Write_Energy(object):
    def __init__(self,data,outfile):
        self.outfile = outfile
        self.data = data

    def JSON(self):
        with open(self.outfile, 'w') as outfile:
            inputs = self.data
            json.dump(inputs,outfile,indent=2)     
