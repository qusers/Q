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
    def __init__(self,ener):
        data = Energy()
        self.data = data.data
        
        self.ener = ener
        
    def QDYN(self):    
        with open (self.ener) as infile:
            for line in infile:
                line = line.strip()
                line = line.split('=')
                line[0] = line[0].strip()
                if len(line) == 1:
                    self.data['frames'].append(int(line[0]))
                    
                else:
                    if line[0] == 'Temp':
                        self.data['Temp'].append(float(line[1]))
                        
                    else:
                        self.data['energies'][line[0]].append(float(line[1]))
                        
        self.data['frames'] = self.data['frames'][1:]
        
        return(self.data)
    
class Write_Energy(object):
    def __init__(self,data,wd):
        self.wd = wd
        self.data = data

    def JSON(self):
        with open(self.wd + '/ener.json', 'w') as outfile:
            inputs = self.data
            json.dump(inputs,outfile,indent=2)     
