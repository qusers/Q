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
                'temperature'   : {
                    'Temp' : None
                },
        
                'bonded'        : {
                    'p' : [None,None,None,None],
                    'w' : [None,None,None,None],
                    'qp': [None,None,None,None], # should be q? 
                },
            
                'nonbonded'     : {
                    'pp' : [None,None],
                    'pw' : [None,None],
                    'ww' : [None,None],
                    'qx' : [None,None],
                },
        
                'restraint'     : {
                    'Uradx'  :   None,
                    'Upolx'  :   None,
                    'Ushell' :   None,
                    'Ufix'   :   None,
                    'Upres'  :   None,
                    'Total'  :   None,
                },
        
                'q-energies'    : {
                    'lambda'    :  [],
                    'SUM'       :  [],
                    'Ubond'	    :  [],
                    'Uangle'	:  [],
                    'Utor'      :  [],
                    'Uimp'      :  [],	
                    'Uvdw'      :  [],	
                    'Ucoul'     :  [],	
                    'Urestr'    :  [],
                },
        
                'total'         : {
                    'Ukin'    :  None,
                    'Upot'    :  None,
                    'Utot'    :  None,
                },   
               }
        
class Read_Energy(object):
    def __init__(self,ener,states):
        data = Energy()
        self.data = data.data
        self.ener = ener
        self.states = states
        
    def QDYN(self):
        total = []
        block = -1
        # Get the total blocks in the energy file
        with open (self.ener) as infile:
            for line in infile:
                if 'interval' in line:
                    interval = line.split()[1]

        interval = int(interval)
        #populate
        for i in range(0,interval+1):
            data = Energy()
            data = data.data
            total.append(data)

        with open (self.ener) as infile:
            #print("reading energies in file {}".format(self.ener))
            i = -1
            for line in infile:
                if len(line) < 2:
                    continue

                if 'lambdas' in line:
                    block = -1
                    continue

                # reset the block
                if 'interval' in line:
                    i += 1
                    block = 0

                if 'type' in line:
                    continue

                if block == -1:
                    # construct lambdas
                    states = int(line[0])
                    block = -2
                
                # I believe this was here for the summary steps
                #if i % 100 != 0:
                #    continue
                    
                #if i % 100 == 0:
                #    i = int(i/100)
                    
                # if block == 0:
                #     data = Energy()
                #     data = data.data
                #     total.append(data)    
                    
                # Find header
                if '[temperature]' in line:
                    block = 1
                    continue

                if '[bonded]' in line:
                    block = 2
                    continue

                if '[nonbonded]' in line:
                    block = 3
                    continue

                if '[restraint]' in line:
                    block = 4
                    continue

                if '[q-energies]' in line:
                    l = -1
                    block = 5
                    for key in total[i]['q-energies']:
                        uh = []
                        for j in range(0,states):
                            uh.append(None)
                        total[i]['q-energies'][key] = uh
                    continue

                if '[total]' in line:
                    block = 6
                    continue

                if block == 1:
                    line = line.split()
                    total[i]['temperature'][line[0]] = float(line[1])

                if block == 2:
                    line = line.split()
                    total[i]['bonded'][line[0]][0] = float(line[1])
                    total[i]['bonded'][line[0]][1] = float(line[2])
                    total[i]['bonded'][line[0]][2] = float(line[3])
                    total[i]['bonded'][line[0]][3] = float(line[4])

                if block == 3:
                    line = line.split()
                    total[i]['nonbonded'][line[0]][0] = float(line[1])
                    total[i]['nonbonded'][line[0]][1] = float(line[2])

                if block == 4:
                    line = line.split()
                    total[i]['restraint'][line[0]] = float(line[1])

                if block == 5:
                    # skip header line
                    if 'lambda' in line:
                        continue
                        
                    line = line.split()
                        
                    l += 1
                    total[i]['q-energies']['lambda'][l] = float(line[0])
                    total[i]['q-energies']['SUM'][l]    = float(line[1])
                    total[i]['q-energies']['Ubond'][l]  = float(line[2])
                    total[i]['q-energies']['Uangle'][l] = float(line[3])
                    total[i]['q-energies']['Utor'][l]   = float(line[4])
                    total[i]['q-energies']['Uimp'][l]   = float(line[5])
                    total[i]['q-energies']['Uvdw'][l]   = float(line[6])
                    total[i]['q-energies']['Ucoul'][l]  = float(line[7])
                    total[i]['q-energies']['Urestr'][l] = float(line[8])

                if block == 6:
                    line = line.split()
                    total[i]['total'][line[0]] = float(line[1])

        return(total)
    
    def JSON(self,outfile):
        with open(outfile) as json_file:
            data = json.load(json_file)
        return(data)
    
class Write_Energy(object):
    def __init__(self,data,outfile):
        self.outfile = outfile
        self.data = data

    def JSON(self):
        with open(self.outfile, 'w') as outfile:
            inputs = self.data
            json.dump(inputs,outfile,indent=2)     
