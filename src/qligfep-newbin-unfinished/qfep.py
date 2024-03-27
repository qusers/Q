import glob
import numpy as np
import argparse
import os
import sys
import itertools
from os import path
import shutil
import math

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../share/')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../env/')))

# Q-GPU packages
import IO
import defaults     as DEFAULTS
import md           as MD
import settings     as SETTINGS
import fep          as FEP
import topology     as TOPOLOGY
import mask         as MASK
import trajectory   as TRAJECTORY
import calc         as CALC
import energy       as ENERGY

class Create_Environment(object):
    """
        Creates the workdirectory environment.
    """
    def __init__(self,wd):
        if not os.path.exists(wd):
            os.mkdir(wd)
            
        else:
            shutil.rmtree(wd)
            os.mkdir(wd)

class Get_Energy(object):
    """
        Creates a topology object.
        Reads: 
                .json and .top topology files.
        Writes:
                .csv and .json topology files.
    """    
    def __init__(self,ener,wd,states):
        self.wd = wd
        for energyfile in ener:       
            infile  = energyfile[0]
            outfile = energyfile[1]
            
            #if outfile != 'out/md_1000_0000.json':
            #    continue                 
            read_ener  = ENERGY.Read_Energy(infile,states)

            # Read the energy from QDYN
            ener_data = read_ener.QDYN()

            # Initiate the write class
            write_ener = ENERGY.Write_Energy(ener_data,self.wd + '/' + outfile)

            # Write the topology in csv and json format
            write_ener.JSON()
            
class Calc_FE(object):
    def __init__(self,ener,wd,states,kT,skip):
        self.ener = ener
        self.wd = wd
        self.states = states
        self.kT = kT
        self.energies = {}
        self.lambdas = {}
        #points to skip
        skip = int(skip)
        dGflist = []
        dGrlist = []
        dGOSlist = []
        dGBARlist = []


        # construct the energy lookup
        for energyfile in ener:
            infile  = energyfile[0]
            outfile = energyfile[1]

            read_ener  = ENERGY.Read_Energy(infile,states)
            ener_data = read_ener.JSON(self.wd + '/' + outfile)

            for frame in range(0,len(ener_data)):
                if not outfile in self.energies:
                    self.energies[outfile] = []

                self.energies[outfile].append(ener_data[frame]['q-energies']['SUM'])

                self.lambdas[outfile] = ener_data[frame]['q-energies']['lambda']

        for i, ifile in enumerate(self.ener):
            print(ifile)
            if ifile == self.ener[-1]:
                continue

            MA1 = self.energies[self.ener[i][1]]
            l_file1 = self.lambdas[self.ener[i][1]]
            l_file2 = self.lambdas[self.ener[i+1][1]]

            # Throw the Q energies to calc module
            dGf = CALC.EXPfwd(MA1,l_file1,l_file2,self.kT,skip)
            dGr = CALC.EXPbwd(MA1,l_file1,l_file2,self.kT,skip)

            dGOS = CALC.overlap_sampling(MA1,l_file1,l_file2,self.kT,skip)
            dGBAR = CALC.BAR(MA1,l_file1,l_file2,self.kT,skip,dGOS)

            dGflist.append(dGf)
            dGrlist.append(dGr)
            dGOSlist.append(dGOS)
            dGBARlist.append(dGBAR)            

        dGflist = np.array(dGflist)
        dGfsum = np.sum(dGflist)
        #print(dGflist)

        dGrlist = np.array(dGrlist)
        dGrlist = dGrlist[::-1]
        dGrsum = np.sum(dGrlist)

        #print(dGrlist)

        dGOSsum = np.sum(np.array(dGOSlist))
        #print(dGOSsum)

        dGBARsum = np.sum(np.array(dGBARlist))
        print(dGBARlist,dGBARsum)        

        dGlist = []
        for i in range(0, len(dGflist)):
            a = (dGflist[i] - dGrlist[i])/2
            dGlist.append(a)
        
        dGavg = np.sum(dGlist)
        dG = np.sum(np.array(dGlist))
        #dG = np.sum(np.array(dGflist))

        with open(self.wd + '/qfep.out', 'w') as outfile:
            outstring = """
Qfep outfile 

dG = {dG:.2f}


""".format(dG=dG)
            outfile.write(outstring)
                
class Init(object):
    def __init__(self, data):
        """ Retrieves a dictionary of user input from qfep.py:
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
        """
        self.environment = data
        # Create user specified work environment
        self.environment['energy_files'] = self.environment['energy_files']
        
        Create_Environment(self.environment['workdir'])
        
    
        # All the qcalc modules need a mask to map coordinate system
        # to the topology/pdb in formation,
        Get_Energy(self.environment['energy_files'],#[0:3],  #Temporarily looking at only 43 files as something goes wrong in the matrix population
                   self.environment['workdir'],
                   self.environment['states'],
                  )

        Calc_FE(self.environment['energy_files'],#[0:3],  #Temporarily looking at only 43 files as something goes wrong in the matrix population
                self.environment['workdir'],
                self.environment['states'],
                self.environment['kT'],
                self.environment['points_to_skip'],
                 )
