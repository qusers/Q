import glob
import numpy as np
import argparse
import os
import sys
import itertools
from os import path
import shutil
import math

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src/share/')))

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
        dGfsum = 0
        dGflist = [0.0]

        # construct the energy lookup
        for energyfile in ener:
            infile  = energyfile[0]
            outfile = energyfile[1]
            
            for frame in range(0,len(self.ener)):
                read_ener  = ENERGY.Read_Energy(infile,states)
                ener_data = read_ener.JSON(self.wd + '/' + outfile)
                if not outfile in self.energies:
                    self.energies[outfile] = [ener_data[frame]['q-energies']['SUM']]
                else:
                    self.energies[outfile].append(ener_data[frame]['q-energies']['SUM'])

                self.lambdas[outfile] = ener_data[frame]['q-energies']['lambda']
                #else:
                #    self.lambdas[outfile].append(ener_data[frame]['q-energies']['lambda'])
            #self.lambdas[outfile]  = [ener_data[0]['q-energies']['lambda'],ener_data[1]['q-energies']['lambda']]
        for i, ifile in enumerate(self.ener):
            if ifile == self.ener[-1]:
                continue
            MA1 = self.energies[self.ener[i][1]]
            l_file1 = self.lambdas[self.ener[i][1]]
            l_file2 = self.lambdas[self.ener[i+1][1]]
            
            # Throw the Q energies to calc module
            dGf = CALC.EXP(MA1,l_file1,l_file2,self.kT,skip)
            dGflist.append(dGf)
            
        dGflist = np.array(dGflist)
        dGfsum = np.sum(dGflist)
        
        with open(self.wd + '/qfep.out', 'w') as outfile:
            outfile.write('{}       {}\n'.format('lambda','dGf'))
            for i in range(0,len(dGflist)):
                outfile.write('{}       {:.3f}\n'.format(self.lambdas[self.ener[i][1]][0],
                                                   dGflist[i]
                                                  ))
            outfile.write('dGfsum = {:.3f}'.format(dGfsum))
                
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
