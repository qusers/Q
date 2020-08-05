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
    def __init__(self,ener,wd):
        self.ener = ener
        self.wd = wd
        
        read_ener  = ENERGY.Read_Energy(self.ener)
        
        # Read the energy from QDYN
        ener_data = read_ener.QDYN()   
        
        # Initiate the write class
        write_ener = ENERGY.Write_Energy(ener_data,self.wd)
        
        # Write the topology in csv and json format
        write_ener.JSON()
        
        avg_Upot = np.average(ener_data['energies']['Upot'])
        sem_Upot = np.std(ener_data['energies']['Upot'])
        
        print(avg_Upot, sem_Upot)
        
            
class Init(object):
    def __init__(self, data):
        """ Retrieves a dictionary of user input from qalc.py:
               {'ener' : ener,
                'wd' : wd,
               }
        """
        self.environment = data
        # Create user specified work environment
        Create_Environment(self.environment['wd'])
    
        # All the qcalc modules need a mask to map coordinate system
        # to the topology/pdb in formation,
        Get_Energy(self.environment['ener'],
                    self.environment['wd'],
                     )
