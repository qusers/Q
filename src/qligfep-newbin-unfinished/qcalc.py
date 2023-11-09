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

class Initiate_Mask(object): 
    def __init__(self, top, wd, ilib):
        self.top = top
        self.wd = wd
        self.ilib = ilib
        
        read_top  = TOPOLOGY.Read_Topology(self.top)

        # Get the extension and read data
        if self.top.split('.')[-1] == 'json':
            self.top_data = read_top.JSON()

        else:
            self.top_data = read_top.Q()
                
        # Read the library
        read_lib = MASK.Read_Library(ilib)
        self.prm2lib = read_lib.Q()
        
        # Now get the mask and store it
        create_mask = MASK.Construct_Mask(self.top_data,
                                          self.prm2lib
                                         )
        
        self.mask = create_mask.Make() 
        
        # Also save to .json file?
        create_mask.write_JSON(self.wd)
        
class Initiate_Trajectory(object):
    def __init__(self,itrj,wd,otrj):
        self.itrj = itrj
        self.otrj = otrj
        self.wd   = wd
        self.mask = self.wd + '/mask.json'
    
        read_traj = TRAJECTORY.Read_Trajectory(self.itrj)
        
        # Maybe add .json parsing?
        traj = read_traj.QDYN()
        
        read_mask = MASK.Read_Mask(self.mask)
        mask = read_mask.JSON()['mask']
        
        if not len(mask) == int(traj['natoms']):
            print(">>> FATAL: number of atoms {} in mask do not match {} " \
                  "atoms in trajectory.".format(len(mask),traj['natoms']))
            sys.exit()
        
        write_traj = TRAJECTORY.Write_Trajectory(self.wd,self.otrj,mask,traj)
        if self.otrj == '.pdb':
            write_traj.PDB()
        write_traj.JSON()
        
class Calculate(object):
    def __init__(self,calc,wd):
        self.calc = calc
        self.wd = wd
        
        read_mask = MASK.Read_Mask(self.wd + '/mask.json')
        mask = read_mask.JSON()
        
        read_traj = TRAJECTORY.Read_Trajectory(self.wd + '/traj.json')
        traj =read_traj.JSON()
        
        if calc == 'number_density':
            calc = CALC.number_density(mask,traj,wd)
        
class Init(object):
    def __init__(self, data):
        """ Retrieves a dictionary of user input from qalc.py:
               {'top' : top,
                 'wd' : wd,
                 'itrj' : itrj,
                 'ilib' : ilib,
                 'otrj' : otrj,
                 'calc' : calc,
               }
        """   
        self.environment = data
        
        # Check if multiple inpit libraries
        self.environment['ilib'] = self.environment['ilib'].strip()
        self.environment['ilib'] = self.environment['ilib'].split(',')
        
        # Create user specified work environment
        Create_Environment(self.environment['wd'])
    
        # All the qcalc modules need a mask to map coordinate system
        # to the topology/pdb in formation,
        Initiate_Mask(self.environment['top'],
                      self.environment['wd'],
                      self.environment['ilib'],
                     )

        Initiate_Trajectory(self.environment['itrj'],
                            self.environment['wd'],
                            self.environment['otrj'],
                           )    

        if self.environment['calc'] != None:
            calcs = self.environment['calc'].split(',')
            for calc in calcs:
                if not os.path.exists(calc + '.calc'):
                    print(">>> FATAL: could not find input file for {}.calc".format(calc))
                    sys.exit()

        Calculate(self.environment['calc'],
                  self.environment['wd'],
                 )
