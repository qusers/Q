import glob
import numpy as np
import argparse
import os
import sys
import itertools
from os import path
import shutil
import json

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src/share/')))

import IO
import functions as f
import potential_energy as Upot
import geometries
import topology as TOPOLOGY        
import defaults as DEFAULTS
import md as MD

class Create_Environment(object):
    """
        Creates the workdirectory environment.
    """
    def __init__(self,top,wd):
        self.top = top
        self.wd = wd
        
        if not os.path.exists(wd):
            os.mkdir(wd)
            
        else:
            shutil.rmtree(wd)
            os.mkdir(wd)        
        
class Prepare_Topology(object):
    """
        Creates a topology object and writes .csv files for qdyn.
    """    
    def __init__(self,top,wd):
        self.top = top
        self.wd = wd
        
        read_top  = TOPOLOGY.Read_Topology(self.top)
        
        # Get the extension and read data
        if self.top.split('.')[-1] == 'json':
            top_data = read_top.JSON()
            
        else:
            top_data = read_top.Q()     
        
        # Initiate the write class
        write_top = TOPOLOGY.Write_Topology(top_data)
        
        # Write the topology in csv and json format
        write_top.CSV(self.wd + '/' + self.top.split('.')[0] + '/')
        write_top.JSON(self.wd + '/' + self.top.split('.')[0] + '/')


class Prepare_MD(object):
    def __init__(self,top,md,wd):
        self.top = top
        self.wd  = wd
        self.md  = md
        read_md  = MD.Read_MD(self.md)
        
        # Get the extension and read data
        if self.md.split('.')[-1] == 'json':
            md_data = read_md.JSON()
            
        else:
            md_data = read_md.Q()     
        
        # Initiate the write class
        write_md = MD.Write_MD(md_data)
        
        # Write md data files (both csv and json file)
        write_md.CSV(self.wd + '/' + self.top.split('.')[0] + '/')
        write_md.JSON(self.wd + '/' + self.top.split('.')[0]+ '/')

class Run_Dynamics(object):
    def __init__(self,rundir):
        #inputs = Prepare_Inputs()
        self.rundir = rundir
        for i, stage in enumerate (MD.MD['stages']):
            Prepare_Inputs.write_csv()(i)


class Init(object):
    def __init__(self, data):
        """ Retrieves a dictionary of user input from qdyn:
               {'top'   :   top,
                'fep'   :   fep,
                'md'    :   md,
                're'    :   re,
                'wd'    :   wd
               }
        """
        self.environment = data
        
        # check extension:
        extensions = ['json','inp','fep','re','top']
        
        if data['top'].split('.')[-1] not in extensions:
            print(data['top'].split('.')[-1])
            print("FATAL: unrecognized extension for {}".format(data['top']))
            sys.exit()
                
        if data['md'].split('.')[-1] not in extensions:
            print("FATAL: unrecognized extension for {}".format(data['md']))
            sys.exit()
            
        if data['re'] != None:
            if data['re'].split('.')[-1] not in extensions:
                print("FATAL: unrecognized extension for {}".format(data['re']))
                sys.exit()
        
        if data['fep'] != None:
            if data['fep'].split('.')[-1] not in extensions:
                print("FATAL: unrecognized extension for {}".format(data['fep']))
                sys.exit()
    
        # running specified stuff
        Prepare_Topology(top = self.environment['top'],
                         wd  = self.environment['wd'],
                        )
        
        Prepare_MD(top = self.environment['top'],
                   wd  = self.environment['wd'],
                   md  = self.environment['md'],
                  )
    #Run_Dynamics(rundir = args.rundir)