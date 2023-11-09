#!/usr/bin/end python3
import glob
import numpy as np
import argparse
import os
import sys
import itertools
from os import path
import shutil
import json

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../env/')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../share/')))

import IO
import topology as TOPOLOGY        
import defaults as DEFAULTS
import md       as MD
import settings as SETTINGS
import fep      as FEP
import restart  as RESTART

class Create_Environment(object):
    """
        Creates the workdirectory environment.
    """
    def __init__(self,top,wd):
        if not os.path.exists(wd.split('/')[0]):
            os.mkdir(wd.split('/')[0])
            
        else:
            shutil.rmtree(wd.split('/')[0])
            os.mkdir(wd.split('/')[0])
        
        if not os.path.exists(wd):
            os.mkdir(wd)
            
        else:
            shutil.rmtree(wd)
            os.mkdir(wd)
            
        # create the output folder for qdyn
        os.mkdir(wd + '/output')
        
class Prepare_Topology(object):
    """
        Creates a topology object.
        Reads: 
                .json and .top topology files.
        Writes:
                .csv and .json topology files.
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
        write_top.CSV(self.wd + '/')
        
        out_json = self.wd + '/' + self.wd.split('/')[1] + '.json'
        write_top.JSON(out_json)

class Prepare_MD(object):
    """
        Creates a topology object.
        Reads: 
                .json and .inp md input files.
        Writes:
                .csv and .json md input files.
    """        
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
        write_md.CSV(self.wd + '/')
        
        out_json = self.wd + '/' + self.md.split('.')[1] + '.json'
        write_md.JSON(out_json)        

class Prepare_FEP(object):       
    """
                     verbose = self.environment['verbose'])
        Creates an FEP object.
        Reads: 
                .json and .fep FEP files.
        Writes:
                .csv and .json FEP files.
    """  
    def __init__(self,fepfile,wd,top):
        self.fepfile = fepfile
        self.wd = wd
        self.top = top
        
        read_fep  = FEP.Read_Fep(self.fepfile)
        
        # Get the extension and read data
        if self.fepfile != None:
            if self.fepfile.split('.')[-1] == 'json':
                fep_data = read_fep.JSON()

            else:
                fep_data = read_fep.Q()     

        else:
            fep_data = read_fep.EMPTY()
            
        # Initiate the write class
        write_fep = FEP.Write_Fep(fep_data)

        # Write the topology in csv and json format
        write_fep.CSV(self.wd + '/')
        
        if self.fepfile != None:
            out_json = self.wd + '/' + self.fepfile.split('.')[1] + '.json'
            write_fep.JSON(out_json)

class Read_Restart(object):       
    """
        Creates restart files.
        Reads: 
                .json and .csv restart files (TO DO .re)
        Writes:
                .json and .csv restart files.
    """  
    def __init__(self,restart,wd,top):
        self.restart = restart
        self.wd = wd
        self.top = top

        read_restart = RESTART.Read_Restart(self.restart)
        
        # Get the extension and read data
        if self.restart != None:
            if self.restart.split('.')[-1] == 'json':
                self.velocities,self.coordinates = read_restart.JSON()

            else:
                self.velocities,self.coordinates = read_restart.CSV()     
                
        # construct empty
        else:
            self.velocities = []
            self.coordinates = []
        
        # Initiate the write class
        write_re = RESTART.Write_Restart(self.velocities, self.coordinates)
            
        # Write md data file
        write_re.CSV(self.wd + '/')
        #write_re.JSON()
        
class Run_Dynamics(object):
    """
        Runs the main dynamics loop.

        Run_Dynamics(wd  = self.environment['wd'],
                     top = self.environment['top'],
                     core = self.environment['core'],
                     clean = self.environment['clean'],
                     md = self.environment['clean'],
                     fep = self.environment['fep'],
                     verbose = self.environment['verbose'])

    """           
    def __init__(self,wd,top,core,clean,md,fep,verbose):
        executables = {
            'q6'        : '{}/bin/q6/qdyn '.format(SETTINGS.ROOT),
            'q6-mpi'    : '{}/bin/q6/qdynp '.format(SETTINGS.ROOT),
            'q7'        : '{}/bin/qdyn '.format(SETTINGS.ROOT),
            'q7-gpu'    : '{}/bin/qdyn --gpu '.format(SETTINGS.ROOT)
        }

        options = wd + '/ ' 

        if 'q7' in core:
            out = IO.run_command(executables[core],options)

        if 'q6' in core:
            q_input = '{}.inp'.format(md.split('.')[0])
            IO.json2Q(md,q_input)

        if verbose == True:
            print(out.decode("utf-8"))
            
        if clean == True:
            for csvfile in glob.glob(wd + '/*csv'):
                os.remove(csvfile)
            
class Init(object):
    def __init__(self, data):
        """ Retrieves a dictionary of user input from qdyn:
               {'top'       :   top,
                'fep'       :   fep,
                'md'        :   md,
                're'        :   re,
                'wd'        :   wd,
                'verbose'   :   verbose
                'clean'   :   clean
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

        if data['fep'] != None:
            if data['fep'].split('.')[-1] not in extensions:
                print("FATAL: unrecognized extension for {}".format(data['fep']))
                sys.exit()
        
        # Create wd
        if '/' in self.environment['top']:
            toproot = self.environment['top'].split('/')[-1]
            
        else:
            toproot = self.environment['top']
        self.environment['wd'] = self.environment['wd']  + '/' + toproot.split('.')[0]   
        
        # INIT
        Create_Environment(top = self.environment['top'],
                           wd  = self.environment['wd'],
                        )

        Prepare_Topology(top = self.environment['top'],
                         wd  = self.environment['wd'],
                        )
        
        Prepare_MD(top = self.environment['top'],
                   wd  = self.environment['wd'],
                   md  = self.environment['md'],
                  )
        
        # Only to run when fep file is specified?
        Prepare_FEP(fepfile = self.environment['fep'],
                    wd  = self.environment['wd'],
                    top = self.environment['top']                 
                   )
        
        Read_Restart(restart = self.environment['re'],
                         wd  = self.environment['wd'],
                        top  = self.environment['top'],
                   )
        
        Run_Dynamics(wd  = self.environment['wd'],
                     top = self.environment['top'],
                     core = self.environment['core'],
                     clean = self.environment['clean'],
                     md = self.environment['md'],
                     fep = self.environment['fep'],
                     verbose = self.environment['verbose'])
