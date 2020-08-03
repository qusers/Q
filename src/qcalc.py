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

import IO
import topology as TOPOLOGY        
import defaults as DEFAULTS
import md       as MD
import settings as SETTINGS
import fep      as FEP

class Mask(object):
    data = {'trajectory': {},
            'mapping': None,
            'topology': None,
           }
        
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

class Prepare_Mapping(object):
    def __init__(self,top,ilib):
        self.top = top
        self.ilib = ilib
                
        # Read topology        
        read_top = TOPOLOGY.Read_Topology(self.top)
        Mask.data['topology'] = read_top.JSON()
                
        mapping  = TOPOLOGY.Mapping(self.ilib,Mask.data['topology'])
        Mask.data['mapping'] = mapping.create_mask()
        
        Mask.data['trajectory']['volume'] = float(Mask.data['topology']['radii'])
        Mask.data['trajectory']['center'] = Mask.data['topology']['solvcenter']
        
class Init(object):
    def __init__(self, data):
        """ Retrieves a dictionary of user input from qalc.py:
               { 'top'  : top,
                 'wd'   : wd,
                 'itrj' : itrj,
                 'ilib' : ilib,
                 'wraj' : wtraj,
                 'calc' : calc,
               }     
        """     
        self.environment = data
        
        # Create user specified work directory
        Create_Environment(wd = self.environment['wd'])
        
        # Create the topology object for reference
        Prepare_Mapping(top = self.environment['top'],
                        ilib  = self.environment['ilib'],
                       )
        
        print(Mask.data['mapping'])
    # Maybe we want specific inputs here (definitely some input file later)
    #if args.itrj != None and args.ilib != None:
    #    Mapping(ilib = args.ilib,
    #            top  = args.top)
    #    Read_Trajectory(itrj = args.itrj)    
    
    #if args.wtraj == True and args.otrj != None:
    #    Write_Trajectory(otrj = args.otrj,
    #                     top = args.top)
        
    #if args.calc != None:
    #    calcs = args.calc.split(',')
    #    for calc in calcs:
    #        if not os.path.exists(calc + '.calc'):
    #            print("FATAL: could not find input file for {}.calc".format(calc))
    #            sys.exit()
                
    #    Calculations(calcs,top = args.top)
