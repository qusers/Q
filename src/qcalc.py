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

import defaults as DEFAULTS
import md       as MD
import settings as SETTINGS
import fep      as FEP
            
class Init(object):
    def __init__(self, data):
        """ Retrieves a dictionary of user input from qalc.py:
               {'top' : top,
                 'otrj' : otrj,
                 'itrj' : itrj,
                 'ilib' : ilib,
                 'wraj' : wtraj,
                 'calc' : calc,
               }
        """
        data = { 'top' : top,
                 'otrj' : otrj,
                 'itrj' : itrj,
                 'ilib' : ilib,
                 'wraj' : wtraj,
                 'calc' : calc,
               }        
        self.environment = data
        
        print(self.environment)
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
