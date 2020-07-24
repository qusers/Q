import glob
import numpy as np
import argparse
import os
import sys
import itertools
from os import path
import shutil

import IO
import functions as f
import potential_energy as Upot
import geometries
import topology         
import defaults

class MD(object):
    def __init__(self):
        """ This class contains properties for the MD run """
        MD.stages       = []         # 'filenames' of regular Q files
        MD.MD           = {
                            'steps':[],
                            'stepsize':[],
                            'temperature':[],
                            'random_seed':[],
                            'initial_temperature':[],
                            'shake_solvent':[],
                            'shake_hydrogens':[],
                            'shake_solute':[],
                            'lrf':[],
                            'bath_coupling':[],
                           }
        
        MD.cutoffs      = {
                            'solute_solvent':[],
                            'solute_solute':[],
                            'solvent_solvent':[],
                            'q_atom':[],
                            'lrf':[],
                          }
        
        MD.sphere       = {
                            'shell_force':[],
                            'shell_radius':[],            
                          }
        
        MD.solvent      = {
                            'radial_force':[],
                            'polarisation':[],
                            'polarisation_force':[],
                          }
        
        MD.intervals    = {
                            'output':[],
                            'trajectory':[],
                            'non_bond':[],
                            'energy':[],
        }
        
        MD.files        = {
                            'topology':[],
                            'trajectory':[],
                            'restart':[],
                            'final':[],
                            'fep':[],
                            'energy':[],            
        }
        
        MD.trajectory   = {
                            'trajectory_atoms':[],
                          }
        
        MD.lambdas      = {
                            'lambdas':[],
                          }
        
        MD.seqrest      = {
                            'seqrest':[],
                          }
        
        MD.posrest      = {
                            'posrest':[],
                          }
        
        MD.distrest      = {
                            'distrest':[],
                          }
    
class Prepare_Topology(object):
    def __init__(self,top):
        self.top = top
        
        # Read the stuff
        self.read_topology()
        self.write_topology()
        
    def read_topology(self):
        top = topology.Read_Topology(self.top)
        self.topology = top.parse_topology()
        
    def write_topology(self):
        write_top = topology.Write_Topology(self.top)
        write_top.write_csv()

class Run_Dynamics(object):
    def __init__(self,top,inp,qdir):
        self.top = top
        self.inp = inp
        self.qdir = qdir
        self.md = MD()
        
        # Run stuff
        if self.qdir != None:
            print("reading Q input files in {}".format(self.qdir))
            self.read_q_inputs()
            
        # write the inputs for the C code, based on the MD class
        self.construct_inputs()
        
    def read_q_inputs(self):
        """ This function reads Q inputfiles from Q-FEP modules"""
        qfiles = []
        block = 0
        stage_ref = -1
        
        # get EQ files
        eqfiles = glob.glob(self.qdir + '/eq*.inp')            
        eqfiles = sorted(eqfiles)
        
        # get MD files
        mdfiles = glob.glob(self.qdir + '/md*.inp')
        mdfiles = list(reversed(sorted(mdfiles)))
        
        # create one list
        qfiles = eqfiles + mdfiles
        
        # get the number of files, and translate to stages
        for qfile in qfiles:
            stage = qfile.split('/')[-1][:-4]
            self.md.stages.append(stage)
            
        # populate the values in the MD topology with empty values
        # this makes sure that undefinied values are skipped but
        # length of array maintained
        # WOULD BE GREAT TO JUST LOOP OVER THESE OBJECTS....!!!
        for i, stage in enumerate(self.md.stages):
            # MD class
            for key in self.md.MD:
                if key in defaults.MD:
                    self.md.MD[key].append(defaults.MD[key])
                    
                else:
                    self.md.MD[key].append(None)
            
            # cutoffs
            for key in self.md.cutoffs:
                if key in defaults.cutoffs:
                    self.md.cutoffs[key].append(defaults.cutoffs[key])
                    
                else:
                    self.md.cutoffs[key].append(None)
            
            # sphere
            for key in self.md.sphere:
                if key in defaults.sphere:
                    self.md.sphere[key].append(defaults.sphere[key])
                    
                else:
                    self.md.sphere[key].append(None)
        
            # solvent
            for key in self.md.solvent:
                if key in defaults.solvent:                
                    self.md.solvent[key].append(defaults.solvent[key])
                    
                else:
                    self.md.solvent[key].append(None)                         
                            
            # intervals
            for key in self.md.intervals:
                if key in defaults.intervals:                
                    self.md.intervals[key].append(defaults.intervals[key])
                    
                else:
                    self.md.intervals[key].append(None)                         

            # md files !! No default parameters
            for key in self.md.files:
                self.md.files[key].append(None)
                
            # trajectory atoms !! No default parameters
            for key in self.md.files:
                self.md.files[key].append(None)
                            
            # trajectory atoms !! No default parameters
            for key in self.md.trajectory:
                self.md.trajectory[key].append(None)
            
            # trajectory atoms !! No default parameters
            for key in self.md.lambdas:
                self.md.lambdas[key].append(None)
                            
            # trajectory atoms !! No default parameters
            for key in self.md.seqrest:
                self.md.seqrest[key].append([])
                                            
            # trajectory atoms !! No default parameters
            for key in self.md.posrest:
                self.md.posrest[key].append([])
                                                            
            # trajectory atoms !! No default parameters
            for key in self.md.distrest:
                self.md.distrest[key].append([])
                
        # now iterate over all files and populate the MD object
        for qfile in qfiles:
            stage_ref += 1
            with open(qfile) as infile:
                for line in infile:
                    if len(line) < 2:
                        continue
                        
                    if '[MD]' in line:
                        block = 1
                        continue
                        
                    if '[cut-offs]' in line:
                        block = 2
                        continue
                            
                    if '[sphere]' in line:
                        block = 3
                        continue
                            
                    if '[solvent]' in line:
                        block = 4
                        continue
                            
                    if '[intervals]' in line:
                        block = 5
                        continue
                            
                    if '[files]' in line:
                        block = 6
                        continue
                                
                    if '[trajectory_atoms]' in line:
                        block = 7
                        continue
                                
                    if '[lambdas]' in line:
                        block = 8
                        continue
                                
                    if '[sequence_restraints]' in line:
                        block = 9
                        continue

                    if '[positional_restraints]' in line:
                        block = 10
                        continue
                        
                    if '[distance_restraints]' in line:
                        block = 11
                        continue                        
                        
                    if block == 1:
                        line = line.split()
                        self.md.MD[line[0]][stage_ref] = line[1]
                            
                    if block == 2:
                        line = line.split()
                        self.md.cutoffs[line[0]][stage_ref] = line[1]
                                
                    if block == 3:
                        line = line.split()
                        self.md.sphere[line[0]][stage_ref] = line[1]
                                    
                    if block == 4:
                        line = line.split()
                        self.md.solvent[line[0]][stage_ref] = line[1]
                                    
                    if block == 5:
                        line = line.split()
                        self.md.intervals[line[0]][stage_ref] = line[1]
                                        
                    if block == 6:
                        line = line.split()
                        self.md.files[line[0]][stage_ref] = line[1]
                                            
                    if block == 7:
                        self.md.trajectory['trajectory_atoms'][stage_ref] = line.strip()
                                                
                    if block == 8:
                        self.md.lambdas['lambdas'][stage_ref] = (line.split())
                                                    
                    if block == 9:
                        self.md.seqrest['seqrest'][stage_ref].append((line.split()))
                                                        
                    if block == 10:
                        self.md.posrest['posrest'][stage_ref].append((line.split()))
                                                        
                    if block == 11:
                        self.md.distrest['distrest'][stage_ref].append((line.split()))
    
    def construct_inputs(self):
        for key in MD.distrest:
            print(key, MD.distrest[key])
        return None
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='Py-MD',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '       == Python based MD engine == ')

    
    parser.add_argument('-t', '--top',
                        dest = "top",
                        default = None,
                        required = True,
                        help = "Q topology file")
    
    parser.add_argument('-i', '--inp',
                        dest = "inp",
                        default = None,
                        required = False,                        
                        help = "Python type MD input file")
        
    parser.add_argument('-d', '--qdir',
                        dest = "qdir",
                        default = None,
                        required = False,                        
                        help = "Directory with standard Q inputfiles")
    
    args = parser.parse_args()
    Prepare_Topology(top = args.top)
    Run_Dynamics(top = args.top,
                 inp = args.inp,
                 qdir = args.qdir
                )