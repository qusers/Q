import glob
import numpy as np
import argparse
import os
import sys
import itertools
from os import path
import shutil
import json

import IO
import functions as f
import potential_energy as Upot
import geometries
import topology         
import defaults

class MD(object):
    def __init__(self):
        """ This class contains properties for the MD run """
        MD.MD = {  'stages':[],
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
                'solute_solvent':[],
                'solute_solute':[],
                'solvent_solvent':[],
                'q_atom':[],
                'lrf':[],
                'shell_force':[],
                'shell_radius':[],
                'radial_force':[],
                'polarisation':[],
                'polarisation_force':[],
                'output':[],
                'trajectory':[],
                'non_bond':[],
                'energy':[],
                'topology':[],
                'trajectory':[],
                'restart':[],
                'final':[],
                'fep':[],
                'energy':[],
                'trajectory_atoms':[],
                'lambdas':[],
                'seqrest':[],
                'posrest':[],
                'distrest':[],
                'thermostat':[],
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
    def __init__(self,top,inp,qdir,json):
        self.top = top
        self.inp = inp
        self.qdir = qdir
        self.md = MD()

        self.write_json = json
        
        # Run stuff
        if self.qdir != None:
            print("reading Q input files in {}".format(self.qdir))
            self.read_q_inputs()
            
        if self.write_json == True:
            self.write_MD()
            
        # write the inputs for the C code, based on the MD class
        self.construct_inputs()
        self.write_csv(0)
        
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
            self.md.MD['stages'].append(stage)
            
        # populate the values in the MD topology with empty values
        # this makes sure that undefinied values are skipped but
        # length of array maintained
        # WOULD BE GREAT TO JUST LOOP OVER THESE OBJECTS....!!!
        restraints = ['seqrest',
                      'posrest',
                      'distrest'
                     ]
        
        for key in self.md.MD:
            for i in range(0, len(self.md.MD['stages'])):
                if key == 'stages':
                    continue

                if key in restraints:
                    self.md.MD[key].append([])
                    continue    

                if key in defaults.MD:
                    self.md.MD[key].append(defaults.MD[key])

                else:
                    self.md.MD[key].append(None)

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
                        self.md.MD[line[0]][stage_ref] = line[1]
                                
                    if block == 3:
                        line = line.split()
                        self.md.MD[line[0]][stage_ref] = line[1]
                                    
                    if block == 4:
                        line = line.split()
                        self.md.MD[line[0]][stage_ref] = line[1]
                                    
                    if block == 5:
                        line = line.split()
                        self.md.MD[line[0]][stage_ref] = line[1]
                                        
                    if block == 6:
                        continue
                        #line = line.split()
                        #self.md.MD[line[0]][stage_ref] = line[1]
                                            
                    if block == 7:
                        self.md.MD['trajectory_atoms'][stage_ref] = line.strip()
                                                
                    if block == 8:
                        self.md.MD['lambdas'][stage_ref] = (line.split())
                                                    
                    if block == 9:
                        if len(self.md.MD['seqrest']) == 0:
                            self.md.MD['seqrest'][stage_ref] = [line.split()]
                        
                        else:
                            self.md.MD['seqrest'][stage_ref].append((line.split()))
                                                                                     
                    if block == 10:
                        if len(self.md.MD['posrest']) == 0:
                            self.md.MD['posrest'][stage_ref] = [line.split()]
                        
                        else:
                            self.md.MD['posrest'][stage_ref].append((line.split()))
                                                                                                                      
                    if block == 11:
                        if len(self.md.MD['distrest']) == 0:
                            self.md.MD['distrest'][stage_ref] = [line.split()]
                        
                        else:
                            self.md.MD['distrest'][stage_ref].append((line.split()))

    def construct_inputs(self):
        #for key in MD.distrest:
        #    continue
            #print(key, MD.distrest[key])
        return None
    
    def write_MD(self):
        print("Writing out json files of the MD inputs")
        with open('test.txt', 'w') as outfile:
            json.dump(self.md.__dict__,outfile)
            
    def write_csv(self,i):
        """
        needs to loop over stages in MD object (for now gets integer 
        of stage)
        
        """
        j = 0
        
        for key in self.md.MD:
            if type(self.md.MD[key][i]) != list:
                j += 1
                
            else:
                j += len(self.md.MD[key][i])
        
        
        # semi-hardcoded, needs fix!
        with open(self.top[:-4] + '/md.csv', 'w') as outfile:
            outfile.write('{}\n'.format(j))
            outfile.write('steps;{}\n'.format(self.md.MD['steps'][i]))        
            outfile.write('stepsize;{}\n'.format(self.md.MD['stepsize'][i]))        
            outfile.write('temperature;{}\n'.format(self.md.MD['temperature'][i]))        
            outfile.write('thermostat;{}\n'.format(self.md.MD['thermostat'][i]))        
            outfile.write('bath_coupling;{}\n'.format(self.md.MD['bath_coupling'][i]))        
            outfile.write('random_seed;{}\n'.format(self.md.MD['random_seed'][i]))        
            outfile.write('initial_temperature;{}\n'.format(self.md.MD['initial_temperature'][i]))        
            outfile.write('shake_solvent;{}\n'.format(self.md.MD['shake_solvent'][i]))        
            outfile.write('shake_hydrogens;{}\n'.format(self.md.MD['shake_hydrogens'][i]))        
            outfile.write('lrf;{}\n'.format(self.md.MD['lrf'][i]))        
            outfile.write('solute_solute;{}\n'.format(self.md.MD['solute_solute'][i]))        
            outfile.write('solvent_solvent;{}\n'.format(self.md.MD['solvent_solvent'][i]))        
            outfile.write('solute_solvent;{}\n'.format(self.md.MD['solute_solvent'][i]))        
            outfile.write('q_atom;{}\n'.format(self.md.MD['q_atom'][i]))        
            outfile.write('shell_radius;{}\n'.format(self.md.MD['shell_radius'][i]))        
            outfile.write('shell_force;{}\n'.format(self.md.MD['shell_force'][i]))        
            outfile.write('radial_force;{}\n'.format(self.md.MD['radial_force'][i]))        
            outfile.write('polarisation;{}\n'.format(self.md.MD['polarisation'][i]))        
            outfile.write('polarisation_force;{}\n'.format(self.md.MD['polarisation_force'][i]))        
            outfile.write('non_bond;{}\n'.format(self.md.MD['non_bond'][i]))        
            outfile.write('output;{}\n'.format(self.md.MD['output'][i]))        
            outfile.write('energy;{}\n'.format(self.md.MD['energy'][i]))        
            outfile.write('trajectory;{}\n'.format(self.md.MD['trajectory'][i]))        
    
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
            
    parser.add_argument('-wp', '--wpython',
                        dest = "wpython",
                        default = False,
                        action = 'store_true',                        
                        help = "Toggle to write out python readable json files of the md")
    
    args = parser.parse_args()
    Prepare_Topology(top = args.top)
    Run_Dynamics(top = args.top,
                 inp = args.inp,
                 qdir = args.qdir,
                 json = args.wpython
                )