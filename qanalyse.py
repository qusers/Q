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
import qdyn
import topology

class Trajectory(object):
    """ Trajectory object class, one nparray per timeframe"""
    def __init__(self):
        Trajectory.frames = []
        Trajectory.coords = []
        Trajectory.natoms = []
        Trajectory.maskPDB = {}
        Trajectory.volume = None
        Trajectory.center = None

class Mapping(object):
    """ Mapping of atom names to atom types"""
    def __init__(self,ilib,top):
        self.ilib = ilib
        self.top = top
        Mapping.prm2lib = {}
        self.trajectory = Trajectory()
        
        # run stuff
        self.read_lib()
        self.read_topology()
        
    def read_lib(self):
        self.lib = IO.read_lib(self.ilib)
    
        # here, generate the lib/prm atom mapping, add the charges to the atom properties in prm
        for resname in self.lib:
            if not resname in Mapping.prm2lib:
                Mapping.prm2lib[resname] = {}
                
            for line in self.lib[resname]['[atoms]']:
                Mapping.prm2lib[resname][int(line[0])] = line[1]
                
    def read_topology(self):
        top = topology.Read_Topology(self.top)
        self.topology = top.parse_topology()
        Trajectory.volume = float(self.topology.radii)
        Trajectory.center = self.topology.solvcenter
        
        resdic = {}
        attypes = {}
        resno = 0
        
        for i in range(0,len(self.topology.residues)):
            resdic[int(self.topology.residues[i])] = self.topology.sequence[i]
            
        for i in range(0,len(self.topology.atypes)):
            attypes[self.topology.atypes[i][0]] = self.topology.atypes[i][1]
            
        for ai in range(0,len(self.topology.atypes)):
            if ai + 1 in resdic:
                resname = resdic[ai + 1]
                resno += 1
                idex = ai
            #atname = self.topology.anames[attypes[ai + 1] - 1]
            atname = Mapping.prm2lib[resname][ai - idex + 1]

            # construct the .pdb matrix
            pdb = ['HETATM',
                   ai + 1,
                   atname,' ',
                   resname,' ',
                   resno,' ',
                   0.0,
                   0.0,
                   0.0,
                   0.0,
                   0.0,
                   ' ',
                   ' '
                  ]

            Trajectory.maskPDB[ai] = pdb

class Read_Trajectory(object):
    def __init__(self,itrj):
        print("Reading trajectory file from Qdyn")
        self.itrj = itrj
        
        # Run some stuff
        self.readfile()
    
    def readfile(self):
        with open (self.itrj) as infile:
            Trajectory.natoms = infile.readline().strip()
            for line in infile:
                #if line == self.natoms
                line = line.strip()
                line = line.split(';')
                # Frame bookkeeping
                if len(line) == 1:
                    Trajectory.frames.append(line)
                
                # Add
                if len(line) == 3:
                    Trajectory.coords.append([float(line[0]),
                                             float(line[1]),
                                             float(line[2])]                                            
                                            )

class Write_Trajectory(object):
    def __init__(self,otrj,top):
        print("Writing trajectory file from Qdyn, for now in .pdb format")
        self.otrj = otrj
        self.top = top
    
        # Run stuff
        self.create_environment()
        self.read_topology()
        self.write_trajectory_pdb()
        
    def create_environment(self):
        self.workdir = os.getcwd()

        self.trajdir = self.workdir + '/' + self.otrj
        if path.exists(self.trajdir) == True:
            os.system('rm -r ' + self.trajdir)
        
        os.mkdir(self.trajdir)
        
    def read_topology(self):
        top = topology.Read_Topology(self.top)
        self.topology = top.parse_topology()
        
    def write_trajectory_pdb(self):
        files = []
        for frame in Trajectory.frames:
            frame = int(frame[0])
            filename = '{:010d}.pdb'.format(frame)
            files.append(filename)
            
            with open(self.trajdir + '/' + filename,'w') as outfile:
                for ai in range(0,int(Trajectory.natoms)):
                    #Trajectory.maskPDB[ai][8] = 0.0
                    #Trajectory.maskPDB[ai][9] =  0.0
                    #Trajectory.maskPDB[ai][10] =  0.0                    
                    # Update the coordinates
                    coord_index = (frame) * int(Trajectory.natoms) + ai
                    coords = Trajectory.coords[coord_index]
                    Trajectory.maskPDB[ai][8] = coords[0]
                    Trajectory.maskPDB[ai][9] = coords[1]
                    Trajectory.maskPDB[ai][10] = coords[2]

                    line = IO.pdb_parse_out(Trajectory.maskPDB[ai])
                    outfile.write(line + '\n')
                    
        with open(self.trajdir + '/load_traj.pml', 'w') as outfile:
            init = '000000.pdb'
            outfile.write('load ' + init + '\n')
            
            for line in files:
                if line == init:
                    continue
                
                outfile.write('load ' + line + ', ' + init + '\n')
                
                
class Calculations(object):
    def __init__(self, calculations, top):
        self.top = top
        self.calcuations = calculations

        # perform calculations
        for calculation in calculations:
            if calculation == 'number_density':
                self.number_density()
            
    def number_density(self):
        print("calculating number densities")
        parameters = {
                      'bins'    : None,
                      'residue' : []
                     }

        center = [float(Trajectory.center[0]),
                  float(Trajectory.center[1]),
                  float(Trajectory.center[2])
                 ]
        # read in parameters from file
        # probably needs some dir logic
        with open('number_density.calc') as infile:
            for line in infile:
                line = line.split()
                if len(line) < 1:
                    continue
                if not line[0] in parameters:
                    print("FATAL: parameter {} not found in inputfile".format(line[0]))
                    sys.exit()
                    
                else:
                    if type(parameters[line[0]]) == list:
                        parameters[line[0]].append(line[1])
                        
                    else:
                        parameters[line[0]] = int(line[1])
        
        # Now we need calculate whether the atom falls within a bin.
        # First construct the bins and calculate the volumes
        steps = Trajectory.volume/float(parameters['bins'])
        maxbin = int(round(round(steps, 2)))
        bins = {}
        empty = []
        
        for i in range(0,len(Trajectory.frames)):
            empty.append(0)

        for i in range(0,maxbin):
            r = (i + 1) * parameters['bins']
            r0 = i * parameters['bins']
            
            # the last bin can be smaller
            if i + 1 == maxbin:
                if r > float(Trajectory.volume):
                    r = r - (r -Trajectory.volume)
            
            bins[i] = [r0,r,empty]
        
        # Now loop over the frames
        for frame in Trajectory.frames:
            frame = int(frame[0])
            for b in bins:    
                for ai in Trajectory.maskPDB:
                    if Trajectory.maskPDB[ai][4] == 'HOH':# and Trajectory.maskPDB[ai][2] == 'O':
                        coord_index = (frame) * int(Trajectory.natoms) + ai
                        if f.euclidian_range(center,Trajectory.coords[coord_index],bins[b][0],bins[b][1]) == True:
                            bins[b][2][frame] += 1
        
        print(bins)
                    
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
        
    parser.add_argument('-ot', '--otrj',
                        dest = "otrj",
                        default = None,
                        help = " Output trajectory file")
            
    parser.add_argument('-it', '--itrj',
                        dest = "itrj",
                        default = None,
                        help = " Input trajectory file")
                
    parser.add_argument('-il', '--ilib',
                        dest = "ilib",
                        default = None,
                        help = "Library files for used topology")
    
    parser.add_argument("-wt", "--wtraj",
                        default = False,
                        action = 'store_true',
                        help="Write trajecotry")
    
    parser.add_argument('-c', '--calc',
                        dest = "calc",
                        required = False,
                        #choices = ['number_density','test'],
                        help = "Comma seperated list of calculations to be performed, requires a *.calc inputfile")
    
    args = parser.parse_args()
    
    # Maybe we want specific inputs here (definitely some input file later)
    if args.itrj != None and args.ilib != None:
        Mapping(ilib = args.ilib,
                top  = args.top)
        Read_Trajectory(itrj = args.itrj)    
    
    if args.wtraj == True and args.otrj != None:
        Write_Trajectory(otrj = args.otrj,
                         top = args.top)
        
    if args.calc != None:
        calcs = args.calc.split(',')
        for calc in calcs:
            if not os.path.exists(calc + '.calc'):
                print("FATAL: could not find input file for {}.calc".format(calc))
                sys.exit()
                
        Calculations(calcs,top = args.top)
