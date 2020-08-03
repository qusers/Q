# Standard Python libraries
import os
import itertools
from os import path
import json

# Q-GPU libraries
import IO

class Trajectory(object):
    """ Trajectory object class, one nparray per timeframe"""
    def __init__(self):
        Trajectory.frames = []
        Trajectory.coords = []
        Trajectory.natoms = []
        Trajectory.maskPDB = {}
        Trajectory.volume = None
        Trajectory.center = None
        
        
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
                