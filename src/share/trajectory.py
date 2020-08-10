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
        self.data = {
                      'frames' : [],
                      'coords' : [],
                      'natoms' : None,            
        }        
        
class Read_Trajectory(object):
    def __init__(self,itrj):
        data = Trajectory()
        self.data = data.data
        self.itrj = itrj
        
    def QDYN(self):    
        with open (self.itrj) as infile:
            self.data['natoms'] = infile.readline().strip()
            for line in infile:
                #if line == self.natoms
                line = line.strip()
                line = line.split(';')
                # Frame bookkeeping
                if len(line) == 1:
                    self.data['frames'].append(line)
                
                # Add
                if len(line) == 3:
                    self.data['coords'].append([float(line[0]),
                                                float(line[1]),
                                                float(line[2])]                                            
                                              )
                    
        return(self.data)
    
    def JSON(self):
        with open(self.itrj) as json_file:
            self.data = json.load(json_file)
        
        return(self.data)
    
class Write_Trajectory(object):
    def __init__(self,wd,otrj,mask,traj):
        self.wd = wd
        self.mask = mask
        self.traj = traj
        
        if otrj != None:
            self.otrj = self.wd + '/traj' + otrj
        
    def JSON(self):
        with open(self.wd + '/traj.json', 'w') as outfile:
            inputs = self.traj
            json.dump(inputs,outfile,indent=2)     
            
    def PDB(self):
        files = []
        for frame in self.traj['frames']:
            frame = int(frame[0])
            filename = '{:010d}.pdb'.format(frame)
            files.append(filename)
            
            with open(self.wd + '/' + filename,'w') as outfile:
                for ai in range(0,int(self.traj['natoms'])):                
                    # Update the coordinates
                    coord_index = (frame) * int(self.traj['natoms']) + ai
                    coords = self.traj['coords'][coord_index]
                    self.mask[ai][8] = coords[0]
                    self.mask[ai][9] = coords[1]
                    self.mask[ai][10] = coords[2]

                    line = IO.pdb_parse_out(self.mask[ai])
                    outfile.write(line + '\n')
                    
        with open(self.wd + '/load_traj.pml', 'w') as outfile:
            init = '000000.pdb'
            outfile.write('load ' + init + '\n')
            
            for line in files:
                if line == init:
                    continue
                
                outfile.write('load ' + line + ', ' + init + '\n')
