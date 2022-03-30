# Standard Python libraries
import os
import itertools
from os import path
import json

# Q-GPU libraries
import IO

class Restart(object):
    """ Trajectory object class, one nparray per timeframe"""
    def __init__(self):
        self.data = { 'frames'     : [],
                      'coords'     : [],
                      'velocities' : [],
                      'natoms'     : None
        }        
        
class Read_Restart(object):
    def __init__(self,restart):
        data = Restart()
        self.data = data.data
        self.restart = restart
        
        if self.restart != None:
            if not(self.restart)[-1] == '/':
                self.restart = self.restart + '/'
        
    def CSV(self):
        # read velocities
        with open (self.restart + 'velocities.csv') as infile:
            self.data['natoms'] = infile.readline().strip()
            for line in infile:
                line = line.strip()
                line = line.split(';')
                # Frame bookkeeping
                if len(line) == 1:
                    frame = int(line[0])
                    self.data['frames'].append(frame)
                
                # Add
                if len(line) == 3:
                    self.data['velocities'].append([float(line[0]),
                                                    float(line[1]),
                                                    float(line[2])]                                            
                                                   )
        # read coordinates
        with open (self.restart + '/coords.csv') as infile:
            for line in infile:
                line = line.strip()
                line = line.split(';')
                if len(line) == 3:
                    self.data['coords'].append([float(line[0]),
                                                float(line[1]),
                                                float(line[2])]                                            
                                               )
        
        last_velo = []
        last_coord = []
        offset = (self.data['frames'][-1]) * \
                  int(self.data['natoms']) 
                
        for ai in range(0,int(self.data['natoms'])):
            ai = offset + ai 
            last_velo.append((self.data['velocities'][ai]))
            last_coord.append((self.data['coords'][ai]))
        
        return(last_velo, last_coord)
    
    def JSON(self):
        with open(self.itrj) as json_file:
            self.data = json.load(json_file)
        
        return(self.data)
    
class Write_Restart(object):
    def __init__(self,velocities,coordinates):
        self.velocities = velocities
        self.coordinates = coordinates
        
    def JSON(self):
        with open(self.wd + '/' + self.top + '/traj.json', 'w') as outfile:
            inputs = self.traj
            json.dump(inputs,outfile,indent=2)     
            
    def CSV(self,wd):
        with open(wd + 'i_coords.csv','w') as outfile:
            outfile.write('{}\n'.format(len(self.coordinates)))
            for line in self.coordinates:
                outfile.write('{};{};{}\n'.format(line[0],line[1],line[2]))
        
        with open(wd + '/i_velocities.csv','w') as outfile:
            outfile.write('{}\n'.format(len(self.velocities)))
            for line in self.velocities:
                outfile.write('{};{};{}\n'.format(line[0],line[1],line[2]))
        