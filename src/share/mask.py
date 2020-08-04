# Standard Python libraries
import os
import itertools
from os import path
import json

# Q-GPU libraries
import IO

class Mask(object):
    def __init__(self):
        self.data = {
                     'mask'     :   {},
                     'volume'   :   {},
                     'center'   :   {}
                    }
        
class Read_Library(object):
    def __init__(self,ilib):
        self.prm2lib = {}
        self.lib = ilib
        
    def Q(self):
        self.lib = IO.read_lib(self.lib)
        for resname in self.lib:
            if not resname in self.prm2lib:
                self.prm2lib[resname] = {}
                
            for line in self.lib[resname]['[atoms]']:
                self.prm2lib[resname][int(line[0])] = line[1]
                
        return(self.prm2lib)

class Construct_Mask(object):
    def __init__(self,top,prm2lib):
        self.top = top
        self.prm2lib = prm2lib
        
        data = Mask()
        self.data = data.data
        
    def Make(self):
        self.data['volume'] = float(self.top['radii'])
        self.data['center'] = self.top['solvcenter']
        
        resdic = {}
        attypes = {}
        resno = 0
        
        for i in range(0,len(self.top['residues'])):
            resdic[int(self.top['residues'][i])] = self.top['sequence'][i]
            
        for i in range(0,len(self.top['atypes'])):
            attypes[self.top['atypes'][i][0]] = self.top['atypes'][i][1]
            
        for ai in range(0,len(self.top['atypes'])):
            if ai + 1 in resdic:
                resname = resdic[ai + 1]
                resno += 1
                idex = ai
            #atname = self.topology.anames[attypes[ai + 1] - 1]
            atname = self.prm2lib[resname][ai - idex + 1]

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

            self.data['mask'][ai] = pdb
        
        return(self.data)
    
    def write_JSON(self,wd):
        """
        .json MD input file
        """
        with open(wd + '/mask.json', 'w') as outfile:
            inputs = self.data
            json.dump(inputs,outfile,indent=2)
            
class Read_Mask(object):
    def __init__(self, mask):
        self.mask = mask
        
    def JSON(self):
        with open(self.mask) as json_file:
            self.data = json.load(json_file)
        
        # convert keys to integers
        self.data['mask'] = {int(key):self.data['mask'][key] for key in self.data['mask']}

        return(self.data)