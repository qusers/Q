# Standard Python libraries
import os
import itertools
from os import path
import json

# Q-GPU libraries
import IO

class Fep():
    def __init__(self):
        self.data = {
                        'q_atoms'                   : None,
                        'q_atypes'                  : [],
                        'q_catypes'                 : None,
                        'q_charges'                 : [],
                        'q_softpairs'               : [],
                        'q_exclpairs'               : [],
                        'q_elscales'                : None,
                        'q_softcores'               : [],
                        'q_bonds'                   : [],
                        'q_cbonds'                  : None,
                        'q_angles'                  : [],
                        'q_cangles'                 : None,
                        'q_torsions'                : [],
                        'q_ctorsions'               : None,
                        'q_impropers'               : [],
                        'q_cimpropers'              : None,
                        'q_angcouples'              : None,
                        'q_torcouples'              : None,
                        'q_imprcouples'             : None,
                        'q_shakes'                  : [],
                        'q_offdiags'                : None,
            
                         # Not implemented in qdyn 
                        'q_states'                  : None,
                        'qquselibrarycharges'       : None,
                        'softcoreusemaxpotential'   : None,
                        'q_monitor_groups'          : None,
                        'q_monitor_group_pairss'    : None,
                        'q_atypemap'                : {},
                    }

class Read_Fep(object):
    """
    Read Q topology file as an input, parse to topology class
    """
    def __init__(self, fepfile, *args, **kwargs):    
        self.fepfile = fepfile

        data = Fep()
        self.data = data.data
    
    def EMPTY(self):
        return(self.data)
        
    def JSON(self):
        with open(self.fepfile) as json_file:
            self.data = json.load(json_file)
        
        return(self.data)
        
    def Q(self):
        block = 0
        # get the number of states
        with open(self.fepfile) as infile:
            for line in infile:
                if len(line.split()) < 1:
                    continue
                line = line.split()    
                if line[0] == 'states':
                        self.data['states'] = int(line[1])    
                        
        # construct states for keywords:
        for key in self.data:
            if type(self.data[key]) == list:
                for i in range(0,self.data['states']):
                    self.data[key].append([])
        
        with open(self.fepfile) as infile:
            for line in infile:
                if len(line.split()) < 1:
                    continue
                
                if '[atoms]' in line:
                    self.data['q_atoms'] = []
                    block = 1
                    continue
            
                if '[FEP]' in line:
                    block = 2
                    continue
                                
                if '[change_charges]' in line:
                    block = 3
                    continue
                                
                if '[atom_types]' in line:
                    self.data['q_catypes'] = []
                    idex = 0
                    block = 4
                    continue
                                
                if '[change_atoms]' in line:
                    block = 5
                    continue
                                
                if '[soft_pairs]' in line:
                    block = 6
                    continue
                                
                if '[excluded_pairs]' in line:
                    block = 7
                    continue
                                
                if '[el_scale]' in line:
                    self.data['q_elscales'] = []                 
                    block = 8
                    continue
                                                    
                if '[softcore]' in line:
                    block = 9
                    continue
                                                    
                if '[monitor_groups]' in line:
                    self.data['q_monitor_groups'] = []                 
                    block = 10
                    continue

                if '[monitor_groups_pairs]' in line:
                    self.data['q_monitor_group_pairs'] = []                         
                    block = 11
                    continue
                                                    
                if '[bond_types]' in line:
                    self.data['q_cbonds'] = [['0','0.0','0.0']]                   
                    block = 12
                    continue
                                                    
                if '[change_bonds]' in line:
                    block = 13
                    continue
                                                    
                if '[angle_types]' in line:
                    self.data['q_cangles'] = [['0','0.0','0.0']]                   
                    block = 14
                    continue
                                                                        
                if '[change_types]' in line:
                    block = 15
                    continue
                                                                        
                if '[torsion_types]' in line:
                    self.data['q_ctorsions'] = [['0','0.0','0.0','0.0']]                   
                    block = 16
                    continue
                                                                        
                if '[change_torsions]' in line:
                    block = 17
                    continue
                                                                        
                if '[improper_types]' in line:
                    self.data['q_cimpropers'] = [['0','0.0','0.0']]                                       
                    block = 18
                    continue
                                                                                          
                if '[change_impropers]' in line:
                    block = 19
                    continue
                                                                                                            
                if '[angle_couplings]' in line:
                    self.data['q_angcouples'] = []                                                       
                    block = 20
                    continue
                                                                                                            
                if '[torsion_couplings]' in line:
                    self.data['q_torcouples'] = []                                                       
                    block = 21
                    continue
                                                                                                            
                if '[improper_couplings]' in line:
                    self.data['q_imprcouples'] = []                                                       
                    block = 22
                    continue
                                                                                                            
                if '[shake_constraints]' in line:
                    block = 23
                    continue
                                                                                                                                
                if '[off-diagonals]' in line:
                    block = 24
                    continue
                    
                    
                # Read stuff
                if block == 0:
                    continue

                if block == 1:
                    line = line.split()
                    self.data['q_atoms'].append(int(line[1]))

                if block == 2:
                    line = line.split()
                    if line[0] == 'states':
                        self.data['q_states'] = int(line[1])                    
                        
                    if line[0] == 'offset':
                        offset = int(line[1])
                        for line in self.data['q_atoms']:
                            self.data['q_atoms'] = self.data['q_atoms'] + offset
                        
                    if line[0] == 'qq_use_library_charges':
                        self.data['qquselibrarycharges'] = line[1]
                                                
                    if line[0] == 'softcore_use_maxpotential':
                        self.data['softcoreusemaxpotential'] = line[1]
                        
                if block == 3:
                    line = line.split()
                    for i in range(0,self.data['states']):
                        self.data['q_charges'][i].append(line[i+1])
                                        
                if block == 4:
                    idex += 1
                    line = line.split()
                    self.data['q_catypes'].append(line)
                    self.data['q_atypemap'][line[0]] = idex
                                                            
                if block == 5:
                    line = line.split()
                    for i in range(0,self.data['states']):
                        atype = self.data['q_atypemap'][line[i+1]]
                        self.data['q_atypes'][i].append(atype)
                    
                if block == 6:
                    line = line.split()   
                    for i in range(0,self.data['states']):
                        self.data['q_softpairs'][i].append(line[i+1])
                                        
                if block == 7:
                    line = line.split()   
                    for i in range(0,self.data['states']):
                        self.data['q_softpairs'][i].append(line[i+1]) 
                        
                if block == 8:
                    line = line.split()   
                    self.data['q_elscales'].append([line])
                                            
                if block == 9:
                    line = line.split() 
                    for i in range(0,self.data['states']):
                        self.data['q_softcores'][i].append(line[i+1]) 
                                                            
                if block == 10:
                    line = line.split() 
                    self.data['q_monitor_groups'].append([line])
                                                                           
                if block == 11:
                    line = line.split() 
                    self.data['q_monitor_group_pairs'].append([line])
                    
                if block == 12:
                    line = line.split()
                    self.data['q_cbonds']    
                    
                if block == 13:
                    line = line.split()
                    for i in range(0,self.data['states']):
                        self.data['q_bonds'][i].append(line[i+1])                     
                    
                if block == 14:
                    line = line.split()
                    self.data['q_cangles'].append(line)  
                    
                if block == 15:
                    line = line.split()
                    for i in range(0,self.data['states']):
                        self.data['q_angles'][i].append(line[i+1])       
                        
                if block == 16:
                    line = line.split()
                    self.data['q_ctorsions'].append(line)  
                
                if block == 17:
                    line = line.split()
                    for i in range(0,self.data['states']):
                        self.data['q_torsions'][i].append(line[i+1])      
                        
                if block == 18:
                    line = line.split()
                    self.data['q_cimpropers'].append(line)  

                if block == 19:
                    line = line.split()
                    for i in range(0,self.data['states']):
                        self.data['q_torsions'][i].append(line[i+1])    
                        
                if block == 20:
                    line = line.split()
                    self.data['q_angcouples'].append(line)  
                                
                if block == 21:
                    line = line.split()
                    self.data['q_torcouples'].append(line)  
                                
                if block == 22:
                    line = line.split()
                    self.data['q_imprcouples'].append(line)  
        
                if block == 23:
                    line = line.split()
                    for i in range(0,self.data['states']):
                        self.data['q_shakes'][i].append(line[i+1])
                                
                if block == 24:
                    line = line.split()
                    for i in range(0,self.data['states']):
                        print("TO DO")
                        
        return(self.data)
                
class Write_Fep(object):        
    """
    Write Python topology object to file
    """
    def __init__(self, data, *args, **kwargs):    
        self.data = data

    def JSON(self,out_json):
        """
        .json MD input file
        """
        with open(out_json, 'w') as outfile:
            inputs = self.data
            json.dump(inputs,outfile,indent=2)        
        
    def CSV(self,wd):
        self.wd = wd
        
        with open(self.wd + 'q_atoms.csv','w') as outfile:
            if self.data['q_atoms'] == None:
                    outfile.write('0\n')
                    
            else:
                outfile.write('{}\n'.format(len(self.data['q_atoms'])))
                for line in self.data['q_atoms']:
                    outfile.write('{}\n'.format(line))
                    
        with open(self.wd + 'q_atypes.csv','w') as outfile:
            if len(self.data['q_atypes']) == 0:
                    outfile.write('0\n')
                    
            else:
                outfile.write('{}\n'.format(len(self.data['q_atypes'][0])*self.data['q_states']))
                for i in range(0,self.data['q_states']):
                    for j in range(0, len(self.data['q_atypes'][i])):
                        outfile.write('{}\n'.format(self.data['q_atypes'][i][j]))
                    
        with open(self.wd + 'q_catypes.csv','w') as outfile:
            if self.data['q_catypes'] == None:
                    outfile.write('0\n')
                    
            else:
                outfile.write('{}\n'.format(len(self.data['q_catypes'])))
                for line in self.data['q_catypes']:
                    outfile.write('{};{};{};{};{};{};{};{}\n'.format(line[0],
                                                                     line[1],
                                                                     line[2],
                                                                     line[3],
                                                                     line[4],
                                                                     line[5],
                                                                     line[6],
                                                                     line[7],
                                                                    ))
        
        with open(self.wd + 'q_charges.csv','w') as outfile:
            if len(self.data['q_charges']) == 0:
                    outfile.write('0\n')
                    
            else:            
                outfile.write('{}\n'.format(len(self.data['q_charges'][0])*self.data['q_states']))
                for i in range(0,self.data['q_states']):
                    for j in range(0, len(self.data['q_charges'][i])):
                        outfile.write('{}\n'.format(self.data['q_charges'][i][j]))
                
        with open(self.wd + 'q_softpairs.csv','w') as outfile:
            if len(self.data['q_softpairs']) == 0:
                    outfile.write('0\n')
                    
            else:                    
                outfile.write('{}\n'.format(len(self.data['q_softpairs'][0])*self.data['q_states']))
                for i in range(0,self.data['q_states']):
                    for j in range(0, len(self.data['q_softpairs'][i])):
                        outfile.write('{}\n'.format(self.data['q_softpairs'][i][j]))
                        
        with open(self.wd + 'q_exclpairs.csv','w') as outfile:
            if len(self.data['q_exclpairs']) == 0:
                    outfile.write('0\n')
                    
            else:                    
                outfile.write('{}\n'.format(len(self.data['q_exclpairs'][0])*self.data['q_states']))
                for i in range(0,self.data['q_states']):
                    for j in range(0, len(self.data['q_exclpairs'][i])):
                        outfile.write('{}\n'.format(self.data['q_exclpairs'][i][j]))
        
        with open(self.wd + 'q_elscales.csv','w') as outfile:
            if self.data['q_elscales'] == None:
                    outfile.write('0\n')
                    
            else:
                outfile.write('{}\n'.format(len(self.data['q_elscales'])))
                for line in self.data['q_elscales']:
                    outfile.write('{}\n'.format(line))
                            
        with open(self.wd + 'q_softcores.csv','w') as outfile:
            if len(self.data['q_softcores']) == 0:
                    outfile.write('0\n')
                    
            else:                        
                outfile.write('{}\n'.format(len(self.data['q_softcores'][0])*self.data['q_states']))
                for i in range(0,self.data['q_states']):
                    for j in range(0, len(self.data['q_softcores'][i])):
                        outfile.write('{}\n'.format(self.data['q_softcores'][i][j]))
                                                    
        with open(self.wd + 'q_bonds.csv','w') as outfile:
            if len(self.data['q_bonds']) == 0:
                    outfile.write('0\n')
                    
            else:                 
                outfile.write('{}\n'.format(len(self.data['q_bonds'][0])*self.data['q_states']))
                for i in range(0,self.data['q_states']):
                    for j in range(0, len(self.data['q_bonds'][i])):
                        outfile.write('{}\n'.format(self.data['q_bonds'][i][j]))
                    
        with open(self.wd + 'q_cbonds.csv','w') as outfile:
            if self.data['q_cbonds'] == None:
                    outfile.write('0\n')
                    
            else:
                outfile.write('{}\n'.format(len(self.data['q_cbonds'])))
                for line in self.data['q_cbonds']:
                    outfile.write('{}\n'.format(line))
                                                                            
        with open(self.wd + 'q_angles.csv','w') as outfile:
            if len(self.data['q_angles']) == 0:
                    outfile.write('0\n')
                    
            else:                 
                outfile.write('{}\n'.format(len(self.data['q_angles'][0])*self.data['q_states']))
                for i in range(0,self.data['q_states']):
                    for j in range(0, len(self.data['q_angles'][i])):
                        outfile.write('{}\n'.format(self.data['q_angles'][i][j]))
                    
        with open(self.wd + 'q_cangles.csv','w') as outfile:
            if self.data['q_cangles'] == None:
                    outfile.write('0\n')
                    
            else:
                outfile.write('{}\n'.format(len(self.data['q_cangles'])))
                for line in self.data['q_cangles']:
                    outfile.write('{}\n'.format(line))
                                                                                                    
        with open(self.wd + 'q_torsions.csv','w') as outfile:
            if len(self.data['q_torsions']) == 0:
                    outfile.write('0\n')
                    
            else:                          
                outfile.write('{}\n'.format(len(self.data['q_torsions'][0])*self.data['q_states']))
                for i in range(0,self.data['q_states']):
                    for j in range(0, len(self.data['q_torsions'][i])):
                        outfile.write('{}\n'.format(self.data['q_torsions'][i][j]))
                    
        with open(self.wd + 'q_ctorsions.csv','w') as outfile:
            if self.data['q_ctorsions'] == None:
                    outfile.write('0\n')
                    
            else:
                outfile.write('{}\n'.format(len(self.data['q_ctorsions'])))
                for line in self.data['q_ctorsions']:
                    outfile.write('{}\n'.format(line))
                                                                                                                            
        with open(self.wd + 'q_impropers.csv','w') as outfile:
            if len(self.data['q_impropers']) == 0:
                    outfile.write('0\n')
                    
            else:                      
                outfile.write('{}\n'.format(len(self.data['q_impropers'][0])*self.data['q_states']))
                for i in range(0,self.data['q_states']):
                    for j in range(0, len(self.data['q_impropers'][i])):
                        outfile.write('{}\n'.format(self.data['q_impropers'][i][j]))
                    
        with open(self.wd + 'q_cimpropers.csv','w') as outfile:
            if self.data['q_cimpropers'] == None:
                    outfile.write('0\n')
                    
            else:
                outfile.write('{}\n'.format(len(self.data['q_cimpropers'])))
                for line in self.data['q_cimpropers']:
                    outfile.write('{}\n'.format(line))
                                            
        with open(self.wd + 'q_angcouples.csv','w') as outfile:
            if self.data['q_angcouples'] == None:
                    outfile.write('0\n')
                    
            else:
                outfile.write('{}\n'.format(len(self.data['q_angcouples'])))
                for line in self.data['q_angcouples']:
                    outfile.write('{}\n'.format(line))
                                            
        with open(self.wd + 'q_torcouples.csv','w') as outfile:
            if self.data['q_torcouples'] == None:
                    outfile.write('0\n')
                    
            else:
                outfile.write('{}\n'.format(len(self.data['q_torcouples'])))
                for line in self.data['q_torcouples']:
                    outfile.write('{}\n'.format(line))
                                            
        with open(self.wd + 'q_imprcouples.csv','w') as outfile:
            if self.data['q_imprcouples'] == None:
                    outfile.write('0\n')
                    
            else:
                outfile.write('{}\n'.format(len(self.data['q_imprcouples'])))
                for line in self.data['q_imprcouples']:
                    outfile.write('{}\n'.format(line))
                                                                                                                                    
        with open(self.wd + 'q_shakes.csv','w') as outfile:
            if len(self.data['q_shakes']) == 0:
                    outfile.write('0\n')
                    
            else:                             
                outfile.write('{}\n'.format(len(self.data['q_shakes'][0])*self.data['q_states']))
                for i in range(0,self.data['q_states']):
                    for j in range(0, len(self.data['q_shakes'][i])):
                        outfile.write('{}\n'.format(self.data['q_shakes'][i][j]))
                    
                    
        with open(self.wd + 'q_offdiags.csv','w') as outfile:
            if self.data['q_offdiags'] == None:
                    outfile.write('0\n')
                    
            else:
                outfile.write('{}\n'.format(len(self.data['q_offdiags'])))
                for line in self.data['q_offdiags']:
                    outfile.write('{}\n'.format(line))