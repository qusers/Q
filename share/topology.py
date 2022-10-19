# Standard Python libraries
import os
import itertools
from os import path
import json
import sys

# Q-GPU libraries
import IO

class Topology():
    def __init__(self):
        self.data = {'header'           : {},
                     'coords'           : [],
                     'atypes'           : [],
                     'catypes'          : {},
                     'nbonds_solute'    : None,
                     'bonds'            : [],
                     'cbonds'           : {},
                     'nangles_solute'   : None,
                     'angles'           : [],
                     'cangles'          : {},
                     'ntorsions_solute' : None,
                     'torsions'         : [],
                     'ctorsions'        : {},
                     'nimpropers_solute': None,
                     'impropers'        : [],
                     'cimpropers'       : {},
                     'charges'          : [],
                     'ccharges'         : {},
                     'ngbr14'           : [],
                     'ngbr14long'       : [],
                     'ngbr23'           : [],
                     'ngbr23long'       : [],
                     'scaling'          : None,
                     'residues'         : [],
                     'anames'           : [],
                     'sequence'         : [],
                     'solvcenter'       : [],
                     'solucenter'       : [],
                     'radii'            : None,
                     'exclusion'        : None,
                     'solvtype'         : None,
                     'excluded'         : [],
                     'topdir'           : None,
                     'coulomb'          : None,
                     '14scaling'        : None,
                    }

class Read_Topology(object):
    """
    Read Q topology file as an input, parse to topology class
    """
    def __init__(self, top, *args, **kwargs):    
        self.top = top
        
        data = Topology()
        self.data = data.data
        
    def JSON(self):
        constants = ['catypes','cbonds','cangles','ctorsions','cimpropers','ccharges']
        to_delete = []
        tmp = {}
        
        with open(self.top) as json_file:
            self.data = json.load(json_file)
            
        for key in self.data:
            if key in constants:
                tmp[key] = {}
                for i in self.data[key]:
                    tmp[key][int(i)] = self.data[key][i]
                    
                    if not key in to_delete:
                        to_delete.append(key)
        
        for key in to_delete:
            del self.data[key]

        self.data = {**tmp, **self.data}
        
        return(self.data)
        
    def Q(self):
        block = 0
        charges_tmp = []
        charges_type_tmp = {}
        header = {}
        cnt = 0
        switch = 0

        with open(self.top) as infile:
            for line in infile:
                    
                if 'Total no. of atoms' in line:
                    self.data['natoms_solute'] = line.split()[1]                    
                    block = 1
                    continue
                    
                if 'No. of integer atom codes' in line:
                    block = 2
                    continue
                    
                if 'No. of bonds' in line:
                    self.data['nbonds_solute'] = line.split()[1]
                    block = 3
                    continue

                if 'No. of bond codes' in line:
                    block = 4
                    continue

                if 'No. of angles' in line:
                    self.data['nangles_solute'] = line.split()[1]
                    block = 5
                    continue
                    
                if 'No. of angle codes' in line:
                    block = 6
                    continue

                if 'No. of torsions' in line:
                    self.data['ntorsions_solute'] = line.split()[1]                    
                    block = 7
                    continue

                if 'No. of torsion codes' in line:
                    block = 8
                    continue
                    
                if 'No. of impropers' in line:
                    self.data['nimpropers_solute'] = line.split()[1]                                        
                    block = 9
                    continue
                                    
                if 'No. of improper codes' in line:
                    block = 10
                    continue
                                                   
                if 'No. of atomic charges' in line:
                    block = 11
                    continue
                                                                   
                if 'No. of charge groups' in line:
                    block = 12
                    charge_groups = []
                    line = line.split()
                    solute_cgp = int(line[1])
                    total_cgp = int(line[0])
                    solvent_cgp = total_cgp - solute_cgp
                    self.data['charge_group_total'] = ['{}'.format(solute_cgp),'{}'.format(solvent_cgp)]
                    switch = 1
                    continue
                                                                                   
                if 'vdW combination rule' in line:
                    block = 13
                    continue
                                                                                   
                if 'Electrostatic 1-4 scaling factor' in line:
                    block = 14
                                                                                                   
                if 'Masses' in line:
                    block = 15
                    Masses = []
                    continue
                                                                                                                   
                if 'sqrt (Aii) normal' in line:
                    Aii_normal = []
                    block = 16
                    continue
                                                                                                                   
                if 'sqrt (Bii) normal' in line:
                    Bii_normal = []
                    block = 17
                    continue
                                                                                                                   
                if 'sqrt (Aii) polar' in line:
                    Aii_polar = []
                    block = 18
                    continue
                                                                                                                   
                if 'sqrt (Bii) polar' in line:
                    Bii_polar = []
                    block = 19
                    continue
                    
                if 'sqrt (Aii) 1-4' in line:
                    Aii_14 = []
                    block = 20
                    continue
                    
                if 'sqrt (Bii) 1-4' in line:
                    Bii_14 = []
                    block = 21
                    continue
                    
                if 'No. of type-2 vdW interactions' in line:
                    block = 22
                    continue
                    
                if 'No. of 1-4 neighbours' in line:
                    ngbr14 = []
                    block = 23
                    continue
                    
                if 'No. of long 1-4 nbrs' in line:
                    ngbr14long = []                                        
                    block = 24
                    continue

                if 'No. of exclusions' in line:
                    ngbr23 = []
                    block = 25
                    continue

                if 'No. of long exclusions' in line:
                    ngbr23long = []
                    block = 26
                    continue

                if 'No. of residues' in line:
                    residues = []
                    block = 27
                    continue

                if 'Sequence' in line:
                    sequence = []
                    block = 28
                    continue

                if 'No. of separate molecules' in line:
                    molecules = []
                    block = 29
                    continue

                if 'No. of atom types' in line:
                    anames = []
                    block = 30   
                    continue

                if 'No. of SYBYL atom types' in line:
                    block = 31
                    continue

                if 'solvent type (0=SPC,1=3-atom,2=general)' in line:
                    line = line.split()
                    self.data['solvtype']=(line[0])   
                    if line[0] != '0':
                         print('FATAL: solvent type other than SPC (0) not supported')
                         sys.exit()             
                    block = 32
                    continue

                if 'No. of excluded atoms' in line:
                    block = 33
                    continue
                    
                # Read stuff
                if block == 0:
                    line = line.split()
                    header[line[0]] = line[1:]
                    
                if block == 1:
                    line = line.split()
                    self.data['coords'].append(line)
                    
                if block == 2:
                    line = line.split()
                    for atype in line:
                        cnt += 1
                        self.data['atypes'].append([cnt,int(atype)])
                        
                if block == 3:
                    line = line.split()
                    self.data['bonds'].append(line)
                    
                if block == 4:
                    line = line.split()
                    self.data['cbonds'][int(line[0])] = [line[1],line[2]]
                    
                if block == 5:
                    line = line.split()
                    self.data['angles'].append(line)
                    
                if block == 6:
                    line = line.split()
                    self.data['cangles'][int(line[0])] = [line[1],line[2]]
                    
                if block == 7:
                    line = line.split()
                    self.data['torsions'].append(line)
                    
                if block == 8:
                    line = line.split()
                    self.data['ctorsions'][int(line[0])] = [line[1],line[2],line[3],line[4]]
                    
                if block == 9:
                    line = line.split()
                    self.data['impropers'].append(line)
                    
                if block == 10:
                    line = line.split()
                    self.data['cimpropers'][int(line[0])] = [line[1],line[2]]
                    
                if block == 11:
                    line = line.split()
                    for charge in line:
                        charges_tmp.append(charge)
                        
                if block == 12:
                    line = line.split()

                    # first get the info line
                    if switch == 1:
                        tmp = []
                        tmp.append(line)
                        switch = 2
                        n_atoms = int(line[0])
                        at_tmp = []
                        continue
                    
                    # even numbers are the atom numbers
                    if switch == 2:
                        if n_atoms > len(at_tmp):
                            for at in line:
                                at_tmp.append(at)

                            if n_atoms == len(at_tmp):
                                # put the atom list in the chatge_group thing
                                tmp.append(at_tmp)
                                charge_groups.append(tmp)
                                switch = 1
                           
                if block == 13:
                    continue
                           
                if block == 14:
                    line = line.split()
                    self.data['14scaling'] = line[0]
                    self.data['coulomb'] = line[1]
                    continue
                           
                if block == 15:
                    line = line.split()
                    for mass in line:
                        Masses.append(mass)
                                   
                if block == 16:
                    line = line.split()
                    for Aii in line:
                        Aii_normal.append(Aii)
                                   
                if block == 17:
                    line = line.split()
                    for Bii in line:
                        Bii_normal.append(Bii)
                                   
                if block == 18:
                    line = line.split()
                    for Aii in line:
                        Aii_polar.append(Aii)
                                   
                if block == 19:
                    line = line.split()
                    for Bii in line:
                        Bii_polar.append(Bii)
                                   
                if block == 20:
                    line = line.split()
                    for Aii in line:
                        Aii_14.append(Aii)
                                           
                if block == 21:
                    line = line.split()
                    for Bii in line:
                        Bii_14.append(Bii)
                        
                if block == 22:
                    continue
                    
                if block == 23:
                    ngbr14.append(line.strip())
                
                if block == 24:
                    ngbr14long.append(line.strip().split())
                    
                if block == 25:
                    ngbr23.append(line.strip())
                
                if block == 26:
                    ngbr23long.append(line.strip().split())
                        
                if block == 27:
                    line = line.split()
                    for residue in line:
                        residues.append(residue)  
                        
                if block == 28:
                    line = line.split()
                    for seq in line:
                        sequence.append(seq)
                        
                if block == 29:
                    line = line.split()
                    for molecule in line:
                        molecules.append(molecule)
                                                
                if block == 30:
                    line = line.split()
                    for aname in line:
                        anames.append(aname)
                                                        
                if block == 31:
                    # for now do nothing with SYBYL atom types
                    continue
                    
                if block == 32:
                    linesplit = line.split()
                    if 'Exclusion' in line:
                        if float(linesplit[0]) > 30.0:
                            print("Sphere sizes exceding 30A are currently not supported")
                            sys.exit()
                        self.data['exclusion'] = linesplit[0]
                        self.data['radii'] = linesplit[1]
                    
                    # capture ugly bug in Q where center = centre
                    if 'Solute centre' in line:
                        self.data['solucenter'] = [linesplit[0],linesplit[1],linesplit[2]]
                                                               
                    if 'Solute center' in line:
                        self.data['solucenter'] = [linesplit[0],linesplit[1],linesplit[2]]
                        
                    # capture ugly bug in Q
                    if 'Solvent center' in line:
                        self.data['solvcenter'] = [linesplit[0],linesplit[1],linesplit[2]]
                                                                
                    if 'Solvent centre' in line:
                        self.data['solvcenter'] = [linesplit[0],linesplit[1],linesplit[2]]
                        
                if block == 33:
                    line = line.strip()
                    for l in line:
                        l = l.strip()
                        if l == 'F':
                            l = '0'
                        else:
                            l = '1'

                        self.data['excluded'].append(l)
        
        # header construct
        self.data['header'] = header
        
        # split coordinates             
        tmp = list(itertools.chain.from_iterable(self.data['coords']))
        self.data['coords'] = IO.split_list(tmp,3)
        
        # split bonds
        tmp = list(itertools.chain.from_iterable(self.data['bonds']))
        self.data['bonds'] = IO.split_list(tmp,3)
        
        # split angles
        tmp = list(itertools.chain.from_iterable(self.data['angles']))
        self.data['angles'] = IO.split_list(tmp,4)        
                
        # split torsions
        tmp = list(itertools.chain.from_iterable(self.data['torsions']))
        self.data['torsions'] = IO.split_list(tmp,5)        
                        
        # split impropers
        tmp = list(itertools.chain.from_iterable(self.data['impropers']))
        self.data['impropers'] = IO.split_list(tmp,5)        
        
        # construct charges
        ctype = 0
        for i in range(0, len(charges_tmp)):
            charge = charges_tmp[i]
            if not charge in charges_type_tmp:
                ctype += 1
                charges_type_tmp[charge] = ctype
                        
            self.data['charges'].append([i+1, charges_type_tmp[charge]])
        
        for key in charges_type_tmp:
            self.data['ccharges'][charges_type_tmp[key]] = key
            
        # construct atom types
        for i in range(0, len(Masses)):
            self.data['catypes'][i+1] = [Masses[i],
                                          Aii_normal[i],
                                          Bii_normal[i],
                                          Aii_polar[i],
                                          Bii_polar[i],
                                          Aii_14[i],
                                          Bii_14[i]
                                         ]
            
        # construct 1-4
        ngbr14 = ''.join(ngbr14)
        self.data['ngbr14'] = [ngbr14[i:i+25] for i in range(0, len(ngbr14), 25)]
        
        #construct 1-4 long
        tmp = list(itertools.chain.from_iterable(ngbr14long))
        self.data['ngbr14long'] = IO.split_list(tmp,2)
                
        # construct 2-3
        ngbr23 = ''.join(ngbr23)
        self.data['ngbr23'] = [ngbr23[i:i+25] for i in range(0, len(ngbr23), 25)]
        
        # construct 2-3 long
        tmp = list(itertools.chain.from_iterable(ngbr23long))
        self.data['ngbr23long'] = IO.split_list(tmp,2)
        
        # Starting atom of every residue 
        self.data['residues'] = residues
        
        # Residue list
        self.data['sequence'] = sequence
                
        # Starting atom of molecules in topology
        self.data['molecules'] = molecules
                        
        # Atomtype names, matches with atype
        self.data['anames'] = anames
        
        # Join exclusions
        self.data['excluded'] = ''.join(self.data['excluded'])
        
        # Charge groups`
        cgp_cnt = 1
        for group in charge_groups:
            cgp_cnt += 1
            for ele in group[1]:
                cgp_cnt += 1
            
        self.data['charge_groups'] = [cgp_cnt, charge_groups]

        return(self.data)
        
class Write_Topology(object):        
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
        #Topology.coords = []
        with open(self.wd + 'coords.csv','w') as outfile:
            outfile.write('{}\n'.format(len(self.data['coords'])))
            outfile.write('{}\n'.format(self.data['natoms_solute']))            
            for line in self.data['coords']:
                outfile.write('{};{};{}\n'.format(line[0],line[1],line[2]))
        
        #Topology.atypes = []
        with open(self.wd + '/atypes.csv','w') as outfile:
            outfile.write('{}\n'.format(len(self.data['atypes'])))
            for line in self.data['atypes']:
                outfile.write('{};{}\n'.format(line[0],line[1]))
                
        #Topology.catypes = {}
        keys = sorted(self.data['catypes'].keys())
        with open(self.wd + '/catypes.csv','w') as outfile:
            outfile.write('{}\n'.format(len(self.data['catypes'])))
            for key in keys:
                outfile.write('{};{};{};{};{};{};{};{}\n'.format(key,
                                                                 self.data['catypes'][key][0],
                                                                 self.data['catypes'][key][1],
                                                                 self.data['catypes'][key][2],
                                                                 self.data['catypes'][key][3],
                                                                 self.data['catypes'][key][4],
                                                                 self.data['catypes'][key][5],
                                                                 self.data['catypes'][key][6],
                                                                ))
                
        #Topology.bonds = []
        with open(self.wd + '/bonds.csv','w') as outfile:
            outfile.write('{}\n'.format(len(self.data['bonds'])))
            outfile.write('{}\n'.format(self.data['nbonds_solute']))            
            for line in self.data['bonds']:
                outfile.write('{};{};{}\n'.format(line[0],line[1],line[2]))      
                
        #Topology.cbonds = {}
        keys = sorted(self.data['cbonds'].keys())
        with open(self.wd + '/cbonds.csv','w') as outfile:
            outfile.write('{}\n'.format(len(self.data['cbonds'])))
            for key in keys:
                outfile.write('{};{};{}\n'.format(key,
                                                  self.data['cbonds'][key][0],
                                                  self.data['cbonds'][key][1]
                                                 ))
        
        #Topology.angles = []
        with open(self.wd + '/angles.csv','w') as outfile:
            outfile.write('{}\n'.format(len(self.data['angles'])))
            outfile.write('{}\n'.format(self.data['nangles_solute']))                        
            for line in self.data['angles']:
                outfile.write('{};{};{};{}\n'.format(line[0],
                                                     line[1],
                                                     line[2],
                                                     line[3]))  
                
        #Topology.cangles = {}
        keys = sorted(self.data['cangles'].keys())        
        with open(self.wd + '/cangles.csv','w') as outfile:
            outfile.write('{}\n'.format(len(self.data['cangles'])))
            for key in keys:
                outfile.write('{};{};{}\n'.format(key,
                                                  self.data['cangles'][key][0],
                                                  self.data['cangles'][key][1]
                                                 ))
                
        #Topology.torsions = []
        with open(self.wd + '/torsions.csv','w') as outfile:
            outfile.write('{}\n'.format(len(self.data['torsions'])))
            outfile.write('{}\n'.format(self.data['ntorsions_solute']))                                    
            for line in self.data['torsions']:
                outfile.write('{};{};{};{};{}\n'.format(line[0],
                                                        line[1],
                                                        line[2],
                                                        line[3],
                                                        line[4],
                                                    ))
                
        #Topology.ctorsions = {}
        keys = sorted(self.data['ctorsions'].keys())                
        with open(self.wd + '/ctorsions.csv','w') as outfile:
            outfile.write('{}\n'.format(len(self.data['ctorsions'])))
            for key in keys:
                outfile.write('{};{};{};{}\n'.format(key,
                                                     self.data['ctorsions'][key][0],
                                                     self.data['ctorsions'][key][1],
                                                     self.data['ctorsions'][key][2],
                                                     self.data['ctorsions'][key][3],
                                                 ))
        #Topology.impropers = []
        with open(self.wd + '/impropers.csv','w') as outfile:
            outfile.write('{}\n'.format(len(self.data['impropers'])))
            outfile.write('{}\n'.format(self.data['nimpropers_solute']))                                                
            for line in self.data['impropers']:
                outfile.write('{};{};{};{};{}\n'.format(line[0],
                                                        line[1],
                                                        line[2],
                                                        line[3],
                                                        line[4],
                                                    ))
        #Topology.cimpropers = {}
        keys = sorted(self.data['cimpropers'].keys())                        
        with open(self.wd + '/cimpropers.csv','w') as outfile:
            outfile.write('{}\n'.format(len(self.data['cimpropers'])))
            for key in keys:
                outfile.write('{};{};{}\n'.format(key,
                                                  self.data['cimpropers'][key][0],
                                                  self.data['cimpropers'][key][1],
                                                 ))
        
        #Topology.charges = []
        with open(self.wd + '/charges.csv','w') as outfile:
            outfile.write('{}\n'.format(len(self.data['charges'])))
            for line in self.data['charges']:
                outfile.write('{};{}\n'.format(line[0],line[1]))
                
        #Topology.ccharges = {}
        keys = sorted(self.data['ccharges'].keys())                                
        with open(self.wd + '/ccharges.csv','w') as outfile:
            outfile.write('{}\n'.format(len(self.data['ccharges'])))
            for key in keys:
                outfile.write('{};{}\n'.format(key,self.data['ccharges'][key]))
                
        #Topology.ngbr14 = []
        with open(self.wd + '/ngbrs14.csv','w') as outfile:
            outfile.write('{}\n'.format(len(self.data['ngbr14'])))
            for line in self.data['ngbr14']:
                outfile.write('{}\n'.format(line))
                
        #Topology.ngbr14long = []
        with open(self.wd + '/ngbrs14long.csv','w') as outfile:
            outfile.write('{}\n'.format(len(self.data['ngbr14long'])))
            for line in self.data['ngbr14long']:
                outfile.write('{};{}\n'.format(line[0],line[1]))
                
        #Topology.ngbr23 = []
        with open(self.wd + '/ngbrs23.csv','w') as outfile:
            outfile.write('{}\n'.format(len(self.data['ngbr23'])))
            for line in self.data['ngbr23']:
                outfile.write('{}\n'.format(line))
                
        #Topology.ngbr23long = []
        with open(self.wd + '/ngbrs23long.csv','w') as outfile:
            outfile.write('{}\n'.format(len(self.data['ngbr23long'])))
            for line in self.data['ngbr23long']:
                outfile.write('{};{}\n'.format(line[0],line[1]))
                                
        #Topology.ngbr23long = []
        with open(self.wd + '/molecules.csv','w') as outfile:
            outfile.write('{}\n'.format(len(self.data['molecules'])))
            for line in self.data['molecules']:
                outfile.write('{}\n'.format(line[0]))
                
        #Topology.excluded = []
        with open(self.wd + '/excluded.csv','w') as outfile:
            outfile.write(self.data['excluded'] + '\n')
            
        #Topo.csv
        with open(self.wd + '/topo.csv','w') as outfile:
            outfile.write('7\n')
            outfile.write(self.data['solvtype'] + '\n')
            outfile.write(self.data['exclusion'] + '\n')
            outfile.write(self.data['radii'] + '\n')
            outfile.write('{};{};{}\n'.format(self.data['solucenter'][0],
                                              self.data['solucenter'][1],
                                              self.data['solucenter'][2],))
            
            outfile.write('{};{};{}\n'.format(self.data['solvcenter'][0],
                                              self.data['solvcenter'][1],
                                              self.data['solvcenter'][2],))
            
            outfile.write(self.data['14scaling'] + '\n')
            outfile.write(self.data['coulomb'] + '\n')

        # Charge groups
        with open(self.wd + '/charge_groups.csv','w') as outfile:
            # TO DO ADD LINE TOTALS
            outfile.write('{}\n'.format(self.data['charge_groups'][0]))
            outfile.write('{};{}\n'.format(self.data['charge_group_total'][0],
                                           self.data['charge_group_total'][1]))

            for charge_group in self.data['charge_groups'][1]:
                outfile.write('{};{}\n'.format(charge_group[0][0],charge_group[0][1]))
                for atom in charge_group[1]:
                    outfile.write('{}\n'.format(atom))       