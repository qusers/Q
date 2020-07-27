# Standard Python libraries
import os
import itertools
from os import path

# Q-GPU libraries
import IO

class Topology():
    def __init__(self):
        Topology.header = {}
        Topology.coords = []            # block 01  ## TODO for reference
        Topology.atypes = []        
        Topology.catypes = {}        
        Topology.nbonds_solute = None
        Topology.bonds = []
        Topology.cbonds = {}
        Topology.nangles_solute = None        
        Topology.angles = []
        Topology.cangles = {}
        Topology.ntorsions_solute = None                
        Topology.torsions = []
        Topology.ctorsions = {}
        Topology.nimpropers_solute = None                
        Topology.impropers = []
        Topology.cimpropers = {}
        Topology.charges = []
        Topology.ccharges = {}
        Topology.ngbr14 = []
        Topology.ngbr14long = []
        Topology.ngbr23 = []
        Topology.ngbr23long = []
        Topology.scaling = None
        Topology.residues = []
        Topology.anames = []
        Topology.sequence = []
        Topology.solvcenter = []
        Topology.solucenter = []
        Topology.radii = None
        Topology.exclusion = None
        Topology.excluded = []

class Read_Topology(object):
    """
    Read Q topology file as an input, parse to topology class
    """
    def __init__(self, top, *args, **kwargs):    
        self.topology = Topology()
        self.topology_file = top
        
    def parse_topology(self):
        block = 0
        charges_tmp = []
        charges_type_tmp = {}
        header = {}
        cnt = 0

        with open(self.topology_file) as infile:
            for line in infile:
                    
                if 'Total no. of atoms' in line:
                    Topology.natoms_solute = line.split()[1]                    
                    block = 1
                    continue
                    
                if 'No. of integer atom codes' in line:
                    block = 2
                    continue
                    
                if 'No. of bonds' in line:
                    Topology.nbonds_solute = line.split()[1]
                    block = 3
                    continue

                if 'No. of bond codes' in line:
                    block = 4
                    continue

                if 'No. of angles' in line:
                    Topology.nangles_solute = line.split()[1]
                    block = 5
                    continue
                    
                if 'No. of angle codes' in line:
                    block = 6
                    continue

                if 'No. of torsions' in line:
                    Topology.ntorsions_solute = line.split()[1]                    
                    block = 7
                    continue

                if 'No. of torsion codes' in line:
                    block = 8
                    continue
                    
                if 'No. of impropers' in line:
                    Topology.nimpropers_solute = line.split()[1]                                        
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
                    continue
                                                                                   
                if 'vdW combination rule' in line:
                    block = 13
                    continue
                                                                                   
                if 'Electrostatic 1-4 scaling factor' in line:
                    block = 14
                    continue
                                                                                                   
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
                    self.topology.coords.append(line)
                    
                if block == 2:
                    line = line.split()
                    for atype in line:
                        cnt += 1
                        self.topology.atypes.append([cnt,int(atype)])
                        
                if block == 3:
                    line = line.split()
                    self.topology.bonds.append(line)
                    
                if block == 4:
                    line = line.split()
                    self.topology.cbonds[line[0]] = [line[1],line[2]]
                    
                if block == 5:
                    line = line.split()
                    self.topology.angles.append(line)
                    
                if block == 6:
                    line = line.split()
                    self.topology.cangles[line[0]] = [line[1],line[2]]
                    
                if block == 7:
                    line = line.split()
                    self.topology.torsions.append(line)
                    
                if block == 8:
                    line = line.split()
                    self.topology.ctorsions[line[0]] = [line[1],line[2],line[3],line[4]]
                    
                if block == 9:
                    line = line.split()
                    self.topology.impropers.append(line)
                    
                if block == 10:
                    line = line.split()
                    self.topology.cimpropers[line[0]] = [line[1],'-2.000',line[2],'1']
                    
                if block == 11:
                    line = line.split()
                    for charge in line:
                        charges_tmp.append(charge)
                        
                if block == 12:
                    continue
                           
                if block == 13:
                    continue
                           
                if block == 14:
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
                    ngbr14long.append(line.strip())
                    
                if block == 25:
                    ngbr23.append(line.strip())
                
                if block == 26:
                    ngbr23long.append(line.strip())
                        
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
                        Topology.exclusion = linesplit[0]
                        Topology.radii = linesplit[1]
                        
                    if 'Solute center' in line:
                        Topology.solucenter = [linesplit[0],linesplit[1],linesplit]
                                        
                    if 'Solvent center' in line:
                        Topology.solvcenter = [linesplit[0],linesplit[1],linesplit]
                        
                if block == 33:
                    line = line.strip()
                    for l in line:
                        l = l.strip()
                        if l == 'F':
                            l = '0'
                        else:
                            l = '1'

                        Topology.excluded.append(l)
        
        # header construct
        Topology.header = header
        
        # split coordinates             
        self.topology.coords = list(itertools.chain.from_iterable(self.topology.coords))
        Topology.coords = IO.split_list(self.topology.coords,3)
        
        # split bonds
        self.topology.bonds = list(itertools.chain.from_iterable(self.topology.bonds))
        Topology.bonds = IO.split_list(self.topology.bonds,3)
        
        # split angles
        self.topology.angles = list(itertools.chain.from_iterable(self.topology.angles))
        Topology.angles = IO.split_list(self.topology.angles,4)        
                
        # split torsions
        self.topology.torsions = list(itertools.chain.from_iterable(self.topology.torsions))
        Topology.torsions = IO.split_list(self.topology.torsions,5)        
                        
        # split impropers
        self.topology.impropers = list(itertools.chain.from_iterable(self.topology.impropers))
        Topology.impropers = IO.split_list(self.topology.impropers,5)        
        
        # construct charges
        ctype = 0
        for i in range(0, len(charges_tmp)):
            charge = charges_tmp[i]
            if not charge in charges_type_tmp:
                ctype += 1
                charges_type_tmp[charge] = ctype
                        
            Topology.charges.append([i+1, charges_type_tmp[charge]])
        
        for key in charges_type_tmp:
            Topology.ccharges[charges_type_tmp[key]] = key
            
        # construct atom types
        for i in range(0, len(Masses)):
            Topology.catypes[i+1] = [Masses[i],
                                          Aii_normal[i],
                                          Bii_normal[i],
                                          Aii_polar[i],
                                          Bii_polar[i],
                                          Aii_14[i],
                                          Bii_14[i]
                                         ]
            
        # construct 1-4
        ngbr14 = ''.join(ngbr14)
        Topology.ngbr14 = [ngbr14[i:i+25] for i in range(0, len(ngbr14), 25)]
        
        #construct 1-4 long
        self.topology.ngbr14long = list(itertools.chain.from_iterable(self.topology.ngbr14long))
        Topology.ngbr14long = IO.split_list(self.topology.ngbr14long,2)
        
        # construct 2-3
        ngbr23 = ''.join(ngbr23)
        Topology.ngbr23 = [ngbr23[i:i+25] for i in range(0, len(ngbr23), 25)]
        
        # construct 2-3 long
        self.topology.ngbr23long = list(itertools.chain.from_iterable(self.topology.ngbr23long))
        Topology.ngbr23long = IO.split_list(self.topology.ngbr23long,2)
        
        # Starting atom of every residue 
        Topology.residues = residues
        
        # Residue list
        Topology.sequence = sequence
                
        # Starting atom of molecules in topology
        Topology.molecules = molecules
                        
        # Atomtype names, matches with atype
        Topology.anames = anames
        
        # Join exclusions
        Topology.excluded = ''.join(Topology.excluded)
        
        return(Topology)
        
class Write_Topology(object):        
    """
    Write Python topology object to file
    """
    def __init__(self, top, *args, **kwargs):    
        self.topology_file = top
        self.header = []
        self.workdir = os.getcwd()
        
    def write_csv(self):
        csvdir = self.workdir + '/' + self.topology_file[:-4]
        if path.exists(csvdir) == True:
            os.system('rm -r ' + csvdir)
        
        os.mkdir(csvdir)
            
        #Topology.coords = []
        with open(csvdir+'/coords.csv','w') as outfile:
            outfile.write('{}\n'.format(len(Topology.coords)))
            outfile.write('{}\n'.format(Topology.natoms_solute))            
            for line in Topology.coords:
                outfile.write('{};{};{}\n'.format(line[0],line[1],line[2]))
        
        #Topology.atypes = []
        with open(csvdir+'/atypes.csv','w') as outfile:
            outfile.write('{}\n'.format(len(Topology.atypes)))
            for line in Topology.atypes:
                outfile.write('{};{}\n'.format(line[0],line[1]))
                
        #Topology.catypes = {}
        keys = sorted(Topology.catypes.keys())
        with open(csvdir+'/catypes.csv','w') as outfile:
            outfile.write('{}\n'.format(len(Topology.catypes)))
            for key in keys:
                outfile.write('{};{};{};{};{};{};{};{}\n'.format(key,
                                                                 Topology.catypes[key][0],
                                                                 Topology.catypes[key][1],
                                                                 Topology.catypes[key][2],
                                                                 Topology.catypes[key][3],
                                                                 Topology.catypes[key][4],
                                                                 Topology.catypes[key][5],
                                                                 Topology.catypes[key][6],
                                                                ))
                
        #Topology.bonds = []
        with open(csvdir+'/bonds.csv','w') as outfile:
            outfile.write('{}\n'.format(len(Topology.bonds)))
            outfile.write('{}\n'.format(Topology.nbonds_solute))            
            for line in Topology.bonds:
                outfile.write('{};{};{}\n'.format(line[0],line[1],line[2]))      
                
        #Topology.cbonds = {}
        keys = sorted(Topology.cbonds.keys())
        with open(csvdir+'/cbonds.csv','w') as outfile:
            outfile.write('{}\n'.format(len(Topology.cbonds)))
            for key in keys:
                outfile.write('{};{};{}\n'.format(key,
                                                  Topology.cbonds[key][0],
                                                  Topology.cbonds[key][1]
                                                 ))
        
        #Topology.angles = []
        with open(csvdir+'/angles.csv','w') as outfile:
            outfile.write('{}\n'.format(len(Topology.angles)))
            outfile.write('{}\n'.format(Topology.nangles_solute))                        
            for line in Topology.angles:
                outfile.write('{};{};{};{}\n'.format(line[0],
                                                     line[1],
                                                     line[2],
                                                     line[3]))  
                
        #Topology.cangles = {}
        keys = sorted(Topology.cangles.keys())        
        with open(csvdir+'/cangles.csv','w') as outfile:
            outfile.write('{}\n'.format(len(Topology.cangles)))
            for key in keys:
                outfile.write('{};{};{}\n'.format(key,
                                                  Topology.cangles[key][0],
                                                  Topology.cangles[key][1]
                                                 ))
                
        #Topology.torsions = []
        with open(csvdir+'/torsions.csv','w') as outfile:
            outfile.write('{}\n'.format(len(Topology.torsions)))
            outfile.write('{}\n'.format(Topology.ntorsions_solute))                                    
            for line in Topology.torsions:
                outfile.write('{};{};{};{};{}\n'.format(line[0],
                                                        line[1],
                                                        line[2],
                                                        line[3],
                                                        line[4],
                                                    ))
                
        #Topology.ctorsions = {}
        keys = sorted(Topology.ctorsions.keys())                
        with open(csvdir+'/ctorsions.csv','w') as outfile:
            outfile.write('{}\n'.format(len(Topology.ctorsions)))
            for key in keys:
                outfile.write('{};{};{};{}\n'.format(key,
                                                     Topology.ctorsions[key][0],
                                                     Topology.ctorsions[key][1],
                                                     Topology.ctorsions[key][2],
                                                     Topology.ctorsions[key][3],
                                                 ))
        #Topology.impropers = []
        with open(csvdir+'/impropers.csv','w') as outfile:
            outfile.write('{}\n'.format(len(Topology.impropers)))
            outfile.write('{}\n'.format(Topology.nimpropers_solute))                                                
            for line in Topology.impropers:
                outfile.write('{};{};{};{};{}\n'.format(line[0],
                                                        line[1],
                                                        line[2],
                                                        line[3],
                                                        line[4],
                                                    ))
        #Topology.cimpropers = {}
        keys = sorted(Topology.cimpropers.keys())                        
        with open(csvdir+'/cimpropers.csv','w') as outfile:
            outfile.write('{}\n'.format(len(Topology.cimpropers)))
            for key in keys:
                outfile.write('{};{};{};{}\n'.format(key,
                                                     Topology.cimpropers[key][0],
                                                     Topology.cimpropers[key][1],
                                                     Topology.cimpropers[key][2],
                                                     Topology.cimpropers[key][3],
                                                 ))
        
        #Topology.charges = []
        with open(csvdir+'/charges.csv','w') as outfile:
            outfile.write('{}\n'.format(len(Topology.charges)))
            for line in Topology.charges:
                outfile.write('{};{}\n'.format(line[0],line[1]))
                
        #Topology.ccharges = {}
        keys = sorted(Topology.ccharges.keys())                                
        with open(csvdir+'/ccharges.csv','w') as outfile:
            outfile.write('{}\n'.format(len(Topology.ccharges)))
            for key in keys:
                outfile.write('{};{}\n'.format(key,Topology.ccharges[key]))
                
        #Topology.ngbr14 = []
        with open(csvdir+'/ngbrs14.csv','w') as outfile:
            outfile.write('{}\n'.format(len(Topology.ngbr14)))
            for line in Topology.ngbr14:
                outfile.write('{}\n'.format(line))
                
        #Topology.ngbr14long = []
        with open(csvdir+'/ngbrs14long.csv','w') as outfile:
            outfile.write('{}\n'.format(len(Topology.ngbr14long)))
            for line in Topology.ngbr14long:
                outfile.write('{};{}\n'.format(line[0],line[1]))
                
        #Topology.ngbr23 = []
        with open(csvdir+'/ngbrs23.csv','w') as outfile:
            outfile.write('{}\n'.format(len(Topology.ngbr23)))
            for line in Topology.ngbr23:
                outfile.write('{}\n'.format(line))
                
        #Topology.ngbr23long = []
        with open(csvdir+'/ngbrs23long.csv','w') as outfile:
            outfile.write('{}\n'.format(len(Topology.ngbr23long)))
            for line in Topology.ngbr23long:
                outfile.write('{};{}\n'.format(line[0],line[1]))