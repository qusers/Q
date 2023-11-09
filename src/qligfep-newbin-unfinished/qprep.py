# NOT USED IN PRODUCTION CODE YET

import glob
import numpy as np
import argparse
import os
import sys
import networkx as nx

import IO
import functions as f
import potential_energy as Upot
import geometries

class Topology():
    def __init__(self):
        Topology.atoms = {}
        Topology.atom_types = {}
        Topology.bonds = {}
        Topology.angles = {}
        Topology.torsions = {}
        Topology.impropers = {}
        Topology.coordinates = []

class Prepare_Topology(object):
    """
    Read in coordinate files and associated forcefield files to construct a 
    topology file.
    """
    def __init__(self, lig, *args, **kwargs):    
        self.topology = Topology()
        self.lig = lig
        self.required = ['.lib','.prm','.pdb']
        self.connectivitiy = nx.Graph()
        
        # some useful mappings
        self.lib2prm = {}
        
    def check_files(self):
        OK = True
        for extension in self.required:
            file = self.lig + extension
            if os.path.isfile((file)) == False:
                print('Could not find {}'.format(file))
                
        if OK == False:
            print('missing required files, quiting now')
            sys.exit()            
        
    def read_prm(self):
        self.prm = IO.read_prm([self.lig + '.prm'])
        
    def read_lib(self):
        self.lib = IO.read_lib(self.lig + '.lib')
        
        # here, generate the lib/prm atom mapping, add the charges to the atom properties in prm
        for resname in self.lib:
            for line in self.lib[resname]['[atoms]']:
                self.lib2prm[line[1]] = line[2]
                
                if line[2] in self.prm['[atom_types]']:
                    # Add the partial charge to the system
                    self.prm['[atom_types]'][line[2]].append(line)
                    
                else:
                    # If the atom is not there, something is wrong
                    print("FATAL, lib atom not found in .prm file, check .prm file")
                    print("FATAL: overlapping atom names: {}".format(line[2]))
                    sys.exit()
                    
    def generate_connectivity(self):
        # Needs a function for CONNECTs
        for resname in self.lib:
            self.connectivitiy.add_edges_from(self.lib[resname]['[bonds]'])
        # construct atom mapping
        for ai in list(self.connectivitiy.nodes):
            for aj in list(self.connectivitiy.nodes):
                if ai == aj:
                    continue
                        
                # construct all possible bonds
                for path in nx.all_simple_paths(self.connectivitiy, source=ai, target=aj, cutoff = 4):
                    # Bonds
                    if len(path) == 2:
                        bond_ref = tuple([self.lib2prm[path[0]], self.lib2prm[path[1]]])
                        bond = tuple(bond_ref)
                        
                        if bond_ref[::-1] in self.topology.bonds:
                            continue
                        
                        # Check if the lib matches the prm
                        if bond in self.prm['[bonds]']:
                            self.topology.bonds[bond_ref] = self.prm['[bonds]'][bond]
                            
                        else:
                            print("FATAL: bond in library not in .prm")
                            print("FATAL: bond {} {}".format(bond[0],bond[1]))
                        
                    # Angles
                    if len(path) == 3:
                        angle_ref = tuple(([self.lib2prm[path[0]],
                                      self.lib2prm[path[1]],
                                      self.lib2prm[path[2]]
                                      ]))                        
                        
                        angle = tuple(sorted(angle_ref))
                        
                        if angle_ref[::-1] in self.topology.angles:
                            continue
                            
                        if angle in self.prm['[angles]']:
                            self.topology.angles[angle_ref] = self.prm['[angles]'][angle]
                            
                        else:
                            print("FATAL: angle in library not in .prm")
                            print("FATAL: angle {} {} {}".format(bond[0],bond[1],bond[2]))
                            
                    # Torsions
                    if len(path) == 4:
                        torsion_ref = tuple([self.lib2prm[path[0]],
                                             self.lib2prm[path[1]],
                                             self.lib2prm[path[2]],
                                             self.lib2prm[path[3]]
                                             ])
                        torsion = tuple(sorted(torsion_ref))
                        
                        if torsion_ref[::-1] in self.topology.torsions:
                            continue
                        
                        if torsion in self.prm['[torsions]']:
                            self.topology.torsions[torsion_ref] = [self.prm['[torsions]'][torsion]]
                        
                        else:
                            print("FATAL: torsion in library not in .prm")
                            print("FATAL: torsion {} {} {} {}".format(bond[0],bond[1],bond[2],bond[3]))
                            
        for improper in self.lib['LIG']['[impropers]']:
            print("TO DO")
            print(improper)
            
    def read_pdb(self):
        # make a function?
        self.pdb = []
        print(self.topology)
        with open(self.lig + '.pdb') as infile:
            # Note, the atom index will be based on the pdb line index
            # Might want to make some mapping if we are to write back to
            # PDB or SDF file
            i = -1
            for line in infile:
                line=IO.pdb_parse_in(line)
                if not line[2] in self.lib2prm[line[2]]:
                    print("FATAL: atom type not found")
                    sys.exit()
                
                at_type = self.lib2prm[line[2]]    
                i += 1
                # Let's see if we want to store more information later
                self.pdb.append(line)
                
                # Construct the atom mapping (1D array should do it)
                self.topology.atoms[i] = at_type
                self.topology.atom_types[at_type] = self.prm['[atom_types]'][at_type]
                        
        # construct a np array of the coordinates in the system
        self.np_pdb = np.array(self.pdb)
        x = self.np_pdb[:,8].astype(float)
        y = self.np_pdb[:,9].astype(float)
        z = self.np_pdb[:,10].astype(float)
        
        # This is the initial coordinate matrix
        self.topology.coordinates.append(np.transpose(np.array([x,y,z])))

    def checktop(self):
        print("Found {} atom types".format(len(self.topology.atoms)))
        print("Found {} bond types".format(len(self.topology.bonds)))
        print("Found {} angle types".format(len(self.topology.angles)))
        print("Found {} torsion types".format(len(self.topology.torsions)))
        print("Found {} improper types".format(len(self.topology.impropers)))
        
    def writetop(self):
        self.topfile = self.lig + '.top'
        with open (self.topfile, 'w') as outfile:
            #print(self.topology.coordinates)
            #print(self.topology.atoms)
            #print(self.topology.atom_types)
            #print(self.topology.bonds)
            #print(self.topology.angles)
            print(self.topology.torsions)
                        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='Py-MD',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '       == Python based MD engine == ')

    
    parser.add_argument('-l', '--lig',
                        dest = "lig",
                        required = True,
                        help = "test")
    
    args = parser.parse_args()
    prepare_top = Prepare_Topology(lig = args.lig)
    
    # Run preperation
    prepare_top.check_files()               # 00
    prepare_top.read_prm()                  # 01
    prepare_top.read_lib()                  # 02
    prepare_top.generate_connectivity()     # 03
    prepare_top.read_pdb()                  # 04
    prepare_top.checktop()                  # 05
    prepare_top.writetop()                  # 06