import re
import shlex
from subprocess import check_output
import os
import stat
import numpy as np
import sys
import shutil
import time
import json

## Global objects
# Translational dictionary for charged residues
charged_res = {'HIS': {'HD1' : 'HID',
                       'HE2' : 'HIE'},
               'GLU': {'HE2' : 'GLH'},
               'ASP': {'HD2' : 'ASH'}
              }    

def replace(string, replacements):
    pattern = re.compile(r'\b(' + '|'.join(replacements.keys()) + r')\b')
    replaced_string = pattern.sub(lambda x: replacements[x.group()], string)
    return replaced_string

def JSONdump(data, outfile):
        """
        Dump data to a json file
        """
        print(outfile)
        with open(outfile, 'w') as outfile:
            json.dump(data,outfile,indent=2) 

def json2Q(json_in, q_out):
    """
        Translate Q7 json file to Q6 input file
    """  
    with open(json_in) as infile:
        data = json.load(infile)

        print(data)

def pdb_parse_in(line, include=('ATOM','HETATM')):
    """
    Takes a pdb file line and parses it into a list, according to Atomic Coordinate Entry Format 
    v3.3
    """
    at_entry = []
    line = line.strip('\n')
    if line.startswith(include):
        at_entry.append(line[0:6])              #  0 ATOM/HETATM
        at_entry.append(int(line[6:11]))        #  1 ATOM serial number
        at_entry.append(line[12:16].strip())    #  2 ATOM name
        at_entry.append(line[16:17])            #  3 Alternate location indicator
        at_entry.append(line[17:21].strip())    #  4 Residue name
        at_entry.append(line[21:22])            #  5 Chain identifier
        at_entry.append(int(line[22:26]))       #  6 Residue sequence number
        at_entry.append(line[26:27])            #  7 Code for insertion of residue
        at_entry.append(float(line[30:38]))     #  8 Orthogonal coordinates for X
        at_entry.append(float(line[38:46]))     #  9 Orthogonal coordinates for Y
        at_entry.append(float(line[46:54]))     # 10 Orthogonal coordinates for Z
        # These entries can be empty
        try:
            at_entry.append(float(line[54:60])) # 11 Occupancy
            
        except:
            at_entry.append(0.0)                # 11 Empty Occupancy
            
        try:
            at_entry.append(float(line[60:66])) # 12 Temperature factor
            
        except:
            at_entry.append(0.0)                # 12 Empty Temperature factor
            
        try:
            at_entry.append(line[76:78])        # 13 Element symbol
            
        except:
            at_entry.append('  ')               # 13 Empty Element symbol
            
        try:
            at_entry.append(line[78:80])        # 14 Charge on atom
            
        except:
            at_entry.append('  ')               # 14 Empty charge
        
    else:
        at_entry = line
    
    return at_entry
    
def pdb_parse_out(line):
    """
    Takes a list and parses it into a pdb writeable line
    """
    #if len(line[2]) <= 3: 
    #    line = '{:6s}{:5d} {:4s}{:1s}{:4s}{:1s}{:4d}{:1s}   '\
    #           '{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}'.format(*line)
    line = "{:6s}{:5d} {:<4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(*line)        
    #elif len(line[2]) == 4: 
    #    line = '{:6s}{:5d} {:4s}{:1s}{:4s}{:1s}{:4d}{:1s}   '\
    #           '{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}'.format(*line)
    return line

def read_prm(prmfiles):
    """
    Takes a list of Q .prm files and merges them, first file is the referene .prm file
    Returns a dicitonary of the merged .prm files
    """    
    block = 0
    prm = {'[options]':[],
            '[atom_types]':{},
            '[bonds]':{},
            '[angles]':{},
            '[torsions]':{},        # NOTE NEEDS FIX
            '[impropers]':{}}
    
    for filename in prmfiles:
        with open(filename) as infile:
            for line in infile:
                if len(line.strip()) == 0:
                    continue
                if line[0][0] == '!':
                    continue
                    
                line = line.strip().split()
    
                if line[0] == '[options]':
                    block = 1
                    continue                
                elif line[0] == '[atom_types]':
                    block = 2
                    continue
                elif line[0] == '[bonds]':
                    block = 3
                    continue
                elif line[0] == '[angles]':
                    block = 4
                    continue
                elif line[0] == '[torsions]':
                    block = 5
                    continue
                if line[0] == '[impropers]':
                    block = 6
                    continue

                if block == 1:
                    prm['[options]'].append(line)
                
                ## Probably want to fix this so we make a function or smth.
                if block == 2:
                    ats = line[0]
                    if ats not in prm['[atom_types]']:
                        prm['[atom_types]'][ats] = line[1:]
                    
                    else:
                        print("FATAL: overlapping atom names, check your .prm file")
                        print("FATAL: overlapping atom names: {}".format(ats))
                        sys.exit()

                elif block == 3:
                    ats = tuple(sorted((line[0],line[1])))
                    if ats not in prm['[bonds]']:
                        prm['[bonds]'][ats] = line[2:]
                        
                    else:
                        print("FATAL: overlapping atom names, check your .prm file")
                        print("FATAL: overlapping atom names: {}".format(ats))
                        sys.exit()

                elif block == 4:
                    ats = tuple(sorted([line[0],line[1],line[2]]))
                    if ats not in prm['[angles]']:
                        prm['[angles]'][ats] = line[3:] 
                        
                    else:
                        print("FATAL: overlapping atom names, check your .prm file")
                        print("FATAL: overlapping atom names: {}".format(ats))
                        sys.exit()
                        
                elif block == 5:
                    ats = tuple(sorted([line[0],line[1],line[2],line[3]]))
                    
                    # Torsions are Fourier with maximum four terms
                    if ats not in prm['[torsions]']:
                        prm['[torsions]'][ats] = [line[4:]]
                    
                    elif len(prm['[torsions]'][ats]) < 5:
                        prm['[torsions]'][ats].append([line[4:]])
                        
                    else:
                        print("FATAL: overlapping atom names, check your .prm file")
                        print("FATAL: overlapping atom names: {}".format(ats))
                        sys.exit()

                elif block == 6:
                    ats = tuple(sorted([line[0],line[1],line[2],line[3]]))
                    if ats not in prm['[impropers]']:
                        prm['[impropers]'][ats] = line[4:]
                    
                    else:
                        print("FATAL: overlapping atom names, check your .prm file")
                        print("FATAL: overlapping atom names: {}".format(ats))
                        sys.exit()                
    return prm

def read_lib(libfiles):
    """
    Takes a list of Q .libraries gives back the dictionary of residues in that library
    """
    libdic = {}
    block = 0    

    for libfile in libfiles:
        with open(libfile) as infile:
            for line in infile:
                line = line.strip()
                if len(line) < 1:
                    continue

                if line[0] == '{':
                    resname = line.split()
                    resname = resname[0]
                    resname = resname.strip('{')
                    resname = resname.strip('}')

                    if resname not in libdic:
                        libdic[resname] = {'[info]':[],
                                           '[atoms]':[],
                                           '[bonds]':[],
                                           '[impropers]':[],
                                           '[charge_groups]':[],
                                           'lib2prm':{},
                                            }
                        continue
                    else:
                        print('Warning, duplicate residue name, did not add parameters')

                if line[0] == "*":
                    continue

                if line[0] == '!':
                    continue
                if line == '[info]':
                    block = 1
                    continue                
                elif line == '[atoms]':
                    block = 2
                    continue
                elif line == '[bonds]':
                    block = 3
                    continue
                elif line == '[impropers]':
                    block = 4
                    continue
                elif line == '[charge_groups]':
                    block = 5
                    continue

                elif line == '[connections]':
                    block = 6
                    continue

                if block == 1:
                    libdic[resname]['[info]'].append(line.split())

                if block == 2:
                    libdic[resname]['[atoms]'].append(line.split())
                    libdic[resname]['lib2prm'][line.split()[1]] = line.split()[2]

                elif block == 3:
                    # the bonded part is a connectiviy graph
                    line = line.split()

                    libdic[resname]['[bonds]'].append((line[0],line[1]))

                elif block == 4:
                    libdic[resname]['[impropers]'].append(line.split())                

                elif block == 5:
                    libdic[resname]['[charge_groups]'].append(line.split())
                
    return(libdic)

def split_list(lst,n):
    return [lst[i:i + n] for i in range(0, len(lst), n)]

def run_command(executable, options, runtime=False, string=False):
    """
    Takes two variables, the executable location and its options (as string), and 
    runs the program.
    
    Returns the output of that program as an unformatted string.
    """
    out = None
    time1 = time.time()
    if string == False:
        args = shlex.split(executable + options)
        out = check_output(args)
        print(' '.join(args))

    else:
        os.system(executable + options)
        out = None
    
    time2 = time.time()    
    if runtime == True:
        print('runtime {:.3f} s'.format((time2-time1)))

    return out

def create_dir(d, overwrite):
    if not os.path.exists(d):
        os.mkdir(d)
        
    elif overwrite == False:
        return

    else:
        if overwrite == True:
            print("Overwriting folder {}".format(d))
        
        shutil.rmtree(d)
        os.mkdir(d)