"""Module containing functions for parsing pdb files."""

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