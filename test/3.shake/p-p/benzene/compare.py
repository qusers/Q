import json

class bcolors:
    OKGREEN = '\033[92m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'    

refdata = {}
passed = True

def parse_energy(energyfile):
    """
        Parses Qdyn energy file 
    """
    total = []
    block = None
    
    energies = {
                'temperature'   : {
                    'Temp' : None
                },
        
                'bonded'        : {
                    'p' : [None,None,None,None],
                    'w' : [None,None,None,None],
                    'qp': [None,None,None,None], # should be q? 
                },
            
                'nonbonded'     : {
                    'pp' : [None,None],
                    'pw' : [None,None],
                    'ww' : [None,None],
                    'qx' : [None,None],
                },
        
                'restraint'     : {
                    'Uradx'  :   None,
                    'Upolx'  :   None,
                    'Ushell' :   None,
                    'Ufix'   :   None,
                    'Upres'  :   None,
                    'Total'  :   None,
                },
        
                'q-energies'    : {
                    'lambda'    :  None,
                    'SUM'       :  None,
                    'Ubond'	    :  None,
                    'Uangle'	:  None,
                    'Utor'      :  None,
                    'Uimp'      :  None,	
                    'Uvdw'      :  None,	
                    'Ucoul'     :  None,	
                    'Urestr'    :  None,
                },
        
                'total'         : {
                    'Ukin'    :  None,
                    'Upot'    :  None,
                    'Utot'    :  None,
                },   
               }
    
    with open(energyfile) as infile:
        i = -1
        for line in infile:
            if len(line) < 2:
                continue
                
            if 'lambdas' in line:
                block = -1
                continue
                
            # reset the block
            if 'interval' in line:
                i += 1
                block = 0
                total.append(energies)
            
            if 'type' in line:
                continue
                
            if block == -1:
                # construct lambdas
                states = int(line[0])
                for i in range(0,states):
                    for key in total[i]['q-energies']:
                        total[i]['q-energies'].append([])
            
            # Find header
            if '[temperature]' in line:
                block = 1
                continue
                
            if '[bonded]' in line:
                block = 2
                continue
                
            if '[nonbonded]' in line:
                block = 3
                continue
                
            if '[restraint]' in line:
                block = 4
                continue
                
            if '[q-energies]' in line:
                l = -1
                block = 5
                continue
                
            if '[total]' in line:
                block = 6
                continue
                
            if block == 1:
                line = line.split()
                total[i]['temperature'][line[0]] = '{:.2f}'.format(float(line[1]))
                                
            if block == 2:
                line = line.split()
                total[i]['bonded'][line[0]][0] = '{:.2f}'.format(float(line[1]))
                total[i]['bonded'][line[0]][1] = '{:.2f}'.format(float(line[2]))
                total[i]['bonded'][line[0]][2] = '{:.2f}'.format(float(line[3]))
                total[i]['bonded'][line[0]][3] = '{:.2f}'.format(float(line[4]))                
                
            if block == 3:
                line = line.split()
                total[i]['nonbonded'][line[0]][0] = '{:.2f}'.format(float(line[1]))
                total[i]['nonbonded'][line[0]][1] = '{:.2f}'.format(float(line[2]))
                                                
            if block == 4:
                line = line.split()
                total[i]['restraint'][line[0]] = '{:.2f}'.format(float(line[1]))
                                                
            if block == 5:
                l += 1
                line = line.split()
                if total[i]['q-energies'][line[0]] != None:
                    total[i]['q-energies'][line[0]][l] = line[1]
                                                
            if block == 6:
                line = line.split()
                total[i]['total'][line[0]] = '{:.2f}'.format(float(line[1]))
                
    return(total)

with open('Q5_data/Q_data.json') as infile:
    Q_data = json.load(infile)

# Get the data, we are only interested in zeroth frame    
QGPU_data = parse_energy('TEST/benzene-vacuum/output/energies.csv')[0]

# nonbonded interactions
if Q_data['solute'][0] != QGPU_data['nonbonded']['pp'][0]:
    print("Energy not matching for {} {} {}, Q5 {} Q7 {}".format('nonbonded',
                                                              'pp',
                                                              'el',
                                                              Q_data['solute'][0],
                                                              QGPU_data['nonbonded']['pp'][0],
                                                             ))
    passed = False
    
if Q_data['solute'][1] != QGPU_data['nonbonded']['pp'][1]:
    print("Energy not matching for {} {} {}, Q5 {} Q7 {}".format('nonbonded',
                                                              'pp',
                                                              'vdw',
                                                              Q_data['solute'][1],
                                                              QGPU_data['nonbonded']['pp'][1],
                                                             )) 
    passed = False
    
if Q_data['solute-solvent'][1] != QGPU_data['nonbonded']['pw'][1]:
    print("Energy not matching for {} {} {}, Q5 {} Q7 {}".format('nonbonded',
                                                              'pw',
                                                              'vdw',
                                                              Q_data['solvent'][1],
                                                              QGPU_data['nonbonded']['ww'][1],
                                                             ))            
    passed = False

# bonded interactions solute
if Q_data['solute'][2] != QGPU_data['bonded']['p'][0]:
    print("Energy not matching for {} {} {}, Q5 {} Q7 {}".format('bonded',
                                                              'p',
                                                              'bond',
                                                              Q_data['solute'][2],
                                                              QGPU_data['bonded']['p'][0],
                                                             ))        
    passed = False
    
if Q_data['solute'][3] != QGPU_data['bonded']['p'][1]:
    print("Energy not matching for {} {} {}, Q5 {} Q7 {}".format('bonded',
                                                              'p',
                                                              'angle',
                                                              Q_data['solute'][3],
                                                              QGPU_data['bonded']['p'][1],
                                                             ))        
    passed = False
    
if Q_data['solute'][4] != QGPU_data['bonded']['p'][2]:
    print("Energy not matching for {} {} {}, Q5 {} Q7 {}".format('bonded',
                                                              'p',
                                                              'torsion',
                                                              Q_data['solute'][4],
                                                              QGPU_data['bonded']['p'][2],
                                                             ))       
    passed = False
    
if Q_data['solute'][5] != QGPU_data['bonded']['p'][3]:
    print("Energy not matching for {} {} {}, Q5 {} Q7 {}".format('bonded',
                                                              'p',
                                                              'improper',
                                                              Q_data['solute'][5],
                                                              QGPU_data['bonded']['p'][3],
                                                             ))           
    passed = False
        
# restraint data
#  total       fix slvnt_rad slvnt_pol     shell    solute
if Q_data['restraints'][0] != QGPU_data['restraint']['Total']:
    print("Energy not matching for {} {}, Q5 {} Q7 {}".format('restraint',
                                                              'Total',
                                                              Q_data['restraints'][0],
                                                              QGPU_data['restraints']['Total'],
                                                             ))        
    passed = False
    
if Q_data['restraints'][1] != QGPU_data['restraint']['Ufix']:
    print("Energy not matching for {} {}, Q5 {} Q7 {}".format('restraint',
                                                              'Ufix',
                                                              Q_data['restraints'][1],
                                                              QGPU_data['restraints']['Ufix'],
                                                             ))      
    passed = False    
    
if Q_data['restraints'][2] != QGPU_data['restraint']['Uradx']:
    print("Energy not matching for {} {}, Q5 {} Q7 {}".format('restraint',
                                                              'Uradx',
                                                              Q_data['restraints'][2],
                                                              QGPU_data['restraints']['Uradx'],
                                                             ))        
    passed = False
    
if Q_data['restraints'][3] != QGPU_data['restraint']['Upolx']:
    print("Energy not matching for {} {}, Q5 {} Q7 {}".format('restraint',
                                                              'Upolx',
                                                              Q_data['restraints'][3],
                                                              QGPU_data['restraints']['Upolx'],
                                                             ))            
    
    passed = False
    
if Q_data['restraints'][4] != QGPU_data['restraint']['Ushell']:
    print("Energy not matching for {} {}, Q5 {} Q7 {}".format('restraint',
                                                              'Ushell',
                                                              Q_data['restraints'][4],
                                                              QGPU_data['restraints']['Ushell'],
                                                             ))       
    passed = False
        
if Q_data['restraints'][5] != QGPU_data['restraint']['Upres']:
    print("Energy not matching for {} {}, Q5 {} Q7 {}".format('restraint',
                                                              'Upres',
                                                              Q_data['restraints'][5],
                                                              QGPU_data['restraints']['Upres'],
                                                             ))         
    passed = False

# Total energies
if Q_data['SUM'][0] != QGPU_data['total']['Utot']:
    print("Energy not matching for {} {}, Q5 {} Q7 {}".format('SUM',
                                                              'total',
                                                              Q_data['SUM'][0],
                                                              QGPU_data['total']['Utot'],
                                                             ))
    passed = False

    
if Q_data['SUM'][1] != QGPU_data['total']['Upot']:
    print("Energy not matching for {} {}, Q5 {} Q7 {}".format('SUM',
                                                              'Upot',
                                                              Q_data['SUM'][1],
                                                              QGPU_data['total']['Upot'],
                                                             ))        
    passed = False
    
if Q_data['SUM'][2] != QGPU_data['total']['Ukin']:
    print("Energy not matching for {} {}, Q5 {} Q7 {}".format('SUM',
                                                              'Ukin',
                                                              Q_data['SUM'][2],
                                                              QGPU_data['total']['Ukin'],
                                                             ))        
    passed = False

        
    
## PRINT PASS ##    
if passed == False:
    print('Passed test? ' + f"{bcolors.FAIL} FALSE {bcolors.ENDC}")
    
if passed == True:
    print('Passed test? ' + f"{bcolors.OKGREEN} TRUE {bcolors.ENDC}")    
