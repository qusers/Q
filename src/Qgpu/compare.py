import json
import numpy as np

# Q6 QGPU variable matching
header = [
    'nonbonded-pp-el',      # 00 QGPU_data['nonbonded']['pp'][0] = Q_data['solute'][0]
    'nonbonded-pp-vdw',     # 01 QGPU_data['nonbonded']['pp'][1] = Q_data['solute'][1]
    'nonbonded-qx-el',      # 02 QGPU_data['nonbonded']['qx'][0] = Q_data['Q-atom'][0]
    'nonbonded-qx-vdw',     # 03 QGPU_data['nonbonded']['qx'][1] = Q_data['Q-atom'][1]
    'bonded-qp-bond',       # 04 QGPU_data['bonded']['qp'][0]    = Q_data['Q-atom'][2]
    'bonded-qp-angle',      # 05 QGPU_data['bonded']['qp'][1]    = Q_data['Q-atom'][3] 
    'bonded-qp-torsion',    # 06 QGPU_data['bonded']['qp'][2]    = Q_data['Q-atom'][4]
    'bonded-qp-improper',   # 07 QGPU_data['bonded']['qp'][3]    = Q_data['Q-atom'][5]
    'nonbonded-pw-el',      # 08 QGPU_data['nonbonded']['pw'][0] = Q_data['solute-solvent'][0]
    'nonbonded-pw-vdw',     # 09 QGPU_data['nonbonded']['pw'][1] = Q_data['solute-solvent'][1]
    'nonbonded-ww-el',      # 10 QGPU_data['nonbonded']['ww'][0] = Q_data['solvent'][0]
    'nonbonded-ww-vdw',     # 11 QGPU_data['nonbonded']['ww'][1] = Q_data['solvent'][1]
    'bonded-w-bond',        # 12 QGPU_data['bonded']['w'][0]     = Q_data['solvent'][2]
    'bonded-w-angle',       # 13 QGPU_data['bonded']['w'][1]     = Q_data['solvent'][3]
    'bonded-w-torsion',     # 14 QGPU_data['bonded']['w'][2]     = Q_data['solvent'][4]
    'bonded-w-improper',    # 15 QGPU_data['bonded']['w'][3]     = Q_data['solvent'][5]
    'bonded-p-bond',        # 16 QGPU_data['bonded']['p'][0]     = Q_data['solute'][2]
    'bonded-p-angle',       # 17 QGPU_data['bonded']['p'][1]     = Q_data['solute'][3]
    'bonded-p-torsion',     # 18 QGPU_data['bonded']['p'][2]     = Q_data['solute'][4]
    'bonded-p-improper',    # 19 QGPU_data['bonded']['p'][3]     = Q_data['solute'][5]
    'restraint-Utot',       # 20 QGPU_data['restraint']['Total'] = Q_data['restraints'][0]
    'restraint-Ufix',       # 21 QGPU_data['restraint']['Ufix']  = Q_data['restraints'][1]
    'restraint-Urad',       # 22 QGPU_data['restraint']['Uradx'] = Q_data['restraints'][2]
    'restraint-Upol',       # 23 QGPU_data['restraint']['Upolx'] = Q_data['restraints'][3]
    'restraint-Ushell',     # 24 QGPU_data['restraint']['Ushell']= Q_data['restraints'][4]
    'restraint-Upres',      # 25 QGPU_data['restraint']['Upres'] = Q_data['restraints'][5]
    'total-Utot',           # 26 QGPU_data['total']['Utot']      = Q_data['SUM'][0]
    'total-Upot',           # 27 QGPU_data['total']['Upot']      = Q_data['SUM'][1]
    'total-Ukin',           # 28 QGPU_data['total']['Ukin']      = QGPU_data['total']['Ukin']
]

class bcolors:
    OKGREEN = '\033[92m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'    

refdata = {}
passed = True

def roundup(number):
    stringnumber = ('{:.2f}'.format(number))
    return(stringnumber)

def compare_energies(Q_data, QGPU_data):
    energies_Q6 = np.full((len(header)), np.nan)
    energies_QGPU = np.full((len(header)), np.nan)

    # 00
    energies_QGPU[0]    = QGPU_data['nonbonded']['pp'][0]
    energies_Q6[0]      = Q_data['solute'][0]
    QGPU_data['nonbonded']['pp'][0] = roundup(QGPU_data['nonbonded']['pp'][0])
    
    # 01
    energies_QGPU[1]    = QGPU_data['nonbonded']['pp'][1]
    energies_Q6[1]      = Q_data['solute'][1]
    QGPU_data['nonbonded']['pp'][1] = roundup(QGPU_data['nonbonded']['pp'][1])

    # 02
    energies_QGPU[2]    = QGPU_data['nonbonded']['qx'][0]
    if 'Q-atom' in Q_data:
        energies_Q6[2]      = Q_data['Q-atom'][0]
    QGPU_data['nonbonded']['qx'][0] = roundup(QGPU_data['nonbonded']['qx'][0])

    # 03
    energies_QGPU[3]    = QGPU_data['nonbonded']['qx'][1]
    if 'Q-atom' in Q_data:
        energies_Q6[3]      = Q_data['Q-atom'][1]
    QGPU_data['nonbonded']['qx'][1] = roundup(QGPU_data['nonbonded']['qx'][1])

    # 04
    energies_QGPU[4]    = QGPU_data['bonded']['qp'][0]
    if 'Q-atom' in Q_data:
        energies_Q6[4]      = Q_data['Q-atom'][2]   
    QGPU_data['bonded']['qp'][0]    = roundup(QGPU_data['bonded']['qp'][0])

    # 05
    energies_QGPU[5]    = QGPU_data['bonded']['qp'][1]
    if 'Q-atom' in Q_data:
        energies_Q6[5]      = Q_data['Q-atom'][3]   
    QGPU_data['bonded']['qp'][1]    = roundup(QGPU_data['bonded']['qp'][1])

    # 06
    energies_QGPU[6]    = QGPU_data['bonded']['qp'][2]
    if 'Q-atom' in Q_data:
        energies_Q6[6]      = Q_data['Q-atom'][4]   
    QGPU_data['bonded']['qp'][2]    = roundup(QGPU_data['bonded']['qp'][2])

    # 07
    energies_QGPU[7]    = QGPU_data['bonded']['qp'][3]
    if 'Q-atom' in Q_data:
        energies_Q6[7]      = Q_data['Q-atom'][5]   
    QGPU_data['bonded']['qp'][3]    = roundup(QGPU_data['bonded']['qp'][3])

    # 08
    energies_QGPU[8]    = QGPU_data['nonbonded']['pw'][0]
    energies_Q6[8]      = Q_data['solute-solvent'][0]   
    QGPU_data['nonbonded']['pw'][0] = roundup(QGPU_data['nonbonded']['pw'][0])

    # 09
    energies_QGPU[9]    = QGPU_data['nonbonded']['pw'][1]
    energies_Q6[9]      = Q_data['solute-solvent'][1]
    QGPU_data['nonbonded']['pw'][1] = roundup(QGPU_data['nonbonded']['pw'][1])

    # 10
    energies_QGPU[10]    = QGPU_data['nonbonded']['ww'][0]
    if 'solvent' in Q_data:
        energies_Q6[10]      = Q_data['solvent'][0]
    QGPU_data['nonbonded']['ww'][0] = roundup(QGPU_data['nonbonded']['ww'][0])

    # 11
    energies_QGPU[11]    = QGPU_data['nonbonded']['ww'][1]
    if 'solvent' in Q_data:
        energies_Q6[11]      = Q_data['solvent'][1]   
    QGPU_data['nonbonded']['ww'][1] = roundup(QGPU_data['nonbonded']['ww'][1])

    # 12
    energies_QGPU[12]    = QGPU_data['bonded']['w'][0]
    if 'solvent' in Q_data:
        energies_Q6[12]      = Q_data['solvent'][2]   
    QGPU_data['bonded']['w'][0]     = roundup(QGPU_data['bonded']['w'][0])

    # 13
    energies_QGPU[13]    = QGPU_data['bonded']['w'][1]
    if 'solvent' in Q_data:
        energies_Q6[13]      = Q_data['solvent'][3]
    QGPU_data['bonded']['w'][1]     = roundup(QGPU_data['bonded']['w'][1])

    # 14
    energies_QGPU[14]    = QGPU_data['bonded']['w'][2]
    if 'solvent' in Q_data:
        energies_Q6[14]      = Q_data['solvent'][4]
    QGPU_data['bonded']['w'][2]     = roundup(QGPU_data['bonded']['w'][2])

    # 15
    energies_QGPU[15]    = QGPU_data['bonded']['w'][3]
    if 'solvent' in Q_data:
        energies_Q6[15]      = Q_data['solvent'][5]
    QGPU_data['bonded']['w'][3]     = roundup(QGPU_data['bonded']['w'][3])

    # 16
    energies_QGPU[16]    = QGPU_data['bonded']['p'][0]
    energies_Q6[16]      = Q_data['solute'][2]   
    QGPU_data['bonded']['p'][0]     = roundup(QGPU_data['bonded']['p'][0])

    # 
    energies_QGPU[17]    = QGPU_data['bonded']['p'][1]
    energies_Q6[17]      = Q_data['solute'][3]   
    QGPU_data['bonded']['p'][1]     = roundup(QGPU_data['bonded']['p'][1])

    # 18
    energies_QGPU[18]    = QGPU_data['bonded']['p'][2]
    energies_Q6[18]      = Q_data['solute'][4]
    QGPU_data['bonded']['p'][2]     = roundup(QGPU_data['bonded']['p'][2])

    # 19
    energies_QGPU[19]    = QGPU_data['bonded']['p'][3] 
    energies_Q6[19]      = Q_data['solute'][5]   
    QGPU_data['bonded']['p'][3]     = roundup(QGPU_data['bonded']['p'][3])

    # 20
    energies_QGPU[20]    = QGPU_data['restraint']['Total']
    energies_Q6[20]      = Q_data['restraints'][0]   
    QGPU_data['restraint']['Total'] = roundup(QGPU_data['restraint']['Total'])

    # 21
    energies_QGPU[21]    = QGPU_data['restraint']['Ufix']
    energies_Q6[21]      = Q_data['restraints'][1]
    QGPU_data['restraint']['Ufix']  = roundup(QGPU_data['restraint']['Ufix'])

    # 22
    energies_QGPU[22]    = QGPU_data['restraint']['Uradx']
    energies_Q6[22]      = Q_data['restraints'][2]   
    QGPU_data['restraint']['Uradx'] = roundup(QGPU_data['restraint']['Uradx'])

    # 23
    energies_QGPU[23]    = QGPU_data['restraint']['Upolx']
    energies_Q6[23]      = Q_data['restraints'][3]   
    QGPU_data['restraint']['Upolx'] = roundup(QGPU_data['restraint']['Upolx'])

    # 24
    energies_QGPU[24]    = QGPU_data['restraint']['Ushell']
    energies_Q6[24]      = Q_data['restraints'][4]   
    QGPU_data['restraint']['Ushell']= roundup(QGPU_data['restraint']['Ushell'])

    # 25
    energies_QGPU[25]    = QGPU_data['restraint']['Upres']
    energies_Q6[25]      = Q_data['restraints'][5]   
    QGPU_data['restraint']['Upres'] = roundup(QGPU_data['restraint']['Upres'])

    # 26
    energies_QGPU[26]    = QGPU_data['total']['Utot']
    energies_Q6[26]      = Q_data['SUM'][0]   
    QGPU_data['total']['Utot']      = roundup(QGPU_data['total']['Utot'])

    # 27
    energies_QGPU[27]    = QGPU_data['total']['Upot']
    energies_Q6[27]      = Q_data['SUM'][1]   
    QGPU_data['total']['Upot']      = roundup(QGPU_data['total']['Upot'])

    #28
    energies_QGPU[28]    = QGPU_data['total']['Ukin']
    energies_Q6[28]      = QGPU_data['total']['Ukin']   
    QGPU_data['total']['Ukin']      = roundup(QGPU_data['total']['Ukin'])

    passed = True
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

    if 'Q-atom' in Q_data:
        if Q_data['Q-atom'][0] != QGPU_data['nonbonded']['qx'][0]:
            print("Energy not matching for {} {} {}, Q5 {} Q7 {}".format('nonbonded',
                                                                    'qx',
                                                                    'el',
                                                                    Q_data['Q-atom'][0],
                                                                    QGPU_data['nonbonded']['qx'][0],
                                                                    ))
            passed = False
            
        if Q_data['Q-atom'][1] != QGPU_data['nonbonded']['qx'][1]:
            print("Energy not matching for {} {} {}, Q5 {} Q7 {}".format('nonbonded',
                                                                    'qx',
                                                                    'vdw',
                                                                    Q_data['Q-atom'][1],
                                                                    QGPU_data['nonbonded']['qx'][1],
                                                                    ))
            passed = False

        # bonded interactions solute
        if Q_data['Q-atom'][2] != QGPU_data['bonded']['qp'][0]:
            print("Energy not matching for {} {} {}, Q5 {} Q7 {}".format('bonded',
                                                                    'qp',
                                                                    'bond',
                                                                    Q_data['Q-atom'][2],
                                                                    QGPU_data['bonded']['qp'][0],
                                                                    ))        
            passed = False
            
        if Q_data['Q-atom'][3] != QGPU_data['bonded']['qp'][1]:
            print("Energy not matching for {} {} {}, Q5 {} Q7 {}".format('bonded',
                                                                    'qp',
                                                                    'angle',
                                                                    Q_data['Q-atom'][3],
                                                                    QGPU_data['bonded']['qp'][1],
                                                                    ))        
            passed = False
            
        if Q_data['Q-atom'][4] != QGPU_data['bonded']['qp'][2]:
            print("Energy not matching for {} {} {}, Q5 {} Q7 {}".format('bonded',
                                                                    'qp',
                                                                    'torsion',
                                                                    Q_data['Q-atom'][4],
                                                                    QGPU_data['bonded']['qp'][2],
                                                                    ))       
            passed = False
            
        if Q_data['Q-atom'][5] != QGPU_data['bonded']['qp'][3]:
            print("Energy not matching for {} {} {}, Q5 {} Q7 {}".format('bonded',
                                                                    'qp',
                                                                    'improper',
                                                                    Q_data['Q-atom'][5],
                                                                    QGPU_data['bonded']['qp'][3],
                                                                    ))           
            passed = False 
        
    if Q_data['solute-solvent'][0] != QGPU_data['nonbonded']['pw'][0]:
        print("Energy not matching for {} {} {}, Q5 {} Q7 {}".format('nonbonded',
                                                                'pw',
                                                                'vdw',
                                                                Q_data['solute-solvent'][0],
                                                                QGPU_data['nonbonded']['pw'][0],
                                                                ))            
        passed = False
        
    if Q_data['solute-solvent'][1] != QGPU_data['nonbonded']['pw'][1]:
        print("Energy not matching for {} {} {}, Q5 {} Q7 {}".format('nonbonded',
                                                                'pw',
                                                                'vdw',
                                                                Q_data['solute-solvent'][1],
                                                                QGPU_data['nonbonded']['pw'][1],
                                                                ))            
        passed = False

    if 'solvent' in Q_data:    
        if Q_data['solvent'][0] != QGPU_data['nonbonded']['ww'][0]:
            print("Energy not matching for {} {} {}, Q5 {} Q7 {}".format('nonbonded',
                                                                    'ww',
                                                                    'el',
                                                                    Q_data['solvent'][0],
                                                                    QGPU_data['nonbonded']['ww'][0],
                                                                    ))     
            passed = False
            
        if Q_data['solvent'][1] != QGPU_data['nonbonded']['ww'][1]:
            print("Energy not matching for {} {} {}, Q5 {} Q7 {}".format('nonbonded',
                                                                    'ww',
                                                                    'vdw',
                                                                    Q_data['solvent'][1],
                                                                    QGPU_data['nonbonded']['ww'][1],
                                                                    ))         
            passed = False   

        # bonded interactions solvent
        if Q_data['solvent'][2] != QGPU_data['bonded']['w'][0]:
            print("Energy not matching for {} {} {}, Q5 {} Q7 {}".format('bonded',
                                                                    'w',
                                                                    'bond',
                                                                    Q_data['solvent'][2],
                                                                    QGPU_data['bonded']['w'][0],
                                                                    ))        
            passed = False
            
        if Q_data['solvent'][3] != QGPU_data['bonded']['w'][1]:
            print("Energy not matching for {} {} {}, Q5 {} Q7 {}".format('bonded',
                                                                    'w',
                                                                    'angle',
                                                                    Q_data['solvent'][3],
                                                                    QGPU_data['bonded']['w'][1],
                                                                    ))      
            passed = False
            
        if Q_data['solvent'][4] != QGPU_data['bonded']['w'][2]:
            print("Energy not matching for {} {} {}, Q5 {} Q7 {}".format('bonded',
                                                                    'w',
                                                                    'torsion',
                                                                    Q_data['solvent'][4],
                                                                    QGPU_data['bonded']['w'][2],
                                                                    ))          
            passed = False
            
        if Q_data['solvent'][5] != QGPU_data['bonded']['w'][3]:
            print("Energy not matching for {} {} {}, Q5 {} Q7 {}".format('bonded',
                                                                    'w',
                                                                    'improper',
                                                                    Q_data['solvent'][5],
                                                                    QGPU_data['bonded']['w'][3],
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
                                                                QGPU_data['restraint']['Total'],
                                                                ))        
        passed = False
        
    if Q_data['restraints'][1] != QGPU_data['restraint']['Ufix']:
        print("Energy not matching for {} {}, Q5 {} Q7 {}".format('restraint',
                                                                'Ufix',
                                                                Q_data['restraints'][1],
                                                                QGPU_data['restraint']['Ufix'],
                                                                ))      
        passed = False    
        
    if Q_data['restraints'][2] != QGPU_data['restraint']['Uradx']:
        print("Energy not matching for {} {}, Q5 {} Q7 {}".format('restraint',
                                                                'Uradx',
                                                                Q_data['restraints'][2],
                                                                QGPU_data['restraint']['Uradx'],
                                                                ))        
        passed = False
        
    if Q_data['restraints'][3] != QGPU_data['restraint']['Upolx']:
        print("Energy not matching for {} {}, Q5 {} Q7 {}".format('restraint',
                                                                'Upolx',
                                                                Q_data['restraints'][3],
                                                                QGPU_data['restraint']['Upolx'],
                                                                ))            
        
        passed = False
        
    if Q_data['restraints'][4] != QGPU_data['restraint']['Ushell']:
        print("Energy not matching for {} {}, Q5 {} Q7 {}".format('restraint',
                                                                'Ushell',
                                                                Q_data['restraints'][4],
                                                                QGPU_data['restraint']['Ushell'],
                                                                ))       
        passed = False
            
    if Q_data['restraints'][5] != QGPU_data['restraint']['Upres']:
        print("Energy not matching for {} {}, Q5 {} Q7 {}".format('restraint',
                                                                'Upres',
                                                                Q_data['restraints'][5],
                                                                QGPU_data['restraint']['Upres'],
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

        
    return passed, energies_Q6, energies_QGPU
