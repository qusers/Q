import json

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
    QGPU_data['nonbonded']['pp'][0] = roundup(QGPU_data['nonbonded']['pp'][0])
    QGPU_data['nonbonded']['pp'][1] = roundup(QGPU_data['nonbonded']['pp'][1])
    QGPU_data['nonbonded']['qx'][0] = roundup(QGPU_data['nonbonded']['qx'][0])
    QGPU_data['nonbonded']['qx'][1] = roundup(QGPU_data['nonbonded']['qx'][1])
    QGPU_data['bonded']['qp'][0]    = roundup(QGPU_data['bonded']['qp'][0])
    QGPU_data['bonded']['qp'][1]    = roundup(QGPU_data['bonded']['qp'][1])
    QGPU_data['bonded']['qp'][2]    = roundup(QGPU_data['bonded']['qp'][2])
    QGPU_data['bonded']['qp'][3]    = roundup(QGPU_data['bonded']['qp'][3])
    QGPU_data['nonbonded']['pw'][0] = roundup(QGPU_data['nonbonded']['pw'][0])
    QGPU_data['nonbonded']['pw'][1] = roundup(QGPU_data['nonbonded']['pw'][1])
    QGPU_data['nonbonded']['ww'][0] = roundup(QGPU_data['nonbonded']['ww'][0])
    QGPU_data['nonbonded']['ww'][1] = roundup(QGPU_data['nonbonded']['ww'][1])
    QGPU_data['bonded']['w'][0]     = roundup(QGPU_data['bonded']['w'][0])
    QGPU_data['bonded']['w'][1]     = roundup(QGPU_data['bonded']['w'][1])
    QGPU_data['bonded']['w'][2]     = roundup(QGPU_data['bonded']['w'][2])
    QGPU_data['bonded']['w'][3]     = roundup(QGPU_data['bonded']['w'][3])
    QGPU_data['bonded']['p'][0]     = roundup(QGPU_data['bonded']['p'][0])
    QGPU_data['bonded']['p'][1]     = roundup(QGPU_data['bonded']['p'][1])
    QGPU_data['bonded']['p'][2]     = roundup(QGPU_data['bonded']['p'][2])
    QGPU_data['bonded']['p'][3]     = roundup(QGPU_data['bonded']['p'][3])
    QGPU_data['restraint']['Total'] = roundup(QGPU_data['restraint']['Total'])
    QGPU_data['restraint']['Ufix']  = roundup(QGPU_data['restraint']['Ufix'])
    QGPU_data['restraint']['Uradx'] = roundup(QGPU_data['restraint']['Uradx'])
    QGPU_data['restraint']['Upolx'] = roundup(QGPU_data['restraint']['Upolx'])
    QGPU_data['restraint']['Ushell']= roundup(QGPU_data['restraint']['Ushell'])
    QGPU_data['restraint']['Upres'] = roundup(QGPU_data['restraint']['Upres'])
    QGPU_data['total']['Utot']      = roundup(QGPU_data['total']['Utot'])
    QGPU_data['total']['Upot']      = roundup(QGPU_data['total']['Upot'])
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

        
    return passed
