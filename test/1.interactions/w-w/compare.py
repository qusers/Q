import json

refdata = {}

#el       vdW      bond     angle   torsion  improper

#total       fix slvnt_rad slvnt_pol     shell    solute

QtoGPUmask = {'solvent': {  'Ubond' :   2,
                            'Uangle':   3,
                            'Utor'  :   4,
                            'Uimp'  :   5,
                            'Ucoul' :   0,
                            'Uvdw'  :   1
                         },
              'restraints': { 'Urestr' :   0,
                              'Ufix'   :   1,
                              'Uradx'  :   2,
                              'Upolx'  :   3,
                              'Ushell' :   4,
                              'Upres'  :   5
                         },
              
              'SUM':     {  'Utot'  :   0,
                            'Upot'  :   1,
                            'Ukin'  :   2,
                         }              
             }

with open('Q5_data/Q_data.json') as infile:
    Q_data = json.load(infile)
    
with open('ref_data.txt') as infile:
    for line in infile:
        line = line.strip()
        line = line.split('=')
        if len(line) > 1:
            refdata[line[0].strip()] = line[1]

for etype in QtoGPUmask:
    for e in QtoGPUmask[etype]:
        value = round(float(refdata[e]),2)
        value = '{:.2f}'.format(value)
        Q_value = Q_data[etype][QtoGPUmask[etype][e]]
        if value != Q_value:
            print("\nEnergies not matching")
            print("{:5s} {:5s}:    (GPU: {:>10s}   Q5: {:>10s})".format(etype, e, value, Q_value))
