import os
import shutil
import glob
import json

jsonfile = '/home/willem/projects/CCR2/Q/ligprep/denovo/2022-09-19_CCR2_denovo_run2/2022-09-21_CCR2_denovo_starmap.json'
ligandref = '/home/willem/projects/CCR2/Q/ligprep/denovo/aligned_ligands'

systems = ['protein', 'water']
cnt = 0
# Change this to where you installed QligFEP
setupFEP = 'python /home/willem/software/Q/qligfep-old/QligFEP.py'
cysbond = '1:SG_246:SG, 82:SG_159:SG, 267:SG_268:SG'

with open(jsonfile) as infile:
    data = json.load(infile)

for system in systems:
    cnt += 1
    directory = str(cnt) + '.' + system
    os.makedirs(directory, exist_ok = True)
    for pair in data['edges']:
        mol1 = pair['from']
        mol2 = pair['to']

        if system == 'water':
            call = setupFEP + ' -l1 ' + mol1 + ' -l2 ' + mol2 + ' -FF OPLS2015 -s water -c SNELLIUS -r 22 -l 0.5'

        if system == 'protein':
            call = setupFEP + ' -l1 ' + mol1 + ' -l2 ' + mol2 + ' -FF OPLS2015 -s protein -c SNELLIUS -r 22 -l 0.5' 
            print(call)

        src = 'FEP_' + mol1 + '-' + mol2
        dst = directory + '/FEP_' + mol1 + '-' + mol2
        os.system(call)
        shutil.move(src, dst)
