import os
import shlex
import math
from subprocess import check_output
import glob

executable = 'python /home/diazalej/private/software/qligfep_Nicolas/analyze_FEP.py '
options = {'-C':'RACKHAM',
           '-F':None}

curdir = os.getcwd()
failed = []
legdirs = glob.glob('*/')
results = {}
fepfiles = []

for i, topdir in enumerate(legdirs):
    results[topdir] = {}
    os.chdir(topdir)
    with open(topdir.split('/')[0] + '.txt', 'w') as outfile:
        outfile.write('{:<50}{:<10}{:<10}{:<10}{:<10}\n'.format('FEP','dG','sem', 'std', 'crashes'))
        for fep in glob.glob('FEP_*'):
            print('Analysing {}'.format(fep))
            options['-F']=fep
            option_list = ' '.join(['{} {}'.format(k,v) for k,v in options.items()])
            args = shlex.split(executable + option_list)
            try:
                out = check_output(args)
                out = out.splitlines()
                out = [line for line in out if '>>' in line.decode('utf-8')]
                out = out[-1].strip().decode("utf-8")
                out = out.split()
                outfile.write('{:<50}{:<10}{:<10}{:<10}{:<10}\n'.format(out[1],out[-4],out[-3], out[-2], out[-1]))
                fep = out[1]
                if i == 0:
                    fepfiles.append(fep)
                results[topdir][fep] = (float(out[-4]), float(out[-3]))
            except:
                print('Failed to read {}'.format(fep))
                failed.append('{} {}\n'.format(fep,topdir.split('.')[1]))

    os.chdir(curdir)
with open('failed.txt', 'w') as outfile:
    for line in failed:
         outfile.write(line)
print('\nddG summary...\n')
print('{:<50}{:<20}{:<30}'.format('FEP', 'ddG(kcal/mol)', 'SEM'))
for fep in fepfiles:
    rp, rw = results['1.protein/'], results['2.water/']
    sem_var = rp[fep][1]**2 + rw[fep][1]**2
    print('{:<50}{:<20.4}{:<30.4}'.format(fep, rp[fep][0] - rw[fep][0], math.sqrt(sem_var)))
