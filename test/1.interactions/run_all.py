import os
import shlex
import subprocess
import sys

# Check version
print("Checking Python version")
print("================================")
print("version is {}.{}.{}".format(sys.version_info[0],
                                   sys.version_info[1],
                                   sys.version_info[2],
                                  ))
print("================================")
    
if sys.version_info[0] < 3: 
    print("Q-GPU only works with Python 3.6 or higher")
    print("Exiting now")
    sys.exit()

if sys.version_info[0] == 3 and sys.version_info[1] < 6:
    print("Q-GPU only works with Python 3.6 or higher")
    print("Exiting now")    
    sys.exit()
    
executable = sys.executable

tests = [['p-p/benzene/',executable, '../../../../bin/qdyn.py -t Q5_data/benzene-vacuum.top -m Q5_data/eq1.inp -d TEST -r Q5_data/'],
         ['q-p/benzene-Na/FEP_benzene',executable, '../../../../../bin/qdyn.py -t Q5_data/Na-benzene-vacuum.top -m Q5_data/eq1.inp -d TEST -f Q5_data/FEP1.fep -r Q5_data/'],
         ['q-p/benzene-Na/FEP_Na',executable, '../../../../../bin/qdyn.py -t Q5_data/Na-benzene-vacuum.top -m Q5_data/eq1.inp -d TEST -f Q5_data/FEP1.fep -r Q5_data/'],
         ['q-p-w/benzene-Na/FEP_benzene',executable, '../../../../../bin/qdyn.py -t Q5_data/Na-benzene-water.top -m Q5_data/eq1.inp -d TEST -f Q5_data/FEP1.fep -r Q5_data/'],
         ['q-p-w/benzene-Na/FEP_Na',executable, '../../../../../bin/qdyn.py -t Q5_data/Na-benzene-water.top -m Q5_data/eq1.inp -d TEST -f Q5_data/FEP1.fep -r Q5_data/'],
         ['q-q/benzene/',executable, '../../../../bin/qdyn.py -t Q5_data/benzene-vacuum.top -m Q5_data/eq1.inp -d TEST -r Q5_data -f Q5_data/FEP1.fep'],
         ['w-p/benzene/',executable, '../../../../bin/qdyn.py -t Q5_data/benzene-water.top -m Q5_data/eq1.inp -d TEST -r Q5_data/'],
         ['w-q/benzene/',executable, '../../../../bin/qdyn.py -t Q5_data/benzene-water.top -m Q5_data/eq1.inp -f Q5_data/FEP1.fep -d TEST -r Q5_data/'],
         ['w-w/water/',executable, '../../../../bin/qdyn.py -t Q5_data/water.top -m Q5_data/eq1.inp -d TEST -r Q5_data/']]

curdir = os.getcwd()

for test in tests:
    print("======================================================")
    print("Running test: {}".format(test[0].split('/')[0]))
    os.chdir(test[0])
    subprocess.check_output(['tar', '-xvf', 'testfiles.tar.gz'])
    args = shlex.split(test[2])
    args = [test[1]] + args
    print(' '.join(args))
    out = subprocess.check_output(args)
    out2 = subprocess.check_output([executable,'compare.py'])
    print(out2.decode('utf-8'))
    print("======================================================")    
    os.system('rm -r TEST Q5_data compare.py')
    os.chdir(curdir)
