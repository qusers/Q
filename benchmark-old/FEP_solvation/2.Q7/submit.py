import os
curdir = os.getcwd()
os.chdir('2.vacuum/FEP_CYS-ALA')
os.system('sbatch submit.sh')
os.chdir('../../')
