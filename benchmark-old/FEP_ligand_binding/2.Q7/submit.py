import os
curdir = os.getcwd()
os.chdir('1.protein/FEP_1h1q-28')
os.system('sbatch submit.sh')
os.chdir('../../')
os.chdir('2.water/FEP_1h1q-28')
os.system('sbatch submit.sh')
os.chdir('../../')
