import os
curdir = os.getcwd()
os.chdir('2.vacuum/FEP_SER-ALA')
os.system('sh submit.sh')
os.chdir('../../')
os.chdir('1.water/FEP_SER-ALA')
os.system('sh submit.sh')
os.chdir('../../')
