import os
curdir = os.getcwd()
os.chdir('1.protein/FEP_1h1q-28')
os.system('sh submit.sh')
os.chdir('../../')
