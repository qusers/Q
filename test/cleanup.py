import os

curdir = os.getcwd()

os.chdir('1.interactions/w-w/')

os.system('rm -r compare.py md0.json Q5_data topo.json TEST')

os.chdir(curdir)
