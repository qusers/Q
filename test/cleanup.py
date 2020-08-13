import os

curdir = os.getcwd()

print('cleaning w-w interaction folder')
os.chdir('1.interactions/w-w/')
os.system('rm -r compare.py md0.json Q5_data topo.json TEST')
os.chdir(curdir)

print('cleaning q-q interaction folder')
os.chdir('1.interactions/q-q/benzene')
os.system('rm -r compare.py eq1.inp FEP1.fep Q5_data/ TEST/')
os.chdir(curdir)
