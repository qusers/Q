import os

curdir = os.getcwd()

print('cleaning w-w interaction folder')
os.chdir('1.interactions/w-w/')
os.system('rm -r compare.py md0.json Q5_data topo.json TEST')
os.chdir(curdir)

print('cleaning q-q interaction folder')
os.chdir('1.interactions/q-q/benzene')
os.system('rm -r compare.py benzene-vacuum.json eq1.json FEP1.json Q5_data/ TEST/')
os.chdir(curdir)

print('cleaning p-p interaction folder')
os.chdir('1.interactions/p-p/benzene')
os.system('rm -r compare.py eq1.json benzene-vacuum.json Q5_data/ TEST/')
os.chdir(curdir)

print('cleaning q-p interaction folder')
os.chdir('1.interactions/q-p/benzene-Na/FEP_benzene/')
os.system('rm -r Na-benzene-vacuum.json compare.py eq1.json Q5_data/ FEP1.json TEST')
os.chdir('../')
os.chdir('FEP_Na')
os.system('rm -r Na-benzene-vacuum.json compare.py eq1.json Q5_data/ FEP1.json TEST')
os.chdir(curdir)

print('cleaning q-p-w interaction folder')
os.chdir('1.interactions/q-p-w/benzene-Na/FEP_benzene/')
os.system('rm -r Na-benzene-vacuum.json compare.py eq1.json Q5_data/ FEP1.json TEST')
os.chdir('../')
os.chdir('FEP_Na')
os.system('rm -r Na-benzene-vacuum.json compare.py eq1.json Q5_data/ FEP1.json TEST')
os.chdir(curdir)
