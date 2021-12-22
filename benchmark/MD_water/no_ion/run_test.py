import glob
import os

spheres=['10A', '15A', '20A', '25A', '30A']
curdir = os.getcwd()
cleanup = True

if cleanup == True:
    os.system('rm -r no_ion/*/water')
    os.system('rm no_ion/*/md.log')

for sphere in spheres:
    inputs = curdir + '/no_ion/' + sphere + '/water'
    os.chdir('no_ion/' + sphere)
    os.system('python ../../../../qdyn.py -t water.top')
    os.mkdir('water/output')
    os.system('../../../../qdyn ' + inputs + ' > md.log')
    os.chdir(curdir)