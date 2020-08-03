import glob
import os

spheres=['10A', '15A', '20A', '25A', '30A']
curdir = os.getcwd()
cleanup = True
# fix to proper handling later
#numb_dens = 'python ../../../../qanalyse.py -t water.top -it water/output/coords.csv -ot MD -il ../../../../FF/OPLS2015.lib -c number_density'

def write_calc(sphere):
    with open('number_density.calc', 'w') as outfile:
        outfile.write('bins    5\n')
        outfile.write('residue HOH\n')

if cleanup == True:
    os.system('rm -r no_ion/*/water')
    os.system('rm no_ion/*/md.log')
    os.system('rm no_ion/*/*.out')

for sphere in spheres:
    inputs = curdir + '/no_ion/' + sphere + '/water'
    os.chdir('no_ion/' + sphere)
 #   write_calc(sphere)    
   os.system('../../../../bin/qdyn.py -t water.top')
    print('Running test {} '.format(sphere))
#    os.system(numb_dens + ' > numb_dens.out')
    os.chdir(curdir)
