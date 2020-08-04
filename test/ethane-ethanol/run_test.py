import os
import glob

eq = glob.glob('FEP_ethane-ethanol/inputfiles/eq*inp')
md = glob.glob('FEP_ethane-ethanol/inputfiles/md*inp')

eq = sorted(eq)
md = sorted(md)[::-1]

#python ../../bin/qdyn.py -t dualtop.top -m eq1.inp -f FEP1.fep -d test
