import glob
from subprocess import check_output
import shlex
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))

import settings

curdir = os.getcwd()
os.chdir(curdir)
opls2Q = settings.ROOT_DIR + '/opls2Q.py'
generate = 'python3 ' + opls2Q


for pdb in glob.glob('*.pdb'):
    name = pdb.split('.')[0]
    options = ' -l ' + name + ' -FF OPLS2015 -m'
    args = shlex.split(generate + options)
    out = check_output(args)
