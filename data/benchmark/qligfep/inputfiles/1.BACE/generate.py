import glob
from subprocess import check_output
import shlex

for pdb in glob.glob('*.pdb'):
    name = pdb.split('.')[0]
    generate = 'python /home/jespers/software/setupFEP/generate_prms.py'
    options = ' -l ' + name + ' -FF OPLS2015 -o Q -m'
    args = shlex.split(generate + options)
    out = check_output(args)
    print out
