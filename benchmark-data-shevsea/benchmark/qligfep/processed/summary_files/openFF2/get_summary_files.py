import glob
import os
import shlex
import subprocess

wd = os.getcwd()
for target in glob.glob('*/'):
    os.chdir(os.path.join(wd, target))
    for leg in glob.glob('*/'):
        os.chdir(os.path.join(wd, target, leg))
        executable = 'rsync -av kebne:/proj/nobackup/uucompbiochem/diazalej/qligfep-benchmarking/openFF2/{}*{}*/analysis/*.txt .'.format(target, leg)
        #print(executable)
        #print(os.getcwd())
        run = subprocess.run(shlex.split(executable))
        os.chdir(wd)
