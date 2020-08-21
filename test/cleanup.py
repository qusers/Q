import os
import shutil

tokeep = ['testfiles.tar.gz',
          'README.md',
          'TODO-ligand.binding.tar.gz',
          'oldQrunfiles.tar.gz',
          'setup.py'
         ]

topdirs = [ 
            '2.physical-properties/1.solvation/1.Q-cpu/',
            '2.physical-properties/1.solvation/2.Q-gpu/2.vacuum/',
            '2.physical-properties/1.solvation/2.Q-gpu/1.water'
          ]

for topdir in topdirs:
    if os.path.exists(topdir):
        shutil.rmtree(topdir)

rootdir = os.getcwd()
rootdirs = [rootdir + '/1.interactions', rootdir + '/2.physical-properties/']

for rootdir in rootdirs:
    for subdir, dirs, files in os.walk(rootdir):
        for file in files:
            if not file in tokeep:
                os.remove(os.path.join(subdir, file))
