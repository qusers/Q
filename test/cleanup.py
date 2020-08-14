import os

tokeep = ['testfiles.tar.gz',
          'README.md'
         ]


rootdir = os.getcwd()
rootdir = rootdir + '/1.interactions'

for subdir, dirs, files in os.walk(rootdir):
    for file in files:
        if not file in tokeep:
            os.remove(os.path.join(subdir, file))