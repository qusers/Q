import glob
import os

QligFEP = '../../../app/QligFEP.py'
# construct FEPfiles
for fep in glob.glob('../1.Q-cpu/*/FEP_*'):
    of = fep
    wd = fep.split('/')[2]
    c  = 'KEBNE'

    print("Setting up FEPs for {}".format(fep))
    
    os.system('python {} -of {} -wd {} -c {}'.format(QligFEP,of,wd,c))
        
# write the submit file
with open('submit.py','w') as outfile:
    outfile.write("import os\n")
    outfile.write("curdir = os.getcwd()\n")
    for fepsubmit in glob.glob('*/FEP_*'):
        print(fepsubmit)
        outfile.write("os.chdir('{}')\n".format(fepsubmit))
        outfile.write("os.system('sbatch submit.sh')\n")
        outfile.write("os.chdir('../../')\n")
