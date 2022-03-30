import glob
import os

QligFEP = '../../../bin/QligFEP.py'
# construct FEPfiles
FEPs = glob.glob('../1.Q5/*/*')
#FEPs = [FEPs[0]]
for fep in FEPs:
    of = fep
    wd = fep.split('/')[2]
    c  = 'KEBNE'

    print("Setting up FEPs for {}".format(fep))
    print('python3 {} -of {} -wd {} -c {}'.format(QligFEP,of,wd,c))
    
    os.system('python3 {} -of {} -wd {} -c {}'.format(QligFEP,of,wd,c))
        
# write the submit file
with open('submit.py','w') as outfile:
    outfile.write("import os\n")
    outfile.write("curdir = os.getcwd()\n")
    for fepsubmit in glob.glob('*/*'):
        outfile.write("os.chdir('{}')\n".format(fepsubmit))
        outfile.write("os.system('sbatch submit.sh')\n")
        outfile.write("os.chdir('../../')\n")
