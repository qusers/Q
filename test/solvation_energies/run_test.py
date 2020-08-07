import os
import glob
import shutil

clean = False
previous = None
replicates = range(0,10)

def cleanup():
    print("Cleaning up")
    if os.path.isdir('results') == True:
        os.system('rm -r results')

if clean == True:
    cleanup()            

def runFEP():    
    os.mkdir('results')
    os.chdir('results')

    EQs = glob.glob('../FEP_ethane-ethanol/inputfiles/eq*inp')
    MDs = glob.glob('../FEP_ethane-ethanol/inputfiles/md*inp')

    EQs = sorted(EQs)
    MDs = sorted(MDs)[::-1]

    for eq in EQs:
        root = eq.split('/')[-1].split('.')[0]
        #os.mkdir(root)
        # Do not use restart in first file
        shutil.copy(eq,root + '.inp')
        if root == 'eq1':
            call = 'python ../../../../bin/qdyn.py -t ../dualtop.top -m {}.inp -f ../FEP1.fep -d {}'.format(root,root)

        else:
            call = 'python ../../../../bin/qdyn.py -t ../dualtop.top -m {}.inp -f ../FEP1.fep -d {} -r {}/dualtop/output/'.format(root,root,previous)

        os.system(call)

        previous = root

    for md in MDs:
        root = md.split('/')[-1].split('.')[0]
        #os.mkdir(root)
        # Do not use restart in first file
        shutil.copy(md,root + '.inp')
        call = 'python ../../../../bin/qdyn.py -t ../dualtop.top -m {}.inp -f ../FEP1.fep -d {} -r {}/dualtop/output/'.format(root,root,previous)

        os.system(call)

        previous = root

def gather_FEPs():
    FEPs = glob.glob('*/FEP_*')
    for FEP in FEPs:
        print(FEP)
        
    for i in replicates:
        print(i)
        
        
gather_FEPs()