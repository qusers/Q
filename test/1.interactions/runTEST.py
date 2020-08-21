import os
import shlex
import subprocess

class Startup(object):
    def __init__(self,arch):
        self.arch =   arch
        
        tests = [['p-p/benzene/','python ../../../../bin/qdyn.py -t Q5_data/benzene-vacuum.top -m Q5_data/eq1.inp -d TEST -r Q5_data/'],
                 ['q-p/benzene-Na/FEP_benzene','python ../../../../../bin/qdyn.py -t Q5_data/Na-benzene-vacuum.top -m Q5_data/eq1.inp -d TEST -f Q5_data/FEP1.fep -r Q5_data/'],
                 ['q-p/benzene-Na/FEP_Na','python ../../../../../bin/qdyn.py -t Q5_data/Na-benzene-vacuum.top -m Q5_data/eq1.inp -d TEST -f Q5_data/FEP1.fep -r Q5_data/'],
                 ['q-p-w/benzene-Na/FEP_benzene','python ../../../../../bin/qdyn.py -t Q5_data/Na-benzene-water.top -m Q5_data/eq1.inp -d TEST -f Q5_data/FEP1.fep -r Q5_data/'],
                 ['q-p-w/benzene-Na/FEP_Na','python ../../../../../bin/qdyn.py -t Q5_data/Na-benzene-water.top -m Q5_data/eq1.inp -d TEST -f Q5_data/FEP1.fep -r Q5_data/'],
                 ['q-q/benzene/','python ../../../../bin/qdyn.py -t Q5_data/benzene-vacuum.top -m Q5_data/eq1.inp -d TEST -r Q5_data -f Q5_data/FEP1.fep'],
                 ['w-p/benzene/','python ../../../../bin/qdyn.py -t Q5_data/benzene-water.top -m Q5_data/eq1.inp -d TEST -r Q5_data/'],
                 ['w-q/benzene/','python ../../../../bin/qdyn.py -t Q5_data/benzene-water.top -m Q5_data/eq1.inp -f Q5_data/FEP1.fep -d TEST -r Q5_data/'],
                 ['w-w/water/','python ../../../../bin/qdyn.py -t Q5_data/water.top -m Q5_data/eq1.inp -d TEST -r Q5_data/']]
        
        self.curdir = os.getcwd()

        for test in tests:
            print("======================================================")
            print("Running test: {}".format(test[0].split('/')[0]))
            os.chdir(test[0])
            subprocess.check_output(['tar', '-xvf', 'testfiles.tar.gz'])
            args = shlex.split(test[1])
            if self.arch == 'gpu':
                args.append('--gpu')
            print(' '.join(args))
            out = subprocess.check_output(args)
            out2 = subprocess.check_output(['python','compare.py'])
            print(out2.decode('utf-8'))
            print("======================================================")    
            os.system('rm -r TEST Q5_data compare.py')
            os.chdir(curdir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='Test',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '       == Test == ')
    
    parser.add_argument('--version', 
                        action='version', 
                        version='%(prog)s 0.1.0')

    parser.add_argument('-a', '--architecture',
                        dest = "arch",
                        required = True,                        
                        choices = ['cpu','gpu'],                        
                        help = "Run tests with either GPU or CPU architecture"
                        )

    args = parser.parse_args()
    
    Startup(arch = args.arch,
           )