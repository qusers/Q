import os
import shutil
import argparse
import sys
import time

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../share/')))

import settings
import IO
import md as MD
import defaults

class Create_Environment(object):
    """
        Creates the workdirectory environment.
    """
    def __init__(self,wd):        
        if not os.path.exists(wd):
            os.mkdir(wd)
            
        else:
            shutil.rmtree(wd)
            os.mkdir(wd)

class Run(object):
    """
        Prepares a Q topology using qprep (Q6.0.7 Fortran code)
        
    """    
    def __init__(self,wd,run):
        self.run = run
        self.spheres = ['10A','15A','20A','25A','30A']

        if run is not None:
            self.spheres = [x for x in self.spheres if any(r for r in self.run if x == r)]

        self.FF = 'OPLS2015'
        self.wd = wd
        self.qprep_inp = """ 
        rl {ROOT}data/ff/{FORCEFIELD}.lib
        rprm {ROOT}data/ff/{FORCEFIELD}.prm
        boundary 1 0.0 0.0 0.0 {SPHERERADIUS}
        solvate 0.0 0.0 0.0 {SPHERERADIUS} 1 HOH
        maketop water_{SPHERERADIUS}
        writetop water_{SPHERERADIUS}.top
        writepdb water_{SPHERERADIUS}.pdb y
        q
        """

        for sphere in self.spheres:
            self.sphere = sphere
            self.writeqprep()
            self.runqprep()
            self.create_MD_json()
            self.run_MD()

    def writeqprep(self):
        radius = self.sphere[:-1]
        formatted_qprep = self.qprep_inp.format(
                                                ROOT         = settings.ROOT,
                                                SPHERERADIUS = radius,
                                                FORCEFIELD   = self.FF)
        workdir = self.wd + '/' + self.sphere
        os.mkdir(workdir)
        qprepfile = workdir + '/qprep.inp'
        with open(qprepfile, 'w') as outfile:
            outfile.write(formatted_qprep)

    def runqprep(self):
        workdir = self.wd + '/' + self.sphere
        qprep = settings.ROOT + 'bin/qprep'
        curdir = os.getcwd()
        os.chdir(workdir)
        os.system(qprep + ' < qprep.inp > qprep.out')
        os.chdir(curdir)

    def create_MD_json(self):
        workdir = self.wd + '/' + self.sphere
        md_data = MD.MD().data
        default_values = defaults.MD
        for key in default_values:
            md_data[key] = default_values[key]

        md_data['steps'] = 1000
        md_data['temperature'] = 298
        md_data['stepsize'] = 1
        md_data['bath_coupling'] = 1
        md_data['initial_temperature'] = 298

        self.data = md_data
        MD.Write_MD.JSON(self,workdir + '/md.json')

    def run_MD(self):
        workdir = self.wd + '/' + self.sphere
        curdir = os.getcwd()
        os.chdir(workdir)
        qdyn_exec = 'python {}bin/qdyn.py '.format(settings.ROOT)

        ## CPU ##
        t = time.localtime()
        cur_time = time.strftime("%H:%M:%S", t)
        print("({}) Running {} size sphere on cpu:".format(cur_time,self.sphere))
        IO.run_command(qdyn_exec, '-t water_{}.top -m md.json -d test'.format(self.sphere[:-1]),runtime=True)

        ## GPU ##
        t = time.localtime()
        cur_time = time.strftime("%H:%M:%S", t)
        print("({}) Running {} size sphere on gpu:".format(cur_time,self.sphere))
        IO.run_command(qdyn_exec, '-t water_{}.top -m md.json -d test --gpu'.format(self.sphere[:-1]),runtime=True)        
        
        ## GO BACK
        os.chdir(curdir)

class Init(object):
    def __init__(self, data):
        """ Retrieves a dictionary of user input from qdyn:
               {'top'       :   top,
                'fep'       :   fep,
                'md'        :   md,
                're'        :   re,
                'wd'        :   wd,
                'verbose'   :   verbose
                'clean'   :   clean
               }
        """

        # INIT
        Create_Environment(
                           wd  = data['wd'],
                        )

        Run(
                           wd  = data['wd'],
                           run = data['run']
                        )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='runBenchmark',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '       == runBenchmark == ')
    
    parser.add_argument('--version', 
                        action='version', 
                        version='%(prog)s 0.1.0')

    parser.add_argument('-wd', '--workdir',
                        dest = "wd",
                        default = None,
                        required = True,
                        help = " Output folder")

    parser.add_argument('-r', '--run',
                        dest = "run",
                        required = False,
                        nargs = "+",
                        help = "Specify which tests to run")                        

    args = parser.parse_args() 

    START = Init(vars(args))