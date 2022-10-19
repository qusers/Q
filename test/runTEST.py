import os
import shlex
import subprocess
import argparse
import sys
import shutil
import json
import numpy as np
from matplotlib import pyplot as plt

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../share/')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../env/')))

import settings
import IO
import topology as TOPOLOGY
import compare
import energy as ENERGY

class Create_Environment(object):
    """
        Creates the workdirectory environment.
    """
    def __init__(self,data):
        wd = data['curtest']        
        if not os.path.exists(wd):
            os.mkdir(wd)
            
        else:
            shutil.rmtree(wd)
            os.mkdir(wd)

class create_MD_input(object):
    def __init__(self,data):
        print('Generating MD input file')
        test = data['test']
        if data['shake'] == True:
            shake = 'on'

        else:
            shake = 'off'
        md_content = \
"""[MD]
steps                     {}
stepsize                  1
temperature               1
bath_coupling             1
random_seed               112
initial_temperature       1
shake_solvent             {}
shake_hydrogens           {}
shake_solute              {}
lrf                       off

[cut-offs]
solute_solvent            99.0
solute_solute             99.0
solvent_solvent           99.0
q_atom                    99.0
lrf                       99.0

[sphere]
shell_force               10.0
shell_radius              {}

[solvent]
radial_force              60.0
polarisation              on
polarisation_force        20.0
charge_correction         off

[intervals]
output                    1
non_bond                  1

[files]
topology                  {}{}
final                     eq1.re
""".format(data['timestep'],
           shake, shake, shake,
           data['testinfo'][data['test']][1],
           data['topdir'],
           data['testinfo'][data['test']][0])
        # Check if we are dealing with a FEP file
        if len(data['testinfo'][test]) == 3:
            fep_part = """fep                       {}{}

[lambdas]
1.000 0.000
""".format(data['inputdir'],data['testinfo'][data['test']][2])
            md_content = md_content + fep_part 

        with open('eq1.inp', 'w') as outfile:
            outfile.write(md_content)

class Run_Q6(object):
    def __init__(self,data):
        print("Running Q6")
        q_command = '{}bin/q6/qdyn_test eq1.inp > eq1.log'.format(settings.ROOT)
        os.system(q_command)

class Parse_Q6_data(object):
    def __init__(self,data):
        print("Parsing Q6 data")
        block = 0
        Q_energies = {}
        velocities = []

        with open('eq1.log') as infile:
            for line in infile:
                if len(line.strip()) < 2:
                    continue        
                
                if 'At_ID' in line and block == 0:
                    line = line.split()
                    velocities.append(line[3])
                    
                if 'Energy summary at step' in line:
                    step = line.split()[-2]
                    Q_energies[step] = {}
                    block = 1
                    continue

                if 'FINAL' in line:
                    # The last step in Q is our last step
                    step = data['timestep']
                    Q_energies[data['timestep']] = {}
                    block = 1
                    continue

                if 'Main loop starts here' in line:
                    block = 2
                    continue
                    
                if block == 1:
                    line = line.split()
                    if line[0] == 'solute':
                        Q_energies[step]['solute'] = line[1:]
                        
                    if line[0] == 'solvent':
                        Q_energies[step]['solvent'] = line[1:]
                                        
                    if line[0] == 'solute-solvent':
                        Q_energies[step]['solute-solvent'] = line[1:]        
                                        
                    if line[0] == 'restraints':
                        Q_energies[step]['restraints'] = line[1:]

                    if line[0] == 'Q-atom':
                        Q_energies[step]['Q-atom'] = line[1:]                     
                                                        
                    if line[0] == 'SUM':
                        Q_energies[step]['SUM'] = line[1:]

        jsondata = json.dumps(Q_energies, indent=1)
        with open('Q_data.json', 'w') as outfile:
            outfile.write(jsondata)

        with open('velocities.csv','w') as outfile:
            outfile.write('{}\n'.format(int(len(velocities)/3)))
            outfile.write('0\n')
            for i, v in enumerate(velocities):
                if (i + 1) % 3 != 0:
                    outfile.write('{};'.format(v))
                
                else:
                    outfile.write('{}\n'.format(v))

        # Parse the topology
        Qtopology = '{}{}'.format(data['topdir'],
                                  data['testinfo'][data['test']][0])
        read_top = TOPOLOGY.Read_Topology(Qtopology)
        top_data = read_top.Q()
        with open('coords.csv','w') as outfile:
            outfile.write('{}\n'.format(len(top_data['coords'])))
            outfile.write('{}\n'.format(top_data['natoms_solute']))            
            for line in top_data['coords']:
                outfile.write('{};{};{}\n'.format(line[0],line[1],line[2]))

class Run_QGPU(object):
    def __init__(self,data):
        print('Running QGPU')
        os.mkdir('tmp')
        shutil.copy('velocities.csv', 'tmp/velocities.csv')
        shutil.copy('coords.csv', 'tmp/coords.csv')
        args = [
                ' {}bin/qdyn.py'.format(settings.ROOT),
                '-t', '{}{}'.format(data['topdir'],
                                   data['testinfo'][data['test']][0]),
                '-m', 'eq1.inp',
                '-d', 'TEST',
                '-r', 'tmp'  
               ]

        # FEP file?
        if len(data['testinfo'][data['test']]) == 3:
            args.append('-f')
            args.append('{}{}'.format(data['inputdir'],data['testinfo'][data['test']][2]))

        if data['verbose'] == True:
            args.append('--verbose')

        if data['arch'] == 'gpu':
            args.append('--gpu')

        args_string = ' '.join(args)

        out = IO.run_command(data['executable'],args_string)
        if data['verbose'] == True:
            print(out.decode("utf-8"))
class bcolors:
    OKGREEN = '\033[92m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'    

class Compare(object):
    def __init__(self,data):
        total_energies_Q6 = []
        total_energies_QGPU = []
        top = '{}'.format(data['testinfo'][data['test']][0][:-4])
        energyfile = '{}/TEST/{}/output/energies.csv'.format(data['curtest'],top)

        #QGPU data
        read_ener  = ENERGY.Read_Energy(energyfile,0)
        QGPU_data = read_ener.QDYN()

        write_ener = ENERGY.Write_Energy(QGPU_data, 'test.json')
        write_ener.JSON()

        #Q6 data
        Q6_data_file = '{}/Q_data.json'.format(data['curtest'])
        with open(Q6_data_file) as infile:
            Q6_data = json.load(infile)

        QGPU_file = 'test.json'
        with open(QGPU_file) as infile:
            QGPU_data = json.load(infile)

        # loop over each step and run the comparison
        for key in Q6_data:
            if int(key) == len(Q6_data) - 1:
                continue
            #print('Comparing energies for frame {}'.format(key))
            Q6_tmp = Q6_data[key]
            QGPU_tmp = QGPU_data[int(key)]

            #print(Q6_tmp,QGPU_tmp)
            passed, energies_Q6, energies_QGPU = compare.compare_energies(Q6_tmp,QGPU_tmp)

            total_energies_Q6.append(energies_Q6)
            total_energies_QGPU.append(energies_QGPU)

            ## PRINT PASS ##    
            if passed == False:
                print('Compared energies for frame {}'.format(key))
                print('Passed test? ' + f"{bcolors.FAIL} FALSE {bcolors.ENDC}")
        
        if passed == True:
            print('Passed test? ' + f"{bcolors.OKGREEN} TRUE {bcolors.ENDC}")

        if data["avg"] == True:
            Q6_mean = np.mean(total_energies_Q6,axis=0)
            Q6_stdev = np.std(total_energies_Q6,axis=0)
            QGPU_mean = np.mean(total_energies_QGPU,axis=0)
            QGPU_stdev = np.std(total_energies_QGPU,axis=0)
            # formatted printing
            for i, headername in enumerate(compare.header):
                print('{} {:.2f} {:.2f} {:.2f} {:.2f}'.format(headername,
                                                              Q6_mean[i],
                                                              Q6_stdev[i],
                                                              QGPU_mean[i],
                                                              QGPU_stdev[i]))

        if data["plot"] == True:
            x = np.arange(0,len(total_energies_Q6))
            y1 = np.asarray(total_energies_Q6)[:, 26]
            y2 = np.asarray(total_energies_QGPU)[:, 26]
            plt.plot(x, y1,label='Q6')
            plt.plot(x, y2,label='QGPU')
            plt.xlabel('time (fs)')
            plt.ylabel('Utotal (kcal/mol)')
            plt.legend()
            plt.savefig(data['wd'] + '/Utot.png')

class Cleanup(object):
    def __init__(self,data):
        wd = data['curtest']  
        shutil.rmtree(wd)

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
        self.data = data
        self.data['curdir'] = os.getcwd()
        self.data['executable'] = sys.executable
        self.data['topdir'] = '{}test/data/topology/'.format(settings.ROOT)
        self.data['inputdir']   = '{}test/data/inputs/'.format(settings.ROOT)
        # Step = step + 1
        self.data['timestep'] = '{}'.format(int(self.data['timestep'])+1)

        if self.data['wd'] == None:
            self.data['wd'] = self.data['curdir'] + '/'
        if self.data['wd'][-1] != '/':
            self.data['wd'] = self.data['wd'] + '/'

        self.data['testinfo'] = {
                    'p-p'               : [
                                            'benzene-vacuum.top',
                                            '20'
                                          ],
                    'q-p_benzene'       : [
                                           'Na-benzene-vacuum.top',
                                           '20',
                                           'FEP_benzene.fep'
                                          ],
                    'q-p_Na'            : [
                                           'Na-benzene-vacuum.top',
                                           '20',
                                           'FEP_Na.fep'
                                          ],
                    'q-p-w_benzene'     : [
                                           'Na-benzene-water.top',
                                           '20',
                                           'FEP_benzene.fep'
                                          ],
                    'q-p-w_Na'          : [
                                           'Na-benzene-water.top',
                                           '20',
                                           'FEP_Na.fep'
                                          ],
                    'q-q'               : [
                                           'benzene-vacuum.top',
                                           '20',
                                           'FEP_benzene.fep'
                                          ],
                    'w-p'               : [
                                           'benzene-water.top',
                                           '20'
                                          ],
                    'w-q'               : [
                                           'benzene-water.top',
                                           '20',
                                           'FEP_benzene.fep'                                            
                                          ],
                    'w-w'               : [
                                           'water.top',
                                           '20'
                                          ],
                    'boundary'          : [
                                           'ala_wat.top',
                                           '14'
                                          ],
                    'polypeptide'       : [
                                           'ala_wat.top',
                                           '15'                       
                                          ],
                    'polypeptide25'     : [
                                           'ala_wat25.top',
                                           '25'
                                          ],

                    'q-q-large_vac'     : [
                                           'dualtop_vacuum.top',
                                           '22',
                                           'dualtop.fep'
                                          ]                                          
                }

        tests = data['testinfo'].keys()
        if self.data['run'] != None:
            tests = [x for x in tests if any(r for r in self.data['run'] if x == r)]
        for test in tests:
            print("\nRunning {}".format(test))
            self.data['test'] = test
            self.data['curtest'] = self.data['wd'] + test
            # INIT
            Create_Environment(self.data)
            
            # Running the actual code
            os.chdir(self.data['curtest'])
            create_MD_input(self.data)
            Run_Q6(self.data)
            Parse_Q6_data(self.data)
            Run_QGPU(self.data)
            failed = Compare(self.data)
            os.chdir(self.data['curdir'])
            print('\n')

            # Cleanup
            if self.data['tokeep'] == 'None':
                Cleanup(self.data)
            if self.data['tokeep'] == 'Failed' and failed == True:
                Cleanup(self.data)

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
                        
    parser.add_argument('-k', '--keep',
                        dest = "tokeep",
                        default = 'None',                        
                        choices = ['All','Failed', 'None'],                        
                        help = "Specify which files will be cleaned up after the test"
                        )

    parser.add_argument('--verbose',
                        action = 'store_true'
    )

    parser.add_argument('-r', '--run',
                        dest = "run",
                        required = False,
                        nargs = "+",
                        help = "Specify which tests to run")

    parser.add_argument('-w', '--wd',
                        dest = "wd",
                        required = False,
                        help = "Specify a working directory")


    parser.add_argument('-t', '--timestep',
                        dest = "timestep",
                        required = True,
                        help = "Specify the number of timesteps for each test")


    parser.add_argument('-s', '--shake',
                        dest = "shake",
                        default = False,
                        action = 'store_true',
                        help = "Turn shake on")

    parser.add_argument('--avg',
                        dest = "avg",
                        default = False,
                        action = 'store_true',
                        help = "Calculate average energy properties")


    parser.add_argument('--plot',
                        dest = "plot",
                        default = False,
                        action = 'store_true',
                        help = "Make a plot of the energies")

    args = parser.parse_args()
    
    START = Init(vars(args))
