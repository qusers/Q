import os
import argparse
import sys
import shutil
import json
import subprocess
import glob
import numpy as np
from matplotlib import pyplot as plt

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src/Qgpu')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../env/')))

import settings
import topology as TOPOLOGY
import md as MD
import fep as FEP
import restart as RESTART
import compare

lambdas = ['eq5', '0744_0256', '0998_0002']


def get_default_testinfo():
    return {
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
                                      ],
                'cdk2'              : [
                                       'cdk2.top',
                                       '22',
                                       'FEPm_cdk2.fep',
                                       'restraints_cdk2.inp'
                                      ],
                'thrombin'          : [
                                       'thrombin.top',
                                       '25',
                                       'FEPm_thrombin.fep',
                                       'restraints_thrombin.inp'
                                      ],
            }


def resolve_path(path, base_dir=None):
    if path is None:
        return None
    if os.path.isabs(path):
        return path
    if base_dir is not None:
        return os.path.abspath(os.path.join(base_dir, path))
    return os.path.abspath(path)

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

def clean_case_name(path, index):
    stem = os.path.splitext(os.path.basename(path))[0]
    safe = ''.join(c if c.isalnum() or c in ('-', '_') else '_' for c in stem)
    return '{}_{}'.format(index + 1, safe)

def parse_input_steps(path):
    section = None
    with open(path) as infile:
        for raw in infile:
            line = raw.split('!')[0].strip()
            if not line:
                continue
            if line.startswith('[') and line.endswith(']'):
                section = line[1:-1].strip().lower().replace('_', '-')
                continue
            if section == 'md':
                fields = line.split()
                if len(fields) >= 2 and fields[0].lower().replace('_', '-') == 'steps':
                    return int(fields[1])
    return 0

def read_direct_config(input_dir):
    config_file = os.path.join(input_dir, 'fep_config.json')
    if not os.path.exists(config_file):
        return {}
    with open(config_file) as infile:
        return json.load(infile)

def choose_fep_file(input_dir, requested):
    if requested is not None:
        return requested
    fep_files = sorted(os.path.basename(path) for path in glob.glob(os.path.join(input_dir, 'FEP*.fep')))
    if len(fep_files) == 1:
        return fep_files[0]
    if 'FEP1.fep' in fep_files:
        return 'FEP1.fep'
    return None

def prepare_direct_case_input(source_inp, case_dir, data):
    source_dir = os.path.dirname(source_inp) or '.'
    for filename in os.listdir(source_dir):
        src = os.path.join(source_dir, filename)
        if os.path.isfile(src):
            shutil.copy2(src, os.path.join(case_dir, filename))

    config = read_direct_config(source_dir)
    fep_file = choose_fep_file(source_dir, data.get('fep_file'))
    seed = str(data.get('seed') if data.get('seed') is not None else 1)
    temperature = str(data.get('temperature') if data.get('temperature') is not None else config.get('temperature', 298))

    case_inp = os.path.join(case_dir, os.path.basename(source_inp))
    with open(case_inp) as infile:
        content = infile.read()

    if 'FEP_VAR' in content:
        if fep_file is None:
            raise RuntimeError('Input uses FEP_VAR but no FEP file was specified or found.')
        content = content.replace('FEP_VAR', fep_file)
    content = content.replace('SEED_VAR', seed)
    content = content.replace('T_VAR', temperature)

    with open(case_inp, 'w') as outfile:
        outfile.write(content)
    return case_inp

class create_MD_input(object):
    def __init__(self,data):
        print('Generating MD input file')
        test = data['test']
        if data['shake']:
            shake = 'on'

        else:
            shake = 'off'

        _lambda = None
        _inv_lambda = None
        # Check if a lambda has been specified

        if data.get('fep_path') is not None and data['lambda'] is not None:
            if not data['lambda'].startswith('eq'):
                str_lambda = data['lambda'].split("_")[0]
                str_inv_lambda = data['lambda'].split("_")[1]
                _lambda = str_lambda[0] + "." + str_lambda[1:4]
                _inv_lambda = str_inv_lambda[0] + "." + str_inv_lambda[1:4]

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
energy                    0
trajectory                0

[files]
topology                  {}
final                     eq1.re
""".format(data['timestep'],
           shake, shake, shake,
           data['testinfo'][data['test']][1],
           data['topology_path'])
        if data.get('fep_path') is not None:
            filename = data['fep_path']
            fep_name = os.path.basename(filename)

            fep_part = """fep                       {}{}

[lambdas]
""".format("" if os.path.isabs(filename) else "",
           filename)
            if _lambda is not None:
                fep_part += _lambda + " " + _inv_lambda + "\n"
            else:
                if fep_name.startswith("FEPm"):
                    fep_part += "0.500 0.500\n"
                else:
                    fep_part += "1.000 0.000\n"
            md_content = md_content + fep_part
        # Check if there are boundary conditions
        if data.get('restraints_path') is not None:
            filename = data['restraints_path']
            with open(filename, 'r') as f:
                restraint_part = f.read()
                md_content = md_content + restraint_part

        with open('eq1.inp', 'w') as outfile:
            outfile.write(md_content)

class Run_Q6(object):
    def __init__(self,data):
        print("Running Q6")
        qdyn = '{}src/q6/bin/q6/qdyn_test'.format(settings.ROOT)
        args = [qdyn, data.get('q_input_file', 'eq1.inp')]
        run_dir = data.get('q_input_dir', os.getcwd())
        log_path = data.get('q6_log', 'eq1.log')
        with open(log_path, 'w') as outfile:
            result = subprocess.run(args, cwd=run_dir, stdout=outfile, stderr=subprocess.STDOUT)
        if result.returncode != 0:
            raise RuntimeError('Q6 failed with exit code {}'.format(result.returncode))
        with open(log_path) as infile:
            q6_log = infile.read()
        if 'ABNORMAL TERMINATION' in q6_log or 'STOP qdyn terminated abnormally' in q6_log:
            reason = 'unknown reason'
            for line in q6_log.splitlines():
                if '>>> Error:' in line:
                    reason = line.strip()
                    break
            raise RuntimeError('Q6 terminated abnormally: {}'.format(reason))

class Parse_Q6_data(object):
    def __init__(self,data):
        print("Parsing Q6 data")
        block = 0
        Q_energies = {}
        velocities = []
        data['q6_final_key'] = None

        with open('eq1.log') as infile:
            for line in infile:
                if len(line.strip()) < 2:
                    continue

                if 'At_ID' in line and block == 0:
                    line = line.split()
                    velocities.append(float(line[3]))

                if 'Energy summary at step' in line:
                    step = line.split()[-2]
                    Q_energies[step] = {}
                    block = 1
                    continue

                if 'FINAL  Energy summary' in line:
                    # The last step in Q is our last step. Match only the
                    # final energy-summary header; FINAL Q-atom energies must
                    # not clear the parsed summary frame.
                    step = data['timestep']
                    Q_energies[step] = {}
                    data['q6_final_key'] = step
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
        if velocities:
            if len(velocities) % 3 != 0:
                raise RuntimeError('Q6 velocity count is not divisible by 3.')
            data['q6_velocities'] = [
                [velocities[i], velocities[i + 1], velocities[i + 2]]
                for i in range(0, len(velocities), 3)
            ]
        else:
            data['q6_velocities'] = []

def prepare_qgpu_csv_input(data):
    csv_dir = 'qgpu_csv'
    if os.path.exists(csv_dir):
        shutil.rmtree(csv_dir)
    os.mkdir(csv_dir)
    os.mkdir(os.path.join(csv_dir, 'output'))

    top_data = TOPOLOGY.Read_Topology(data['topology_path']).Q()
    TOPOLOGY.Write_Topology(top_data).CSV(csv_dir + '/')

    md_data = MD.Read_MD('eq1.inp').Q()
    MD.Write_MD(md_data).CSV(csv_dir + '/')

    if data.get('fep_path') is not None:
        fep_data = FEP.Read_Fep(data['fep_path']).Q()
    else:
        fep_data = FEP.Read_Fep(None).EMPTY()
    FEP.Write_Fep(fep_data).CSV(csv_dir + '/')

    if len(data.get('q6_velocities', [])) != len(top_data['coords']):
        raise RuntimeError('Q6 velocity count does not match topology atom count.')
    RESTART.Write_Restart(data['q6_velocities'], top_data['coords']).CSV(csv_dir + '/')
    return csv_dir

class Run_QGPU(object):
    def __init__(self,data):
        print('Running QGPU ({})'.format(data['qgpu_input']))
        if data.get('direct_inp') and data['qgpu_input'] != 'inp':
            raise RuntimeError('Direct .inp tests support --qgpu-input inp only.')
        qdyn = '{}bin/qdyn'.format(settings.ROOT)
        args = [qdyn]
        if data['arch'] == 'gpu':
            args.append('--gpu')
        if data['qgpu_input'] == 'csv':
            args.extend(['--csv-dir', prepare_qgpu_csv_input(data)])
        else:
            args.append(data.get('q_input_file', 'eq1.inp'))

        run_dir = data.get('q_input_dir', os.getcwd())
        log_path = data.get('qgpu_log', 'qgpu.log')
        with open(log_path, 'w') as outfile:
            result = subprocess.run(args, cwd=run_dir, stdout=outfile, stderr=subprocess.STDOUT)
        if result.returncode != 0:
            raise RuntimeError('QGPU failed with exit code {}'.format(result.returncode))
        if data['verbose']:
            with open(log_path) as infile:
                print(infile.read())

class bcolors:
    OKGREEN = '\033[92m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'    

def empty_qgpu_energy():
    return {
        'temperature': {'Temp': None},
        'bonded': {'p': [None, None, None, None],
                   'w': [None, None, None, None],
                   'qp': [None, None, None, None]},
        'nonbonded': {'pp': [None, None],
                      'pw': [None, None],
                      'ww': [None, None],
                      'qx': [None, None]},
        'restraint': {'Uradx': None,
                      'Upolx': None,
                      'Ushell': None,
                      'Ufix': None,
                      'Upres': None,
                      'Total': None},
        'q-energies': {'lambda': [],
                       'SUM': [],
                       'Ubond': [],
                       'Uangle': [],
                       'Utor': [],
                       'Uimp': [],
                       'Uvdw': [],
                       'Ucoul': [],
                       'Urestr': []},
        'total': {'Ukin': None,
                  'Upot': None,
                  'Utot': None},
    }

def parse_qgpu_log(filename):
    frames = []
    current = None
    block = None

    with open(filename) as infile:
        for raw in infile:
            line = raw.strip()
            if not line:
                continue

            if (
                line == '== INITIAL ENERGIES'
                or line.startswith('== STEP ')
                or line.startswith('== FINAL ENERGIES')
            ):
                if current is not None:
                    frames.append(current)
                current = empty_qgpu_energy()
                block = None
                continue

            if current is None:
                continue

            if line in ('[temperature]', '[bonded]', '[nonbonded]', '[restraint]', '[q-energies]', '[total]'):
                block = line.strip('[]')
                continue

            if line.startswith('type') or line.startswith('lambda'):
                continue

            parts = line.split()
            if block == 'temperature' and len(parts) >= 2 and parts[0] == 'Temp':
                current['temperature'][parts[0]] = float(parts[1])
            elif block == 'bonded' and len(parts) >= 5 and parts[0] in current['bonded']:
                current['bonded'][parts[0]] = [float(x) for x in parts[1:5]]
            elif block == 'nonbonded' and len(parts) >= 3 and parts[0] in current['nonbonded']:
                current['nonbonded'][parts[0]] = [float(x) for x in parts[1:3]]
            elif block == 'restraint' and len(parts) >= 2 and parts[0] in current['restraint']:
                current['restraint'][parts[0]] = float(parts[1])
            elif block == 'q-energies' and len(parts) >= 9 and parts[0].replace('.', '', 1).isdigit():
                current['q-energies']['lambda'].append(float(parts[0]))
                current['q-energies']['SUM'].append(float(parts[1]))
                current['q-energies']['Ubond'].append(float(parts[2]))
                current['q-energies']['Uangle'].append(float(parts[3]))
                current['q-energies']['Utor'].append(float(parts[4]))
                current['q-energies']['Uimp'].append(float(parts[5]))
                current['q-energies']['Ucoul'].append(float(parts[6]))
                current['q-energies']['Uvdw'].append(float(parts[7]))
                current['q-energies']['Urestr'].append(float(parts[8]))
            elif block == 'total' and len(parts) >= 2 and parts[0] in current['total']:
                current['total'][parts[0]] = float(parts[1])

    if current is not None:
        frames.append(current)
    return frames

class Compare(object):
    def __init__(self,data):
        compare.ENERGY_TOLERANCE = data.get('tolerance', 0.0)
        self.failed = False
        total_energies_Q6 = []
        total_energies_QGPU = []

        #Q6 data
        Q6_data_file = '{}/Q_data.json'.format(data['curtest'])
        with open(Q6_data_file) as infile:
            Q6_data = json.load(infile)

        QGPU_data = parse_qgpu_log('{}/qgpu.log'.format(data['curtest']))
        with open('QGPU_data.json', 'w') as outfile:
            json.dump(QGPU_data, outfile, indent=2)

        final_key = data.get('q6_final_key')
        ordered_keys = sorted(Q6_data.keys(), key=lambda key: int(key) if str(key).isdigit() else 10**12)
        comparable_keys = [key for key in ordered_keys if key != final_key]
        if not comparable_keys:
            raise RuntimeError('No Q6 energy frames were parsed from {}'.format(Q6_data_file))

        # loop over each step and run the comparison
        compared = False
        for key in comparable_keys:
            compared = True
            #print('Comparing energies for frame {}'.format(key))
            Q6_tmp = Q6_data[key]
            frame_index = int(key) if str(key).isdigit() else len(total_energies_Q6)
            if frame_index >= len(QGPU_data):
                frame_index = len(total_energies_Q6)
            if frame_index >= len(QGPU_data):
                self.failed = True
                print('Missing QGPU frame for Q6 frame {}'.format(key))
                print('Passed test? ' + f"{bcolors.FAIL} FALSE {bcolors.ENDC}")
                continue
            QGPU_tmp = QGPU_data[frame_index]

            #print(Q6_tmp,QGPU_tmp)
            passed, energies_Q6, energies_QGPU = compare.compare_energies(Q6_tmp,QGPU_tmp)

            total_energies_Q6.append(energies_Q6)
            total_energies_QGPU.append(energies_QGPU)

            ## PRINT PASS ##    
            if passed is False:
                self.failed = True
                print('Compared energies for frame {}'.format(key))
                print('Passed test? ' + f"{bcolors.FAIL} FALSE {bcolors.ENDC}")
        
        if compared and not self.failed:
            print('Passed test? ' + f"{bcolors.OKGREEN} TRUE {bcolors.ENDC}")

        if data["avg"]:
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

        if data["plot"]:
            x = np.arange(0,len(total_energies_Q6))
            y1 = np.asarray(total_energies_Q6)[:, 26]
            y2 = np.asarray(total_energies_QGPU)[:, 26]
            plt.plot(x, y1,label='Q6')
            plt.plot(x, y2,label='QGPU')
            plt.xlabel('time (fs)')
            plt.ylabel('Utotal (kcal/mol)')
            plt.legend()
            plt.savefig(data['wd'] + '/Utot.png')

    def __bool__(self):
        return self.failed

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
        self.data['topdir'] = os.path.join(settings.ROOT, 'test/data/topology')
        self.data['inputdir'] = os.path.join(settings.ROOT, 'test/data/inputs')

        if self.data['wd'] is None:
            self.data['wd'] = self.data['curdir'] + '/'
        if self.data['wd'][-1] != '/':
            self.data['wd'] = self.data['wd'] + '/'

        self.data['testinfo'] = get_default_testinfo()

        if self.data['custom_top'] is not None:
            custom_info = [
                os.path.basename(self.data['custom_top']),
                self.data['custom_shell_radius']
            ]
            if self.data['custom_fep'] is not None:
                custom_info.append(os.path.basename(self.data['custom_fep']))
            if self.data['custom_restraints'] is not None:
                while len(custom_info) < 3:
                    custom_info.append(None)
                custom_info.append(os.path.basename(self.data['custom_restraints']))
            self.data['testinfo'][self.data['custom_name']] = custom_info

        if self.data['inp'] is not None:
            self.run_direct_inputs()
            return

        # Step = step + 1
        self.data['timestep'] = '{}'.format(int(self.data['timestep'])+1)

        tests = data['testinfo'].keys()
        if self.data['run'] is not None:
            tests = [x for x in tests if any(r for r in self.data['run'] if x == r)]
        for test in tests:
            print("\nRunning {}".format(test))
            self.data['test'] = test
            self.data['curtest'] = self.data['wd'] + test
            _topfile = data['testinfo'][data['test']][0]
            if len(data['testinfo'][test]) >= 3 and data['lambda'] is not None and test != self.data['custom_name']:
                _topfile = _topfile.split(".")[0] + "_" + data['lambda'] + "." + _topfile.split(".")[1]
            self.data['topfile'] = _topfile
            if test == self.data['custom_name'] and self.data['custom_top'] is not None:
                self.data['topology_path'] = self.data['custom_top']
                self.data['fep_path'] = self.data['custom_fep']
                self.data['restraints_path'] = self.data['custom_restraints']
            else:
                self.data['topology_path'] = os.path.join(self.data['topdir'], self.data['topfile'])
                self.data['fep_path'] = None
                self.data['restraints_path'] = None
                if len(data['testinfo'][test]) >= 3:
                    self.data['fep_path'] = os.path.join(self.data['inputdir'], data['testinfo'][test][2])
                if len(data['testinfo'][test]) >= 4:
                    self.data['restraints_path'] = os.path.join(self.data['inputdir'], data['testinfo'][test][3])
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
            if self.data['tokeep'] == 'Failed' and failed:
                Cleanup(self.data)

    def run_direct_inputs(self):
        if self.data['qgpu_input'] != 'inp':
            raise RuntimeError('Direct .inp tests support --qgpu-input inp only.')

        if self.data['wd'] is None:
            self.data['wd'] = self.data['curdir'] + '/'
        if self.data['wd'][-1] != '/':
            self.data['wd'] = self.data['wd'] + '/'

        for index, inp_file in enumerate(self.data['inp']):
            inp_path = os.path.abspath(inp_file)
            if not os.path.exists(inp_path):
                raise RuntimeError('Input file does not exist: {}'.format(inp_file))

            print("\nRunning {}".format(inp_path))
            self.data['direct_inp'] = True
            self.data['curtest'] = self.data['wd'] + clean_case_name(inp_path, index)

            Create_Environment(self.data)
            case_inp = prepare_direct_case_input(inp_path, self.data['curtest'], self.data)
            self.data['q_input_dir'] = self.data['curtest']
            self.data['q_input_file'] = os.path.basename(case_inp)
            steps = parse_input_steps(case_inp)
            self.data['timestep'] = str(steps) if steps > 0 else 'FINAL'
            self.data['q6_log'] = os.path.join(self.data['curtest'], 'eq1.log')
            self.data['qgpu_log'] = os.path.join(self.data['curtest'], 'qgpu.log')

            Run_Q6(self.data)
            Run_QGPU(self.data)

            os.chdir(self.data['curtest'])
            Parse_Q6_data(self.data)
            failed = Compare(self.data)
            os.chdir(self.data['curdir'])
            print('\n')

            if self.data['tokeep'] == 'None':
                Cleanup(self.data)
            if self.data['tokeep'] == 'Failed' and failed:
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
                        required = False,
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
    
    parser.add_argument('--lambda',
                        dest = "lambda",
                        default = None,
                        required = False,
                        help = "Specify a particular phase of the perturbation")

    parser.add_argument('--custom-top',
                        dest = "custom_top",
                        default = None,
                        required = False,
                        help = "Path to a custom topology file to add as a test")

    parser.add_argument('--custom-shell-radius',
                        dest = "custom_shell_radius",
                        default = '25',
                        required = False,
                        help = "Shell radius to use with --custom-top")

    parser.add_argument('--custom-fep',
                        dest = "custom_fep",
                        default = None,
                        required = False,
                        help = "Optional FEP file for --custom-top")

    parser.add_argument('--custom-restraints',
                        dest = "custom_restraints",
                        default = None,
                        required = False,
                        help = "Optional restraints file for --custom-top")

    parser.add_argument('--custom-name',
                        dest = "custom_name",
                        default = 'custom',
                        required = False,
                        help = "Test name to use with --custom-top")

    parser.add_argument('--tolerance',
                        dest = "tolerance",
                        type = float,
                        default = 0.0,
                        required = False,
                        help = "Energy comparison tolerance (default: 0.0 = exact match)")

    parser.add_argument('--qgpu-input',
                        dest = "qgpu_input",
                        default = 'inp',
                        choices = ['inp', 'csv'],
                        required = False,
                        help = "Choose QGPU input path: native .inp parser or generated CSV directory")

    parser.add_argument('--inp',
                        dest = "inp",
                        nargs = '+',
                        default = None,
                        required = False,
                        help = "Run and compare explicit Q .inp files directly")

    parser.add_argument('--fep-file',
                        dest = "fep_file",
                        default = None,
                        required = False,
                        help = "Replacement for FEP_VAR in direct --inp mode")

    parser.add_argument('--seed',
                        dest = "seed",
                        default = None,
                        required = False,
                        help = "Replacement for SEED_VAR in direct --inp mode")

    parser.add_argument('--temperature',
                        dest = "temperature",
                        default = None,
                        required = False,
                        help = "Replacement for T_VAR in direct --inp mode")

    args = parser.parse_args()
    if args.inp is None and args.timestep is None:
        parser.error('-t/--timestep is required unless --inp is used')

    data = vars(args)
    data['custom_top'] = resolve_path(data['custom_top'])
    data['custom_fep'] = resolve_path(data['custom_fep'])
    data['custom_restraints'] = resolve_path(data['custom_restraints'])

    START = Init(data)
