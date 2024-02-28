import argparse
import glob
import os
import re
from pathlib import Path
from itertools import product

import numpy as np

from .functions import avg_sem
from .IO import read_qfep, read_qfep_verbose, run_command
from .logger import logger
from .settings.settings import CLUSTER_DICT


def info_from_run_file(file_path: Path):
    """Extract the FEP temperature from a run file."""
    info = {'temperature': None, 'replicates': None}
    run_files = sorted(list(file_path.glob('run*.sh')))
    if len(run_files) == 0:
        logger.error(f'No run files found in {file_path}')
    elif len(run_files) > 1:
        logger.warning(f'Multiple run files found in {file_path}!! Using the first one.')
    run_file = run_files[0]
    temp_pattern = re.compile(r'temperature=(\d+)')
    replic_pattern = re.compile(r'run=(\d+)')
    with run_file.open('r', encoding='utf-8') as _file:
        for line in _file:
            temp_match = temp_pattern.search(line)
            replicate_match = replic_pattern.search(line)
            if temp_match:
                info['temperature'] = temp_match.group(1)
            if replicate_match:
                info['replicates'] = replicate_match.group(1)
    if any([v is None for v in info.values()]):
        logger.error(f'Could not extract temperature and/or replicates from {run_file}')
    return info

try:
    import matplotlib
    matplotlib.use('Agg')    
    import matplotlib.pyplot as plt
    plot = True
except ModuleNotFoundError:
    print('cannot import matplotlib, skipping plot generation')
    plot = False

class FepReader(object):
    """Class to analyze FEP output files. This wrapper class will read, structure the files
    and enable input/output operation with the analyzed data.
    """    
    def __init__(self, system) -> None:
        self.cwd = Path.cwd()
        self.data = {}
        self.system = None
        self.load_system(system)
        
    def load_system(self, system:str):
        """This method will load the FEP directories within a system directory. The reason
        for this is the way setupFEP is structured. When running it, it will create two
        different diretories for storing the perturbation; `1.water` and `2.protein`."""
        fep_dirs = sorted(list((self.cwd / system).glob('FEP_*')))
        self.data.update({system: {}})
        self.data[system].update({_dir.name : {'root': str(_dir.absolute())} for _dir in fep_dirs})
        self.system = system
        
    def load_new_system(self, system: str):
        """This method will load a new system into the data dictionary."""
        self.load_system(system)
        self.read_fep_inputs()
    
    def read_fep_inputs(self):
        for fep in self.data[self.system].keys():
            _dir = Path(self.data[self.system][fep]['root'])
            run_info = info_from_run_file(_dir / 'inputfiles')
            temperature = run_info['temperature']
            replicates = run_info['replicates']
            inputs = sorted(list(_dir.glob('inputfiles/md*.inp')))
            fep_files = sorted(list(_dir.glob('inputfiles/FEP*.fep')))
            fep_stages = []
            for fep_file in fep_files:
                fep_stages.append(fep_file.stem)
            if len(fep_stages) > 1: # Are there cases where we'll have more than 1?
                logger.warning(f'Multiple FEP files found in {fep}: {fep_files}!! Using the first...')
            fep_stage = fep_stages[0]
            lambda_sum = len(fep_files) * (len(inputs)-1)
            self.data[self.system][fep].update(
                {
                    'lambda_sum': lambda_sum,
                    'fep_stage': fep_stage,
                    'temperature': str(temperature),
                    'replicates': replicates
                }
            )
            
    def read_perturbations(self):
        methods_list = ['dG', 'dGf', 'dGr', 'dGos', 'dGbar']
        feps = [k for k in self.data[self.system].keys()]
        
        for fep in feps:
            fep_dict = self.data[self.system][fep]
            _dir = Path(fep_dict['root'])
            replicate_root = _dir / fep_dict['fep_stage'] / fep_dict['temperature']
            replicate_qfep_files = sorted( # here we use the int to sort the replicates
                list(replicate_root.glob('*/qfep.out')), key=lambda x: int(x.parent.name)
            )
            energies = {}
            method_results = {method: {} for method in methods_list}
            failed_replicates = []
            all_replicates = [i for i in range(1, int(fep_dict['replicates']) + 1)]
            stage = self.data[self.system][fep]['fep_stage']
            for rep in replicate_qfep_files:
                repID = int(rep.parent.name)
                try:
                    energies[repID] = read_qfep(rep)
                except OSError as e:
                    logger.warning(
                        f"Failed to retrieve energies for: {fep}, {stage} - rep.{repID}. Error: \n{e}"
                    )
                    failed_replicates.append(repID)
                    energies[repID] = np.array([np.nan] * len(methods_list))  # Assuming 5 methods
                    
            # per different type of energy, populate the methods dictionary
            for mname in methods_list:
                method_idx = methods_list.index(mname)
                method_energies = np.array([energies[repID][method_idx] for repID in all_replicates])
                print('energies', method_energies)
                
                method_results[mname] = {
                    'energies': method_energies.tolist(),
                    'avg': np.nanmean(method_energies),
                    'sem': np.nanstd(method_energies) / np.sqrt(method_energies.shape)
                }
            
            self.data[self.system].update({'CrashedReplicates': failed_replicates})
            self.data[self.system].update({'FEP_result': method_results})
            
    # now write me a function that will take the systems within self.data.keys(), try to find `protein` and `water` and for each of the nodes, it will calculate the ddG according to the correct perturbation id (self.data[self.system][fep][FEP_result])
    def calculate_ddG(self):
        self.data.update({'ddG': {}})
        systems = ['1.water', '2.protein']
        # assert both systems have the same FEPs
        prot_feps = sorted([k for k in self.data['2.protein'].keys()])
        water_feps = sorted([k for k in self.data['1.water'].keys()])
        assert prot_feps == water_feps, 'FEPs do not match between protein and water!!'
        for fep in water_feps:
            w_fep = self.data['1.water'][fep]
            p_fep = self.data['2.protein'][fep]
            w_result = w_fep['FEP_result']
            p_result = p_fep['FEP_result']
            
            for _sys in systems:
                # check for inconsistencies
                if w_fep['fep_stage'] != p_fep['fep_stage']:
                    logger.error(f'FEP stages do not match between {fep} and {node}!!')
                    continue
                if w_fep['temperature'] != p_fep['temperature']:
                    logger.error(f'Temperatures do not match between {fep} and {node}!!')
                    continue
                if w_fep['lambda_sum'] != p_fep['lambda_sum']:
                    logger.error(f'Lambda sums do not match between {fep} and {node}!!')
                    continue
            
            ddG = p_result['dG']['avg'] - w_result['dG']['avg']
            ddG_sem = np.sqrt(p_result['dG']['sem']**2 + w_result['dG']['sem']**2)
            self.data['ddG'].update({fep: {'ddG': ddG, 'ddG_sem': ddG_sem}})

class Run(object):
    """
    """
    def __init__(self, FEP, color, PDB, cluster, *args, **kwargs):
        self.cluster=cluster
        self.FEP = FEP.strip('/')
        self.energies = {}
        self.FEPstages = []
        FEPfiles = glob.glob(self.FEP + '/inputfiles/FEP*.fep')
        inputs = glob.glob(self.FEP + '/inputfiles/md*.inp')
        FEPfiles.sort()
        self.failed = []
        for FEPfile in FEPfiles:
            FEPstage = FEPfile.split('/')[-1]
            FEPstage = FEPstage.split('.')[0]
            self.FEPstages.append(FEPstage)
            
        self.lambda_sum = len(FEPfiles) * (len(inputs)-1)
        
        colors = {'blue':['navy','lightblue'],
                  'red' :['darkred','mistyrose']
                 }
        
        self.color = colors[color]
        
    def create_environment(self):
        self.analysisdir = self.FEP + '/analysis'
        # Add overwrite function?
        if os.path.isdir(self.analysisdir) is not True:
            os.mkdir(self.analysisdir)
    
    def read_FEPs(self):
        methods_list = ['dG', 'dGf', 'dGr', 'dGos', 'dGbar']
        methods = {'dG'     : {},
                   'dGf'    : {},
                   'dGr'    : {},
                   'dGos'   : {},
                   'dGbar' :  {}
                  }
        results = {}
        out = []
        
        FEPs = sorted(glob.glob(self.FEP + '/*/*/*/qfep.out'))
        for filename in FEPs:
            i = -1
            file_parse = filename.split('/')
            FEP = file_parse[1]
            temperature = file_parse[2]
            replicate = file_parse[3]
            
            try:
                energies = read_qfep(filename)
            except:
                print("Could not retrieve energies for: " + filename)
                energies = [np.nan, np.nan, np.nan, np.nan, np.nan]
                self.failed.append(replicate)
            #try:
            #    energies = read_qfep(filename)
            #except:
            #    print "Could not retrieve energies for: " + filename
            #    energies = [np.nan, np.nan, np.nan, np.nan, np.nan]

            for key in methods_list:
                i += 1
                try:
                    methods[key][FEP].append(energies[i])
                except:
                    methods[key][FEP] = [energies[i]]
                    
            # Construct for the energy figure
            if not replicate in self.energies:
                self.energies[replicate] = {}
                
            self.energies[replicate][FEP] = read_qfep_verbose(filename)
        for method in methods:
            dG_array = []
            for key in methods[method]:
                print(method, key, methods[method][key])
                dG_array.append(methods[method][key])
            dG_array = np.array(dG_array)
            dG_array = dG_array.astype(np.float)
            dG = avg_sem(dG_array)
            results[method]='{:6.2f}{:6.2f}'.format(*dG)
            
        for method in methods_list:
            out.append(results[method])

        print(self.FEP, '{} {} {} {} {}'.format(*out))
        
    def read_mdlog(self):
        mapping = {}
        cnt = -1
        # Add temperature variable later
        md_files = glob.glob(self.FEP + '/FEP*/*/*/md*.log')        
        md_files.sort()
        md_ref = glob.glob(self.FEP + '/inputfiles/md*.inp')
        windows = len(glob.glob(self.FEP + '/inputfiles/md*.inp')) - 1
        stages = len(glob.glob(self.FEP + '/inputfiles/FEP*.fep'))
        for ref in md_ref:
            w = ref.split('/')[-1].split('_')[2].split('.')[0]
            cnt += 1
            mapping[w]=cnt

        for md_file in md_files:
            stage = md_file.split('/')[1][-1]
            l = md_file.split('/')[-1].split('_')[2].split('.')[0]
            offset = (int(stage) - 1) * int(windows)
            cumulative_l = (mapping[l] + offset)
            (cumulative_l)
            
    def plot_data(self):
        y_axis = {}
        x_axis = range(0,self.lambda_sum+1)
        avg = []
        for replicate in self.failed:
            del self.energies[replicate]
        for replicate in self.energies:
            y_axis[replicate] = [0]
            dG = 0
            for FEPstage in self.FEPstages:
                for energy in self.energies[replicate][FEPstage][0][1:]:
                    energy = dG + energy
                    y_axis[replicate].append(energy)
                dG+=self.energies[replicate][FEPstage][0][-1]
        
        for y in y_axis:
            for i,energy in enumerate(y_axis[y]):
                if len(avg) < self.lambda_sum + 1:
                    avg.append(energy)
                else:
                    avg[i] += energy

            plt.plot(x_axis,y_axis[y],color=self.color[1])
        y_avg = [x / len(y_axis) for x in avg]
        plt.plot(x_axis,y_avg,color=self.color[0])
        axes = plt.gca()
        axes.set_xlim([0,self.lambda_sum])
        plt.xlabel(r'cumulative $\lambda$', fontsize=18)
        plt.ylabel(r'$\Delta$G (kcal/mol)', fontsize=16)        
        plt.savefig(self.analysisdir+'/dG.png',dpi=300,transparent=True)
        
    def write_re2pdb(self):
        curdir = os.getcwd()
        os.chdir(self.FEP + '/analysis')
        if not os.path.exists('pdbs'):
            os.mkdir('pdbs')

        libfiles = glob.glob('../inputfiles/*.lib')
        re_files = glob.glob('../FEP*/*/*/*.re')
        topology = glob.glob('../inputfiles/*.top')[0]
        
        with open('../inputfiles/qprep.inp') as f:
            protlib = f.readline()

        with open('re2pdb.inp', 'w') as outfile:
            outfile.write('{}'.format(protlib))
            
            for libfile in libfiles:
                outfile.write('rl {}\n'.format(libfile))

            outfile.write('rt {}\n'.format(topology))

            for re_file in re_files:
                pdb_out = re_file.split('/')[-1][:-3]
                repeat = '{:02d}'.format(int(re_file.split('/')[3]))
                pdb_out = 'pdbs/{}_{}'.format(repeat, pdb_out)
                outfile.write('rx {}\n'.format(re_file))
                outfile.write('wp {}.pdb\n'.format(pdb_out))
                outfile.write('y\n')
            
            outfile.write('mask none\n')
            outfile.write('mask not excluded\n')
            outfile.write('wp pdbs/complexnotexcluded.pdb\n')
            outfile.write('y\n')
            
            outfile.write('q\n')
            
        qprep = CLUSTER_DICT[self.cluster]['QPREP']
        options = ' < re2pdb.inp > re2pdb.out'
        # Somehow Q is very annoying with this < > input style so had to implement
        # another function that just calls os.system instead of using the preferred
        # subprocess module....
        run_command(qprep, options, string = True)
            
        os.chdir(curdir)
        
            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='protPREP',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '       == Analyse FEP == ')

    
    parser.add_argument('-F', '--FEP',
                        dest = "FEP",
                        required = True,
                        help = "name of FEP directory (FEP_$)")
    
    parser.add_argument('-pdb', '--PDB',
                        dest = "PDB",
                        required = False,
                        default = False,
                        action = 'store_true',
                        help = "Add this argument if you want .pdb files of the trajectory")
    
    parser.add_argument('-c', '--color',
                        dest = "color",
                        required = False,
                        default = 'blue',
                        choices = ['blue', 'red'],
                        help = "color for the plot")
    
    parser.add_argument('-C', '--cluster',
                        dest = "cluster",
                        required = True,
                        help = "cluster information")
    
    
    args = parser.parse_args()
    run = Run(FEP = args.FEP,
              color = args.color,
              cluster = args.cluster,
              PDB = args.PDB
             )
    
    run.create_environment()
    run.read_FEPs()
    run.read_mdlog()
    
    if plot:
        run.plot_data()
        
    if args.PDB:
        run.write_re2pdb()
