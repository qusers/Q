import argparse
import os
import re
from pathlib import Path
import matplotlib.pyplot as plt

import numpy as np

from .IO import read_qfep, read_qfep_verbose, run_command
from .logger import logger, setup_logger
from .settings.settings import Q_PATHS


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

class FepReader(object):
    """Class to analyze FEP output files. This wrapper class will read, structure the files
    and enable input/output operation with the analyzed data.
    """    
    def __init__(self, system:str, target_name:str) -> None:
        """Initialize the FEP reader class. This class will store the FEP information inside
        `self.data` and will be used to analyze the results & generate plots.

        Args:
            system: which system to be loaded first. This should be a directory containing
                the FEP directories, named by default with the format `FEP_*`.
            target_name: name of the target protein so we can load the correct FEP directories.
        """        
        self.cwd = Path.cwd()
        self.data = {}
        self.system = None
        self.target_name = target_name
        self.load_system(system)
        self.read_fep_inputs()
        
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
        """Load basic FEP information from the input files, including temperature,
        number of replicates, and lambda sum. Information is stored in`self.data`."""
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
        """Read the ran perturbations. Running this method will populate the `self.data` dictionary
        with the FEPs and their respective delta-G's for the loaded system (self.system)."""
        methods_list = ['dG', 'dGf', 'dGr', 'dGos', 'dGbar']
        feps = [k for k in self.data[self.system].keys()]
        
        for fep in feps:
            logger.debug(f'Reading FEP: {fep}')
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
                if rep.stat().st_size == 0: # if the file is empty, try runnign qfep again
                    logger.warning(f'Empty qfep.out file: {rep}. Trying to run qfep again...')
                    qfep = Q_PATHS['QFEP']
                    os.chdir(str(rep.parent.absolute()))
                    options = ' < qfep.inp > qfep.out'
                    run_command(qfep, options, string = True)
                    os.chdir(str(self.cwd.absolute()))
                logger.debug(f'    Reading qfep.out file: {rep}')
                repID = int(rep.parent.name)
                try:
                    # TODO: shall we also support the verbose output? -> see IO.read_qfep_verbose
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
            
    def calculate_ddG(self, water_sys:str = '1.water', protein_sys:str = '2.protein'):
        """After running `read_perturbations`, for both the water and the protein systems,
        running this method will calculate the ddG for each FEP and store it in `self.data`
        under the key `result`.

        Args:
            water_sys: name of the water system that was read. Defaults to '1.water'.
            protein_sys: name of the protein system that was read. Defaults to '2.protein'.
        """        
        self.data.update({'result': {}})
        systems = [water_sys, protein_sys]
        # assert both systems have the same FEPs
        prot_feps = sorted([k for k in self.data[protein_sys].keys()])
        water_feps = sorted([k for k in self.data[water_sys].keys()])
        assert prot_feps == water_feps, 'FEPs do not match between protein and water!!'
        for fep in water_feps:
            w_fep = self.data[water_sys][fep]
            p_fep = self.data[protein_sys][fep]
            w_result = w_fep['FEP_result']
            p_result = p_fep['FEP_result']
            
            for _sys in systems:
                # check for inconsistencies
                if w_fep['fep_stage'] != p_fep['fep_stage']:
                    logger.error(f'FEP stages do not match between water/{fep} and protein/{fep}.')
                    continue
                if w_fep['temperature'] != p_fep['temperature']:
                    logger.error(f'Temperatures do not match between water/{fep} and protein/{fep}.')
                    continue
                if w_fep['lambda_sum'] != p_fep['lambda_sum']:
                    logger.error(f'Lambda sums do not match between water/{fep} and protein/{fep}.')
                    continue
            
            for method in w_result.keys():
                new_key = f'g{method}'
                ddG = p_result[method]['avg'] - w_result[method]['avg']
                ddG_sem = np.sqrt(p_result[method]['sem']**2 + w_result[method]['sem']**2)
                self.data['result'][new_key].update({fep: {new_key: ddG, f'{new_key}_sem': ddG_sem}})
                
    def create_ddG_plot(
            self,
            method,
            margin: float = 1.0,
            xylims: tuple | None = None,
            output_path: str | None = None,
        ):
        """Creates the ddG plot for the FEP that has already been analyzed. The plot will
        show the experimental (X axis) vs mean predicted values (Y axis), with error bars
        representing the standard error of the mean (SEM).

        Args:
            method: the energy method to be used for the plot. Must be one of the keys in
                the result dictionary.
            margin: margin value to be added/subtracted to the max/min values obtained. Defaults to 1.0.
            xylims: if values are passed, x&y min will be xylims[0] and max will be [1]. Defaults to None.

        Returns:
            the matplotlib figure and axis objects (fig, ax).
        """        
        if 'result' not in self.data.keys():
            logger.error('No results to plot. Run calculate_ddG first.')
            return
        elif method not in self.data['result'].keys():
            logger.error(
                f"{method} not in result dictionary. Pick one of the following: "
                f"{', '.join(self.data['result'].keys())}"
                )
            return
        data_dict = self.data['result'][method]
        avg_values = [data['avg'] for data in data_dict.values()]
        sem_values = [data['sem'] for data in data_dict.values()]
        exp_values = [data['exp_data'] for data in data_dict.values()]
        
        # Calculate RMSE & correlation coefficient
        rmse = np.sqrt(np.mean((np.array(avg_values) - np.array(exp_values))**2))
        correlation_coef = np.corrcoef(exp_values, avg_values)[0, 1]
        
        if xylims is not None:
            assert len(xylims) == 2, 'xylims must be a tuple with 2 elements.'
            assert xylims[0] < xylims[1], 'xylims[0] must be smaller than xylims[1].'
            min_val = xylims[0]
            max_val = xylims[1]
        else:
            all_values = avg_values + exp_values
            min_val = min(all_values) - margin
            max_val = max(all_values) + margin
        
        fig, ax = plt.subplots()
        
        plt.errorbar(exp_values, avg_values, yerr=sem_values, fmt='o', color='black', ecolor='lightgray', elinewidth=3, capsize=0)
        plt.plot([min_val, max_val], [min_val, max_val], 'k-', linewidth=2) # Black identity line
        
        # Highlight predictions within 1 and 2 kcal/mol of the experimental affinity
        ax.fill_between([min_val, max_val], [min_val - 1, max_val - 1], [min_val + 1, max_val + 1], color='darkgray', alpha=0.5)
        ax.fill_between([min_val, max_val], [min_val - 2, max_val - 2], [min_val + 2, max_val + 2], color='lightgray', alpha=0.5)
        
        # Annotating the plot with τ and RMSE # TODO: figure out how to place it...
        plt.text(min_val, max_val, f'$\\tau = {correlation_coef:.2f}$, RMSE = {rmse:.2f} kcal/mol', fontsize=12, verticalalignment='top')
        
        # set labels, make it square and add legend
        plt.xlabel('$\Delta\Delta G_{exp} [kcal/mol]$')
        plt.ylabel('$\Delta\Delta G_{pred} [kcal/mol]$')
        plt.xlim(min_val, max_val)
        plt.ylim(min_val, max_val)
        ax.set_aspect('equal', adjustable='box')
        ax.legend(['Identity line', 'Within 1 kcal/mol', 'Within 2 kcal/mol'], loc='upper left')
        if output_path is not None:
            if isinstance(output_path, str):
                output_path = Path(output_path)
            else: 
                assert isinstance(output_path, Path), 'output_path must be a string or a Path object.'
            if output_path.isdir():
                logger.info(f'Using default name to save the plot file under the dir: {output_path}')
                fig.savefig(f'{self.target_name}_{method}_ddG_plot.png', dpi=300, bbox_inches='tight')
            elif output_path.exits():
                logger.warning(f'File {output_path} already exists. Overwriting...')
            else:
                logger.info(f'Saving the plot to {output_path}')
            assert output_path.stem == '.png', 'output_path must be a .png file.'
            fig.savefig(str(output_path), dpi=300, bbox_inches='tight')
        return fig, ax


def parse_arguments() -> argparse.Namespace:
    """Method to parse the arguments."""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = 'Analyze FEP output files and generate plots.')
    parser.add_argument('-p', '--protein-dir',
                        dest = "protein_dir",
                        required = False,
                        default = '2.protein',
                        help = (
                            "Path to the directory containing the protein system FEPs. "
                            "Will default to `2.protein` in the current working directory."
                            ))

    parser.add_argument('-w', '--water-dir',
                        dest = "water_dir",
                        required = False,
                        default = '2.protein',
                        help = (
                            "Path to the directory containing the water system FEPs. "
                            "Will default to `1.water` in the current working directory."
                            ))
    
    parser.add_argument('-j', '--json-file',
                        dest = "json_file",
                        required = True,
                        help = (
                            "Path to the .json file containing the mapping of the perturbations. "
                            "This should be the same file used to run the setupFEP script."
                        ))
    
    parser.add_argument('-t', '--target',
                        dest = "target",
                        required = True,
                        help = "Name of the protein target; used to save the plot."
                        )
    
    parser.add_argument('-m', '--method',
                        required = False,
                        default = 'ddG',
                        choices = ['ddG', 'ddGf', 'ddGr', 'ddGos', 'ddGbar'],
                        help = 'Energy method to be used for the plot. Defaults to ddG.'
                        )

    parser.add_argument('-l', '--log-level',
                        required = False,
                        default = 'INFO',
                        help = "Set the log level for the logger. Defaults to INFO.",
                        choices = ["TRACE", "DEBUG", "INFO", "SUCCESS", "WARNING", "ERROR", "CRITICAL"],
                        )
                        
    return parser.parse_args()

def main(args):
    setup_logger(level=args.log_level)
    fep_reader = FepReader(system = args.water_dir, target_name = args.target)
    fep_reader.read_perturbations()
    fep_reader.load_new_system(system = args.protein_dir)
    fep_reader.read_perturbations()
    fep_reader.calculate_ddG()
    fep_reader.create_ddG_plot(method = args.method)
    # TODO: add method to populate json file
    

def main_exe():
    args = parse_arguments()
    main(args)

# def write_re2pdb(self): # TODO: port this to the new class
#     curdir = os.getcwd()
#     os.chdir(self.FEP + '/analysis')
#     if not os.path.exists('pdbs'):
#         os.mkdir('pdbs')

#     libfiles = glob.glob('../inputfiles/*.lib')
#     re_files = glob.glob('../FEP*/*/*/*.re')
#     topology = glob.glob('../inputfiles/*.top')[0]
    
#     with open('../inputfiles/qprep.inp') as f:
#         protlib = f.readline()

#     with open('re2pdb.inp', 'w') as outfile:
#         outfile.write('{}'.format(protlib))
        
#         for libfile in libfiles:
#             outfile.write('rl {}\n'.format(libfile))

#         outfile.write('rt {}\n'.format(topology))

#         for re_file in re_files:
#             pdb_out = re_file.split('/')[-1][:-3]
#             repeat = '{:02d}'.format(int(re_file.split('/')[3]))
#             pdb_out = 'pdbs/{}_{}'.format(repeat, pdb_out)
#             outfile.write('rx {}\n'.format(re_file))
#             outfile.write('wp {}.pdb\n'.format(pdb_out))
#             outfile.write('y\n')
        
#         outfile.write('mask none\n')
#         outfile.write('mask not excluded\n')
#         outfile.write('wp pdbs/complexnotexcluded.pdb\n')
#         outfile.write('y\n')
        
#         outfile.write('q\n')
        
#     qprep = CLUSTER_DICT[self.cluster]['QPREP']
#     options = ' < re2pdb.inp > re2pdb.out'
#     run_command(qprep, options, string = True)
        
#     os.chdir(curdir)
        
            
if __name__ == "__main__":
    main_exe()