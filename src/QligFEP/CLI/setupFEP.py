"""Module to generate (setup) all FEP files for the directory you're working on."""

from pathlib import Path
import argparse
import json
import subprocess
import logging
import shutil

logging.basicConfig(
    filename='cli_calls.log',
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s')

def ligpairs_from_json(json_file):
    with open(json_file, 'r') as infile:
        json_dict = json.load(infile)
    try:
        edges = json_dict['edges'] # should be a list of dictionaries
    except KeyError:
        raise KeyError('Could not find "edges" in json file')
    ligpairs = [(e['from'], e['to']) for e in edges]
    return ligpairs

def create_call(**kwargs):
    """Function to dynamically create a call to QligFEP.cli.main_exe() based on the kwargs."""    
    template = (
        'qligfep -l1 {lig1} -l2 {lig2} -FF {FF} -s {system} -c {cluster} -R {replicates} '
        '-S {sampling} -r {sphereradius} -l {start} -w {windows} -T {temperature}'
        )
    if 'cysbond' in kwargs and kwargs['cysbond'] is not None:
        template += ' -b {cysbond}'
    return template.format(**kwargs)

def submit_command(command :str) -> None:
    """Function to submit a command using the subprocess module.

    Args:
        command: string with the command to be submitted.
    """
    try:
        result = subprocess.run(command, shell=True, text=True, capture_output=True, check=True)
        logging.info(f"Command executed successfully: {command}\nOutput: {result.stdout}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed: {command}\nError: {e.stderr}")
        raise
    except Exception as e:
        logging.error(f"Unexpected error: {e}")
        raise

def parse_arguments():
    parser = argparse.ArgumentParser(
        prog='QligFEP',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description =
        "Generate all FEP files for the directory you're working on.")
    parser.add_argument('-FF', '--forcefield',
                        dest = "FF",
                        required = True,
                        choices = ['OPLS2005', 'OPLS2015', 'AMBER14sb', 'CHARMM36', 'CHARMM22', 'CHARMM_TEST'],
                        help = "Forcefield to be used")

    parser.add_argument('-c', '--cluster',
                        dest = "cluster",
                        required = True,
                        help = "cluster you want to submit to, cluster specific parameters added to settings"
                       )

    parser.add_argument('-r', '--sphereradius',
                        dest = "sphereradius",
                        required = False,
                        default = '25',
                        help = "size of the simulation sphere"
                       )

    parser.add_argument('-b', '--cysbond',
                        dest = "cysbond",
                        default = None,
                        help = "Temporary function to add cysbonds at1:at2,at3:at4 etc."
                       )

    parser.add_argument('-l', '--start',
                        dest = "start",
                        default = '0.5',
                        choices = ['1', '0.5'],
                        help = "Starting FEP in the middle or endpoint"
                       )

    parser.add_argument('-T', '--temperature',
                        dest = "temperature",
                        default = '298',
                        help = "Temperature(s), mutliple tempereratures given as 'T1,T2,...,TN'"
                       )

    parser.add_argument('-R', '--replicates',
                        dest = "replicates",
                        default = '10',
                        help = "How many repeats should be run"
                       )

    parser.add_argument('-S', '--sampling',
                        dest = "sampling",
                        default = 'sigmoidal',
                        choices = ['linear', 'sigmoidal', 'exponential', 'reverse_exponential'],
                        help = "Lambda spacing type to be used"
                       )
    parser.add_argument('-w', '--windows',
                        dest = "windows",
                        default = '100',
                        help = "Total number of windows that will be run"
                       )
    parser.add_argument('-j', '--json_map',
                        dest='json_map',
                        help = (
                            "Path for the '.json' QmapFEP file. If not given, the script will "
                            "look for a single '.json' file in the current directory and raise "
                            "an error if there are more than one."
                        ),
                        default = None,
                       )
    return parser.parse_args()

def main_exe():
    args = parse_arguments()
    
    cwd = Path.cwd()
    systems = ['water', 'protein']
    sys_directories = [cwd / '1.water', cwd / '1.protein']
    # make sure that the default directories for running the FEP calculations are there
    for sys_dir in sys_directories:
        if not sys_dir.exists():
            sys_dir.mkdir()
    
    if args.json_map is None: # Try to load a json file from cwd
        json_files = list(cwd.glob('*.json'))
        if len(json_files) == 1:
            args.json_map = json_files[0]
        else:
            raise FileNotFoundError('No QmapFEP json file found in the current directory')

    lig_pairs = ligpairs_from_json(args.json_map)
    for system, sys_dir in zip(systems, sys_directories):
        for pair in lig_pairs:
            lig1 = pair[0]
            lig2 = pair[1]
            
            temp_dir = cwd  / f'FEP_{lig1}_{lig2}'
            command = create_call(
                lig1 = lig1,
                lig2 = lig2,
                FF = args.FF,
                system=system,
                cluster = args.cluster,
                sphereradius = args.sphereradius,
                cysbond = args.cysbond,
                start = args.start,
                temperature=args.temperature,
                replicates = args.replicates,
                sampling = args.sampling,
                windows = args.windows
            )
            
            dst = sys_dir / ('FEP_' + lig1 + '_' + lig2)
            submit_command(command)
            shutil.move(temp_dir, dst)

if __name__ == '__main__':
    main_exe()
