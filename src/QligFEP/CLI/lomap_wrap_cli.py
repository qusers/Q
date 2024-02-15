"""Module to wrap the lomap package for the QligFEP CLI."""

import json
import lomap
import argparse
from pathlib import Path
from loguru import logger
from typing import Optional
from multiprocessing import cpu_count

# QligFEP imports
from ..chemIO import MoleculeIO

class LomapWrap(object):
    """Class to wrap the lomap package for the QligFEP CLI."""
    def __init__(self, inp: str, out: Optional[str]=None, time=30, verbose='info', **kwargs):
        self.inp = inp
        self.out = self._parse_output(out)
        self.cores = self._setup_cores()
        self.lomap_args = {
            'directory': None,
            'parallel': self.cores,
            'time': time,
            'verbose':verbose,
            **kwargs
        }
        self._check_input()
        
    def _check_input(self):
        """Method to check self.inp for the correct file format."""
        inpath = Path(self.inp)
        if inpath.suffix == '.sdf':
            if any([len(list(inpath.glob('*.sdf'))) < 2, len(list(inpath.glob('*.mol2'))) > 2]):
                logger.warning('You are using a directory as input, certify that you have the desired ligand files.')
            self.lomap_args.update({'directory': self.inp})
        elif inpath.is_dir():
            handler = MoleculeIO(self.inp)
            Path(self.out).mkdir(parents=True, exist_ok=False)
            logger.info(f'Writing {self.inp} to separate `.sdf` files to be stored in {self.out}.')
            handler.write_sdf_separate(self.out)
            self.lomap_args.update({'directory': self.out})
    
    def _parse_output(self, output):
        if output is None:
            output = Path(self.inp).parent / Path(self.inp).stem
        else:
            output = Path(output).absolute()
        return output
            
    def _setup_cores(self):
        """Method to check the number of cores."""
        if cpu_count() < 2:
            logger.warning('You are using only one core. This might take a while.')
            cores = 1
        elif cpu_count() < 8:
            cores = 4
        else:
            cores = 8
        return cores
    
    def format_graph_data(self, nx_graph):
        organized_data = []

        # Iterate over all edges to gather necessary information
        for source, target, data in nx_graph.edges(data=True):
            # Extract edge data
            weight = data.get("similarity", 0)  # Default to 0 if not found
            strict_flag = data.get("strict_flag", None)  # Include if needed

            # Extract node names from nodes data
            source_node_data = nx_graph.nodes[source]
            target_node_data = nx_graph.nodes[target]
            source_name = source_node_data.get("fname_comp").replace(
                ".sdf", ""
            )  # Assuming fname_comp holds the file name
            target_name = target_node_data.get("fname_comp").replace(".sdf", "")

            # Construct the dictionary for the current edge
            edge_dict = {
                "weight": weight,
                "source": source,
                "target": target,
                "strict_flag": strict_flag,
                "from": source_name,
                "to": target_name,
            }
            # Append the constructed dictionary to the list
            organized_data.append(edge_dict)

        return {'edges': organized_data}
    
    def run_lomap(self):
        db_mol = lomap.DBMolecules("haha", output=True)
        # Calculate the similarity matrices
        strict, loose = db_mol.build_matrices()
        # Generate the NetworkX graph and output the results
        nx_graph = db_mol.build_graph()
        result_dict = self.format_graph_data(nx_graph)
        with open(self.out / 'lomap.json', 'w') as f:
            json.dump(result_dict, f, indent=4)
        
def parse_arguments() -> argparse.Namespace:
    """Method to parse the arguments."""
    parser = argparse.ArgumentParser(description='Wrap the lomap package for the QligFEP CLI.')
    parser.add_argument('--input', '-i', type=str, help='Input file or directory.')
    parser.add_argument('--output', '-o', type=str, help='Output directory.')
    parser.add_argument('--time', '-t', type=int, default=30, help='Time to run the lomap package.')
    parser.add_argument('--verbose', '-v', type=str, default='info', help='Verbosity level.')
    return parser.parse_args()
    
def main(args):
    # TODO: could implement other lomap arguemnts here
    lomap = LomapWrap(args.input, args.output, args.time, args.verbose)
    lomap.run_lomap()

def main_exe():
    """Main method to execute the lomap package."""
    args = parse_arguments()
    main(args)

if __name__ == '__main__':
    main_exe()