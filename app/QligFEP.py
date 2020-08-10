import argparse
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src/')))

import QligFEP

class Startup(object):
    """
    Create dual topology FEP files based on two ligands
    """
    def __init__(self, 
                 lig1, 
                 lig2, 
                 FF, 
                 system, 
                 cluster, 
                 sphereradius, 
                 cysbond, 
                 start, 
                 temperature, 
                 replicates,
                 sampling,
                 *args, 
                 **kwargs):
        
        data = {'self'          : None,
                'lig1'          : None,
                'lig2'          : None,
                'FF'            : None,
                'system'        : None,
                'cluster'       : None,
                'sphereradius'  : None,
                #'cysbond' : None,      NEEDS TO GO
                'start'         : None,
                'temperature'   : None,
                'replicates'    : None,
                'sampling'      : None,
               }
        
        self.lig1 = lig1
        self.lig2 = lig2
        self.FF = FF
        self.system = system
        self.rootdir = os.getcwd()
        self.cluster = cluster
        self.sphereradius = sphereradius
        self.cysbond = cysbond
        self.start = start
        self.include = ['ATOM', 'HETATM']
        self.temperature = temperature
        self.replicates = replicates
        self.sampling = sampling

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='QligFEP',
        version='1.0',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '       == Generate FEP files for dual topology ligand FEP == ')

    
    parser.add_argument('-l1', '--lig_1',
                        dest = "lig1",
                        required = True,
                        help = "name of ligand 1")
    
    parser.add_argument('-l2', '--lig_2',
                        dest = "lig2",
                        required = True,
                        help = "name of ligand 2")
    
    parser.add_argument('-FF', '--forcefield',
                        dest = "FF",
                        required = True,
                        choices = ['OPLS2005', 'OPLS2015', 'AMBER14sb', 'CHARMM36', 'CHARMM22', 'CHARMM_TEST'],
                        help = "Forcefield to be used")
    
    parser.add_argument('-s', '--system',
                        dest = "system",
                        required = True,
                        choices = ['water', 'protein', 'vacuum'],
                        help = "what type of system we are setting up")
    
    parser.add_argument('-c', '--cluster',
                        dest = "cluster",
                        required = True,
                        help = "cluster you want to submit to, cluster specific parameters added to settings"
                       )
    
    parser.add_argument('-r', '--sphereradius',
                        dest = "sphereradius",
                        required = False,
                        default = '15',
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
                        default = 'linear',
                        choices = ['linear', 'sigmoidal', 'exponential', 'reverse_exponential'],
                        help = "Lambda spacing type to be used"
                       )
    
    parser.add_argument('-w', '--windows',
                        dest = "windows",
                        default = '50',
                        help = "Total number of windows that will be run"
                       )
    
    args = parser.parse_args()
    Startup(lig1 = args.lig1,
            lig2 = args.lig2,
            FF= args.FF,
            system = args.system,
            cluster = args.cluster,
            sphereradius = args.sphereradius,
            cysbond = args.cysbond,
            start = args.start,
            temperature = args.temperature,
            replicates = args.replicates,
            sampling = args.sampling
           )