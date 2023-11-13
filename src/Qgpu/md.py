# Standard Python libraries
import os
import sys
import itertools
from os import path
import json

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../env/')))

# Q-GPU libraries
import IO
import defaults as DEFAULTS

class MD():
    def __init__(self):
        """ This class contains properties for the MD run """
        self.data = { 'stages'              :   None,
                      'steps'               :   None,
                      'stepsize'            :   None,
                      'temperature'         :   None,
                      'random_seed'         :   None,
                      'initial_temperature' :   None,
                      'shake_solvent'       :   None,
                      'shake_hydrogens'     :   None,
                      'shake_solute'        :   None,
                      'lrf'                 :   None,
                      'bath_coupling'       :   None,
                      'solute_solvent'      :   None,
                      'solute_solute'       :   None,
                      'solvent_solvent'     :   None,
                      'q_atom'              :   None,
                      'lrf_cutoff'          :   None,
                      'shell_force'         :   None,
                      'shell_radius'        :   None,
                      'radial_force'        :   None,
                      'polarisation'        :   None,
                      'polarisation_force'  :   None,
                      'output'              :   None,
                      'trajectory'          :   None,
                      'non_bond'            :   None,
                      'energy'              :   None,
                      'topology'            :   None,
                      'trajectory'          :   None,
                      'restart'             :   None,
                      'final'               :   None,
                      'fep'                 :   None,
                      'trajectory_atoms'    :   None,
                      'lambdas'             :   [],
                      'seqrest'             :   [],
                      'posrest'             :   [],
                      'distrest'            :   [],
                      'anglerest'           :   [],
                      'wallrest'            :   [],
                      'thermostat'          :   None,
                    }

class Read_MD(object):
    """
    Read Q topology file as an input, parse to topology class
    """
    def __init__(self, md, *args, **kwargs):    
        self.md = md            
        data = MD()
        self.data = data.data
    
    def JSON(self):
        with open(self.md) as json_file:
            self.data = json.load(json_file)
        
        return(self.data)
    
    def Q(self):
        """ This function reads Q inputfiles from Q-FEP modules"""
        block = 0

        for key in self.data:
            if key in DEFAULTS.MD:
                self.data[key] = DEFAULTS.MD[key]

        # now iterate over all files and populate the MD object
        with open(self.md) as infile:
            for line in infile:
                if len(line) < 2:
                    continue

                if '[MD]' in line:
                    block = 1
                    continue

                if '[cut-offs]' in line:
                    block = 2
                    continue

                if '[sphere]' in line:
                    block = 3
                    continue

                if '[solvent]' in line:
                    block = 4
                    continue

                if '[intervals]' in line:
                    block = 5
                    continue

                if '[files]' in line:
                    block = 6
                    continue

                if '[trajectory_atoms]' in line:
                    block = 7
                    continue

                if '[lambdas]' in line:
                    block = 8
                    continue

                if '[sequence_restraints]' in line:
                    block = 9
                    continue

                if '[atom_restraints]' in line:
                    block = 10
                    continue

                if '[distance_restraints]' in line:
                    block = 11
                    continue                        

                if '[wall_restraints]' in line:
                    block = 12
                    continue                        

                if '[angle_restraints]' in line:
                    block = 13
                    continue                        

                if block == 1:
                    line = line.split()
                    self.data[line[0]] = line[1]

                if block == 2:
                    line = line.split()
                    self.data[line[0]] = line[1]

                if block == 3:
                    line = line.split()
                    self.data[line[0]] = line[1]

                if block == 4:
                    line = line.split()
                    self.data[line[0]] = line[1]

                if block == 5:
                    line = line.split()
                    self.data[line[0]] = line[1]

                if block == 6:
                    continue

                if block == 7:
                    self.data['trajectory_atoms'] = line.strip()

                if block == 8:
                    self.data['lambdas'] = (line.split())

                if block == 9:
                    if len(self.data['seqrest']) == 0:
                        self.data['seqrest'] = [line.split()]

                    else:
                        self.data['seqrest'].append((line.split()))

                if block == 10:
                    if len(self.data['posrest']) == 0:
                        self.data['posrest'] = [line.split()]

                    else:
                        self.data['posrest'].append((line.split()))

                if block == 11:
                    if len(self.data['distrest']) == 0:
                        self.data['distrest'] = [line.split()]

                    else:
                        self.data['distrest'].append((line.split()))

                if block == 12:
                    if len(self.data['wallrest']) == 0:
                        self.data['wallrest'] = [line.split()]

                    else:
                        self.data['wallrest'].append((line.split()))
                
                # Not sure about this keyword
                if block == 13:
                    if len(self.data['anglerest']) == 0:
                        self.data['anglerest'] = [line.split()]

                    else:
                        self.data['anglerest'].append((line.split()))

        return(self.data)
    
class Write_MD(object):        
    """
    Write Python topology object to file
    """
    def __init__(self, data, *args, **kwargs):
        self.data = data
    
    def CSV(self,wd):
        """
        .csv file format for qdyn GPU version
        
        """
        self.wd = wd
        lists = ['lambdas','seqrest','distrest']
        
        j = 31
        
        # Get length of lists (lambdas, restraints)
        for l in lists:
            if self.data[l] != None:
                j += len(self.data[l])
            
        
        # Add lines for restraints
        # TO DO, restraints not implemented yet
        
        with open(self.wd + 'md.csv', 'w') as outfile:
            outfile.write('{}\n'.format(j))
            outfile.write('steps;{}\n'              .format(self.data['steps']))        
            outfile.write('stepsize;{}\n'           .format(self.data['stepsize']))        
            outfile.write('temperature;{}\n'        .format(self.data['temperature']))        
            outfile.write('thermostat;{}\n'         .format(self.data['thermostat']))        
            outfile.write('bath_coupling;{}\n'      .format(self.data['bath_coupling']))        
            outfile.write('random_seed;{}\n'        .format(self.data['random_seed']))        
            outfile.write('initial_temperature;{}\n'.format(self.data['initial_temperature']))        
            outfile.write('shake_solvent;{}\n'      .format(self.data['shake_solvent']))       
            outfile.write('shake_solute;{}\n'       .format(self.data['shake_solute'])) 
            outfile.write('shake_hydrogens;{}\n'    .format(self.data['shake_hydrogens']))        
            outfile.write('lrf;{}\n'                .format(self.data['lrf']))   
            outfile.write('charge_groups;{}\n'      .format('on')) # charge group, always on        
            outfile.write('solute_solute;{}\n'      .format(self.data['solute_solute']))        
            outfile.write('solvent_solvent;{}\n'    .format(self.data['solvent_solvent']))        
            outfile.write('solute_solvent;{}\n'     .format(self.data['solute_solvent']))        
            outfile.write('q_atom;{}\n'             .format(self.data['q_atom']))        
            outfile.write('shell_radius;{}\n'       .format(self.data['shell_radius']))        
            outfile.write('shell_force;{}\n'        .format(self.data['shell_force']))        
            outfile.write('radial_force;{}\n'       .format(self.data['radial_force']))        
            outfile.write('polarisation;{}\n'       .format(self.data['polarisation']))        
            outfile.write('polarisation_force;{}\n' .format(self.data['polarisation_force']))        
            outfile.write('non_bond;{}\n'           .format(self.data['non_bond']))        
            outfile.write('output;{}\n'             .format(self.data['output']))        
            outfile.write('energy;{}\n'             .format(self.data['energy']))        
            outfile.write('trajectory;{}\n'         .format(self.data['trajectory']))

            # possible to have multiple lambdas, need to loop
            ### LAMBDAS ###
            if self.data['lambdas'] != None:
                outfile.write('{};lambdas\n'.format(len(self.data['lambdas'])))
                for l in self.data['lambdas']:
                    outfile.write('{}\n'.format(l))
                    
            else:
                outfile.write('0;lambdas\n')

            ### SEQREST ###
            if self.data['seqrest'] != None:
                outfile.write('{};{}\n'.format(len(self.data['seqrest']),'seqrest'))
                for r in self.data['seqrest']:
                    outfile.write('{};{};{};{};{}\n'.format(*r))

            else:
                outfile.write('0;{}\n'.format('seqrest'))    
                
            ### POSREST ###
            if self.data['posrest'] != None:
                outfile.write('{};{}\n'.format(len(self.data['posrest']),'posrest'))
                for r in self.data['posrest']:
                    outfile.write('{};{};{};{};{};{}\n'.format(*r))

            else:
                outfile.write('0;{}\n'.format('posrest'))    
                
            ### DISTREST ###
            if self.data['distrest'] != None:
                outfile.write('{};{}\n'.format(len(self.data['distrest']),'distrest'))
                for r in self.data['distrest']:
                    outfile.write('{};{};{};{};{};{}\n'.format(*r))

            else:
                outfile.write('0;{}\n'.format('distrest'))    
                
            ### ANGLEREST ###
            if self.data['anglerest'] != None:
                outfile.write('{};{}\n'.format(len(self.data['anglerest']),'anglerest'))
                for r in self.data['anglerest']:
                    outfile.write('{};{};{};{};{};{}\n'.format(*r))

            else:
                outfile.write('0;{}\n'.format('anglerest'))    
                
            ### WALLREST ###
            if self.data['wallrest'] != None:
                outfile.write('{};{}\n'.format(len(self.data['wallrest']),'wallrest'))
                for r in self.data['wallrest']:
                    outfile.write('{};{};{};{};{};{}\n'.format(*r))

            else:
                outfile.write('0;{}\n'.format('wallrest'))    

    def JSON(self,out_json):
        """
        .json MD input file
        """
        with open(out_json, 'w') as outfile:
            inputs = self.data
            json.dump(inputs,outfile,indent=2)    
