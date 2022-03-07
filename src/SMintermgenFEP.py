## currently working with input = sdf file with molecules for which you want to calculate all the intermediates
## output is sdf file with intermediates for each pair?

## functionalities to add
## able to work with molecule pairs that have more than 1 Rgroup
## protonate/deprotonate molecules at pH 7 & remove charged
## add function to permutate that changes tokens into tokens from smaller fragment (currently only removing tokens)
## weld_r_groups gives index & aromaticity errors & sometimes returns SMILES with fragments & also for different rdkit versions get different number of intermediates back
## probably better to work with something else than dataframe for intermediates from each better

## functions that are in this class that might be useful to use from src/shared instead:
## find_largest()
## get_rgroups()
## also reading of sdf file might be done elsewhere

from cmath import nan
from rdkit import rdBase, Chem
rdBase.DisableLog('rdApp.error')
from rdkit.Chem import rdFMCS
from rdkit.Chem import rdRGroupDecomposition as rdRGD
from rdkit.Chem import PandasTools


import pandas as pd
import itertools
import re
import os

import itertools

from collections import defaultdict
from rdkit.Chem.rdchem import EditableMol


class GenInterm():
    def __init__(self, path, namef, wd):
        self.path = path
        self.namef = namef
        self.wd = wd
        
        with open(self.path + self.namef, 'rb') as reader:
            suppl = Chem.ForwardSDMolSupplier(reader)
            # create list of molecules and remove if value is none
            self.mols = [x for x in suppl if x]

        # create list of all molecule combinations
        self.pairs = []
        for subset in itertools.combinations(self.mols, 2):
            self.pairs.append(subset)

        # Catch if path exist?
        if not os.path.isdir(self.wd):
            os.mkdir(self.wd)
        
        intermdir = self.wd + '/' + 'SMinterm'
        if not os.path.isdir(intermdir):
            os.mkdir(intermdir)


    def generate_intermediates(self):
        # Function that actually runs it all
        df_all = pd.DataFrame(columns=['Intermediate', 'Large', 'Small'])
        
        # for each of the pairs create the intermediates
        ## for loop could probably be prettier
        for self.i, self.pair in enumerate(self.pairs):
            self.find_largest()
            self.get_rgroups()
            for self.column in self.columns:
                self.tokenize()
                self.permutate()
            self.weld()
            # save original large & small molecule for keeping track of where the intermediates came from
            ## not sure if better to save as smiles or as mol
            self.df_interm['Large'] = self.df[self.df['Largest'] == True]['Molecule'].values[0]
            self.df_interm['Small'] = self.df[self.df['Largest'] == False]['Molecule'].values[0]
            self.df_interm['PairNum'] = self.i
            self.df_interm['MultipleRgroups'] = self.multiple
        
            df_all = pd.concat([df_all, self.df_interm], ignore_index=True)

        # save sdf with columns; molecules of intermediate, large & smalle fragment 
        ## could be improved
        ## for testing also saving information on pair number (so that easy to reproduce) 
        ## if final intermediate is fragmented & if pairs had multiple r-groups
        PandasTools.WriteSDF(df_all, self.wd + '/' + 'SMinterm/' + self.namef, molColName='Intermediate', properties=['Large', 'Small', 'PairNum', 'MultipleRgroups']) 

    
    def find_largest(self):
        """ 
        Identify the largest molecule of pair
        """
        ## currently selected based on num atoms
        larger = self.pair[0].GetNumAtoms() > self.pair[1].GetNumAtoms()
        d = {'Molecule': self.pair, 'Largest': [larger, not larger]}
        self.df_largest = pd.DataFrame(data=d)

    
    def get_rgroups(self):
        """ 
        Use maximum common substructure of two molecules to get the differing R-Groups
        """
        self.multiple = False
        # find maximim common substructure & pass if none found
        ## possible to use different comparison functions
        res=rdFMCS.FindMCS(self.pair, matchValences=True, ringMatchesRingOnly=True)

        core = Chem.MolFromSmarts(res.smartsString)
        res,_ = rdRGD.RGroupDecompose([core],self.pair,asSmiles=True,asRows=False)
        self.df_rgroup= pd.DataFrame(res)

        if len(self.df_rgroup.columns) > 2:
            print(f'Multiple rgroups for pair number {self.i}')
            self.multiple = True
            self.columns = self.df_rgroup.columns[1:]
        else:
            self.columns = ['R1']

        # combine largest and rgroup info
        self.df = pd.concat([self.df_largest, self.df_rgroup], axis=1)   


    def tokenize(self):
        """ 
        Tokenize SMILES of the small and large R group
        """
        # set regex for SMILES 
        pattern =  r"(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\\\|\\|\/|:|~|@|\?|>|\*|\$|\%[0-9]{2}|[0-9])"
        regex = re.compile(pattern)

        # get index of the largest molecule
        idx_large = self.df.index[self.df['Largest'] == True]
        idx_small = self.df.index[self.df['Largest'] == False]

        rgroup_large = self.df.at[idx_large[0],self.column]
        rgroup_small = self.df.at[idx_small[0],self.column]
        
        # for each r group delete one token and save as intermediate
        self.tokens = [token for token in regex.findall(rgroup_large)]
        self.tokens_small = [token for token in regex.findall(rgroup_small)]
        

    def permutate(self):
        """ 
        Remove tokens or edit tokens from R-groups of the largest molecule. Currently only able to remove tokens
        """
        # create new dataframe for storing adapted rgroups [has Core, R1, ... Rn]
        self.df_interm = pd.DataFrame(columns = self.df.columns[1:]).drop(columns = ['Largest'])
           
        # for large rgroup go over all the combinations with length in range 1 shorter than largest fragment - to 1 larger than shortest fragment 
        ## ask willem if shortest intermediate fragment can have as much atoms as shortest fragment or should always be 1 bigger
        ## maybe handle connection token differently
        for i in range(len(self.tokens) - 1, len(self.tokens_small) - 1, -1):
            for subset in itertools.combinations(self.tokens, i):
                # in some cases connection token will be removed, discard those cases
                ## does not take into account situation with multiple connections in rgroup, like for pair 7 
                ## C1CC1[*:1].[H][*:1].[H][*:1]
                connection = [item for item in subset if re.match('\\[\\*:.\\]', item)]
                if connection:
                    interm = ''.join(subset)
                    # save fragments with valid SMILES
                    if Chem.MolFromSmiles(interm) is not None:
                        self.df_interm.loc[self.df_interm.shape[0], self.column] = interm
        
        self.df_interm = self.df_interm.drop_duplicates(subset=self.column)   

        self.df_interm['Core'] = self.df.at[0,'Core']
        # in case of multiple rgroups also add unchanged rgroups to df
        if self.multiple == True:
            for rgroup in self.columns:
                if rgroup != self.column:
                    self.df_interm[rgroup] = self.df.at[0,rgroup]


    def weld(self):
        """ 
        Put modified rgroups back on the core, returns intermediate
        """
        self.df_interm['Intermediate'] = nan
        # drop doplicate R groups to save time
        self.df_interm = self.df_interm.drop_duplicates(subset='R1')

        for index, row in self.df_interm.iterrows():
            try:
                combined_smiles = row['Core']
                for column in self.columns:
                    combined_smiles = combined_smiles + '.' + row[column]
                mol_to_weld = Chem.MolFromSmiles(combined_smiles)
                welded_mol = self.weld_r_groups(mol_to_weld)
                self.df_interm.at[index, 'Intermediate'] = welded_mol
            except AttributeError:
                pass
            except IndexError:
                pass
            except Chem.rdchem.AtomKekulizeException:
                pass

        ## not 100% sure this works on mol objects
        self.df_interm = self.df_interm.drop_duplicates(subset='Intermediate')
        self.df_interm = self.df_interm.drop(columns = ['R1', 'Core'])
        self.df_interm = self.df_interm.dropna()
        

    def weld_r_groups(self, input_mol):
        """Adapted from 
        https://sourceforge.net/p/rdkit/mailman/rdkit-discuss/thread/CANPfuvAqWxR%2BosdH6TUT-%2B1Fy85fUXh0poRddrEQDxXmguJJ7Q%40mail.gmail.com/"""
        # First pass loop over atoms and find the atoms with an AtomMapNum
        join_dict = defaultdict(list)
        for atom in input_mol.GetAtoms():
            map_num = atom.GetAtomMapNum()
            if map_num > 0:
                join_dict[map_num].append(atom)

        # Second pass, transfer the atom maps to the neighbor atoms
        for idx, atom_list in join_dict.items():
            if len(atom_list) == 2:
                atm_1, atm_2 = atom_list
                nbr_1 = [x.GetOtherAtom(atm_1) for x in atm_1.GetBonds()][0]
                nbr_1.SetAtomMapNum(idx)
                nbr_2 = [x.GetOtherAtom(atm_2) for x in atm_2.GetBonds()][0]
                nbr_2.SetAtomMapNum(idx)

        # Nuke all of the dummy atoms
        new_mol = Chem.DeleteSubstructs(input_mol, Chem.MolFromSmarts('[#0]'))

        # Third pass - arrange the atoms with AtomMapNum, these will be connected
        bond_join_dict = defaultdict(list)
        for atom in new_mol.GetAtoms():
            map_num = atom.GetAtomMapNum()
            if map_num > 0:
                bond_join_dict[map_num].append(atom.GetIdx())

        # Make an editable molecule and add bonds between atoms with correspoing AtomMapNum
        em = EditableMol(new_mol)
        for idx, atom_list in bond_join_dict.items():
            if len(atom_list) == 2:
                start_atm, end_atm = atom_list
                em.AddBond(start_atm, end_atm,
                            order=Chem.rdchem.BondType.SINGLE)

        final_mol = em.GetMol()

        # remove the AtomMapNum values
        for atom in final_mol.GetAtoms():
            atom.SetAtomMapNum(0)
        final_mol = Chem.RemoveHs(final_mol)

        return final_mol