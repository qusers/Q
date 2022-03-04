import argparse
import os
import sys

from rdkit.Chem import PandasTools

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

from SMintermgenFEP import GenInterm

if __name__ == '__main__':

    generator = GenInterm(path='../data/benchmark/qligfep/sdf/', namef='Tyk2_ligands.sdf', wd='Data')
    generator.generate_intermediates()

    df = PandasTools.LoadSDF('Data/SMinterm/Tyk2_ligands.sdf',smilesName='SMILES',molColName='Intermediate',
           includeFingerprints=False)
    print(df['SMILES'])
