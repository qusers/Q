import argparse
import os
import sys

from rdkit.Chem import PandasTools

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

from SMintermgenFEP import GenInterm

if __name__ == '__main__':

    generator = GenInterm(path='../data/benchmark/qligfep/sdf/', namef='Tyk2_ligands.sdf', wd='Data', insert_small = True)
    generator.generate_intermediates()

    df = PandasTools.LoadSDF('Data/SMinterm/Tyk2_ligands.sdf',smilesName='SMILES',molColName='Intermediate',
           includeFingerprints=False)
    try:   
       print(df['SMILES'])
       df['Fragmented'] = df.apply(lambda row: '.' in row.SMILES, axis=1)
       print(df.Fragmented.value_counts())
       df = df.drop(columns = ['Large', 'Small', 'Intermediate'])
       df.to_csv('Data/SMinterm/test.csv')
    except:
           pass
