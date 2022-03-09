import argparse
import io
import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
from networkx.drawing.nx_agraph import graphviz_layout
from networkx.readwrite import json_graph
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import PandasTools
from rdkit.Chem import AllChem
from rdkit.Chem import rdDepictor
from rdkit import Chem
from rdkit.Chem import rdFMCS
import sys
import os
import json
import shutil

# import Q modules
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src/share/')))

from share import settings as s
from share import Rgroup
from share import plot

class GenRgroup():
    def __init__(self, IO_sdf, sdf, wd, env=True):
        self.sdf = sdf
        self.wd = wd
        #self.suppl = Chem.ForwardSDMolSupplier(in_sdf)
        #self.suppl2 = Chem.SDMolSupplier(in_sdf)
        self.frame = PandasTools.LoadSDF(self.sdf,smilesName='SMILES',molColName='Molecule',
           includeFingerprints=True)

        rdDepictor.SetPreferCoordGen(True)

        if env == True:
            # Catch if path exist?
            if not os.path.isdir(self.wd):
                os.mkdir(self.wd)
        
            imgdir = self.wd + '/' + 'img'
            if not os.path.isdir(imgdir):
                os.mkdir(imgdir)

    def generate_images(self):
        # Function that actually runs it all
        self.compute_coords()
        self.get_core()
        self.generate_Rgroups()
        self.draw()

    def compute_coords(self):
        
        smis = self.frame['SMILES']
        self.cids = list(self.frame.ID)
        self.ms = [Chem.MolFromSmiles(x) for x in smis]
        for m in self.ms:
            rdDepictor.Compute2DCoords(m)

    def get_core(self):
        # Change to ForwardSDMolSupplier at some point
        with Chem.SDMolSupplier(self.sdf) as cdk2mols:
            res=rdFMCS.FindMCS(cdk2mols).smartsString
        #res=rdFMCS.FindMCS(self.sdf).smartsString
        self.core = Chem.MolFromSmarts(res)
        rdDepictor.Compute2DCoords(self.core)

    def generate_Rgroups(self):
        ps = Chem.AdjustQueryParameters.NoAdjustments()
        ps.makeDummiesQueries=True
        self.qcore = Chem.AdjustQueryProperties(self.core,ps)
        mhs = [Chem.AddHs(x,addCoords=True) for x in self.ms]
        self.mms = [x for x in mhs if x.HasSubstructMatch(self.qcore)]
        for m in self.mms:
            for atom in m.GetAtoms():
                atom.SetIntProp("SourceAtomIdx",atom.GetIdx())
        
        self.groups,_ = Chem.rdRGroupDecomposition.RGroupDecompose([self.qcore],self.mms,asSmiles=False,asRows=True)

    def draw(self):
        # Make the label definition non-dependent on user input
        lbls=[*self.groups[0]]
        lbls.remove("Core")
        overview, images = Rgroup.draw_multiple(self.mms,
                                                self.groups,
                                                self.qcore,
                                                lbls=lbls,
                                                #legends=self.cids,
                                                subImageSize=(300,250))

        # save file with unique identifier (?)
        for i, image in enumerate(images):
            image.save(self.wd + '/img/' + self.cids[i] + '.png', format="png")

class MapGen():
    def __init__(self, in_sdf, metric, o, wd):
        self.suppl = Chem.ForwardSDMolSupplier(in_sdf)
        self.metric = metric
        self.otxt = o
        self.lig_dict = {}
        self.simF = None    # The similarity function to calculate distance between ligands
        self.wd = wd
        self._set_similarity_function()

    def make_fp(self, mol):
        if self.metric == 'Tanimoto':
            fp = FingerprintMols.FingerprintMol(mol)
        elif self.metric == 'MFP':
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        elif self.metric == 'SMILES':
            fp = Chem.MolToSmiles(mol, isomericSmiles=True)
        else:
            fp = None
        return fp

    def _set_similarity_function(self):
        """Set the similarity function to be used with selected metric."""
        if self.metric in ['Tanimoto', 'MFP']:
            from rdkit.DataStructs import FingerprintSimilarity as f
        elif self.metric == 'SMILES':
            from Bio.pairwise2 import align
            f = align.globalms
        self.simF = f

    def _ligands_score(self, data, lig_i, lig_j):
        """Return a similarity score between lig_i and lig_j using self.simF."""
        if lig_i == lig_j:
            return 100.0 if self.metric in ['SMILES', 'MCS'] else 1.0
        if self.metric in ['Tanimoto', 'MFP']:
            return self.simF(data['FP'][lig_i], data['FP'][lig_j])
        if self.metric == 'MCS':
            score = self.simF([data['Mol'][lig_i],
                               data['Mol'][lig_j]],
                              atomCompare=rdFMCS.AtomCompare.CompareAny)
            return score.numAtoms + score.numBonds
        if self.metric == 'SMILES':
            return (100 - self.simF(
                data['FP'][lig_i], data['FP'][lig_j], 1, -1, -0.5, -0.05)[0][2]) / 100

    def set_ligdict(self):
        for mol in self.suppl:
            charge = Chem.rdmolops.GetFormalCharge(mol)
            v = self.lig_dict.setdefault(charge, {'Name': [], 'Mol': [], 'FP': [], 'AtNumb':[]})
            v['Name'].append(mol.GetProp('_Name'))
            v['AtNumb'].append(mol.GetNumHeavyAtoms())
            v['Mol'].append(mol)
            if self.metric != 'MCS':
                v['FP'].append(self.make_fp(mol))

    def sim_mx(self):
        # Ensure the self.lig_dict has been created
        if not self.lig_dict:
            self.set_ligdict()
        for charge, ligand in self.lig_dict.items():
            seed = Chem.rdFMCS.FindMCS([mol for mol in self.lig_dict[charge]['Mol']], atomCompare=rdFMCS.AtomCompare.CompareAny, bondCompare=rdFMCS.BondCompare.CompareAny).smartsString
            sim_df = pd.DataFrame()
            dif_df = pd.DataFrame()
            ligands_done = []
            for i, lig1 in enumerate(ligand['Name']):
                for j, lig2 in enumerate(ligand['Name']):
                    if lig2 not in ligands_done:
                        pair = [ligand['Mol'][i], ligand['Mol'][j]]
                        mcs = Chem.rdFMCS.FindMCS(pair, atomCompare=rdFMCS.AtomCompare.CompareAny, bondCompare=rdFMCS.BondCompare.CompareAny, seedSmarts=seed)
                        d1, d2 = abs(mcs.numAtoms - ligand['AtNumb'][i]), abs(mcs.numAtoms - ligand['AtNumb'][j])
                        sim_df.loc[lig2, lig1] = self._ligands_score(ligand, i, j) - 0.001*(d1+d2)
                        dif_df.loc[lig2, lig1] = d1 + d2
                            
                ligands_done.append(lig1)
            self.lig_dict[charge]['sim_df'] = sim_df
            self.lig_dict[charge]['dif_df'] = dif_df

    def set_ligpairs(self):
        for charge, value in self.lig_dict.items():
            sim_dict = {}
            dif_dict = {}
            for i in value['sim_df'].index:
                for j in value['sim_df'].columns:
                    if value['sim_df'].loc[i, j] not in [1.0, 100] and not pd.isnull(value['sim_df'].loc[i, j]):
                        sim_dict['{} {}'.format(i, j)] = round(value['sim_df'].loc[i, j], 3)
                        dif_dict['{} {}'.format(i, j)] = value['dif_df'].loc[i, j]

            if self.metric in ['Tanimoto', 'MFP', 'MCS']:
                sim_dict = {k: v for k, v in sorted(sim_dict.items(), key=lambda item: item[1], reverse=True)}
            elif self.metric == 'SMILES':
                sim_dict = {k: v for k, v in sorted(sim_dict.items(), key=lambda item: item[1], reverse=False)}

            value['sim_dict'] = sim_dict
            value['dif_dict'] = dif_dict

    def process_map(self):
        self.set_ligdict()
        self.sim_mx()
        self.set_ligpairs()
        self.make_map()
        self.savemapJson()
        self.copyhtml()
        #self.savemap()
        #self.make_graph()

    def intersection(self, edge_list, candidate_edge):
        r1, r2 = candidate_edge.split()[0], candidate_edge.split()[1]
        for edge in edge_list:
            if r1 == edge[0] or r1 == edge[1] or r2 == edge[0] or r2 == edge[1]:
                return True  # Shortcut comparing: it's already True
        return False

    def not_ingraph(self, node_list, candidate_edge):
        r1, r2 = candidate_edge.split()[0], candidate_edge.split()[1]
        return r1 not in node_list or r2 not in node_list

    def outer_nodes(self, G):
        node_list = [node for node, val in G.degree() if val == 1]
        return node_list

    def make_map(self):
        for charge, lig in self.lig_dict.items():
            H = nx.Graph()
            lpd = self.lig_dict[charge]['sim_dict']
            dpd = self.lig_dict[charge]['dif_dict']
            if len(lig['Name']) == 1:  # In case one ligand is found alone in a charge group
                raise AttributeError('Ligand found alone in charge group.')

            # 1. Make shortest spannign tree (SPT)
            incomplete = True
            while incomplete:
                for pert, score in lpd.items():
                    if len(H.nodes) == len(lig['Name']):
                        incomplete = False
                        break
                    l1, l2 = pert.split()[0], pert.split()[1]
                    if H.has_edge(l1, l2) or H.has_edge(l2, l1):
                        continue
                    if len(H.nodes) == 0 or self.intersection(H.edges, pert) and self.not_ingraph(H.nodes, pert):
                        H.add_edge(l1, l2, weight=score)
                        break
            
            # 2. Close Cycles
            while len(self.outer_nodes(H)) != 0:
                for pert, score in lpd.items():
                    l1, l2 = pert.split()[0], pert.split()[1]
                    if l1 in self.outer_nodes(H) or l2 in self.outer_nodes(H):
                        if ((l1, l2) not in H.edges or (l2, l1) not in H.edges):
                            H.add_edge(l1, l2, weight=score)
                            break

            # 3. Ensure diameter constraint
            d = nx.diameter(H)
            while d > 8:
                for pert, score in lpd.items():
                    l1, l2 = pert.split()[0], pert.split()[1]
                    if (l1, l2) not in H.edges() and (l2, l1) not in H.edges() and self.intersection(H.edges(), pert):
                        H.add_edge(l1, l2, weight=score)
                        nd = nx.diameter(H)
                        if nd == d:
                            H.remove_edge(l1, l2)
                            continue
                        else:
                            d = nd
                            break

            # 4. k-edge augmentation 
            t = 1.0
            k_edge_comp = False
            no_complement = False
            while not k_edge_comp and not no_complement:
                complement = None
                try:
                    avail = {(k.split()[0], k.split()[1]):v for k, v in lpd.items() if v > t}
                    complement = list(nx.k_edge_augmentation(H, k = 2.0, avail=avail))
                except:
                    pass
                if complement == None:
                    t -= 0.05
                else:
                    ilegal = []
                    for edge in complement:
                        score = lpd['{} {}'.format(edge[0], edge[1])]
                        dif = dpd['{} {}'.format(edge[0], edge[1])]
                        e = (edge[0], edge[1])
                        if (dif >= 6) and score <= 0.70:
                            ilegal.append(edge)
                            continue
                        else:
                            w = lpd['{} {}'.format(edge[0], edge[1])]
                            H.add_edge(*edge, weight=w)
                    if len(ilegal) == len(complement):
                        no_complement = True
                    k_edge_comp = nx.is_k_edge_connected(H, 2.0)

            # 5. Break big cycles
            big_cycles = False
            cycles = nx.cycle_basis(H)
            breaking_edges = []
            for cyc in cycles:
                if len(cyc) >= 8:
                    big_cycles = True

            while big_cycles:
                cycles = nx.cycle_basis(H)
                c = 0
                for cyc in cycles:
                    if len(cyc) >= 8:
                        for pert, score in lpd.items():
                            l1, l2 = pert.split()[0], pert.split()[1]
                            if l1 in cyc and l2 in cyc and (l1, l2) not in H.edges() and (l2, l1) not in H.edges():
                                H.add_edge(l1, l2, weight=score)
                                breaking_edges.extend([(l1, l2), (l2, l1)])
                                break
                    else:
                        c +=1
                if c == len(cycles):
                    big_cycles = False

            # 6. Remove redundant edges
            for edge in sorted(H.edges(data=True), key=lambda t: t[2].get('weight', 1)):
                try:
                    dif = dpd['{} {}'.format(edge[0], edge[1])]
                except:
                    dif = dpd['{} {}'.format(edge[1], edge[0])]
                if ((edge[2]['weight'] <= 0.70) or (dif >= 6)) and ((edge[0], edge[1]) not in breaking_edges):
                    H.remove_edge(*edge[:2])
                    if not nx.is_connected(H):
                        H.add_edge(*edge[:2], weight=edge[2]['weight'])
                    elif nx.diameter(H) > 7:
                        H.add_edge(*edge[:2], weight=edge[2]['weight'])
                    elif not nx.is_k_edge_connected(H, k=2):
                        H.add_edge(*edge[:2], weight=edge[2]['weight'])
                    elif len(self.outer_nodes(H)) > 0:
                        H.add_edge(*edge[:2], weight=edge[2]['weight'])
            print('Number of edges: {}'.format(len(H.edges())))
            lig['Graph'] = H
    
    def make_graph(self):
        """ Plotting of the FEP network graph. """
        for charge, lig in self.lig_dict.items():
            G = nx.DiGraph()
            for e in lig['Graph'].edges():
                G.add_edge(e[0], e[1])
            dG_pos = graphviz_layout(G, prog='neato')
            labels = nx.get_edge_attributes(G, 'weight')
            nx.draw(G, pos=dG_pos, with_labels=True, node_size=1000)
            nx.draw_networkx_edges(G, dG_pos, arrowstyle='-|>', arrowsize=10)
            nx.draw_networkx_edge_labels(G, dG_pos, edge_labels=labels)
            plt.show()
    
    def clean_name(self, str):
        return str.replace('.', '-').replace('_', '-').replace('(', '-').replace(')', '-')

    def savemap(self):
        with open('{}'.format(self.otxt), 'w') as outfile:
            for charge, lig in self.lig_dict.items():
                for e in lig['Graph'].edges():
                    idx, jdx = lig['Name'].index(e[0]), lig['Name'].index(e[1])
                    e0, e1 = self.clean_name(e[0]), self.clean_name(e[1])
                    an0, an1 = lig['AtNumb'][idx], lig['AtNumb'][jdx]
                    if an0 > an1: 
                        outfile.write('{} {}\n'.format(e0, e1))
                    elif an1 > an0:
                        outfile.write('{} {}\n'.format(e1, e0))
                    else:
                        outfile.write('{} {}\n'.format(e0, e1))
        outfile.close()

    def savemapJson(self):
        for charge, lig in self.lig_dict.items():
            graph = lig['Graph']

            # Some data reformatting needs to be done for visjs
            data = json_graph.node_link_data(graph)
            data['edges'] = data.pop('links')
            for edge in data['edges']:
                edge['from'] = edge.pop('source')
                edge['to'] = edge.pop('target')
                edge['payload'] = {"ddG":"Test","ddGexpt":None}

            for node in data['nodes']:
                node['label'] = node["id"]   # maybe need unique identifiers?
                node["shape"] = "image"      # to be changed to img location
                node["image"] = "./img/{}.png".format(node["id"])
                node['payload'] = {"dG":"Test","dGexpt":None}
                #node["size"]  = 40

        with open('{}/{}.json'.format(self.wd, self.otxt), 'w') as outfile:
            outfile.write(json.dumps(data,indent = 4))

    def copyhtml(self):
        shutil.copy(s.ENV['ROOT'] + 'app/QmapFEP.html', '{}/QmapFEP.html'.format(self.wd))

class LoadData(object):
    """
        	Load experimental or user data
    """
    def __init__(self, o, ifile): # ifile needs to come from API
        self.ifile = ifile
        self.mapfile = o

        self.get_extension()
        self.read_file()
        self.populate_json()
        self.write()

    def get_extension(self):
        self.extension = self.ifile.split('.')[-1]

    def read_file(self):
        self.expt = {}
        extensions = ['csv'] # Maybe I will
        if self.extension not in extensions:
            print("can't read file")
            sys.exit()

        if self.extension == 'csv':
            with open(self.ifile) as infile:
                i = -1
                for line in infile:
                    i += 1
                    if i == 0:
                        # skip header
                        continue
                    line = line.split()
                    self.expt[line[0].strip('"')] = float(line[1])

    def populate_json(self):
        with open(self.mapfile) as infile:
            self.QmapFEPdata = json.load(infile)

        for node in self.QmapFEPdata['nodes']:
            node["payload"]['dGexpt'] = self.expt[node['id']]

        for edge in self.QmapFEPdata['edges']:
            edge["payload"]['ddGexpt'] = self.expt[edge['to']] - self.expt[edge['from']]
            edge["payload"]['ddGexpt'] = '{:.2f}'.format(edge["payload"]['ddGexpt'])

    def write(self):
        with open(self.mapfile, 'w') as outfile:
            outfile.write(json.dumps(self.QmapFEPdata,indent=4))

class Analyze(object):
    def __init__(self, o, datadir): # TO DO datadir not in API currently
        self.mapfile = o
        self.datadir = datadir

        self.readmap()
        self.populate_map()
        self.write()

    def readmap(self):
        with open(self.mapfile) as infile:
            self.data = json.load(infile)

    def populate_map(self):
        tmp = {}
        for system in ['protein', 'water']:
            with open(self.datadir + '/' + system + '_whole.txt') as infile:
                for line in infile:
                    if len(line) == 0:
                        continue
                    if 'crashes' in line:
                        continue
                    line = line.split()
                    pert = line[0].split('_')
                    From = pert[1].split('-')[0]
                    To = pert[1].split('-')[1]

                    if system == 'protein':
                        tmp[(From,To)] = {'protein' : float(line[1])}
                        tmp[(To,From)] = {'protein' : -1. * float(line[1])}
                    
                    if system == 'water':
                        tmp[(From,To)]['water'] = float(line[1])
                        tmp[(To,From)]['water'] = -1. * float(line[1])

        for edge in self.data['edges']:
            ddG = (tmp[edge['from'],edge['to']])["protein"] - (tmp[edge['from'],edge['to']])["water"]
            edge["payload"]["ddG"] = "{:.2f}" .format(ddG)
        
    def write(self):
        with open(self.mapfile, 'w') as outfile:
            outfile.write(json.dumps(self.data,indent=4))

class GenPlot(object):
    def __init__(self,o,wd):
        self.mapfile = o
        self.wd = wd

        self.load()
        self.plot()

    def load(self):
        with open(self.mapfile) as infile:
            self.data = json.load(infile)   

    def plot(self):
        x = []
        y = []
        error = []
        # generate ddG plot
        for edge in self.data['edges']:
            x.append(float(edge["payload"]["ddGexpt"]))
            y.append(float(edge["payload"]["ddG"]))
        plot.linplot(x=x,y=y,d='ddG')

class Init(object):
    def __init__(self,data):
        # Put file in memory stream. This allows the server to read uploaded file
        #  into memory and pass it as an io.BytesIO() to MapGen        
        with open(data['isdf'], "rb") as f:
            metric  = data['metric']
            o       = data['o']
            wd   = data['wd']
            with io.BytesIO(f.read()) as fio:
                # Generate Rgroup images
                Rg = GenRgroup(fio,data['isdf'],wd)
                Rg.generate_images()

                # Create the network
                mg = MapGen(fio, metric, o, wd)
                mg.process_map()

                # Show interactive map
                ## TO DO
