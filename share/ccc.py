import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
import numpy as np
import matplotlib.pyplot as plt


class CCC:
    """ Cycle closure correction algorithm.
    Input a list of edges and its associated ddG values.
    Output corrected ddG values without hysteresis."""
    def __init__(self, edges, E, Esem, workdir):
        self.edges = edges
        self.E = E
        self.Esem = Esem
        self.wrkdir = workdir

    def independent_rows(self, A, tol=1e-05):
        """
        Return an array composed of independent columns of A.
        Note the answer may not be unique; this function returns one of many
        possible answers.
        http://stackoverflow.com/q/13312498/190597 (user1812712)
        http://math.stackexchange.com/a/199132/1140 (Gerry Myerson)
        http://mail.scipy.org/pipermail/numpy-discussion/2008-November/038705.html
            (Anne Archibald)
        A = np.array([(2,4,1,3),(-1,-2,1,0),(0,0,2,2),(3,6,2,5)])
        independent_columns(A)
        np.array([[1, 4],
                [2, 5],
                [3, 6]])
        """
        Q, R = np.linalg.qr(A.T)
        independent = np.where(np.abs(R.diagonal()) > tol)[0]
        return A[independent, :]


    def generate_cycles(self):
        input_list = self.edges.copy()
        tmp = {}
        # construct energy dictionary, including invert
        for i in range(0, len(input_list)):
            tmp[input_list[i]] = (self.E[i], self.Esem[i])
            invert_edge = (input_list[i][1],input_list[i][0])
            tmp[invert_edge] = (-self.E[i], self.Esem[i])

        # make the graph
        G = nx.Graph(input_list)

        # add the weights, these are the ddG values
        for edge in G.edges():
            G[edge[0]][edge[1]]['weight'] = tmp[edge][0]

        G = G.to_directed()
        simp_cyc = list(nx.simple_cycles(G))
        all_cycles = [i for i in simp_cyc if len(i) > 2]
        converged_cycles = []
        bad_cycles = []

        # identify bad cycles:
        for cycle in all_cycles:
            ddGsum = 0.
            sem = 0.
            for i in range(0,len(cycle)):
                if i < len(cycle) - 1:
                    query = (cycle[i], cycle[i+1])
                else:
                    query = (cycle[i], cycle[0])
                
                ddGsum += tmp[query][0]
                sem += (tmp[query][1] ** 2)

            sem_pooled = np.sqrt(sem)
            if abs(ddGsum) - sem_pooled  > 0.0:
                bad_cycles.append([cycle, ddGsum, sem_pooled])

            else:
                converged_cycles.append(cycle)

        return all_cycles, tmp, G
    
    def make_cccMatrix(self, all_cycles):
        Mg = np.zeros(shape=(len(all_cycles), len(self.edges)))
        k = -1
        for cyc in all_cycles:
            k += 1
            for j in range(len(cyc)):
                if j == (len(cyc) - 1):
                    edge = (cyc[j], cyc[0])
                else:
                    edge = (cyc[j], cyc[j + 1])
                for i in range(len(self.edges)):
                    if edge == self.edges[i]:
                        Mg[k, i] = 1
                    elif (edge[1], edge[0]) == self.edges[i]:
                        Mg[k, i] = -1
        Sr = []
        for i in range(Mg.shape[0]):
            nz = np.nonzero(Mg[i])[0]
            Sr.append(list(nz))

        Sr_copy = Sr.copy()
        while k:
            i = 0
            k = False
            while i < len(Sr):
                for j in range(len(Sr)):
                    if i == j:
                        continue
                    if np.in1d(Sr[i], Sr[j]).any():
                        Sr_copy[i] = list(set(Sr_copy[i] + Sr[j]))
                        Sr_copy[j] = []
                        k = True

                Sr_copy = [i for i in Sr_copy if i]
                Sr = Sr_copy.copy()
                i = i + 1

        return Sr, Mg
    
    def get_independent(self, Sr, Mg):
        Sf = Sr.copy()
        for i in range(0, len(self.edges)):
            tmp_list = [y for x in Sr for y in x]
            if i not in tmp_list:
                Sf.append([i])
        Ms_dict = {} #Cycle closure connectivity matrix
        k = 0
        for Lf in Sf:
            k += 1
            tmp_M = []
            for i in range(Mg.shape[0]):
                new_row = Mg[i][Lf]
                if np.all((new_row == 0)):
                    continue
                else:
                    tmp_M.append(new_row)

            E_lf = [self.E[i] for i in Lf]
            Ms_dict[k] = (np.array(tmp_M), E_lf, Lf)
            if not tmp_M:
                Ms_dict[k] = (np.array([0]), E_lf, Lf)

        Mr_dict = {} #Reduced cycle closure connectivity matrix
        E_list = [] #Error calculation dict
        for key, v in Ms_dict.items():
            if v[0].shape[0] > 1 and v[0].shape[1] > 1:
                tmp = self.independent_rows(v[0])
                Mr_dict[key] = (tmp, v[1])
            else:
                Mr_dict[key] = (v[0], v[1])

            for i in range(v[0].shape[0]):
                tmp = np.zeros(shape=(1, len(self.edges)))[0]
                non_z = np.nonzero(v[0][i])[0]
                for j in non_z:
                    idx = v[2][j]
                    tmp[idx] = abs((sum([v[1][j] * v[0][i][j] for j in non_z]) / np.sqrt(len(non_z))))
                E_list.append(tmp)
            
        E_n = []
        for edge in range(len(self.edges)):
            E_max = 0
            for l in range(len(E_list)):
                if E_max < E_list[l][edge]:
                    E_max = round(E_list[l][edge], 3)

            E_n.append(E_max)

        return Mr_dict, Sf, E_n

    def make_corrections(self, Mr_dict, Sf):
        F_list = []
        for isg in Mr_dict.keys():
            if not np.any(Mr_dict[isg][0]):
                F_list.append(Mr_dict[isg][1])
            else:
                ME_ = np.dot(-Mr_dict[isg][0], Mr_dict[isg][1])  # -ME
                B = np.dot(Mr_dict[isg][0], Mr_dict[isg][0].T)  # Mt
                C = np.dot(np.linalg.inv(B), ME_)
                F = Mr_dict[isg][1] + np.dot(Mr_dict[isg][0].T, C)
                F_list.append(F)
            

        Srf = [val for sublist in Sf for val in sublist]
        F = [val for sublist in F_list for val in sublist]

        F_n = []
        for i in range(len(self.edges)):
            idx = Srf.index(i)
            F_n.append(round(F[idx], 3))

        return self.edges, F_n