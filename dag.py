## DAG module
import numpy as np
from functools import reduce

scm = {
    "A": ["C"],
    "B": ["E", "Y"],
    "C": [],
    "D": ["C"],
    "E": ["A","D"],
    "X": ["E", "A"],
    "Y": []
}

outcome = "Y"
exposure = "D"

set([1]) == {1}

def build_adj_matrix(scm):
    d = len(scm)
    mat = np.zeros((d,d))
    for i, key_i in enumerate(scm):
        for j, key_j in enumerate(scm):
            if key_j in scm[key_i]:
                mat[i,j] = 1
    return mat

class Dag():
    def __init__(self, scm, outcome, exposure):
        self.scm = scm # list of dictionaries
        self.outcome = outcome # outcome variable
        self.exposure = exposure # exposure variable
        self.adj_matrix = build_adj_matrix(scm) # build adj matrix
        self.matrix_index = {k:i for i,k in enumerate(scm)}

    def find_all_paths(self):
        index_exp = self.matrix_index[self.exposure]
        parents_of_exp = [k for i,k in enumerate(self.scm) if self.adj_matrix[i,index_exp]==1]
        explored = {k:True if k in parents_of_exp else False for k in self.scm.keys()}
        L = [parents_of_exp]
        L_i = L[0]
        counter = 0
        tree = []
        while (len(L_i)!=0 and set(L_i)!={self.outcome}):
            L += [[]]
            tree += [[]]
            vertices_to_explore = [i for i in L[counter] if i!=self.outcome] # explore only non-outcome variables
            counter_vertex = 0
            for vertex in vertices_to_explore:
                # Explose all edges except the one leading to exposure:
                current_pos = self.matrix_index[vertex]
                neighbours = [k for i,k in enumerate(self.scm) if self.adj_matrix[current_pos,i]==1 and k!=self.exposure]
                neighbours += [k for i,k in enumerate(self.scm) if self.adj_matrix[i,current_pos]==1 and k!=self.exposure]
                for neighbour in neighbours:
                    if explored[neighbour] == False:
                        L[counter+1] += neighbour
                        if neighbour!="Y":
                            explored[neighbour] = True
                        tree[counter] += [(vertex,neighbour)]
                counter_vertex += 1
            counter += 1
            L_i = L[counter]
        self.all_paths_through_parents = tree

    def valid_adjustment_sets(self):
        pass

    def is_valid_adjustment_set(self, proposed_set):
        pass

dag = Dag(scm, outcome, exposure)

dag.adj_matrix
dag.find_all_paths()

dag.all_paths_through_parents
