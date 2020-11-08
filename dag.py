## DAG module
import numpy as np
from functools import reduce

scm = {
    "A": ["C"],
    "B": ["E", "Y"],
    "C": ["Y"],
    "D": ["C"],
    "E": ["A","D"],
    "X": ["E", "A"],
    "Y": []
}

outcome = "Y"
exposure = "D"

def build_adj_matrix(scm):
    d = len(scm)
    mat = np.zeros((d,d))
    for i, key_i in enumerate(scm):
        for j, key_j in enumerate(scm):
            if key_j in scm[key_i]:
                mat[i,j] = 1
    return mat

def connect_paths(connected_paths, tree):
    for path in connected_paths:
        counter = 0
        while counter < (len(tree)-1):
            children = tree[counter+1]
            for child in children:
                if path[-1]==child[0]:
                    path += child[1]
            counter += 1
    return connected_paths

# A modified breadth-first-search:
def bfs(starting_nodes, explored, scm, outcome, exposure, adj_matrix, matrix_index):
    L = [starting_nodes]
    L_i = L[0]
    counter = 0
    tree = []
    while (len(L_i)!=0 and set(L_i)!={outcome}):
        L += [[]]
        tree += [[]]
        vertices_to_explore = [i for i in L[counter] if i!=outcome] # explore only non-outcome variables
        counter_vertex = 0
        for vertex in vertices_to_explore:
            # Explose all edges except the one leading to exposure:
            current_pos = matrix_index[vertex]
            neighbours = [k for i,k in enumerate(scm) if adj_matrix[current_pos,i]==1 and k!=exposure]
            neighbours += [k for i,k in enumerate(scm) if adj_matrix[i,current_pos]==1 and k!=exposure]
            for neighbour in neighbours:
                if explored[neighbour] == False:
                    L[counter+1] += neighbour
                    if neighbour!="Y":
                        explored[neighbour] = True
                    tree[counter] += [[vertex,neighbour]]
            counter_vertex += 1
        if len(tree[counter])==0:
            tree = tree[0:counter]
        counter += 1
        L_i = L[counter]
    return tree

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
        # Run modified breadth-first-search starting from parents:
        unconnected_tree = bfs(
            parents_of_exp,
            explored,
            self.scm,
            self.outcome,
            self.exposure,
            self.adj_matrix,
            self.matrix_index
        )
        # Connect related edges in tree:
        tree = connect_paths(unconnected_tree[0],unconnected_tree)
        self.all_paths_through_parents = tree

    def valid_adjustment_sets(self):
        pass

    def is_valid_adjustment_set(self, proposed_set):
        pass

dag = Dag(scm, outcome, exposure)

dag.adj_matrix
dag.find_all_paths()

dag.all_paths_through_parents
