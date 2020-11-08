## DAG module
import numpy as np
from functools import reduce

scm = {
    "A": ["D","Y"],
    "D": [],
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

    def find_parents_of_exposure(self):
        index_exp = self.matrix_index[self.exposure]
        parents_of_exp = [k for i,k in enumerate(self.scm) if self.adj_matrix[i,index_exp]==1]
        return parents_of_exp

    def find_all_paths_to_outcome(self, departure):
        explored = {k:True if k in departure else False for k in self.scm.keys()}
        # Run modified breadth-first-search starting from parents:
        unconnected_tree = bfs(
            departure,
            explored,
            self.scm,
            self.outcome,
            self.exposure,
            self.adj_matrix,
            self.matrix_index
        )
        # Connect related edges in tree:
        tree = connect_paths(unconnected_tree[0],unconnected_tree)
        return tree

    def irrelevant_parent_path(self, path):
        if self.outcome not in path:
            return float('inf')
        else:
            return 0

    def is_valid_adjustment_set(self, proposed_set):
        if self.outcome in proposed_set or self.exposure in proposed_set:
            raise Exception("You supplied the outcome variables. What's wrong with you?")
        # Find all descendants of exposure variable:
        self.find_all_descendants()
        # Find all paths through parents:
        self.find_all_paths_through_parents()
        # Update backdoor paths:
        tree = self.all_paths_through_parents
        # Paths before set was evaluated:
        potential_backdoor_paths = list(map(lambda x: {y:(self.is_collider(path=x, vertex=y) + self.irrelevant_parent_path(x)) for y in x},tree))
        # Evaluates with proposed set:
        final_backdoor_paths = list(map(lambda x: {y:x[y]+(self.is_collider(path=list(x.keys()), vertex=y)*(-2)) + 1 if y in proposed_set else x[y] for y in x},potential_backdoor_paths))
        self.final_backdoor_paths = final_backdoor_paths
        back_door_closed = list(map(lambda x: sum(list(x.values()))>=1,final_backdoor_paths))
        self.back_door_closed = back_door_closed
        if any([node in self.descendants for node in proposed_set]):
            return False
        else:
            return all(back_door_closed)

    def find_all_paths_through_parents(self):
        parents = self.find_parents_of_exposure()
        if len(parents)>0:
            tree = reduce(lambda x,y: x+y, [self.find_all_paths_to_outcome(parent) for parent in parents])
        else:
            tree = []
        self.all_paths_through_parents = tree

    def find_all_descendants(self):
        departure = self.exposure
        descendent_paths = self.find_all_paths_to_outcome(departure)
        # Filter for paths leading into exposure:
        descendent_paths = [path for path in descendent_paths if self.is_descendant(self.exposure, path[1])]
        if len(descendent_paths)==0:
            raise Exception("Exposure has no descendants. No causal effect to measure.")
        descendants = list(set(reduce(lambda x,y: x+y, descendent_paths)))
        if self.outcome not in descendants:
            raise Exception("Outcome variable is not a descendant of exposure variables. No causal effect to measure.")
        descendants = list(filter(lambda x: x not in [self.exposure, self.outcome], descendants))
        self.descendants = descendants

    # Check if node is collider on given path:
    def is_collider(self, path, vertex):
        if vertex == path[0] or vertex==path[-1]:
            return int(False)
        else:
            idx = path.index(vertex)
            parent = path[idx-1]
            child = path[idx+1]
            col = self.matrix_index[vertex]
            row_parent = self.matrix_index[parent]
            row_child = self.matrix_index[child]
            return int(self.adj_matrix[row_parent,col]==1 and self.adj_matrix[row_child,col]==1)

    def is_descendant(self, parent, child):
        return self.adj_matrix[self.matrix_index[parent], self.matrix_index[child]]==1

dag = Dag(scm, outcome, exposure)

dag.is_valid_adjustment_set("A")
