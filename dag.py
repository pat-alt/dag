## DAG module
import numpy as np
from functools import reduce
import networkx as nx # only for charting

# Helper function that builds adjacency matrix:
def build_adj_matrix(scm):
    d = len(scm)
    mat = np.zeros((d,d))
    for i, key_i in enumerate(scm):
        for j, key_j in enumerate(scm):
            if key_j in scm[key_i]:
                mat[i,j] = 1
    return mat

# Helper function that connects realted edges:
def connect_paths(starting_edges, tree):
    final_paths = []
    for edges in starting_edges:
        L = [[edges]]
        counter = 0
        while counter < (len(tree)-1):
            L += [[]]
            children = tree[counter+1]
            for edge in L[counter]:
                for child in children:
                    path_to_here = edge
                    if path_to_here[-1]==child[0]:
                        L[counter+1] += [path_to_here + [child[1]]]
            counter += 1
        paths = [i for i in L if len(i)>0][-1]
        final_paths += paths
    return final_paths

# A modified breadth-first-search from starting nodes to outcome:
def bfs(starting_nodes, explored, scm, outcome, exposure, adj_matrix, matrix_index):
    L = [starting_nodes]
    L_i = L[0]
    counter = 0
    tree = []
    while (len(L_i)!=0 and set(L_i)!={outcome}):
        L += [[]]
        tree += [[]]
        nodes_to_explore = [i for i in L[counter] if i!=outcome] # explore only non-outcome variables
        counter_node = 0
        for node in nodes_to_explore:
            # Explose all edges except the one leading to exposure:
            current_pos = matrix_index[node]
            neighbours = [k for i,k in enumerate(scm) if adj_matrix[current_pos,i]==1 and k!=exposure]
            neighbours += [k for i,k in enumerate(scm) if adj_matrix[i,current_pos]==1 and k!=exposure]
            for neighbour in neighbours:
                if explored[neighbour] == False:
                    L[counter+1] += neighbour
                    if neighbour!=outcome:
                        explored[neighbour] = True
                    tree[counter] += [[node,neighbour]]
            counter_node += 1
        if len(tree[counter])==0:
            tree = tree[0:counter]
        counter += 1
        L_i = L[counter]
    return tree

# Helper function that prepares graph to be plotted:
def create_edges(scm):
    L=[]
    for k, v in scm.items():
        for val in v:
            L.append((k,val))
    return L

## Dag class: ----
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

    # Paths through parents that do not lead to outcome:
    def irrelevant_parent_path(self, path):
        if self.outcome not in path:
            return float('inf')
        else:
            return 0

    # Main method:
    def is_valid_adjustment_set(self, proposed_set):
        ## 1.) Simple checks: ----
        if self.outcome in proposed_set or self.exposure in proposed_set:
            raise Exception("You should not supply the outcome or exposure variable as a proposed set")
        # Find all descendants of exposure variable:
        descendants = self.find_descendants([], self.exposure)
        if len(descendants)==0:
            raise Exception("Exposure has no descendants. No causal effect to measure.")
        if self.outcome not in descendants:
            raise Exception("Outcome not caused by exposure.")
        if self.exposure in descendants:
            raise Exception("Exposure is a descendant of the outcome variable. Are you sure you supplied an acyclical graph?")
        # Check for descendants
        includes_descendant = any([node in descendants for node in proposed_set])
        if includes_descendant:
            print("Proposed set includes descendant of exposure variable.")
            return False
        ## 2.) Backdoor checks: ----
        else:
            # Find all paths through parents:
            self.find_all_paths_through_parents()
            # Update backdoor paths:
            tree = self.all_paths_through_parents
            # Paths before set was evaluated:
            potential_backdoor_paths = list(map(lambda x: {y:(self.is_collider(path=x, node=y) + self.irrelevant_parent_path(x)) for y in x},tree))
            # Evaluates with proposed set:
            final_backdoor_paths = list(
                map(lambda x: {y:x[y]+(self.is_collider(path=list(x.keys()), node=y)*(-2)) + 1 if y in proposed_set else x[y] for y in x},potential_backdoor_paths)
            )
            self.final_backdoor_paths = final_backdoor_paths
            back_door_closed = list(map(lambda x: sum(list(x.values()))>=1,final_backdoor_paths))
            self.back_door_closed = back_door_closed
            return all(back_door_closed)

    def find_all_paths_through_parents(self):
        parents = self.find_parents_of_exposure()
        if len(parents)>0:
            tree = reduce(lambda x,y: x+y, [self.find_all_paths_to_outcome(parent) for parent in parents])
        else:
            tree = []
        self.all_paths_through_parents = tree

    def find_descendants(self, descendants, node):
        scm = self.scm
        exposure = self.exposure
        len_before = len(descendants)
        # add descendants
        for desc in scm[node]:
            if desc not in descendants:
                descendants.append(desc)
        len_after = len(descendants)
        if len_after == len_before:
            return descendants
        else:
            for nested_node in descendants:
                descendants = self.find_descendants(descendants, nested_node)
        return descendants

    # Check if node is collider on given path:
    def is_collider(self, path, node):
        if node == path[0] or node==path[-1]:
            return int(False)
        else:
            idx = path.index(node)
            parent = path[idx-1]
            child = path[idx+1]
            col = self.matrix_index[node]
            row_parent = self.matrix_index[parent]
            row_child = self.matrix_index[child]
            return int(self.adj_matrix[row_parent,col]==1 and self.adj_matrix[row_child,col]==1)

    def is_descendant(self, parent, child):
        return self.adj_matrix[self.matrix_index[parent], self.matrix_index[child]]==1

    # Plot graph:
    def plot(self):
        gr = nx.DiGraph()
        gr.add_nodes_from(self.scm)
        edges = create_edges(self.scm)
        gr.add_edges_from(edges)
        nx.draw(gr, with_labels=True)
