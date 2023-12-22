import networkx as nx
from collections import defaultdict
import sys

class Comigration:
    def __init__(self):
        self.comig2mig = dict()
        self.mig2comig = dict()
        self.n_comig = 0

    def addComig(self, l):
        self.comig2mig[self.n_comig] = l
        for i in l:
            self.mig2comig[i] = self.n_comig
        self.n_comig += 1
        return self.n_comig - 1

    def addMig(self, comig, mig):
        self.comig2mig[comig].append(mig)
        self.mig2comig[mig] = comig

    def findComig(self, mig):
        if mig in self.mig2comig:
            return self.mig2comig[mig]
        else:
            return None
    
    def getMig(self, comig):
        return self.comig2mig[comig]

    def comigs(self):
        return self.comig2mig.keys()
    

def read_from_file(clone_tree_filename, labeling_filename):
    with open(clone_tree_filename, 'r') as clone_tree_file:
        tree = nx.DiGraph()
        for edge in clone_tree_file:
            try:
                a, b = edge.split()
            except:
                raise ValueError('Ill-formatted input clone tree file')
            tree.add_edge(a,b)
        if not nx.is_tree(tree):
            raise ValueError('Input clone tree is not a valid tree')

    with open(labeling_filename, 'r') as labeling_file:
        attrs = dict()
        for label in labeling_file:
            try:
                a, l = label.split()
            except:
                raise ValueError('Ill-formatted input leaf labeling file')
            attrs[a] = {'label': l}
    nx.set_node_attributes(tree, attrs)
    return tree

def get_migrations(tree):
    M = []
    for u, v in tree.edges:
        if tree.nodes[u]['label'] != tree.nodes[v]['label']:
            M.append((u, v))
    return M

def greedyComigration(T, M=None):
    if M is None:
        M = get_migrations(T)
    C = Comigration()
    label2comig = defaultdict(list)
    root = [i for i in T.nodes if T.in_degree(i) == 0][0]
    for u, v in nx.dfs_edges(T, root):
        if (u,v) in M:
            b = True
            for c in label2comig[(T.nodes[u]['label'], T.nodes[v]['label'])]:
                bb = True
                for uu, vv in C.getMig(c):
                    if nx.has_path(T, vv, u):
                        bb = False
                        break
                if bb:
                    cc = C.addMig(c, (u,v))
                    b = False
                    break
            if b:
                cc = C.addComig([(u,v)])
                label2comig[(T.nodes[u]['label'], T.nodes[v]['label'])].append(cc)
    return C


def buildComigrationGraph(T, M, C, u):
    if T.out_degree(u) == 0:
        G = nx.DiGraph()
        G.add_nodes_from(C.comigs())
        return G, set()
    G = nx.DiGraph()
    G.add_nodes_from(C.comigs())
    X = set()
    for v in T.neighbors(u):
        Gtvc, Xv = buildComigrationGraph(T, M, C, v)
        G.add_edges_from(Gtvc.edges)
        if (u,v) in M:
            Cs = C.findComig((u,v))
            for Ct in Xv:
                G.add_edge(Cs, Ct)
            X.add(Cs)
        else:
            X = X.union(Xv)
    return G, X

def TCC(T, M=None):
    root = [i for i in T.nodes if T.in_degree(i) == 0][0]
    if M is None:
        M = get_migrations(T)
    C = greedyComigration(T, M)
    comig_graph, _ = buildComigrationGraph(T, M, C, root)
    return nx.is_directed_acyclic_graph(comig_graph)

import time
t = time.time()
T = read_from_file(sys.argv[1], sys.argv[2])
print(TCC(T), time.time() - t)


