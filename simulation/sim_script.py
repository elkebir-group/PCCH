import networkx as nx
import random
import copy
from collections import Counter
# random.seed(2)

def pick_partition(i):
    partition = []
    part = random.randint(1, i - 1)
    partition.append(part)
    k = i - part
    while k > 0:
        part = random.randint(1,k)
        partition.append(part)
        k = k - part
    return partition

def make_paths(site, partition):
    paths =[]
    init_seq = [chr(i+ord('a')) for i in range(site)]
    init_str = init_seq[:2]
    for p in range(len(partition)):
        pos = sum(partition[:p])
        pstr = init_str + [ init_seq[(i + 2 + pos)%site] for i in range(partition[p])]
        init_str = pstr[-2:]
        paths.append(pstr)
    return paths


def make_biased_tree(paths, site):
    tree = nx.DiGraph()
    ell = {0:'o'}
    j=1
    pot = []
    vpot = []
    for path in paths:
        tree.add_edge(0, j)
        ell[j] = path[0]
        j += 1
        for i in range(1, len(path)):
            tree.add_edge(j-1, j)
            ell[j]=path[i]
            if i > 1 and i < len(path) - 1:
                tree.add_edge(j, 100*j)
                ell[100*j] = path[i-1]
                pot.append([path[i], path[i-1], path[i]])
            if i == 1 and len(path) > 3:
                ell[j] += chr(site + ord('a'))
                site += 1
                vpot.append([ell[j][1], ell[j][0]])
            j += 1
    for path in pot + vpot:
        tree.add_edge(0, j)
        ell[j] = path[0]
        j += 1
        for i in range(1, len(path)):
            tree.add_edge(j-1, j)
            ell[j]=path[i]
            j += 1
    return tree, ell, len(pot) + 1, len(vpot), site

def complete_tree(tree, ell):
    orig_nodes = copy.deepcopy(tree.nodes)
    nleaves = {}
    for node in orig_nodes:
        if tree.out_degree(node) == 0:
            nleaves[node] = 0
        else:
            neighbors = []
            for i in tree.predecessors(node):
                neighbors += list(ell[i])
            for i in tree.neighbors(node):
                neighbors += list(ell[i])
            nleaves[node] = Counter(neighbors).most_common(1)[0][1] + 1
    for node in orig_nodes:
        if node == 0:
            k = -1
        else:
            k = 1
        j = 1
        for ell_ in ell[node]:
            for l in range(nleaves[node]):
                tree.add_edge(node, k * (10000*node+100*j+l))
                ell[k * (10000*node+100*j+l)] = ell_
            j += 1

if __name__ == "__main__":
    site = random.randint(1,8)
    partition = pick_partition(site)
    site = sum(partition)
    paths = make_paths(site, partition)
    tree, ell, comig_diff, label_diff, site = make_biased_tree(paths, site)

    first = set(tree.edges)
    complete_tree(tree, ell)
    second = set(tree.edges)

    leaves = [i for i in tree.nodes if tree.out_degree(i)==0]
    
    partition.sort()
    part_string = ''
    for part in partition:
        part_string += f'_{part}'
    with open(f'sim_data/c_{comig_diff}_l_{label_diff}_s_{site}{part_string}.tree', 'w+') as ft, \
         open(f'sim_data/c_{comig_diff}_l_{label_diff}_s_{site}{part_string}.labeling', 'w+') as fl:
        for i,j in tree.edges:
            ft.write(f'{i}\t{j}\n')
        for i in ell:
            if i in leaves:
                fl.write(f'{i}\t{ell[i]}\n')

