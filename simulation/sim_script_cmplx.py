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

def make_paths(site, partition, paths, arcs):
    init_seq = [chr(i+ord('a')) for i in range(site)]
    temp_paths =[]
    init_str = init_seq[:2]
    for p in range(len(partition)):
        pos = sum(partition[:p])
        pstr = init_str + [ init_seq[(i + 2 + pos)%site] for i in range(partition[p])]
        init_str = pstr[-2:]
        temp_paths.append(pstr)
    paths = paths + temp_paths
    
    k = site
    arc_paths = []
    for arc in arcs:
        arc_paths.append([init_seq[arc[0]], init_seq[(arc[0]+1)%site], chr(k+ord('a')), init_seq[arc[1]], init_seq[(arc[1]+1)%site]])
        k += 1
    return paths, arc_paths, k



def make_biased_tree(paths, arcs, site):
    tree = nx.DiGraph()
    ell = {0:'o'}
    j=1
    vpot = []
    for path in paths:
        tree.add_edge(0, j)
        ell[j] = path[0]
        j += 1
        for i in range(1, len(path)):
            tree.add_edge(j-1, j)
            ell[j]=path[i]
            if i == 1 and len(path) > 3:
                ell[j] += chr(site + ord('a'))
                site += 1
                vpot.append([ell[j][1], ell[j][0]])
            for arc in arcs:
                if arc[0] in ell[j-1] and arc[1] in ell[j]:
                    sp = nx.shortest_path(tree, source=0, target=j)
                    path_cms = set()
                    for i in range(len(sp)-1):
                        path_cms.add((ell[sp[i]], ell[sp[i+1]]))
                    if (arc[3],arc[4]) not in path_cms:
                        tree.add_edge(j, -100*j)
                        ell[-100*j]=arc[2]
                        tree.add_edge(-100*j, -100*j-1)
                        ell[-100*j-1] = arc[3]
                        tree.add_edge(-100*j-1, -100*j-2)
                        ell[-100*j-2] = arc[4]
                        arcs.remove(arc)
                        break
            j += 1
    for path in arcs+vpot:
        tree.add_edge(0, j)
        ell[j] = path[0]
        j += 1
        for i in range(1, len(path)):
            tree.add_edge(j-1, j)
            ell[j]=path[i]
            j += 1
    return tree, ell, len(vpot), site

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
    arcs = [(1,3)]
    site = sum(partition)
    paths = []
    paths, arc_paths, site = make_paths(site, partition, paths, arcs)
    # print(paths, arc_paths)
    tree, ell, label_diff, site = make_biased_tree(paths, arc_paths, site)
    # print(tree.edges)
    complete_tree(tree, ell)
    leaves = [i for i in tree.nodes if tree.out_degree(i)==0]
    partition.sort()
    part_string = ''
    for part in partition:
        part_string += f'_{part}'
    arc_string = ''
    for arc in arcs:
        arc_string += f'_{arc[0]}-{arc[1]}'
    with open(f'sim_data/cmplx_l_{label_diff}_s_{site}{part_string}{arc_string}.tree', 'w+') as ft, \
         open(f'sim_data/cmplx_l_{label_diff}_s_{site}{part_string}{arc_string}.labeling', 'w+') as fl:
        for i,j in tree.edges:
            ft.write(f'{i}\t{j}\n')
        for i in ell:
            if i in leaves:
                fl.write(f'{i}\t{ell[i]}\n')

