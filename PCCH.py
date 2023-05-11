import os
import argparse
import time
import networkx as nx
import gurobipy as gp
from gurobipy import GRB

def process_args():
    """ Parses the command line arguments """

    parser = argparse.ArgumentParser(description='PCCH')

    parser.add_argument('phylogeny', type=str, help='Input phylogeny')
    parser.add_argument('leaf_labeling', type=str, help='Input leaf labeling')

    parser.add_argument('-p', '--primary', type=str, help='Primary location')
    parser.add_argument('-c', '--colormap', metavar='COLORMAP', type=str, help='Color map file', action='store')
    parser.add_argument('--log', action='store_true', default=False, help='Outputs Gurobi logging')
    parser.add_argument('-o', '--output', action='store', default=None, help='Output folder')
    parser.add_argument('-t', '--threads', type=int, help='Number of threads')
    
    return parser.parse_args()


class CloneTree:
    """ Represents a clone tree

        Member variables:
            tree: nx.DiGraph
                Clone Tree
            leaves: list(nx.DiGraph.nodes)
                List of leaves
            root: nx.DiGraph.nodes
                Root
            sites: list(str)
                List of anatomical sites
            primary_site: str
                Primary anatomical site
            nodes: list(str)
                List of nodes
            edges: list(str)
                List of edges
            n_nodes: int
                Number of nodes
            n_edges: int
                Number of edges
            n_sites: int
                Number of sites
            paths: dict(str: list(str))
                Path from root to a leaf. paths[leaf] returns the corresponding path
            max_height: int
                Maximum height of the tree
    """
    
    def __init__(self, clone_tree_filename = None, leaf_labeling_filename = None, edges = None, labels = None, timestamps = None, primary_site = None):
        """
        builds the clone tree from files or list of edges with labels

        :param clone_tree_filename: str
            name of the file containing list of edges
        :param leaf_labeling filename: str
            name of the file containing list of leaf labels
        :param edges: list(str)
            list of edges
        :param labels: dict(str, dict('label', str))
            dict of node labels
        :param primary_site: str
            Primary site
        """

        if clone_tree_filename is not None:
            with open(clone_tree_filename, 'r') as clone_tree_file:
                self.tree = nx.DiGraph()
                for edge in clone_tree_file:
                    try:
                        a, b = edge.split()
                    except:
                        raise ValueError('Ill-formatted input clone tree file')
                    self.tree.add_edge(a,b)
                if not nx.is_tree(self.tree):
                    raise ValueError('Input clone tree is not a valid tree')

            self.leaves = [i for i in self.tree.nodes if self.tree.out_degree(i) == 0]
            self.root = [i for i in self.tree.nodes if self.tree.in_degree(i) == 0][0]
            self.sites = []
            with open(leaf_labeling_filename, 'r') as leaf_labeling_file:
                attrs = dict()
                for label in leaf_labeling_file:
                    try:
                        a, l = label.split()
                    except:
                        raise ValueError('Ill-formatted input leaf labeling file')
                    attrs[a] = {'label': l}
                    if l not in self.sites:
                        self.sites.append(l)
            nx.set_node_attributes(self.tree, attrs)
            for leaf in self.leaves:
                if 'label' not in self.tree.nodes[leaf]:
                    raise ValueError(f'Site label of clone {leaf} is missing')
            if primary_site is not None and primary_site not in self.sites:
                raise ValueError('Primary site not found')
            else:
                self.primary_site = primary_site
            self.n_sites = len(self.sites)
        else:
            self.tree = nx.DiGraph()
            for edge in edges:
                self.tree.add_edge(edge[0], edge[1])
            if not nx.is_tree(self.tree):
                raise RuntimeError('Output is not a tree. Possible bug')
            nx.set_node_attributes(self.tree, labels)
            if timestamps is not None:
                nx.set_edge_attributes(self.tree, timestamps)
            self.primary_site = primary_site
            self.leaves = [i for i in self.tree.nodes if self.tree.out_degree(i) == 0]
            self.root = [i for i in self.tree.nodes if self.tree.in_degree(i) == 0][0]

        self.nodes = [v for v in self.tree.nodes]
        self.edges = [e for e in self.tree.edges]

        self.n_nodes = len(self.nodes)
        self.n_edges = len(self.edges)
        self.n_leaves = len(self.leaves)

    def get_label(self, node):
        """
        returns label of the node

        :param node: str
            query node
        :returns str

        """
        return self.tree.nodes[node]['label']

    def get_timestamp(self, edge):
        """
        returns label of the node

        :param node: str
            query node
        :returns str

        """
        return self.tree.edges[edge]['timestamp']

    def get_parent_arc(self, node):
        """
        returns the parent of the node
        :param node: str
            query node
        :returns str

        """
        return [k for k in self.tree.in_edges(nbunch=node)][0]

    def get_children_arcs(self, node):
        """
        returns the list of children for the input node
        :param node: str
            query node
        :returns str
        """
        return [k for k in self.tree.out_edges(nbunch=node)]

    def infer_paths(self):
        """ 
        sets paths variable and max height of the tree
        """
        self.paths = {}
        self.max_height = 0
        for u in self.leaves:
            self.paths[u] = []
            path_u = nx.shortest_path(self.tree, self.root, u)
            for ii in range(len(path_u) - 1):
                self.paths[u].append( (path_u[ii], path_u[ii+1]) )
            if len(self.paths[u]) > self.max_height:
                self.max_height = len(self.paths[u])

    def write_dot(self, filename, colormap=None, vertex_labels=True, timestamps=False):
        """

        :param filename: 
        :param colormap:  (Default value = None)
        :param vertex_labels:  (Default value = True)

        """
        with open(filename, 'w+') as f:
            f.write('digraph T {\n\t{\n\t\trank=same\n')
            node_edge_index = {}
            ind = 0
            for l in self.leaves:
                f.write(f'\t\t{ind} [penwidth=3,colorscheme=set19,color={colormap[self.get_label(l)]},label="{l}\\n{self.get_label(l)}"]\n')
                node_edge_index[l] = ind
                ind += 1
            f.write('\t}\n')
            for v in self.tree.nodes:
                if v not in self.leaves:
                    if vertex_labels == True:
                        f.write(f'\t{ind} [penwidth=3,colorscheme=set19,color={colormap[self.get_label(v)]},label="{v}"]\n')
                    else:
                        f.write(f'\t{ind} [penwidth=3,colorscheme=set19,color=0,label=""]\n')
                    node_edge_index[v] = ind
                    ind += 1
            for i,j in self.tree.edges:
                if vertex_labels == True:
                    if self.get_timestamp((i,j)) != -1:
                        f.write(f'\t{node_edge_index[i]} -> {node_edge_index[j]} [penwidth=3,colorscheme=set19,' + 
                            f'color=\"{colormap[self.get_label(i)]};0.5:{colormap[self.get_label(j)]}\",' + 
                            f'label=\"{self.get_timestamp((i,j))}\"]\n')
                    else:
                        f.write(f'\t{node_edge_index[i]} -> {node_edge_index[j]} [penwidth=3,colorscheme=set19,' + 
                            f'color=\"{colormap[self.get_label(i)]};0.5:{colormap[self.get_label(j)]}\"]\n')
                else:
                    f.write(f'\t{node_edge_index[i]} -> {node_edge_index[j]} [penwidth=3,colorscheme=set19,' + 
                        f'color=0]\n')
            f.write('}\n')

    def write_tree(self, filename):
        """

        :param filename: 

        """
        with open(filename, 'w+') as f:
            for edge in self.tree.edges:
                f.write(f'{edge[0]} {edge[1]}\n')

    def write_labeling(self, filename):
        """

        :param filename: 

        """
        with open(filename, 'w+') as f:
            for node in self.tree.nodes:
                f.write(f'{node} {self.get_label(node)}\n')
        

class ILPSolver:
    """ """

    def __init__(self, clone_tree, logfile=None, n_threads=None):
        self.clone_tree = clone_tree
        self.nSolutions = 1
        
        self.node_index = { p:i for i,p in enumerate(clone_tree.nodes) }
        self.edge_index = { p:i for i,p in enumerate(clone_tree.edges) }
        self.site_index = { p:i for i,p in enumerate(clone_tree.sites) }
        
        self.m = gp.Model('ILP')
        self.m._logfile = logfile
        self.m.setParam(GRB.param.LogToConsole, 0)
        self.m.setParam(GRB.param.LogFile, logfile)
        if n_threads is not None:
            self.m.setParam(GRB.Param.Threads, n_threads)

        self.add_ILP_vars()
        self.add_vertex_leaf_constraints()
        self.add_part_constraints()
        self.add_constraints_for_g()
        if self.clone_tree.primary_site is not None:
            self.set_primary_site()

        self.set_optimization_function()

        self.m.setParam(GRB.Param.PoolSolutions, 1)
        self.m.setParam(GRB.Param.PoolSearchMode, 2)
        self.m.setParam(GRB.Param.MIPGap, 0)

        start_time = time.time()
        self.m.optimize()
        self.total_time = time.time() - start_time

        self.nActualSolutions = self.m.SolCount
        self.n_migrations, self.n_comigrations = self.compute_summary()


    def add_ILP_vars(self):
        """
        sets up ILP variables
        """
        self.l = self.m.addMVar((self.clone_tree.n_nodes, self.clone_tree.n_sites), vtype=GRB.BINARY, name="l")
        self.g = self.m.addMVar((self.clone_tree.n_edges, self.clone_tree.n_sites, self.clone_tree.n_sites, self.clone_tree.n_edges), \
            vtype=GRB.BINARY, name="g")
        self.pi = self.m.addMVar((self.clone_tree.n_edges, self.clone_tree.n_sites, self.clone_tree.n_sites), vtype=GRB.BINARY, name="pi")


    def add_vertex_leaf_constraints(self):
        """ """
        for v in range(self.clone_tree.n_nodes):
            sum1 = 0
            for s in range(self.clone_tree.n_sites):
                sum1 += self.l[ v, s ]
            self.m.addConstr(sum1 == 1)
            if self.clone_tree.nodes[v] in self.clone_tree.leaves:
                hat_ell_v = self.site_index[self.clone_tree.get_label(self.clone_tree.nodes[v])]
                self.m.addConstr( self.l[ v, hat_ell_v ] == 1 )


    def add_part_constraints(self):
        """ """
        for e in range(self.clone_tree.n_edges):
            sum1 = 0
            sum2 = 0
            sum3 = 0
            for s in range(self.clone_tree.n_sites):
                for t in range(self.clone_tree.n_sites):
                    sum1 += self.pi[ e, s, t ]
                    if e < self.clone_tree.n_edges - 1:
                        sum2 += self.pi[ e+1, s, t ]
                    for uv in range(self.clone_tree.n_edges):
                        sum3 += self.g[ uv, s, t, e ]
            self.m.addConstr( sum1 <= 1)
            self.m.addConstr( sum3 >= sum1 )
            if e < self.clone_tree.n_edges - 1:
                self.m.addConstr(sum1 >= sum2)

        for vw in range(self.clone_tree.n_edges):
            v = self.clone_tree.edges[vw][0]
            if v != self.clone_tree.root:
                uv = self.edge_index[self.clone_tree.get_parent_arc(v)]
                for E in range(self.clone_tree.n_edges):
                    sum1 = 0
                    sum2 = 0
                    for s in range(self.clone_tree.n_sites):
                        for t in range(self.clone_tree.n_sites):
                            for e in range(E + 1):
                                sum1 += self.g[ uv, s, t, e ]
                                sum2 += self.g[ vw, s, t, e ]
                    self.m.addConstr( sum1 >= sum2 )


    def add_constraints_for_g(self):
        for uv in range(self.clone_tree.n_edges):
            u = self.node_index[self.clone_tree.edges[uv][0]]
            v = self.node_index[self.clone_tree.edges[uv][1]]
            sum1 = 0
            for s in range(self.clone_tree.n_sites):
                sum2 = 0
                sum3 = 0
                for t in range(self.clone_tree.n_sites):
                    for e in range(self.clone_tree.n_edges):
                        sum1 += self.g[ uv, s, t, e ]
                        sum2 += self.g[ uv, s, t, e ]
                        sum3 += self.g[ uv, t, s, e ]
                        self.m.addConstr(self.g[ uv, s, t, e ] <= self.pi[ e, s, t ])
                self.m.addConstr(sum2 == self.l[ u, s ])
                self.m.addConstr(sum3 == self.l[ v, s ])
            self.m.addConstr(sum1 == 1)


    def set_optimization_function(self):
        """ """
        sum1 = 0
        for s in range(self.clone_tree.n_sites):
            for t in range(self.clone_tree.n_sites):
                if s != t:
                    for uv in range(self.clone_tree.n_edges):
                        for e in range(self.clone_tree.n_edges):
                            sum1 += self.g[ uv,s,t,e ]
        sum2 = 0
        for s in range(self.clone_tree.n_sites):
            for t in range(self.clone_tree.n_sites):
                if s != t:
                    for e in range(self.clone_tree.n_edges):
                        sum2 += self.pi[ e,s,t ]

        self.m.setObjective( ( (self.clone_tree.n_sites ** 2) * (self.clone_tree.n_sites + 1) ) * sum1 + (self.clone_tree.n_sites + 1) \
            * sum2, GRB.MINIMIZE)
        
    def set_primary_site(self):
        self.m.addConstr( self.l[self.node_index[self.clone_tree.root], self.site_index[self.clone_tree.primary_site]] == 1)


    def compute_summary(self, soln_ind=None):
        """

        :param soln_ind:  (Default value = None)

        """
        if soln_ind == None:
            self.m.setParam(GRB.Param.SolutionNumber, 0)
        else:
            self.m.setParam(GRB.Param.SolutionNumber, soln_ind)
        
        sum1 = 0
        for s in range(self.clone_tree.n_sites):
            for t in range(self.clone_tree.n_sites):
                if s != t:
                    for uv in range(self.clone_tree.n_edges):
                        for e in range(self.clone_tree.n_edges):
                            sum1 += self.g[ uv,s,t,e ].X

        sum2 = 0
        for s in range(self.clone_tree.n_sites):
            for t in range(self.clone_tree.n_sites):
                if s != t:
                    for e in range(self.clone_tree.n_edges):
                        sum2 += self.pi[ e,s,t ].X

        return round(sum1), round(sum2)


    def infer_tree(self):
        """

        :param soln: 

        """
        # self.m.setParam(GRB.Param.SolutionNumber, soln)
        edges = self.clone_tree.edges
        labels = dict()

        ts_map = {}
        ee = 1
        for tt in range(self.clone_tree.n_edges):
            for s in range(self.clone_tree.n_sites):
                for t in range(self.clone_tree.n_sites):
                    if s!= t and self.pi[ tt,s,t ].X > 0.5:
                        ts_map[tt] = ee
                        ee += 1
        for v in range(self.clone_tree.n_nodes):
            for s in range(self.clone_tree.n_sites):
                if self.l[ v,s ].X > 0.5:
                    labels[ self.clone_tree.nodes[v] ] = { 'label' : self.clone_tree.sites[s] }

        timestamps = dict()
        for uv in range(self.clone_tree.n_edges):
            for s in range(self.clone_tree.n_sites):
                for t in range(self.clone_tree.n_sites):
                    for e in range(self.clone_tree.n_edges):
                        if self.g[ uv,s,t,e ].X > 0.5:
                            if e in ts_map:
                                timestamps[ self.clone_tree.edges[uv] ] = { 'timestamp' : ts_map[e] }
                            else:
                                timestamps[ self.clone_tree.edges[uv] ] = { 'timestamp' : -1 }

        return CloneTree(edges=edges, labels=labels, timestamps=timestamps)
    

def process_colormap(colormap_filename = None, clone_tree = None):
    """

    :param colormap_filename:  (Default value = None)
    :param clone_tree:  (Default value = None)

    """
    colormap = {}
    if colormap_filename is not None:
        with open(colormap_filename, 'r') as f:
            for line in f:
                site, color = line.split()
                colormap[site] = color
    else:
        colormap = {s:i for i,s in enumerate(clone_tree.sites)}

    return colormap


if __name__ == '__main__':
    args = process_args()
    clone_tree = CloneTree(clone_tree_filename=args.phylogeny, leaf_labeling_filename=args.leaf_labeling, primary_site=args.primary)

    if args.output is None:
        output_str = '.'
    else:
        output_str = args.output
        if not os.path.exists(output_str):
            os.makedirs(output_str)

    if args.primary is None:
        primary_str = 'ALL'
    else:
        primary_str = args.primary

    if args.log == True:
        logfile = f'{output_str}/{primary_str}-log.txt'
    else:
        logfile = ''

    soln = ILPSolver(clone_tree, logfile=logfile, n_threads=args.threads)


    colormap = process_colormap(colormap_filename=args.colormap, clone_tree=clone_tree)

    T_prime = soln.infer_tree()
    primary_str = T_prime.get_label(T_prime.root)
    T_prime.write_dot(f'{output_str}/T-{primary_str}.dot', colormap=colormap)
    T_prime.write_tree(f'{output_str}/T-{primary_str}.tree')
    T_prime.write_labeling(f'{output_str}/T-{primary_str}.labeling')
    summary = soln.compute_summary()
    print(f'{primary_str}-\t{summary[0]}\t{summary[1]}\tOptimal\t\t{soln.total_time}')
