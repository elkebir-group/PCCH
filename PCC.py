import os
import argparse
import time
import networkx as nx
import gurobipy as gp
from gurobipy import GRB

def process_args():
    """ Parses the command line arguments """

    parser = argparse.ArgumentParser(description='PCC')

    parser.add_argument('phylogeny', type=str, help='Input phylogeny')
    parser.add_argument('vertex_labeling', type=str, help='Input location labeling')

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
    
    def __init__(self, clone_tree_filename = None, vertex_labeling_filename = None, primary_site = None):
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

        self.nodes = self.tree.nodes
        self.leaves = [i for i in self.tree.nodes if self.tree.out_degree(i) == 0]
        self.root = [i for i in self.tree.nodes if self.tree.in_degree(i) == 0][0]
        self.sites = []
        with open(vertex_labeling_filename, 'r') as vertex_labeling_file:
            attrs = dict()
            for label in vertex_labeling_file:
                try:
                    a, l = label.split()
                except:
                    raise ValueError('Ill-formatted input leaf labeling file')
                attrs[a] = {'label': l}
                if l not in self.sites:
                    self.sites.append(l)
        nx.set_node_attributes(self.tree, attrs)
        for node in self.tree.nodes:
            if 'label' not in self.tree.nodes[node]:
                raise ValueError(f'Site label of clone {node} is missing')
        if primary_site is not None and primary_site not in self.sites:
            raise ValueError('Primary site not found')
        else:
            self.primary_site = primary_site
        self.n_sites = len(self.sites)

        self.nodes = [v for v in self.tree.nodes]
        self.edges = [e for e in self.tree.edges]

        self.n_nodes = len(self.nodes)
        self.n_edges = len(self.edges)
        self.n_leaves = len(self.leaves)

        self.mig_edges = [(u,v) for u,v in self.tree.edges 
                             if self.tree.nodes[u]['label'] != self.tree.nodes[v]['label']]
        self.n_migs = len(self.mig_edges)

    def set_timestamps(self, timestamps):
        self.timestamps = {}
        for uv in self.edges:
            if uv in timestamps:
                self.timestamps[uv] = {'timestamp':timestamps[uv]}
            else:
                self.timestamps[uv] = {'timestamp': -1}
        nx.set_edge_attributes(self.tree, self.timestamps)

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
        self.mig_index = { p:i for i,p in enumerate(clone_tree.mig_edges) }
        self.clone_tree.infer_paths()
        self.X = set()
        for leaf in self.clone_tree.paths:
            path = self.clone_tree.paths[leaf]
            current = None
            for uv in path:
                if uv in self.clone_tree.mig_edges:
                    if current is None:
                        current = uv
                    else:
                        self.X.add((current, uv))
                        current = uv
        # print(self.X)
        self.m = gp.Model('ILP')
        self.m._logfile = logfile
        self.m.setParam(GRB.param.LogToConsole, 0)
        self.m.setParam(GRB.param.LogFile, logfile)
        if n_threads is not None:
            self.m.setParam(GRB.Param.Threads, n_threads)

        self.add_ILP_vars()
        # self.add_vertex_leaf_constraints()
        self.add_timestamp_constraints()
        self.add_comigration_constraints()
        self.set_additional_constraints()
        # self.add_constraints_for_g()
        # if self.clone_tree.primary_site is not None:
        #     self.set_primary_site()

        self.set_optimization_function()

        self.m.setParam(GRB.Param.PoolSolutions, 1)
        self.m.setParam(GRB.Param.PoolSearchMode, 2)
        self.m.setParam(GRB.Param.MIPGap, 0)

        start_time = time.time()
        self.m.optimize()
        self.total_time = time.time() - start_time

        # self.nActualSolutions = self.m.SolCount
        self.n_comigrations = self.compute_summary()


    def add_ILP_vars(self):
        """
        sets up ILP variables
        """
        # self.l = self.m.addMVar((self.clone_tree.n_nodes, self.clone_tree.n_sites), vtype=GRB.BINARY, name="l")
        # self.g = self.m.addMVar((self.clone_tree.n_edges, self.clone_tree.n_sites, self.clone_tree.n_sites, self.clone_tree.n_edges), \
        #     vtype=GRB.BINARY, name="g")
        # self.pi = self.m.addMVar((self.clone_tree.n_edges, self.clone_tree.n_sites, self.clone_tree.n_sites), vtype=GRB.BINARY, name="pi")
        self.l = self.m.addMVar((self.clone_tree.n_migs, self.clone_tree.n_migs), vtype=GRB.BINARY, name="l")
        self.pi = self.m.addMVar((self.clone_tree.n_migs, self.clone_tree.n_sites, self.clone_tree.n_sites), 
                                 vtype=GRB.BINARY, name="pi")
        
    def add_timestamp_constraints(self):
        for uv in self.clone_tree.mig_edges:
            sum = 0
            for e in range(self.clone_tree.n_migs):
                sum += self.l[self.mig_index[uv], e]
            self.m.addConstr(sum == 1)
        for E in range(self.clone_tree.n_migs):
            for uv,uv2 in self.X:
                sum1 = 0
                sum2 = 0
                for e in range(E):
                    sum1 += self.l[self.mig_index[uv] , e]
                    sum2 += self.l[self.mig_index[uv2] , e]
                self.m.addConstr(sum1 >= sum2)

    def add_comigration_constraints(self):
        for e in range(self.clone_tree.n_migs):
            sum = 0
            for s in range(self.clone_tree.n_sites):
                self.m.addConstr(self.pi[e,s,s]==0)
                for t in range(self.clone_tree.n_sites):
                    if s != t:
                        sum += self.pi[e, s, t]
            self.m.addConstr(sum <= 1)
        
        for u, v in self.clone_tree.mig_edges:
            for e in range(self.clone_tree.n_migs):
                self.m.addConstr(self.pi[e, self.site_index[self.clone_tree.get_label(u)], self.site_index[self.clone_tree.get_label(v)]] >= 
                    self.l[self.mig_index[(u,v)], e])
                
    def set_additional_constraints(self):
        for e in range(self.clone_tree.n_migs - 1):
            sum1 = 0
            sum2 = 0
            for s in range(self.clone_tree.n_sites):
                for t in range(self.clone_tree.n_sites):
                    sum1 += self.pi[e,s,t]
                    sum2 += self.pi[e+1,s,t]
            self.m.addConstr(sum1 >= sum2)

    def set_optimization_function(self):
        """ """
        sum = 0
        for s in range(self.clone_tree.n_sites):
            for t in range(self.clone_tree.n_sites):
                if s != t:
                    for e in range(self.clone_tree.n_migs):
                        sum += self.pi[ e,s,t ]

        self.m.setObjective(sum, GRB.MINIMIZE)


    def compute_summary(self):
        """

        :param soln_ind:  (Default value = None)

        """
        sum = 0
        for s in range(self.clone_tree.n_sites):
            for t in range(self.clone_tree.n_sites):
                if s != t:
                    for e in range(self.clone_tree.n_migs):
                        sum += self.pi[ e,s,t ].X
        return round(sum)


    def infer_tree(self):
        """

        :param soln: 

        """
        # self.m.setParam(GRB.Param.SolutionNumber, soln)
        edges = self.clone_tree.edges
        labels = dict()

        timestamps = {}
        l = []
        ee = 1
        for tt in range(self.clone_tree.n_migs):
            for e in range(self.clone_tree.n_migs):
                if self.l[ tt, e ].X > 0.5:
                    timestamps[ self.clone_tree.mig_edges[tt] ] = e
                    l.append((self.clone_tree.mig_edges[tt], e))
        # print(l)
        pi = []
        for e in range(self.clone_tree.n_migs):
            for s in range(self.clone_tree.n_sites):
                for t in range(self.clone_tree.n_sites):
                    if self.pi[e,s,t].X>0.5:
                        pi.append((e, self.clone_tree.sites[s], self.clone_tree.sites[t]))
        # print(pi)

        return timestamps
    

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
    clone_tree = CloneTree(clone_tree_filename=args.phylogeny, vertex_labeling_filename=args.vertex_labeling)

    if args.output is None:
        output_str = '.'
    else:
        output_str = args.output
        if not os.path.exists(output_str):
            os.makedirs(output_str)

    if args.log == True:
        logfile = f'{output_str}/{primary_str}-log.txt'
    else:
        logfile = ''

    soln = ILPSolver(clone_tree, logfile=logfile, n_threads=args.threads)


    colormap = process_colormap(colormap_filename=args.colormap, clone_tree=clone_tree)

    clone_tree.set_timestamps(soln.infer_tree())
    # primary_str = T_prime.get_label(T_prime.root)
    clone_tree.write_dot(f'{output_str}/T.dot', colormap=colormap)
    # T_prime.write_tree(f'{output_str}/T-{primary_str}.tree')
    clone_tree.write_labeling(f'{output_str}/T.labeling')
    summary = soln.compute_summary()
    print(f'{summary}\tOptimal\t\t{soln.total_time}')
