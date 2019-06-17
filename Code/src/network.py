"""
Network class is declared in this file. All meaningful elements of simulations
are effectively derived from an instance of a network.
"""

VISUALISE = False   # Set to True if working with Jupyter notebooks 
                    # in order to visualise networks and assignment matrices
import numpy as np
import networkx as nx
from copy import deepcopy
if VISUALISE:
    from visualiser import Vis
np.set_printoptions(linewidth=150)

def mat_conv(M, v):
    """Performs convolution of assignment matrix with flow vector and returns
    result. len(M) zeros are prepended to the flow vector to obtain a result
    with len(v) elements. 
    
    Arguments:
        M {list} -- Multi-step assignemnt matrix, from 0 to tau_max-1.
            List of 2-dimensional arrays. Each matrix in the list
            corresponds to one time step.
        v {list} -- Flow vector, from 0 to n_t-1 for full convolution.
            List of 1-dimensional arrays. Each flow vector in the list
            corresponds to one time sample.
    
    Returns:
        {list} -- List of flow vectors, from 0 to n_t-1. Each element in the
            list corresponds to one time sample.
    """
    n_t = len(v)        # the observation time horizon
    tau_max = len(M)    # the maximum path length
    full_flows = [np.zeros(len(v[0])) for _ in range(tau_max-1)] + v
    res = [np.zeros(M[0].shape[0])] * n_t
    for t in range(n_t):                    # 0 to n_T-1
        for tau in range(tau_max):          # 0 to tau_max-1
            res[t] = res[t] + M[tau] @ full_flows[t-tau+tau_max-1]
    return res

class Network():
    """
    Representation of a network as a set of nodes and links.
    A NetworkX graph object is used to leverage all the built-in functionalities
    and algorithms. https://networkx.github.io/
    """

    def __init__(self, uni_bi='bi', h=3, w=3, nodes=None, links=None, seed=None, verbose=False):
        """Builds network and checks validity of links with respect to nodes
        
        Arguments:
            uni_bi {str} -- Type of network to generate automatically, either
                'uni' for unidirectional networks or 'bi' for bidirectional
                networks. Set to None for passing nodes and links
                manually (default: {'bi'})
            h {int} -- Height of network (default: {3})
            w {int} -- Width of network (default: {3})
            nodes {int list} -- Pass list of nodes in the network manually
                (default: {None})
            links {(int*int) list} -- Pass list of links in the network manually
                (default: {None})
            seed {int} -- Random seed used throughout random generations 
                (default: {None})
        """
        self.G = nx.DiGraph()
        if uni_bi:
            grid = nx.grid_2d_graph(h, w)
            self.G.update(grid)
            self.G = nx.convert_node_labels_to_integers(self.G)
            if uni_bi == 'bi':
                rev_edges = []
                for (n1,n2) in self.G.edges:
                    rev_edges.append((n2, n1))
                self.G.add_edges_from(rev_edges)
        else:
            self.G.add_nodes_from(nodes)    
            self.G.add_edges_from(links)

        self.nodes = list(self.G.nodes)
        self.links = list(self.G.edges)
        if seed:
            np.random.seed(seed)
        else:
            np.random.seed()

        # Verify links contain valid nodes  only
        for link in self.links:
            if link[0] not in self.nodes:
                raise ValueError('Node {n} was not declared in the list of nodes'.format(n=link[0]))
            elif link[1] not in self.nodes:
                raise ValueError('Node {n} was not declared in the list of nodes'.format(n=link[1]))

        self.verbose = verbose

    def assign_link_costs(self, costs='rigid'):
        """Assigns costs to all links, generating them if necessary
        
        Keyword Arguments:
            costs {str or list} -- Either string to say which costs to use, or
            list of costs for all links.
                'rigid' assigns a cost of 1 to each node
                'mult_int' assigns random non-zero positive integers
                'real' assigns random non-zero real valued numbers
                (default: {'rigid'})
        """
        self.costs = costs
        if costs == 'rigid':
            nx.set_edge_attributes(self.G, 1, 'cost')
        elif costs == 'mult_int':
            rand_costs = np.random.randint(1,5,len(self.G.edges))
            c_dict = {e:rand_costs[i]
                for (i,e) in enumerate(self.G.edges)
            }
            nx.set_edge_attributes(self.G, c_dict, 'cost')
        elif costs == 'real':
            rand_costs = np.random.rand(len(self.G.edges))*4
            # Do not want any link cost to be zero
            while np.any(np.isclose(rand_costs, 0)):
                rand_costs = np.random.rand(len(self.G.edges))*4
            c_dict = {e:rand_costs[i]
                for (i,e) in enumerate(self.G.edges)
            }
            nx.set_edge_attributes(self.G, c_dict, 'cost')
        elif type(costs) == list:
            c_dict = {e:costs[i]
                for (i,e) in enumerate(self.G.edges)
            }
            nx.set_edge_attributes(self.G, c_dict, 'cost')
        else:
            raise ValueError('argument to `costs` not recognised')
        self.c = [cost for (_,_,cost) in self.G.edges.data('cost')]

    def compute_spl_matrix(self):
        """Computes shortest path length matrix F, which gives minimal cost to
        reach any node i from any origin o.
        """
        # Initialise to infinity for all nodes, if a node is not reachable
        # f will remain infinity
        self.F = np.ones((len(self.nodes), len(self.nodes))) * np.inf
        for spl in nx.shortest_path_length(self.G, weight='cost'):
            for node in spl[1].keys():
                self.F[node, spl[0]] = spl[1][node]

    def find_all_paths(self, tau_max=4, assignment='random'):
        """Finds all feasible loop-free paths in given network, which are
        shorter than tau_max. For this, shortest path matrix F is also computed.
        
        Keyword Arguments:
            tau_max {int or str} -- Maximum path length. If set to
                'mpl', tau_max is set to the length of the
                longest shortest path (default: {4})
            assignment {str} -- Assignment, either 'random' or 'shortest_path'.
                Both assignments are loop-free. (default: {'random'})
        """
        self.compute_spl_matrix()
        self.tau_max = tau_max
        if self.tau_max == 'mpl':
            self.tau_max = int(np.ceil((np.max(self.F[self.F < np.inf]))))
        if min(self.c) >= self.tau_max:
            print('No paths given current link costs, reassigning costs')
            self.assign_link_costs(costs=self.costs)
        self.assignment = assignment
        self.paths = []
        if self.assignment == 'shortest_path':
            for n1 in self.G.nodes:
                for n2 in self.G.nodes:
                    if n1 != n2:
                        try:
                            s_paths = nx.algorithms.shortest_paths.generic.all_shortest_paths(
                                self.G, n1, n2, weight='cost'
                            )
                            for p in s_paths:
                                self.paths.append(p)
                        except nx.exception.NetworkXNoPath:
                            continue
        elif self.assignment == 'random':
            for n1 in self.G.nodes:
                for n2 in self.G.nodes:
                    if n1 != n2:
                        try:
                            s_paths = nx.algorithms.simple_paths.all_simple_paths(
                                self.G, n1, n2, cutoff=self.tau_max
                            )
                            for p in s_paths:
                                self.paths.append(p)
                        except nx.exception.NetworkXNoPath:
                            continue
        else:
            raise ValueError('Assignment \'{a}\' is not recognised'.format(a=assignment))
        # Store paths as sequence of link indexes as well
        self.paths_links = []
        for path in self.paths:
            tmp_plink = []
            for n_i in range(len(path)-1):
                tmp_plink.append(self.links.index((path[n_i],path[n_i+1])))
            self.paths_links.append(tmp_plink)

        # Remove paths that exceed the maximum path length
        if self.assignment == 'shortest_path':
            # Maximum path length allowed
            mpl = min(np.amax(self.F), self.tau_max)
        else:
            mpl = self.tau_max
        pl_ps = sorted(list(filter(
            lambda pl_p: sum(self.c[l] for l in pl_p[0]) <= mpl,
            zip(self.paths_links, self.paths)
        )))
        self.paths_links, self.paths = [list(t) for t in zip(*pl_ps)]

        # Reduce tau_max if it is greater than longest path
        if self.assignment == 'random':
            mpl = max(map(
                lambda pl: sum(self.c[l] for l in pl),
                self.paths_links
            ))
        if self.tau_max > int(np.ceil(mpl)):
            self.tau_max = int(np.ceil(mpl))
            if self.verbose:
                print('Maximum path length tau_max changed to {tm}'.format(
                    tm=self.tau_max)
                )

    def generate_od_pairs(self):
        """
        Generate set of Origin-Destination pairs for the given network.
        """
        assert(self.paths), 'Paths were not determined, or no paths exist.'

        self.od_pairs = sorted(list(
            set(map(lambda p: (p[0],p[-1]), self.paths))
        ))
        self.origins = sorted(list(
            set(map(lambda odp: odp[0], self.od_pairs))
        ))
        self.destinations = sorted(list(
            set(map(lambda odp: odp[-1], self.od_pairs))
        ))

    def compute_path_assignment_matrix(self):
        """
        Computes deterministic path assignment matrix, which depends on
        topology of the network, assignment strategy and type of link costs.
        NB: Random assignment is currently only supported for rigid costs.
        """
        assert(self.paths), 'Paths were not determined, or no paths exist.'
        # Initialise path assignment matrix and arrival/departure matrices
        self.Delta = np.zeros((len(self.links), len(self.paths)))
        self.Delta_ms = [
            np.zeros((len(self.links), len(self.paths)), np.int8)
            for _ in range(self.tau_max)
        ]
        self.T_minus = np.zeros((len(self.nodes), len(self.nodes)))
        self.T_plus = np.zeros((len(self.nodes), len(self.nodes)))
        # Compute
        if self.assignment == 'random':
            assert(all(cst == 1 for cst in self.c)), 'Random assignment only supported for rigid models'
            self.lc = np.ones((len(self.links), len(self.origins)))
            for l_idx in range(len(self.links)):
                for p_idx in range(len(self.paths_links)):
                    for tau in range(self.tau_max):
                        if tau < len(self.paths_links[p_idx]):
                            self.Delta_ms[tau][l_idx,p_idx] = (
                                l_idx == self.paths_links[p_idx][tau]
                            )
        elif self.assignment == 'shortest_path':
            self.T_minus = np.ceil(self.F)
            #self.T_minus = self.T_minus.astype(int)
            self.T_plus = np.floor(self.F + 1)
            #self.T_plus = self.T_plus.astype(int)
            self.lc = np.zeros((len(self.links), len(self.origins)))
            for (l_i,l) in enumerate(self.links):
                for (o_i,o) in enumerate(self.origins):
                    # Difference may be negative, or zero whenever the link is
                    # not used by the origin
                    if (
                        self.T_plus[l[0],o] != np.inf and
                        self.T_plus[l[1],o] != np.inf
                    ):
                        self.lc[l_i,o_i] = max(
                            self.T_minus[l[1],o] - self.T_plus[l[0],o] + 1,
                            1
                        )
                    else:
                        self.lc[l_i,o_i] = 1


            for (p_idx,pl) in enumerate(self.paths_links):
                o = self.links[pl[0]][0]
                for l_idx in pl:
                    i = self.links[l_idx][0]
                    j = self.links[l_idx][1]
                    for tau in range(self.tau_max):
                        if (
                            tau+1 >= self.T_plus[i,o] and
                            tau+1 <= self.T_minus[j,o] and
                            self.T_minus[j,o] != np.inf
                            ):
                            self.Delta_ms[tau][l_idx,p_idx] = 1
        # Obtain single-step path assignment matrix as sum of the steps
        self.Delta = sum(d for d in self.Delta_ms)

    def generate_random_proportions(
        self,
        like_paper=True,
        min_prop=0,
        max_prop=100
        ):
        """Generates random fixed proportions O to OD flows and OD to path flows.
        
        Keyword Arguments:
            like_paper {bool} -- Whether to generate proportions with same mean
                and variance as paper (default: {True})
            min_prop {int} -- Minimum flow proportion (default: {0})
            max_prop {int} -- Maximum flow proportion (default: {100})
        """

        assert(min_prop >=  0), 'Minimum proportion must be positive'
        assert(max_prop > 0), 'Maximum proportion must be strictly positive'
        assert(min_prop < max_prop+1), 'Maximum proportion must be greater than minimum proportion plus 1'

        self.o_od_proportions = np.zeros((len(self.origins), len(self.od_pairs)))
        self.od_path_proportions = np.zeros((len(self.od_pairs), len(self.paths)))
        if like_paper:
            # Generate path proportions
            path_proportions = 2 + np.random.rand(len(self.paths))
            path_proportions = path_proportions/np.sum(path_proportions)
            od_probs = np.zeros(len(self.od_pairs))
            for od_idx in range(len(self.od_pairs)):
                # Distribution of OD flow over path flows
                p_idxs = list(filter(
                    lambda p_i: 
                        self.paths[p_i][0] == self.od_pairs[od_idx][0] and
                        self.paths[p_i][-1] == self.od_pairs[od_idx][1],
                    range(len(self.paths))
                ))
                od_fractions = [path_proportions[i] for i in p_idxs]
                od_probs[od_idx] = np.sum(od_fractions)
                od_fractions = od_fractions/od_probs[od_idx] # Normalise
                for (i,p_idx) in enumerate(p_idxs):
                    self.od_path_proportions[od_idx,p_idx] = od_fractions[i]
            for o_idx in range(len(self.origins)):
                # Distribution of O flow over OD flows
                od_idxs = list(filter(
                    lambda odp_i: self.od_pairs[odp_i][0] == self.origins[o_idx], 
                    range(len(self.od_pairs))
                ))
                o_fractions = [od_probs[i] for i in od_idxs]
                o_fractions = o_fractions/np.sum(o_fractions) # Normalise
                for (i,od_idx) in enumerate(od_idxs):
                    self.o_od_proportions[o_idx,od_idx] = o_fractions[i]

        else:
            for o_idx in range(len(self.origins)):
                # Distribution of O flow over OD flows
                od_idxs = list(filter(
                    lambda odp_i: self.od_pairs[odp_i][0] == self.origins[o_idx], 
                    range(len(self.od_pairs))
                ))
                o_fractions = np.zeros((1,len(od_idxs)))
                while np.sum(o_fractions) == 0:
                    o_fractions = np.random.randint(
                        min_prop, max_prop,
                        len(od_idxs)
                    )
                # Normalise distribution
                o_fractions = o_fractions/np.sum(o_fractions)
                for (i,od_idx) in enumerate(od_idxs):
                    self.o_od_proportions[o_idx,od_idx] = o_fractions[i]
            for od_idx in range(len(self.od_pairs)):
                # Distribution of OD flow over path flows
                p_idxs = list(filter(
                    lambda p_i: 
                        self.paths[p_i][0] == self.od_pairs[od_idx][0] and
                        self.paths[p_i][-1] == self.od_pairs[od_idx][1],
                    range(len(self.paths))
                ))
                od_fractions = np.zeros((1,len(p_idxs)))
                while np.sum(od_fractions) == 0:
                    od_fractions = np.random.randint(
                        min_prop, max_prop,
                        len(p_idxs)
                    )
                # Normalise distribution
                od_fractions = od_fractions/np.sum(od_fractions)
                for (i,p_idx) in enumerate(p_idxs):
                    self.od_path_proportions[od_idx,p_idx] = od_fractions[i]
    
    def compute_assignment_matrix(self):
        """
        Computes assignment matrix A, and O-flow assignment matrix P based on
        traffic proportions of the network. Assignment matrix will be different
        based on assignment strategy selected when finding all paths, and on
        type of costs selected when generating link costs.
        """
        assert(self.Delta.size is not None), 'Path assignment matrix was not computed.'
        assert(self.o_od_proportions.size is not None), 'O to OD proportions have not been generated.'
        assert(self.od_path_proportions.size is not None), 'OD to path proportions have not been generated.'
        assert(self.assignment is not None), 'Assignment strategy has not been determined.'

        # Initialise traffic assignment matrices
        self.A = np.zeros((len(self.links), len(self.od_pairs)))
        self.P = np.zeros((len(self.links), len(self.origins)))
        self.A_ms = [np.zeros((len(self.links), len(self.od_pairs))) for i in range(self.tau_max)]
        self.P_ms = [np.zeros((len(self.links), len(self.origins))) for i in range(self.tau_max)]

        for (od_i, od) in enumerate(self.od_pairs):
            p_idxs = list(filter(
                lambda p_i: self.paths[p_i][0] == od[0] and self.paths[p_i][-1] == od[1],
                range(len(self.paths))
            ))
            for p_i in p_idxs:
                for tau in range(self.tau_max):
                    self.A_ms[tau][:,od_i] += self.Delta_ms[tau][:,p_i]*self.od_path_proportions[od_i,p_i]
        self.A = sum(a for a in self.A_ms)

        for (o_i, o) in enumerate(self.origins):
            od_idxs = list(filter(
                lambda od_i: self.od_pairs[od_i][0] == o,
                range(len(self.od_pairs))
            ))
            for od_i in od_idxs:
                for tau in range(self.tau_max):
                    self.P_ms[tau][:,o_i] += self.A_ms[tau][:,od_i]*self.o_od_proportions[o_i,od_i]
        self.P = sum(p for p in self.P_ms)

    def duplicate_network(self):
        """
        Returns a copy of the current network with the same topology and
        assignment.
        """
        return deepcopy(self)
