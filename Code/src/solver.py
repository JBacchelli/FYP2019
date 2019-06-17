"""
Solver class is declared in this file. Notably, constraint matrices are
generated here for single-step, multi-step models, and for random, shortest path
assignment strategies.
"""

import numpy as np

class Solver():
    """
    Class containing methods for obtaining constraint matrices to use in
    Matlab CVX optimization routines. The constraints that can be obtained are:

    for all assignment types:
        (c2: probability constraint)
            Fully included in speed constraint
        c3: observability constraint
            The sum of flows out of each origin during first step
            should sum up to 1
        c4: speed constraint
            The speed constraint includes reachability and for 
            shortest path assignment only, allowed time steps 
            for each origin/link pair
        c5: flow constraint, which is made of
            - c5_in_edges: indicates whether link flows in a node
            - c5_out_edges: indicates whether link flows out of node
            - c5_in_check: indicates whether node is distinct from origin
                and for multi-step model whether time-step is valid
            - c5_out_check: indicates whether node is distinct from origin
                and for multi-step model whether time-step is valid
            The total O-flow that enters a node i (!=o) must be
            greater or equal than the total O-flow that leaves 
            that same node. One constraint per node in the network
    for multi-step shortest path assignment only:
        (c6: shortest path constraint)
            Fully included in speed constraint
        c7: duplicate counts constraint, which is made of
            - c7_enter_link: indicates O-flows that enter link at
                that time step
            - c7_end_link: indicates O-flows that are still in link at
                a time step later than that in which they entered it
    """

    net = None

    def __init__(self, net):
        """Constructor of Solver class
        
        Arguments:
            net {Network} -- Network instance, where assignment matrices and
                link costs have already been set.
        """
        self.net = net

    def get_single_step_constraints(self):
        """
        Generate constraint matrices for single-step model of network.
        """
        dims = self.net.P.shape

        #C4
        net2 = self.net.duplicate_network()
        # Generate positive flows through each link where flows are allowed
        net2.generate_random_proportions(1,2)
        net2.compute_assignment_matrix()
        if np.array_equal(net2.P, self.net.P):
            print('Assignment matrix P is the same')
        self.c4 = (net2.P > 0)
        #C3
        self.c3 = np.full(dims, False)
        for (l_i,l) in enumerate(self.net.links):
            for (o_i,o) in enumerate(self.net.origins):
                if l[0] == o:# and self.c4[l_i,o_i]:
                    self.c3[l_i,o_i] = True
        #C5
        if len(self.net.origins) <  len(self.net.nodes) and self.net.verbose:
            print('\n\nWARNING: Less origins than nodes in the network, check this is expected\n\n')

        self.c5_in_edges = np.full(
            (len(self.net.links), len(self.net.nodes)),
            False
        )
        self.c5_out_edges = np.full(
            (len(self.net.links), len(self.net.nodes)),
            False
        )
        self.c5_in_check = np.full(
            (len(self.net.origins), len(self.net.nodes)),
            False
        )
        self.c5_out_check = np.full(
            (len(self.net.origins), len(self.net.nodes)),
            False
        )
        for (l_i, l) in enumerate(self.net.links):
            for (n_i,n) in enumerate(self.net.nodes):
                if l[1] == n:
                    self.c5_in_edges[l_i, n_i] = True
                elif l[0] == n:
                    self.c5_out_edges[l_i, n_i] = True
        for (o_i,o) in enumerate(self.net.origins):
            for (n_i,n) in enumerate(self.net.nodes):
                if o != n:
                    self.c5_in_check[o_i, n_i] = True
                    self.c5_out_check[o_i, n_i] = True

        #C6 is expressed by c4 in the single-step model
        #C7 cannot be enforced in the single-step model

    def get_multi_step_constraints(self):
        """
        Generates constraints for multi-step model of the network.
        """
        P = np.stack(self.net.P_ms)

        #C3
        self.c3 = np.full(P.shape, False)
        for (l_i,l) in enumerate(self.net.links):
            for (o_i,o) in enumerate(self.net.origins):
                if l[0] == o:
                    # We only need to set this for P[0], constraints 4 and 6
                    # will take care of later time steps
                    self.c3[0,l_i,o_i] = True

        #C4
        net2 = self.net.duplicate_network()
        # Generate positive flows through each link where flows are allowed
        net2.generate_random_proportions(1,2)
        net2.compute_assignment_matrix()
        if np.array_equal(net2.P_ms, self.net.P_ms):
            print('Assignment matrix P_ms is the same')
        P2 = np.stack(net2.P_ms)
        self.c4 = (P2 > 0)

        # Initialising constraint 5 matrices
        if len(self.net.origins) <  len(self.net.nodes) and self.net.verbose:
            print('\n\nWARNING: Less origins than nodes in the network, check this is expected\n\n')
        self.c5_in_edges = np.full(
            (len(self.net.links), len(self.net.nodes)),
            False
        )
        self.c5_out_edges = np.full(
            (len(self.net.links), len(self.net.nodes)),
            False
        )
        self.c5_in_check = np.full(
            (self.net.tau_max, len(self.net.origins), len(self.net.nodes)),
            False
        )
        self.c5_out_check = np.full(
            (self.net.tau_max, len(self.net.origins), len(self.net.nodes)),
            False
        )
        
        if self.net.assignment == 'random':# or self.net.assignment == 'shortest_path':
            #C5
            for (l_i, l) in enumerate(self.net.links):
                for (n_i,n) in enumerate(self.net.nodes):
                    if l[1] == n:
                        self.c5_in_edges[l_i, n_i] = True
                    elif l[0] == n:
                        self.c5_out_edges[l_i, n_i] = True
            for (o_i, o) in enumerate(self.net.origins):
                for (n_i,n) in enumerate(self.net.nodes):
                    if o != n:
                        for tau in range(self.net.tau_max-1):
                            self.c5_in_check[tau, o_i, n_i] = True
                        for tau in range(1, self.net.tau_max):
                            self.c5_out_check[tau, o_i, n_i] = True
        elif self.net.assignment == 'shortest_path':
            #C5
            # For shortest path assignment, only need to take into account 
            # two time steps for flow contraint. The rest is set to 0 by C4.
            for (l_i, l) in enumerate(self.net.links):
                for (n_i,n) in enumerate(self.net.nodes):
                    if l[1] == n:
                        self.c5_in_edges[l_i, n_i] = True
                    elif l[0] == n:
                        self.c5_out_edges[l_i, n_i] = True
            for (o_i, o) in enumerate(self.net.origins):
                for (n_i,n) in enumerate(self.net.nodes):
                    if o != n:
                        arr_tau = self.net.T_minus[n_i,o] - 1
                        dep_tau = self.net.T_plus[n_i,o] - 1
                        if (
                                arr_tau < self.net.tau_max and
                                dep_tau < self.net.tau_max and
                                dep_tau >= 0 and
                                arr_tau >= 0 and
                                arr_tau <= dep_tau
                            ):
                            self.c5_in_check[int(arr_tau), o_i, n_i] = True
                            self.c5_out_check[int(dep_tau), o_i, n_i] = True

            #C7
            # time-step when O-flow enters link
            self.c7_enter_link = np.zeros(
                (len(self.net.links), len(self.net.origins)),
                dtype=int
            )
            # time-step when O-flow reaches end of link
            self.c7_end_link = np.zeros(
                (len(self.net.links), len(self.net.origins)),
                dtype=int
            )
            # Elements of P2 that are greater than 0 should remain constant
            # throughout their traversal of the link, by shortest path
            # assumption
            for (l_i,l) in enumerate(self.net.links):
                for (o_i,o) in enumerate(self.net.origins):
                    # Assume l is link ij
                    # enter_tau is time step when O flow leaves i (enters link ij)
                    enter_tau = self.net.T_plus[l[0],o] - 1
                    # leave_tau is time step when O flow reaches j (reaches end of link ij)
                    leave_tau = self.net.T_minus[l[1],o] - 1
                    if leave_tau < self.net.tau_max - 1 and leave_tau > enter_tau:
                        self.c7_enter_link[l_i,o_i] = enter_tau
                        self.c7_end_link[l_i,o_i] = leave_tau
        
        # Correctly reshaping constraint matrices for usage in Matlab routine
        # P is stacked along dimension 1 in order to obtain a
        # n_l by (n_o*tau_max) 2-dimensional matrix
        P = np.concatenate(P, axis=1)
        self.c3 = np.concatenate(self.c3, axis=1)
        self.c4 = np.concatenate(self.c4, axis=1)
        self.c5_in_check = np.concatenate(self.c5_in_check, axis=0)
        self.c5_out_check = np.concatenate(self.c5_out_check, axis=0)