"""Script to quickly test Network's class features.
"""

import numpy as np
from network import Network
from solver import Solver

def network_check(net, verbose=False):
    # ----------------------------
    # Check that all paths are shorter than maximum path length
    for pl in net.paths_links:
        assert(
            sum(net.c[l] for l in pl) <= net.tau_max
        ), 'Path_link {pl} is longer than tau_max'.format(pl=pl)
    else:
        if verbose:
            print('All paths are within max path length')
    # ----------------------------
    # If shortest path assignment, check that paths between nodes are the shortest ones
    if net.assignment == 'shortest_path':
        for (node, os) in enumerate(net.F):
            for (o, spl) in enumerate(os):
                path_o_node_idxs = list(filter(
                    lambda i: net.paths[i][0] == o and net.paths[i][-1] == node,
                    range(len(net.paths))
                ))
                for p_idx in path_o_node_idxs:
                    assert(
                        sum(net.c[l] for l in net.paths_links[p_idx]) == spl
                    ), 'Path link {pl} is not a shortest path'.format(pl=net.paths_links[p_idx])
        else:
            if verbose:
                print('All paths are shortest paths')
    # ----------------------------
    # Check that matrices P and P_ms are valid, as well as constraints by
    # verifying that all constraints are satified
    solver_ss = Solver(net)
    # ----------------------------
    #   Single-step
    solver_ss.get_single_step_constraints()
    # Build link count vector: times a flow is counted over a link in the
    # shortest path case. Otherwise set everything to 1 since random accepts
    # only rigid costs

    for (o_i,o) in enumerate(net.origins):
        c3_cons = np.sum(np.divide(
            net.P[:,o_i][solver_ss.c3[:,o_i]],
            net.lc[:,o_i][solver_ss.c3[:,o_i]]
            ))
        if not np.isclose(c3_cons, 1.0):
            print('FAILED: SS C3 observability constraint NOT satisifed for origin {o}'.format(o=o))
            exit()
    else:
        if verbose:
            print('SS C3 observability constraint satisifed')
    
    if not (
            np.all(net.P[solver_ss.c4] >= 0) and 
            np.all(net.P[solver_ss.c4] <= net.lc[solver_ss.c4]) and
            np.all(net.P[np.logical_not(solver_ss.c4)] == 0)
        ):
        check1 = net.P[solver_ss.c4] - net.lc[solver_ss.c4]
        check2 = net.P[np.logical_not(solver_ss.c4)]
        if max(check1) > 1e-10 or max(np.abs(check2)) > 1e-10:
            print('FAILED: SS C4 speed constraint NOT satisifed')
            exit()
    elif verbose:
        print('SS C4 speed constraint satisifed')

    for (n_i,_) in enumerate(net.nodes):
        rows_in = solver_ss.c5_in_edges[:,n_i]
        cols_in = solver_ss.c5_in_check[:,n_i]
        rows_out = solver_ss.c5_out_edges[:,n_i]
        cols_out = solver_ss.c5_out_check[:,n_i]
        inflow_n = np.sum(np.divide(
            net.P[rows_in,:][:,cols_in],
            net.lc[rows_in,:][:,cols_in]
            ), 0
        )
        outflow_n = np.sum(np.divide(
            net.P[rows_out,:][:,cols_out],
            net.lc[rows_out,:][:,cols_out]
            ), 0
        )
        check = inflow_n - outflow_n
        if not np.all(np.greater_equal(inflow_n, outflow_n)):
            if min(check) < -1e-10:
                print('FAILED: SS C5 flow constraint NOT satisifed')
                exit()
    else:
        if verbose:
            print('SS C5 flow constraint satisifed')
    # ----------------------------
    #   Multi-step
    solver_ms = Solver(net)
    solver_ms.get_multi_step_constraints()
    # c_arr_ms = np.ceil(np.array(net.c))
    P_ms = np.concatenate(net.P_ms, axis=-1)
    for (o_i,o) in enumerate(net.origins):
        c3_cons = np.sum(P_ms[:,o_i][solver_ms.c3[:,o_i]])
        if not np.isclose(c3_cons,1.0):
            print('FAILED: MS C3 observability constraint NOT satisifed for origin {o}'.format(o=o))
            exit()
    else:
        if verbose:
            print('MS C3 observability constraint satisifed')

    if not (
            np.all(P_ms[solver_ms.c4] >= 0) and 
            np.all(P_ms[solver_ms.c4] <= 1) and
            np.all(P_ms[np.logical_not(solver_ms.c4)] == 0)
        ):
        check1 = P_ms[solver_ms.c4]
        check2 = P_ms[np.logical_not(solver_ms.c4)]
        if max(check1) - 1 > 1e-10 or max(np.abs(check2)) > 1e-10:
            print('FAILED: MS C4 speed constraint NOT satisifed')
            exit()
    elif verbose:
        print('MS C4 speed constraint satisifed')

    for (n_i,_) in enumerate(net.nodes):
        rows_in = solver_ms.c5_in_edges[:,n_i]
        cols_in = solver_ms.c5_in_check[:,n_i]
        rows_out = solver_ms.c5_out_edges[:,n_i]
        cols_out = solver_ms.c5_out_check[:,n_i]
        inflow_n = np.sum(P_ms[rows_in,:][:,cols_in], 0)
        outflow_n = np.sum(P_ms[rows_out,:][:,cols_out], 0)
        check = inflow_n - outflow_n
        if not np.all(inflow_n - outflow_n >= 0):
            if min(check) < -1e-10:
                print('FAILED: MS C5 flow constraint NOT satisifed')
                exit()
    else:
        if verbose:
            print('MS C5 flow constraint satisifed')
    
    if net.assignment == 'shortest_path':
        for (l_i,_) in enumerate(net.links):
            for (o_i,o) in enumerate(net.origins):
                enter_tau = solver_ms.c7_enter_link[l_i,o_i]
                end_tau = solver_ms.c7_end_link[l_i,o_i]
                for tau in range(enter_tau+1, end_tau+1):
                    if not np.isclose(P_ms[l_i,o_i+tau*len(net.origins)], P_ms[l_i,o_i+enter_tau*len(net.origins)]):
                        print('FAILED: MS C7 flow constraint NOT satisifed')
                        exit()
        else:
            if verbose:
                print('MS C7 flow constraint satisifed')

def make_test(uni_bi, h, w, tau_max, costs, assment, num_tests=1000):
    for _ in range(num_tests):
        net = Network(uni_bi, h=h, w=w)
        net.assign_link_costs(costs)
        net.find_all_paths(tau_max=tau_max, assignment=assment)
        net.generate_od_pairs()
        net.compute_path_assignment_matrix()
        net.generate_random_proportions(like_paper=True)
        net.compute_assignment_matrix()
        # net.P_ms = [np.ones(net.P.shape) for _ in range(net.tau_max)]
        network_check(net)
    else:
        print(
            '{uni_bi}directional, {h} by {w}, tau_max:{tau_max}, {costs} link costs, {assment} assignment:\n\tAll tests passed successfully\n'.format(
                uni_bi=uni_bi,
                h=h,
                w=w,
                tau_max=tau_max,
                costs=costs,
                assment=assment
            ))

if __name__ == '__main__':
    # Change parameters to test different assignments and topologies
    N_TESTS = 1000
    make_test('uni', 3, 3, 4, 'rigid', 'random', N_TESTS)
    make_test('uni', 3, 3, 4, 'rigid', 'shortest_path', N_TESTS)
    make_test('uni', 3, 3, 4, 'real', 'shortest_path', N_TESTS)
    make_test('uni', 3, 3, 'mpl', 'real', 'shortest_path', N_TESTS)
    make_test('bi', 3, 3, 4, 'rigid', 'random', N_TESTS)
    make_test('bi', 3, 3, 4, 'rigid', 'shortest_path', N_TESTS)
    make_test('bi', 3, 3, 4, 'real', 'shortest_path', N_TESTS)
    make_test('bi', 3, 3, 'mpl', 'real', 'shortest_path', N_TESTS)
    make_test('bi', 3, 3, 2, 'real', 'shortest_path', N_TESTS)
    make_test('bi', 8, 8, 4, 'real', 'shortest_path', 50)