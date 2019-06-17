"""
Script containing useful methods to export data produced in Python in
Matlab readable formats, so that it can safely be loaded and used with the CVX
optimisation routine.
"""

import os
import numpy as np
import scipy.io as io

from network import Network
from solver import Solver

class ToMatlab:
    """
    Class for converting network, assignment and flow structures generated in
    Python, to Matlab matrices so that they can be used for optimization problem
    """

    def __init__(
            self,
            net,
            step='multi',
            trials=5,
            m_dir='../../OflowEstimationFull/tests/'
        ):
        self.net = net
        self.step = step
        self.trials = trials
        self.m_dir = m_dir

    def save_assignment(self):
        sp = self.net.assignment == 'shortest_path'
        io.savemat(self.m_dir+'shortest_path', {'shortest_path':sp})

    def save_trials(self):
        io.savemat(self.m_dir+'trials', {'trials':self.trials})

    def save_tau_max(self):
        io.savemat(self.m_dir+'tau_max', {'tau_max':self.net.tau_max})

    def convert_P(self, ext):
        if self.step == 'multi':
            # 4-dimensional array
            P_py = np.zeros((
                len(self.net.links),    # links
                len(self.net.origins),  # origins
                self.net.tau_max,       # time steps
                self.trials             # trials
            ))
            for tr in range(self.trials):
                self.net.generate_random_proportions()
                self.net.compute_assignment_matrix()
                P_py[:,:,:,tr] = np.dstack(self.net.P_ms)
                # print('P'+ext)
                # print(P_py[:,:,:,0])

        if self.step == 'single':
            # 3-dimensional array
            P_py = np.zeros((
                len(self.net.links),    # links
                len(self.net.origins),  # origins
                self.trials             # trials
            ))
            for tr in range(self.trials):
                self.net.generate_random_proportions()
                self.net.compute_assignment_matrix()
                P_py[:,:,tr] = self.net.P

        io.savemat(self.m_dir+'P_'+ext, {'P_'+ext:P_py})

    def convert_o_list(self):
        o_list_py = np.array(self.net.origins) + 1
        io.savemat(self.m_dir+'o_list', {'o_list':o_list_py})

    def convert_e_list(self):
        e_list_py = np.array(self.net.links) + 1
        io.savemat(self.m_dir+'e_list', {'e_list':e_list_py})

    def convert_od_list(self):
        od_list_py = np.array(self.net.od_pairs) + 1
        io.savemat(self.m_dir+'od_list', {'od_list':od_list_py})

    def convert_lc(self):
        lc_py = self.net.lc
        io.savemat(self.m_dir+'lc', {'lc':lc_py})

    def convert_constraints(self):
        solver = Solver(self.net)
        if self.step == 'multi':
            solver.get_multi_step_constraints()
            if self.net.assignment == 'shortest_path':
                io.savemat(
                    self.m_dir+'constraints',
                    {
                        'c3':solver.c3,
                        'c4':solver.c4,
                        'c5_in_edges':solver.c5_in_edges,
                        'c5_in_check':solver.c5_in_check,
                        'c5_out_edges':solver.c5_out_edges,
                        'c5_out_check':solver.c5_out_check,
                        'c7_enter_link':solver.c7_enter_link,
                        'c7_end_link':solver.c7_end_link
                    }
                )
            elif self.net.assignment == 'random':
                io.savemat(
                    self.m_dir+'constraints',
                    {
                        'c3':solver.c3,
                        'c4':solver.c4,
                        'c5_in_edges':solver.c5_in_edges,
                        'c5_in_check':solver.c5_in_check,
                        'c5_out_edges':solver.c5_out_edges,
                        'c5_out_check':solver.c5_out_check
                    }
                )
        elif self.step == 'single':
            solver.get_single_step_constraints()
            io.savemat(
                    self.m_dir+'constraints',
                    {
                        'c3':solver.c3,
                        'c4':solver.c4,
                        # 'c5_inflows':solver.c5_inflows,
                        # 'c5_outflows':solver.c5_outflows,
                        'c5_in_edges':solver.c5_in_edges,
                        'c5_in_check':solver.c5_in_check,
                        'c5_out_edges':solver.c5_out_edges,
                        'c5_out_check':solver.c5_out_check
                    }
                )

    def convert_data(self):
        
        # One will be P_target, the second will be P_initialise
        self.save_assignment()
        self.save_trials()
        self.save_tau_max()
        self.convert_P('target')
        self.convert_P('initialise')
        self.convert_o_list()
        self.convert_e_list()
        self.convert_od_list()
        self.convert_lc()
        self.convert_constraints()

if __name__ == '__main__':
    # SETTING UP EXPERIMENT
    STEP = 'multi'
    COSTS = 'real'
    ASSMENT = 'shortest_path'
    H = 3
    W = 3
    TAU_MAX = 4
    TRIALS = 3
    DIRECTION = 'bi'

    DATA = 'histogram' # 'histogram' or 'plot'
    PARAM = 'n_T' # n_T or 'size'

    # NO NEED TO MODIFY AFTER THIS POINT
    if DATA == 'histogram':
        m_dir = '../../OflowEstimationFull/tests/{step}_{costs}_{assment}_{h}by{w}_{tau_max}taumax_{trials}trials/'.format(
        # m_dir = '../../temp/tests/RANDCONSC4C5_{step}_{costs}_{assment}_{h}by{w}_{tau_max}taumax_{trials}trials/'.format(
            step=STEP,
            costs=COSTS,
            assment=ASSMENT,
            h=H,
            w=W,
            tau_max=TAU_MAX,
            trials=TRIALS
        )
        if not os.path.exists(m_dir):
            os.mkdir(m_dir)
        else:
            print('WARNING: test data in {dir} will be overwritten'.format(dir=m_dir))
            ans = input('Press any key to continue, or \'n\' to stop execution')
            if ans == 'n':
                print('Stopping execution.')
                exit()
            else:
                print('Overwriting files.')

        net = Network(DIRECTION, h=H, w=W)
        net.assign_link_costs(COSTS)
        net.find_all_paths(tau_max=TAU_MAX, assignment=ASSMENT)
        net.generate_od_pairs()
        net.compute_path_assignment_matrix()
        conv = ToMatlab(
            net, step=STEP, trials=TRIALS,
            m_dir=m_dir
        )
        conv.convert_data()
    elif DATA == 'plot':
        m_dir = '../../OflowEstimationFull/tests/{param}_{step}_{costs}_{assment}_{h}by{w}_{tau_max}taumax_{trials}trials/'.format(
            param=PARAM,
            step=STEP,
            costs=COSTS,
            assment=ASSMENT,
            h=H,
            w=W,
            tau_max=TAU_MAX,
            trials=1
        )
        if not os.path.exists(m_dir):
            os.mkdir(m_dir)
        else:
            print('WARNING: test data in {dir} will be overwritten'.format(dir=m_dir))
            ans = input('Press any key to continue, or \'n\' to stop execution')
            if ans == 'n':
                print('Stopping execution.')
                exit()
            else:
                print('Overwriting files.')
        if PARAM == 'n_T':
            TRIALS = 1
            net = Network(DIRECTION, h=H, w=W)
            net.assign_link_costs(COSTS)
            net.find_all_paths(tau_max=TAU_MAX, assignment=ASSMENT)
            net.generate_od_pairs()
            net.compute_path_assignment_matrix()
            conv = ToMatlab(
                net, step=STEP, trials=1,
                m_dir=m_dir
            )
            conv.convert_data()
        elif PARAM == 'size':
            TRIALS=1
            #size: 4, 6, 9, 12, 16,  
            hw_list = [(2,2), (2,3), (3,3), (3,4), (4,4), (4,5), (5,5), (5,6), (6,6), (7,7), (8,8)]
            for i in range(len(hw_list)):
                (H,W) = hw_list[i]
                net = Network(DIRECTION, h=H, w=W)
                net.assign_link_costs(COSTS)
                net.find_all_paths(tau_max=TAU_MAX, assignment=ASSMENT)
                net.generate_od_pairs()
                net.compute_path_assignment_matrix()
                m_dir = '../../OflowEstimation/tests/{param}_{step}_{costs}_{assment}_{h}by{w}_{tau_max}taumax_{trials}trials/'.format(
                    param=PARAM,
                    step=STEP,
                    costs=COSTS,
                    assment=ASSMENT,
                    h=H,
                    w=W,
                    tau_max=TAU_MAX,
                    trials=1
                )
                if not os.path.exists(m_dir):
                    os.mkdir(m_dir)
                conv = ToMatlab(
                    net, step=STEP, trials=1,
                    m_dir=m_dir
                )
                conv.convert_data()