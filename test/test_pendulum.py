'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2022-12-22
Copyright © Department of Physics, Tsinghua University.  All rights reserved
'''

import py_invariant_manifold as pyIM
import numpy as np
from math import sin, cos, sqrt,pi
import matplotlib.pyplot as plt
from scipy import linalg as LA

'''
    \dot{x} = p/m
    \dot{p} = - m \omega^2 sin(x)
'''
class physical_pendulum(pyIM.InvariantManifoldSolverPy):
    omega:float
    m:float
    def __init__(self, Kmax:int, m=1, omega=1):
        self.m = m
        self.omega = omega
        super().__init__(2, 1, Kmax)

    def set_F(self, x0):
        # expand for x = x_0 + u_0 and p = u_1 
        # F_0 = u_1/m
        # F_1 = - m \omega^2 \sin(u_0 + x_0)
        self.add_term_F(0, pyIM.CIndexVec([0, 1]), 1/self.m)

        curr_taylor_coeff = -self.m*(self.omega**2)
        for curr_k in range(1, self.Kmax):
            curr_taylor_coeff /= curr_k
            if(curr_k % 4 == 1):
                diff_val = np.cos(x0)
            elif(curr_k % 4 == 3):
                diff_val = - np.cos(x0)
            
            if(curr_k %2):
                self.add_term_F(1, pyIM.CIndexVec([curr_k, 0]), curr_taylor_coeff * diff_val)

    def set_J_and_P(self, stable:bool=True):
        J = self.jacobian_by_F()
        eigv, eigvec = LA.eig(J)
        order = np.argsort(eigv.real)
        if(not stable):
            order = order[::-1]

        lam = eigv[order]
        P = eigvec[:,order]

        self.init_without_T(P, lam)


def main_test_F():
    solver = physical_pendulum(103)
    solver.set_Kceil(101)
    xp_matrix = np.vstack([np.linspace(-5*pi, 5*pi, 1000).reshape((1,-1)), np.zeros((1,1000))])
    # print(xp_matrix)
    solver.set_F(pi)
    solver.print_poly_info(pyIM.POLY_F)

    solver.set_J_and_P()
    # for j in range(15):
    #     solver.solve_step()
    #     solver.print_poly_info(pyIM.POLY_W)
    #     # solver.print_poly_info(pyIM.POLY_f)
    #     # solver.print_poly_info(pyIM.POLY_Ek)
    #     # curr_k = solver.get_k()
    #     # for err_k in range(curr_k):
    #     #     solver.calculate_err(err_k)
    #     #     solver.print_poly_info(pyIM.POLY_Ek)

    #     input()

    solver.solve_to_Kceil()
    tol = 1e-3

    solver.calculate_err(101)
    err_coeff = solver.eval(pyIM.POLY_Ek, np.array([[1.0]]))
    s_range = (tol/np.max(np.abs(err_coeff))) ** (1/101)

    # solver.print_poly_info(pyIM.POLY_W)
    s_matrix = np.linspace(-s_range,s_range, 1000).reshape((1,-1))
    W_val = solver.eval(pyIM.POLY_W, s_matrix)
    plt.plot(W_val[0,:] + pi, W_val[1,:], 'b-')
    plt.savefig("invariant saddle.png")





if __name__ == '__main__':
    main_test_F()