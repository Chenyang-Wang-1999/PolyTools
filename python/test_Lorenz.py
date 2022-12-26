'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2022-12-23
Copyright © Department of Physics, Tsinghua University.  All rights reserved
'''

import py_invariant_manifold as pyIM
import sympy as sym
import numpy as np
from scipy import linalg as LA
import matplotlib.pyplot as plt

class LorenzSys(pyIM.PolynomialSystemSolver):
    def __init__(self, Kmax, sigma, rho, beta, manifold_dim):
        x,y,z = sym.symbols('x,y,z')
        polys = [
            sigma*(y - x),
            rho*x - y - x*z, 
            x*y - beta*z
        ]
        super().__init__(polys, [x,y,z], manifold_dim, Kmax)


def calculate_eigenvalues():
    Kmax = 100
    sigma = 10
    rho = 28
    beta = 8/3

    lorenz_sys = LorenzSys(Kmax, 10, 28, 8/3)
    lorenz_sys.set_F_around_x([0,0,0])
    J = lorenz_sys.jacobian_by_F()
    print(J)
    eigv, eigvec = LA.eig(J)
    order = np.argsort(eigv.real)
    eigv = eigv[order]
    eigvec = eigvec[:,order]
    print(eigv)
    print(eigvec)

def calculate_stable_manifold():
    Kmax = 200
    sigma = 10
    rho = 28
    beta = 8/3

    lorenz_sys = LorenzSys(Kmax, 10, 28, 8/3, 2)
    lorenz_sys.set_F_around_x([0,0,0])

    # eigen values
    J = lorenz_sys.jacobian_by_F()
    eigv, eigvec = LA.eig(J)
    order = np.argsort(eigv.real)
    eigv = eigv[order]
    eigvec = eigvec[:,order]   
    lorenz_sys.init_without_T(eigvec, eigv)

    # solve
    lorenz_sys.set_Kceil(100)
    lorenz_sys.solve_to_Kceil(False)

    # retrieve data 
    conv_range = lorenz_sys.estimite_conv_range(20, tol = 1e-5)
    conv_range = min([conv_range, 20])

    s1 = np.linspace(-conv_range, conv_range, 100)
    s2 = np.linspace(-conv_range, conv_range, 100)
    s1_mesh, s2_mesh = np.meshgrid(s1, s2)
    np.savetxt('s1_mesh.csv',s1_mesh, delimiter=',')
    np.savetxt('s2_mesh.csv',s2_mesh, delimiter=',')
    s_pts = np.zeros((2, s1_mesh.size))
    s_pts[0,:] = s1_mesh.flatten()
    s_pts[1,:] = s2_mesh.flatten()

    W_pts = lorenz_sys.eval(pyIM.POLY_W, s_pts)
    print(np.max(np.abs(W_pts.imag)))
    np.savetxt('w1_mesh.csv',W_pts[0,:].reshape(s1_mesh.shape), delimiter=',')
    np.savetxt('w2_mesh.csv',W_pts[1,:].reshape(s1_mesh.shape), delimiter=',')
    np.savetxt('w3_mesh.csv',W_pts[2,:].reshape(s1_mesh.shape), delimiter=',')

    lorenz_sys.print_poly_info(pyIM.POLY_W)
    

    # plot
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.plot(W_pts[0,:], W_pts[1,:], W_pts[2,:], '.')
    fig.savefig("Lorenz_stable.png")
    plt.close(fig)

if __name__ == '__main__':
    # calculate_eigenvalues()
    calculate_stable_manifold()