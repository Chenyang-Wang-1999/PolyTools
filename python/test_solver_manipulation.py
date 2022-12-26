'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2022-12-22
Copyright © Department of Physics, Tsinghua University.  All rights reserved
'''

import py_invariant_manifold as pyIM
import numpy as np
import sympy as sym
from scipy import linalg as LA

def test_jacobian():
    solver = pyIM.InvariantManifoldSolverPy(3, 2, 100)
    solver.add_term_F(0, pyIM.CIndexVec([1,0,0]), 1.0)
    solver.add_term_F(0, pyIM.CIndexVec([0,1,0]), 2.0)
    solver.add_term_F(1, pyIM.CIndexVec([0,0,1]), 1.0)
    solver.add_term_F(1, pyIM.CIndexVec([1,0,0]), -1.0)
    solver.add_term_F(1, pyIM.CIndexVec([2,0,0]), -1.0)
    solver.add_term_F(2, pyIM.CIndexVec([0,0,1]), 1.0)
    solver.add_term_F(2, pyIM.CIndexVec([0,1,0]), 1.0)
    solver.add_term_F(2, pyIM.CIndexVec([1,0,0]), 1.0)

    J = solver.jacobian_by_F()
    print(J)

def test_vec():
    A = [1,2,3,4]
    A_cvec = pyIM.CIndexVec(A)
    A_cvec[2] = 5
    print(A)
    print(A_cvec)

def test_add_terms():
    solver = pyIM.InvariantManifoldSolverPy(3, 2, 100)
    solver.add_term_F(0, pyIM.CIndexVec([1,0,0]), 1.0)
    solver.add_term_F(0, pyIM.CIndexVec([0,1,0]), 2.0)
    solver.add_term_F(1, pyIM.CIndexVec([0,0,1]), 1.0)
    solver.add_term_F(1, pyIM.CIndexVec([1,0,0]), -1.0)
    solver.add_term_F(1, pyIM.CIndexVec([2,0,0]), -1.0)
    solver.add_term_F(2, pyIM.CIndexVec([0,0,1]), 1.0)
    solver.add_term_F(2, pyIM.CIndexVec([0,1,0]), 1.0)
    solver.add_term_F(2, pyIM.CIndexVec([1,0,0]), 1.0)

    J = solver.jacobian_by_F()
    print("J:", J)
    eigv,eigvec = LA.eig(J)
    print("eigv:", eigv)
    print("eigvec:", eigvec)
    solver.init_without_T(eigvec, eigv)

    solver.print_poly_info(pyIM.POLY_F)
    solver.print_poly_info(pyIM.POLY_f)
    solver.print_poly_info(pyIM.POLY_W)
    solver.print_poly_info(pyIM.POLY_s)

def test_get_poly():
    solver = pyIM.InvariantManifoldSolverPy(3, 2, 100)
    solver.add_term_F(0, pyIM.CIndexVec([1,0,0]), 1.0)
    solver.add_term_F(0, pyIM.CIndexVec([0,1,0]), 2.0)
    solver.add_term_F(1, pyIM.CIndexVec([0,0,1]), 1.0)
    solver.add_term_F(1, pyIM.CIndexVec([1,0,0]), -1.0)
    solver.add_term_F(1, pyIM.CIndexVec([2,0,0]), -1.0)
    solver.add_term_F(2, pyIM.CIndexVec([0,0,1]), 1.0)
    solver.add_term_F(2, pyIM.CIndexVec([0,1,0]), 1.0)
    solver.add_term_F(2, pyIM.CIndexVec([1,0,0]), 1.0)
 
    J = solver.jacobian_by_F()
    print("J:", J)
    eigv,eigvec = LA.eig(J)
    print("eigv:", eigv)
    print("eigvec:", eigvec)
    solver.init_without_T(eigvec, eigv)   

    from_k = 0
    # F_data = solver.get_poly_data(pyIM.POLY_W, from_k)
    x,y,z = sym.symbols('x,y,z')
    F_poly = solver.get_poly_sympy(pyIM.POLY_F,[x,y,z],from_k=1 )
    for scalar_data in F_poly:
        print(scalar_data)

    solver.print_poly_info(pyIM.POLY_F)
    print(solver.get_poly_var_dim(pyIM.POLY_F))
    print(solver.get_poly_val_dim(pyIM.POLY_F))
    print(solver.get_poly_var_dim(pyIM.POLY_f))
    print(solver.get_poly_val_dim(pyIM.POLY_f))
    print(solver.get_poly_var_dim(pyIM.POLY_W))
    print(solver.get_poly_val_dim(pyIM.POLY_W))
    print(solver.get_poly_var_dim(pyIM.POLY_s))
    print(solver.get_poly_val_dim(pyIM.POLY_s))



if __name__ == '__main__':
    # test_jacobian()
    # test_add_terms()
    # test_vec()
    test_get_poly()