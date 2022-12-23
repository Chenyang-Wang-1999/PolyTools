'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2022-12-22
Copyright © Department of Physics, Tsinghua University.  All rights reserved
'''

import py_invariant_manifold._invariant_manifold as c_invariant_manifold
# import _invariant_manifold as c_invariant_manifold
import sympy as sym
import numpy as np

# definitions
POLY_F = 0
POLY_f = 1
POLY_W = 2
POLY_s = 3
POLY_Ek = 4

TOL = 1e-15

class CIndexVec(c_invariant_manifold.CIndexVec):
    def __init__(self, vec):
        super().__init__(vec)

class CScalarVec(c_invariant_manifold.CScalarVec):
    def __init__(self, vec):
        super().__init__(vec)

class InvariantManifoldSolverPy:
    c_solver: c_invariant_manifold._CInvariantManifoldSolver
    phys_dim:int 
    manifold_dim:int 
    Kmax:int
    def __init__(self, phys_dim:int, manifold_dim:int, Kmax:int):
        self.c_solver = c_invariant_manifold._CInvariantManifoldSolver(phys_dim, manifold_dim, Kmax)
        self.phys_dim = self.c_solver.phys_dim
        self.manifold_dim = self.c_solver.manifold_dim
        self.Kmax = self.c_solver.Kmax

    def init_without_T(self, P, lam):
        self.c_solver.init_without_T(P, lam)

    def init_with_T(self, P, T, lam):
        self.c_solver.init_with_T(P, T, lam)

    def jacobian_by_F(self):
        return self.c_solver.jacobian_by_F()

    def add_term_F(self, val_id, order_vec, coeff):
        self.c_solver.add_term_F(val_id, order_vec, coeff)

    def get_curr_term(self, order_vec):
        return self.c_solver.get_curr_term(order_vec)

    def set_data_accessor(self, which_val, which_poly, from_k):
        return self.c_solver.set_data_accessor(which_val, which_poly, from_k)

    def load_data_and_move(self):
        return self.c_solver.load_data_and_move()

    def print_poly_info(self, which_poly):
        return self.c_solver.print_poly_info(which_poly)
    
    def get_poly_val_dim(self, which_poly):
        return self.c_solver.get_poly_val_or_var_dim(which_poly, False)

    def get_poly_var_dim(self, which_poly):
        return self.c_solver.get_poly_val_or_var_dim(which_poly, True)

    def set_Kceil(self, Kceil):
        if(Kceil >= self.Kmax):
            raise ValueError("Kceil should not exceed Kmax - 1\n")
        else:
            self.c_solver.Kceil = Kceil

    def get_curr_Kceil(self):
        return self.c_solver.Kceil

    def solve_step(self):
        self.c_solver.solve_step()

    def solve_to_Kceil(self, print_log:bool = True):
        self.c_solver.solve_to_Kceil(print_log)

    def eval(self, which_poly, x_arr):
        return self.c_solver.eval(which_poly, x_arr)

    def get_k(self):
        return self.c_solver.get_k()

    def calculate_err(self, err_k):
        self.c_solver.calculate_err(err_k)

    def clear_all(self):
        self.c_solver.clear_all()

    def estimite_conv_range(self, start_k=2, tol=1e-3):
        if(start_k < 2):
            start_k = 2
        min_range = 1e200
        curr_k = self.get_k()
        for k in range(start_k, min(curr_k+2, self.Kmax-1)):
            self.calculate_err(k)
            curr_err = self.eval(POLY_Ek, np.ones((self.manifold_dim,1)))
            # curr_range = (tol/ np.max(np.abs(curr_err))) ** (1/k)
            curr_range_inv = (np.max(np.abs(curr_err)) / tol) ** (1/k)
            if(curr_range_inv > 1/min_range):
                min_range = 1/curr_range_inv

        return min_range

    # retrieve data
    def get_poly_data(self, which_poly, from_k):
        poly_data = []
        poly_var_dim = self.get_poly_var_dim(which_poly)
        poly_val_dim = self.get_poly_val_dim(which_poly)

        for val_id in range(poly_val_dim):
            poly_scalar_data = []

            continue_flag = self.set_data_accessor(val_id, which_poly, from_k)
            while(continue_flag):
                curr_order = CIndexVec([0]*poly_var_dim)
                curr_coeff = self.get_curr_term(curr_order)
                poly_scalar_data.append((curr_coeff,curr_order))

                continue_flag = self.load_data_and_move()

            poly_data.append(poly_scalar_data)

        return poly_data

    def get_poly_sympy(self, which_poly, vars, from_k=1):
        poly_data = self.get_poly_data(which_poly, from_k)
        var_dim = self.get_poly_var_dim(which_poly)
        polys = []
        for j in range(len(poly_data)):
            curr_poly = sym.Poly(0.0, vars[0])
            for term in poly_data[j]:
                curr_expr = term[0]
                for var_id in range(var_dim):
                    curr_expr *= (vars[var_id] ** term[1][var_id])
                curr_term = sym.Poly(curr_expr)
                curr_poly += curr_term
            polys.append(curr_poly)
        return polys

class PolynomialSystemSolver(InvariantManifoldSolverPy):
    sys_poly: list
    sys_var_list: list

    def __init__(self, sys_poly:list, sys_var_list:list, manifold_dim, Kmax:int):
        self.sys_poly = sys_poly
        self.sys_var_list = sys_var_list
        super().__init__(len(sys_poly), manifold_dim, Kmax)

    def set_F_around_x(self, x0):
        # substitute list
        subs_list = []
        for j, var in enumerate(self.sys_var_list):
            subs_list.append((var, var+x0[j]))

        # substitute
        for val_id in range(len(self.sys_poly)):
            curr_poly = self.sys_poly[val_id]
            curr_poly_around_x0 = (curr_poly.subs(subs_list)).as_poly(*(self.sys_var_list), domain='CC')

            curr_poly_terms = curr_poly_around_x0.monoms()
            curr_poly_coeff = curr_poly_around_x0.coeffs()
            for j, curr_term in enumerate(curr_poly_terms):
                if(not (any(curr_term))):
                    if(np.abs(curr_poly_coeff[j]) > TOL):
                        raise ValueError("x0 is not a fixed point")
                    else:
                        continue
                
                self.add_term_F(val_id, CIndexVec(curr_term), curr_poly_coeff[j])


def test_polynomial_system():
    x,y,z = sym.symbols('x,y,z')
    var_list  = [x,y,z]
    sys_poly = [sym.Poly(x+2*y, *var_list, domain='CC'), sym.Poly(y-2*z + 1, *var_list, domain='CC'), 
        sym.Poly(z**2, *var_list, domain='CC')]
    poly_solver = PolynomialSystemSolver(sys_poly, var_list, 1, 100)
    poly_solver.set_F_around_x([2,-1+1e-16,0])
    J = poly_solver.jacobian_by_F()
    print(J)

if __name__ == '__main__':
    test_polynomial_system()