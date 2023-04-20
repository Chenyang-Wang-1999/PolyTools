'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2022-12-22
Copyright © Department of Physics, Tsinghua University. All rights reserved
'''

MODE = 'cc'
if(MODE == 'cc'):
    import py_invariant_manifold._invariant_manifold_cc as c_invariant_manifold
elif(MODE == 'rc'):
    import py_invariant_manifold._invariant_manifold_rc as c_invariant_manifold
elif(MODE == 'rr'):
    import py_invariant_manifold._invariant_manifold_rr as c_invariant_manifold

# import _invariant_manifold as c_invariant_manifold
import sympy as sym
import numpy as np
import csv
from scipy import io as spio
from scipy import linalg as LA

# definitions
POLY_F = 0
POLY_f = 1
POLY_W = 2
POLY_s = 3
POLY_Ek = 4


TOL = 1e-12

class CIndexVec(c_invariant_manifold.CIndexVec):
    def __init__(self, vec):
        super().__init__(vec)

class CScalarVec(c_invariant_manifold.CScalarVec):
    def __init__(self, vec):
        super().__init__(vec)

if(MODE == 'rc'):
    class CVarScalarVec(c_invariant_manifold.CVarScalarVec):
        def __init__(self, vec):
            super().__init__(vec)
else:
    class CVarScalarVec(c_invariant_manifold.CScalarVec):
        def __init__(self, vec):
            super().__init__(vec)

class CPolyLinkedList(c_invariant_manifold._CPolyLinkedList):
    def __init__(self, dim):
        super().__init__(dim)

    def reinit(self, dim):
        super().reinit(dim)
    
    def remove_zeros(self):
        super().remove_zeros()

    def copy(self, another):
        super().copy(another)

    def print_info(self):
        super().print_info()
    
    def destructive_add_self(self, another):
        super().destructive_add_self(another)

    def destructive_add(self, another):
        result = CPolyLinkedList(self.dim)
        super().destructive_add(another, result)
        return result

    def add_self(self, another):
        super().add_self(another)

    def add(self, another):
        result = CPolyLinkedList(self.dim)
        super().add(another, result)
        return result

    def scalar_mul_self(self, k):
        super().scalar_mul_self(k)
    
    def scalar_mul(self, k):
        result = CPolyLinkedList(self.dim)
        super().scalar_mul(k, result)
        return result

    def scale_var_self(self, k_list):
        if(not isinstance(k_list, CScalarVec)):
            k_list = CScalarVec(k_list)
        super().scale_var_self(k_list)

    def scale_var(self, k_list):
        if(not isinstance(k_list, CScalarVec)):
            k_list = CScalarVec(k_list)
        new_poly = CPolyLinkedList(self.dim)
        super().scale_var(k_list, new_poly)
        return new_poly

    def neg_self(self):
        super().neg_self()

    def neg(self):
        result = CPolyLinkedList(self.dim)
        super().neg(result)
        return result

    def destructive_subs_self(self, another):
        super().destructure_subs_self(another)

    def subs_self(self, another):
        super().subs_self(another)

    def destructive_subs(self, another):
        result = CPolyLinkedList(self.dim)
        super().destructure_subs(another, result)
        return result
    
    def subs(self, another):
        result = CPolyLinkedList(self.dim)
        super().subs(another, result)
        return result

    def derivative(self, var_id:int):
        result = CPolyLinkedList(self.dim)
        super().derivative(var_id, result)
        return result

    def eval(self, x):
        if(not isinstance(x, CVarScalarVec)):
            x = CVarScalarVec(x)
        return super().eval(x)

    def eval_diff(self, diff_order, x):
        if(not isinstance(diff_order, CIndexVec)):
            diff_order = CIndexVec(diff_order)

        if(not isinstance(x, CVarScalarVec)):
            x = CVarScalarVec(x)

        return super().eval_diff(diff_order, x)

    def batch_eval(self, diff_order, x_arr):
        if(not isinstance(diff_order, CIndexVec)):
            diff_order = CIndexVec(diff_order)
        return super().batch_eval(diff_order, x_arr)

    def batch_add_elements(self, coeffs:CScalarVec, orders:CIndexVec):
        if(not isinstance(coeffs, CScalarVec)):
            coeffs = CScalarVec(coeffs)
        if(not isinstance(orders, CIndexVec)):
            orders = CIndexVec(orders)
        super().batch_add_elements(coeffs, orders)

    def batch_get_data(self):
        _coeffs = CScalarVec([])
        _orders = CIndexVec([])
        # print("get_data_begin")
        super().batch_get_data(_coeffs, _orders)
        # print("get_data_finished")
        coeffs = np.array(_coeffs)
        orders = np.array(_orders, dtype=int).reshape((-1,self.dim))
        return (coeffs, orders)

    def init_with_data(self, coeff, orders):
        _coeffs = CScalarVec(coeff)
        _orders = CIndexVec(orders)
        super().init_with_data(_coeffs, _orders)

    # format: coeffs: np.array(complex), orders: np.array(int)
    def save_to_file(self, fname):
        coeffs, orders = self.batch_get_data()
        data = {}
        data['coeffs'] = coeffs
        data['orders'] = orders
        data['dim'] = self.dim
        spio.savemat(fname, data)

    def to_str(self, varname:list):
        coeffs, orders = self.batch_get_data()
        poly_str = ''
        for j in range(len(coeffs)):
            if(isinstance(coeffs[j], complex)):
                # coeff_str = str(coeffs[j])
                # poly_str += "+" + "%s"%(coeff_str.replace('j','i'))
                poly_str += "+" + "(%+.18e%+.18e*i)"%(coeffs[j].real, coeffs[j].imag)
            else:
                poly_str += "+" + "(%s)"%(str(coeffs[j]))
            for k in range(self.dim):
                if(orders[j, k] > 1):
                    poly_str += '*(%s**%d)'%(varname[k], orders[j,k])
                elif(orders[j, k] == 1):
                    poly_str += '*%s'%(varname[k])
        poly_str += ';'
        return poly_str   

class CMonomial(c_invariant_manifold._Monomial):
    def __init__(self, coeff, orders) -> None:
        super().__init__(coeff, CIndexVec(orders))

    def eval(self, x):
        return super().eval(CVarScalarVec(x))

    def derivative(self, var_id):
        return super().derivative(var_id) 


def poly_load_from_file(fname):
    data = spio.loadmat(fname)
    coeffs = data['coeffs'].flatten()
    orders = data['orders'].flatten()
    dim = data['dim'][0,0]

    new_poly = CPolyLinkedList(dim)
    new_poly.init_with_data(coeffs, orders)
    return new_poly
        

class InvariantManifoldSolverPy:
    c_solver: c_invariant_manifold._CInvariantManifoldSolver
    phys_dim:int 
    manifold_dim:int 
    Kmax:int
    scale_factor: complex

    curr_poly:list
    curr_poly_raw:list
    # deprecated
    # curr_sym_poly:list
    # curr_sym_Dpoly:list
    def __init__(self, phys_dim:int, manifold_dim:int, Kmax:int):
        self.c_solver = c_invariant_manifold._CInvariantManifoldSolver(phys_dim, manifold_dim, Kmax)
        self.phys_dim = self.c_solver.phys_dim
        self.manifold_dim = self.c_solver.manifold_dim
        self.Kmax = self.c_solver.Kmax
        self.curr_poly = [[]] * 3
        self.curr_poly_raw = [[]] * 3
        self.scale_factor = 1.0

    def __del__(self):
        del(self.c_solver)
        

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
        self.phys_dim = self.c_solver.phys_dim
        self.manifold_dim = self.c_solver.manifold_dim
        self.Kmax = self.c_solver.Kmax
        self.curr_poly = [[]] * 3
        self.curr_poly_raw = [[]] * 3
        self.scale_factor = 1.0

        self.c_solver.clear_all()

    def reinit_with_scale_factor(self, scale_factor):
        # reinitialize with scale factor
        # u = k * s
        self.scale_factor = scale_factor
        F_poly = self.get_poly(POLY_F)

        # clear and append
        self.clear_all()
        for val_id in range(self.phys_dim):
            F_poly[val_id].scalar_mul_self(self.scale_factor)
            curr_coeffs, curr_orders = F_poly[val_id].batch_get_data()
            for term_id in range(len(curr_coeffs)):
                self.add_term_F(val_id, CIndexVec(curr_orders[term_id, :]), curr_coeffs[term_id])


    def estimate_conv_range(self, start_k=2, tol=1e-3):
        if(start_k < 2):
            start_k = 2
        min_range = 1e200
        curr_k = self.get_k()
        for k in range(start_k, min(curr_k+2, self.Kmax-1)):
            self.calculate_err(k)
            curr_err = self.eval(POLY_Ek, np.ones((self.manifold_dim,1)))/self.scale_factor
            # curr_range = (tol/ np.max(np.abs(curr_err))) ** (1/k)
            curr_range_inv = (np.max(np.abs(curr_err)) / tol) ** (1/k)
            if(curr_range_inv > 1/min_range):
                min_range = 1/curr_range_inv

        return min_range

    # overriding functions
    def tangent_vector(self, P):
        if(len(self.curr_poly[POLY_F]) == 0):
            self.update_poly(POLY_F)
        
        curr_tangent = np.zeros((self.phys_dim,1), dtype = complex)
        for curr_dim in range(self.phys_dim):
            curr_tangent[curr_dim] = self.curr_poly[POLY_F][curr_dim].eval(P)
        return curr_tangent

    def jacobian(self, P):
        if(len(self.curr_poly[POLY_F]) == 0):
            self.update_poly(POLY_F)
        
        jac = np.zeros((self.phys_dim, self.phys_dim), dtype = complex)
        for curr_dim in range(self.phys_dim):
            for curr_diff in range(self.phys_dim):
                diff_order = CIndexVec([0]*self.phys_dim)
                diff_order[curr_diff] = 1
                jac[curr_dim, curr_diff] = self.curr_poly[POLY_F][curr_dim].eval_diff(diff_order, P)
        return jac


    def get_poly_for_val(self, which_poly, val_id):
        curr_poly = CPolyLinkedList(self.get_poly_var_dim(which_poly))
        self.c_solver.get_poly(which_poly,curr_poly , val_id)
        return curr_poly

    def get_poly(self, which_poly):
        all_polys = []
        for val_id in range(self.get_poly_val_dim(which_poly)):
            all_polys.append(self.get_poly_for_val(which_poly, val_id))
        return all_polys

    # retrieve data
    def get_poly_data(self, which_poly, from_k):
        # old
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

    def get_poly_sympy(self, which_poly, poly_vars, from_k=1):
        poly_data = self.get_poly_data(which_poly, from_k)
        var_dim = self.get_poly_var_dim(which_poly)
        polys = []
        for j in range(len(poly_data)):
            curr_poly = sym.Poly(0.0, poly_vars[0])
            for k, term in enumerate(poly_data[j]):
                print(k/len(poly_data[j]) * 100, '%')
                curr_expr = term[0]
                for var_id in range(var_dim):
                    curr_expr *= (poly_vars[var_id] ** term[1][var_id])
                curr_term = sym.Poly(curr_expr)
                curr_poly += curr_term
            polys.append(curr_poly)
        return polys

    # def update_poly_sympy(self, which_poly, poly_vars, from_k=1):
    #     new_poly = self.get_poly_sympy(which_poly, poly_vars, from_k)
    #     self.curr_sym_poly[which_poly] = new_poly
    #     # if(poly_id < len(self.curr_sym_poly)):
    #     #     self.curr_sym_poly[poly_id] = new_poly
    #     # else:
    #     #     self.curr_sym_poly.append(new_poly)

    # def update_Dpoly_sympy(self, which_poly, poly_vars, from_k=1):
    #     if(len(self.curr_sym_poly[which_poly]) == 0):
    #         self.update_poly_sympy(which_poly, poly_vars, from_k)
    #     new_poly = self.curr_sym_poly[which_poly]
    #     new_Dpoly =  []
    #     for j in range(self.get_poly_val_dim(which_poly)):
    #         new_Dpoly_dim =[]
    #         for k in range(self.get_poly_var_dim(which_poly)):
    #             new_Dpoly_dim.append(new_poly[j].diff(k))
    #         new_Dpoly.append(new_Dpoly_dim)

    #     self.curr_sym_Dpoly[which_poly] = new_poly
    #     # if(poly_id < len(self.curr_sym_poly)):
    #     #     self.curr_sym_poly[poly_id] = new_Dpoly
    #     # else:
    #     #     self.curr_sym_poly.append(new_Dpoly)

    def update_poly(self, which_poly):
        self.curr_poly_raw[which_poly] = self.get_poly(which_poly)
        if(which_poly == POLY_F):
            curr_poly = []
            for val_id in range(self.get_poly_val_dim(which_poly)):
                curr_poly.append(self.curr_poly_raw[which_poly][val_id].scalar_mul(1/self.scale_factor))
        else:
            curr_poly = []
            for val_id in range(self.get_poly_val_dim(which_poly)):
                curr_poly.append(self.curr_poly_raw[which_poly][val_id].scale_var(CScalarVec([1/self.scale_factor]*self.get_poly_var_dim(which_poly))))
        self.curr_poly[which_poly] = curr_poly

    # def sympy_eval(self, which_poly, x_arr):
    #     # evaluate the polynomial by sympy
    #     val_dim = self.get_poly_val_dim(which_poly)
    #     var_dim = self.get_poly_var_dim(which_poly)

    #     poly_vars = sym.symbols('x_{0:%d}'%(var_dim))
    #     if(len(self.curr_sym_poly[which_poly]) == 0):
    #         self.update_poly_sympy(which_poly, poly_vars, 1)

    #     vals = np.zeros((val_dim, x_arr.shape[1]), dtype=complex)
    #     for j in range(x_arr.shape[1]):
    #         for k in range(val_dim):
    #             vals[k,j] = self.curr_sym_poly[which_poly][k].eval(list(x_arr[:,j]))
    #     return vals
        # sym_poly = self.get_poly_sympy(which_poly, poly_vars, 1)

    # def load_F_from_file(self, val_id:int, fname:str):
    #     data_fp = open(fname, "r")
    #     reader = csv.reader(data_fp)

    #     for line in reader:
    #         if(len(line) != self.phys_dim + 1):
    #             data_fp.close()
    #             self.clear_all()
    #             raise ValueError("data type do not match phys dim")

    #         coeff = complex(line[0])
    #         orders  = CIndexVec([0]*(len(line) - 1))
    #         for j in range(1, len(line)):
    #             orders[j-1] = int(line[j])

    #         self.add_term_F(val_id, orders, coeff)

    #     data_fp.close()

    def calculate_err_angle(self, s_arr):
        # calculate the angle between the tangent vector and the manifold

        # get W
        if(len(self.curr_poly[POLY_W]) == 0):
            self.update_poly(POLY_W)

        err_vec = np.zeros((s_arr.shape[1],))
        for j in range(s_arr.shape[1]):
            curr_s = list(s_arr[:,j])
            curr_W = np.zeros((self.phys_dim,), dtype = complex)
            curr_DW = np.zeros((self.phys_dim, self.manifold_dim), dtype=complex)
            for k in range(self.phys_dim):
                curr_W[k] = self.curr_poly[POLY_W][k].eval(CVarScalarVec(curr_s))
                for m in range(self.manifold_dim):
                    diff_order = CIndexVec([0]*self.manifold_dim)
                    diff_order[m] = 1
                    curr_DW[k,m] = self.curr_poly[POLY_W][k].eval_diff(diff_order, CVarScalarVec(curr_s))

            # get tangent vector
            tan_vec = self.tangent_vector(curr_W)
            tan_vec /= LA.norm(tan_vec)

            # normalize W
            for m in range(self.manifold_dim):
                curr_DW[:,m] /= LA.norm(curr_DW[:,m])

            # solve least square
            res = LA.lstsq(curr_DW, tan_vec)[1]
            if(len(res)):
                err_vec[j] = res[0]
            else:
                err_vec[j] = 0.0
        return err_vec

    def move_to_converge_range(self, s_points, tol = 1e-12, tol_scale = 1e-4, scale_bound = 1.0):
        new_s_points = np.zeros_like(s_points, dtype=complex)
        for j in range(s_points.shape[1]):
            print("point id:",j)
            curr_s = s_points[:,j].reshape((-1,1))
            curr_err = self.calculate_err_angle(curr_s)[0]
            # print(curr_err)

            scale_max = 1.0
            scale_min = 1.0
            if(curr_err > tol):
                while(curr_err > tol):
                    scale_min /= 2
                    curr_err = self.calculate_err_angle(curr_s*scale_min)[0]
                    # print(scale_min,curr_err)
            else:
                while(curr_err <= tol):
                    scale_max *= 2
                    if(curr_err > scale_bound):
                        break
                    curr_err = self.calculate_err_angle(curr_s*scale_max)[0]
                    # print(scale_max,curr_err)
            
            # bisect to converge
            curr_scale = np.sqrt(scale_max * scale_min)
            while(np.abs(scale_max/scale_min - 1)> tol_scale):
                curr_err = self.calculate_err_angle(curr_s * curr_scale)[0]
                # print("curr err:", curr_err)
                # print("curr scale:", curr_scale)
                # print("scale max:", scale_max)
                # print("scale min:", scale_min)
                if(curr_err > tol):
                    scale_max = curr_scale
                else:
                    scale_min = curr_scale
                curr_scale = np.sqrt(scale_max * scale_min)
            
            new_s_points[:,j] = s_points[:,j] * curr_scale
        return new_s_points

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
                        print(np.abs(curr_poly_coeff[j]))
                        raise ValueError("x0 is not a fixed point")
                    else:
                        continue
                
                self.add_term_F(val_id, CIndexVec(curr_term), curr_poly_coeff[j])

# print_poly_data, name = fname + "%d.csv"
def print_poly_data(poly_data, fname:str):
    for j in range(len(poly_data)):
        fp = open(fname + "%d.csv"%(j), "w")
        writer = csv.writer(fp)
        for item in poly_data[j]:
            curr_row = [item[0]]
            for p in item[1]:
                curr_row.append(p)
            writer.writerow(curr_row)
        fp.close()

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