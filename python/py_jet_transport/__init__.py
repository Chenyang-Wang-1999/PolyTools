'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2022-12-22
Copyright © Department of Physics, Tsinghua University. All rights reserved
'''

MODE = 'rr'
if(MODE == 'cc'):
    import py_jet_transport._jet_transport_cc as c_jet_transport
elif(MODE == 'rc'):
    import py_jet_transport._jet_transport_rc as c_jet_transport
elif(MODE == 'rr'):
    import py_jet_transport._jet_transport_rr as c_jet_transport

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

class CIndexVec(c_jet_transport.CIndexVec):
    def __init__(self, vec):
        super().__init__(vec)

class CScalarVec(c_jet_transport.CScalarVec):
    def __init__(self, vec):
        super().__init__(vec)

if(MODE == 'rc'):
    class CVarScalarVec(c_jet_transport.CVarScalarVec):
        def __init__(self, vec):
            super().__init__(vec)
else:
    class CVarScalarVec(c_jet_transport.CScalarVec):
        def __init__(self, vec):
            super().__init__(vec)

class CPolyLinkedList(c_jet_transport._CPolyLinkedList):
    def __init__(self, dim):
        super().__init__(dim)

    def reinit(self, dim):
        super().reinit(dim)
    
    def remove_zeros(self):
        super().remove_zeros()

    def copy_to(self, another):
        super().copy_to(another)

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

class CMonomial(c_jet_transport._Monomial):
    def __init__(self, coeff, orders) -> None:
        super().__init__(coeff, CIndexVec(orders))

    def eval(self, x):
        return super().eval(CVarScalarVec(x))

    def derivative(self, var_id):
        return super().derivative(var_id) 

class CSeriesVec(c_jet_transport._CSeriesVec):
    def __init__(self, var_dim, val_dim, Kmax):
        super().__init__(var_dim, val_dim, Kmax)
    def reinit(self):
        super().reinit()
    def copy_to(self, another):
        super().copy_to(another)
    def add_term(self, val_id, term:CMonomial):
        super().add_term(val_id, term)
    def destructive_add_poly(self, val_id, poly):
        super().destructive_add_poly(val_id, poly)
    def get_poly(self, val_id, curr_order):
        poly = CPolyLinkedList(self.var_dim)
        super().get_poly(val_id, curr_order, poly)
        return poly


def poly_load_from_file(fname):
    data = spio.loadmat(fname)
    coeffs = data['coeffs'].flatten()
    orders = data['orders'].flatten()
    dim = data['dim'][0,0]

    new_poly = CPolyLinkedList(dim)
    new_poly.init_with_data(coeffs, orders)
    return new_poly
        

class cJetTransport(c_jet_transport._cJetTransport):
    var_dim:int
    param_dim:int
    max_order:int
    def __init__(self, var_dim:int, param_dim:int, max_order:int) -> None:
        super().__init__(var_dim, param_dim, max_order)

    def set_ode_fun(self, dim:int, poly:CPolyLinkedList):
        super().set_ode_fun(dim, poly)
    def set_sol(self, dim:int, poly:CPolyLinkedList):
        super().set_sol(dim, poly)
    def get_ode_fun(self, dim:int, order:int):
        poly = CPolyLinkedList(self.var_dim + self.param_dim)
        super().get_ode_fun(dim, order, poly)
        return poly
    def get_sol(self, dim:int, order:int):
        poly = CPolyLinkedList(self.var_dim + self.param_dim)
        super().get_sol(dim, order, poly)
        return poly
    def get_sol_time(self, dim:int, order:int):
        poly = CPolyLinkedList(self.var_dim + self.param_dim + 1)
        super().get_sol_time(dim, order, poly)
        return poly
    def get_dvars_time(self, dim:int, order:int):
        poly = CPolyLinkedList(self.var_dim + self.param_dim + 1)
        super().get_dvars_time(dim, order, poly)
        return poly
    def set_sol_by_scalars(self, coeffs:CScalarVec, err_order:int):
        super().set_sol_by_scalars(coeffs, err_order)
    def get_sol_as_scalars(self,err_order:int):
        coeffs = CScalarVec([])
        super().get_sol_as_scalars(coeffs, err_order)
        return coeffs
    def get_dvars_as_scalars(self,err_order:int):
        coeffs = CScalarVec([])
        super().get_dvars_as_scalars(coeffs, err_order)
        return coeffs
    def get_series_indices(self, j:int):
        return super().get_series_indices(j)
    def calculate_dvars(self, err_order:int):
        super().calculate_dvars(err_order)
    def set_time_variable(self):
        super().set_time_variable()
    def derivative_iterate(self):
        super().derivative_iterate()
    def Poincare_projection(self, err_order:int, normal_vec:CScalarVec):
        super().Poincare_projection(err_order, normal_vec)

def poly_comp(f:CPolyLinkedList, v:CSeriesVec, k:int, res:CPolyLinkedList):
    c_jet_transport.poly_comp(f, v, k, res)

def poly_multiplication(f:CPolyLinkedList, g: CPolyLinkedList):
    h = CPolyLinkedList(f.dim)
    c_jet_transport.poly_multiplication(f,g,h)
    return h