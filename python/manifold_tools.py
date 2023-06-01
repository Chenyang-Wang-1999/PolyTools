'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2023-05-25
Copyright © Department of Physics, Tsinghua University. All rights reserved
'''

import py_jet_transport as jet
import numpy as np
from scipy import linalg as LA
from scipy import sparse
from scipy.sparse import linalg as spLA
from math import sin,cos,sqrt,pi

from bisect import bisect

'''
    dPx/dt = - Px + omega0 * Py + Px * Pz
    dPy/dt = - Py - omega0 * Px
    dPz/dt = - Gamma1 * Pz - Px^2 + alpha

    param list: [Px, Py, Pz, omega0, Gamma1, alpha]
'''

def get_feedback_model_poly():
    FeedbackModelPolys = [jet.CPolyLinkedList(6) for j in range(3)]
    coeffs = jet.CScalarVec([ -1, 1, 1])
    orders = jet.CIndexVec([
        1, 0, 0, 0, 0, 0,
        0, 1, 0, 1, 0, 0,
        1, 0, 1, 0, 0, 0
    ])
    FeedbackModelPolys[0].batch_add_elements(coeffs, orders)
    coeffs = jet.CScalarVec([-1, -1])
    orders = jet.CIndexVec([
        0, 1, 0, 0, 0, 0,
        1, 0, 0, 1, 0, 0
    ])
    FeedbackModelPolys[1].batch_add_elements(coeffs, orders)
    coeffs = jet.CScalarVec([-1, -1, 1])
    orders = jet.CIndexVec([
        0, 0, 1, 0, 1, 0,
        2, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1
    ])
    FeedbackModelPolys[2].batch_add_elements(coeffs, orders)
    return FeedbackModelPolys

FeedbackModelPolys = get_feedback_model_poly()

def set_expansion_vec(P0, param0, expand_vector:jet.CSeriesVec):
    var_dim = expand_vector.var_dim
    for j in range(3):
        curr_order = jet.CIndexVec([0]*var_dim)
        curr_order[j] = 1
        expand_vector.add_term(j, jet.CMonomial(P0[j], jet.CIndexVec([0]*var_dim)))
        expand_vector.add_term(j, jet.CMonomial(1, curr_order))

    for j in range(3):
        expand_vector.add_term(3+j, jet.CMonomial(param0[j], jet.CIndexVec([0]*4)))



def test_00Z():
    omega0 = 1.2
    Gamma1 = 1
    alpha = 2*Gamma1
    P0 = np.array([0,0,alpha/Gamma1])
    subs_vec = jet.CSeriesVec(4, 6, 5)
    set_expansion_vec(P0, np.array([omega0, Gamma1, alpha]), subs_vec)
    subs_vec.add_term(2, jet.CMonomial(1.0, jet.CIndexVec([0,0,0,1])))
    subs_vec.add_term(5, jet.CMonomial(Gamma1, jet.CIndexVec([0,0,0,1])))
    print(subs_vec.get_poly(0,0).to_str(['x','y','z','mu']))
    print(subs_vec.get_poly(0,1).to_str(['x','y','z','mu']))
    print(subs_vec.get_poly(0,2).to_str(['x','y','z','mu']))
    print(subs_vec.get_poly(1,0).to_str(['x','y','z','mu']))
    print(subs_vec.get_poly(1,1).to_str(['x','y','z','mu']))
    print(subs_vec.get_poly(2,0).to_str(['x','y','z','mu']))
    print(subs_vec.get_poly(2,1).to_str(['x','y','z','mu']))
    print(subs_vec.get_poly(3,0).to_str(['x','y','z','mu']))
    print(subs_vec.get_poly(3,1).to_str(['x','y','z','mu']))
    print(subs_vec.get_poly(4,0).to_str(['x','y','z','mu']))
    print(subs_vec.get_poly(4,1).to_str(['x','y','z','mu']))
    print(subs_vec.get_poly(5,0).to_str(['x','y','z','mu']))
    print(subs_vec.get_poly(5,1).to_str(['x','y','z','mu']))
    model_expansion = [jet.CPolyLinkedList(4) for j in range(3)]
    for val_id in range(3):
        for k in range(5):
            print(val_id, k)
            curr_poly = jet.CPolyLinkedList(4)
            jet.poly_comp(FeedbackModelPolys[val_id], subs_vec, k, curr_poly)
            model_expansion[val_id].destructive_add_self(curr_poly)
            print(model_expansion[val_id].to_str(['x','y','z','mu']))

if __name__ == '__main__':
    test_00Z()
    # test_map_to_mat()