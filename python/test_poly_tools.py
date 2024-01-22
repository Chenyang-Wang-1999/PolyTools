'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2023-07-10
Copyright © Department of Physics, Tsinghua University.  All rights reserved
'''

import numpy as np
import poly_tools as pt
def main_test_partial_eval():
    my_poly = pt.CPolyLinkedList(3)
    my_poly.batch_add_elements(
        pt.CScalarVec([
            1,2,3
        ]),
        pt.CIndexVec([
            1,0,0,
            1,1,1,
            2,0,1
        ])
    )
    print(my_poly.to_str(['x','y','z']))

    # partial evaluation of z = 0
    new_poly = my_poly.partial_eval(
        pt.CVarScalarVec([0.0]),
        pt.CIndexVec([2]), # z
        pt.CIndexVec([0,1]) # dof map: 0:x, 1:y
    )
    print(new_poly.to_str(['x','y']))

    new_poly = my_poly.partial_eval(
        pt.CVarScalarVec([0]),
        pt.CIndexVec([2]), # z
        pt.CIndexVec([1,0]) # dof map: 0:x, 1:y
    )
    print(new_poly.to_str(['y','x']))

    new_poly = my_poly.scale_var(pt.CScalarVec([2,1,-1j]))
    print(new_poly.to_str(['x','y','z']))

def test_max_and_min_orders():
    my_poly = pt.CPolyLinkedList(3)
    my_poly.batch_add_elements(
        pt.CScalarVec([
            1,2,3
        ]),
        pt.CIndexVec([
            1,0,0,
            1,1,1,
            2,0,1
        ])
    )
    print(my_poly.to_str(['x','y','z']))



def test_Laurent_basic():
    # test set by vectors
    coeffs = pt.CScalarVec([
        1, 2, 3, 4, 5, 6,7
    ])
    orders = pt.CLaurentIndexVec([
        0, 0, 1,
        -1, 1, 2, 
        3, -2, 1,
        4, 4, 4, 
        -3, -2, 5,
        2, 3, 2,
        -1, -1, -1 
    ])
    my_Laurent = pt.CLaurent(3)
    my_Laurent.set_Laurent_by_terms(coeffs, orders)
    my_Laurent.reduction()
    print(my_Laurent.dim)
    print(my_Laurent.num_max_orders)
    print(my_Laurent.denom_orders)
    print(my_Laurent.num.to_str(pt.CStrVec(['x', 'y', 'z'])))

    my_Laurent_2 = pt.CLaurent(3)
    my_Laurent_2.set_Laurent(my_Laurent.num, pt.CMonomial(1j, my_Laurent.denom_orders))
    print(my_Laurent_2.dim)
    print(my_Laurent_2.num_max_orders)
    print(my_Laurent_2.denom_orders)
    print(my_Laurent_2.num.to_str(pt.CStrVec(['x', 'y', 'z'])))

def test_Laurent_reduction():
    poly = pt.CPolyLinkedList(4)
    scalars = pt.CScalarVec([
        1,2,3,4,5
    ])
    orders = pt.CIndexVec([
        2, 2, 4, 3,
        1, 2, 1, 2, 
        5, 3, 2, 5,
        4, 4, 3, 2,
        1, 4, 2, 1
    ])
    poly.batch_add_elements(scalars, orders)
    print(poly.to_str(['x', 'y', 'z', 'w']))
    my_Laurent = pt.CLaurent(4)
    my_Laurent.set_Laurent(poly, pt.CMonomial(1, pt.CIndexVec([3,2,2,3])))
    print(my_Laurent.num_max_orders)
    print(my_Laurent.denom_orders)
    print(my_Laurent.num.to_str(pt.CStrVec(['x', 'y', 'z', 'w'])))

def test_Laurent_flip():
    val = pt.CScalarVec(np.random.randn(4) + np.random.randn(4)*1j)
    poly = pt.CPolyLinkedList(4)
    scalars = pt.CScalarVec([
        1,2,3,4,5
    ])
    orders = pt.CIndexVec([
        2, 2, 4, 3,
        1, 2, 1, 2, 
        5, 3, 2, 5,
        4, 4, 3, 2,
        1, 4, 2, 1
    ])
    poly.batch_add_elements(scalars, orders)
    print(poly.to_str(['x', 'y', 'z', 'w']))
    my_Laurent = pt.CLaurent(4)
    my_Laurent.set_Laurent(poly, pt.CMonomial(1, pt.CIndexVec([3,2,2,3])))
    print(my_Laurent.num_max_orders)
    print(my_Laurent.denom_orders)
    print(my_Laurent.num.to_str(pt.CStrVec(['x', 'y', 'z', 'w'])))
    print(my_Laurent.eval(val))
    print(my_Laurent.num.eval(val)/np.prod(np.array(val)**np.array(my_Laurent.denom_orders)))
    print()

    flip_id = 3
    new_Laurent = my_Laurent.flip_variable(flip_id)
    print(my_Laurent.num_max_orders)
    print(my_Laurent.denom_orders)
    print(my_Laurent.num.to_str(pt.CStrVec(['x', 'y', 'z', 'w'])))
    print()
    print(new_Laurent.num_max_orders)
    print(new_Laurent.denom_orders)
    print(new_Laurent.num.to_str(pt.CStrVec(['x', 'y', 'z', 'w'])))
    val[flip_id] = 1/val[flip_id]
    print(new_Laurent.eval(val))
    print()

def test_Laurent_partial_eval():
    # 1. generate a random Laurent polynomial
    dim = 3
    N_terms = 100
    max_deg = 5
    min_deg = -4
    all_degs = np.random.randint(min_deg, max_deg, size=N_terms * dim)
    all_coeffs = np.random.rand(N_terms) + 1j * np.random.rand(N_terms)

    f = pt.CLaurent(dim)
    f.set_Laurent_by_terms(pt.CScalarVec(all_coeffs), pt.CLaurentIndexVec(all_degs))

    # 2. generate x0
    x0 = np.random.rand(dim) + 1j * np.random.rand(dim)

    # 3. calculate partial evaluation
    f1 = f.partial_eval(pt.CVarScalarVec([x0[1]]), pt.CIndexVec([1]), pt.CIndexVec([0,2]))

    # eval
    print(f.eval(pt.CVarScalarVec(x0)))
    print(f1.eval(pt.CVarScalarVec(x0[[0,2]])))

def test_Laurent_derivative():
    # 1. generate a random Laurent polynomial
    dim = 3
    N_terms = 100
    max_deg = 5
    min_deg = -4
    all_degs = np.random.randint(min_deg, max_deg, size=N_terms * dim)
    all_coeffs = np.random.rand(N_terms) + 1j * np.random.rand(N_terms)

    f = pt.CLaurent(dim)
    f.set_Laurent_by_terms(pt.CScalarVec(all_coeffs), pt.CLaurentIndexVec(all_degs))

    diff_var = 1
    df = f.derivative(diff_var)

    # 2. generate x
    x0 = np.random.rand(dim) + 1j * np.random.rand(dim)
    dx = np.random.rand() * 1e-5
    x1 = x0.copy()
    x1[diff_var] += dx

    # 3. evaluate
    fx = f.eval(pt.CVarScalarVec(x0))
    fx2 = f.eval(pt.CVarScalarVec(x1))
    dfx = df.eval(pt.CVarScalarVec(x0))
    print(fx2 - fx)
    print(dfx*dx)
    print(abs(fx2 - fx - dfx*dx)/abs(dfx*dx))

def test_Laurent_to_str():
    coeffs = pt.CScalarVec([
        1, -1.15, -0.3, -0.3, -0.3, -0.3, 0.02, 0.02, 0.02, 0.02, 
        -0.5, -0.4, 0.04, 0.04
    ])
    degs = pt.CLaurentIndexVec([
        2, 0, 0, 
        0, 0, 0,
        1, 1, 0,
        1, -1, 0,
        1, 0, -1,
        1, 0, 1,
        0, 2, 0,
        0, -2, 0,
        0, 0, 2,
        0, 0, -2,
        0, 1, 1,
        0, -1, -1,
        0, 1, -1,
        0, -1, 1 
    ])

    f = pt.CLaurent(3)
    f.set_Laurent_by_terms(coeffs, degs)
    print(f.to_str(pt.CStrVec(['x', 'y', 'z'])))

if __name__ == '__main__':
    # main_test_partial_eval()
    # test_Laurent_basic()
    # test_Laurent_reduction()
    # test_Laurent_flip()
    # test_Laurent_partial_eval()
    # test_Laurent_derivative()
    test_Laurent_to_str()