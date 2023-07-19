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



def test_Laurant_basic():
    # test set by vectors
    coeffs = pt.CScalarVec([
        1, 2, 3, 4, 5, 6,7
    ])
    orders = pt.CLaurantIndexVec([
        0, 0, 1,
        -1, 1, 2, 
        3, -2, 1,
        4, 4, 4, 
        -3, -2, 5,
        2, 3, 2,
        -1, -1, -1 
    ])
    my_laurant = pt.CLaurant(3)
    my_laurant.set_Laurant_by_terms(coeffs, orders)
    my_laurant.reduction()
    print(my_laurant.dim)
    print(my_laurant.num_max_orders)
    print(my_laurant.denom_orders)
    print(my_laurant.num.to_str(pt.CStrVec(['x', 'y', 'z'])))

    my_laurant_2 = pt.CLaurant(3)
    my_laurant_2.set_Laurant(my_laurant.num, pt.CMonomial(1j, my_laurant.denom_orders))
    print(my_laurant_2.dim)
    print(my_laurant_2.num_max_orders)
    print(my_laurant_2.denom_orders)
    print(my_laurant_2.num.to_str(pt.CStrVec(['x', 'y', 'z'])))

def test_laurant_reduction():
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
    my_laurant = pt.CLaurant(4)
    my_laurant.set_Laurant(poly, pt.CMonomial(1, pt.CIndexVec([3,2,2,3])))
    print(my_laurant.num_max_orders)
    print(my_laurant.denom_orders)
    print(my_laurant.num.to_str(pt.CStrVec(['x', 'y', 'z', 'w'])))

def test_laurant_flip():
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
    my_laurant = pt.CLaurant(4)
    my_laurant.set_Laurant(poly, pt.CMonomial(1, pt.CIndexVec([3,2,2,3])))
    print(my_laurant.num_max_orders)
    print(my_laurant.denom_orders)
    print(my_laurant.num.to_str(pt.CStrVec(['x', 'y', 'z', 'w'])))
    print(my_laurant.eval(val))
    print(my_laurant.num.eval(val)/np.prod(np.array(val)**np.array(my_laurant.denom_orders)))
    print()

    flip_id = 3
    new_laurant = my_laurant.flip_variable(flip_id)
    print(my_laurant.num_max_orders)
    print(my_laurant.denom_orders)
    print(my_laurant.num.to_str(pt.CStrVec(['x', 'y', 'z', 'w'])))
    print()
    print(new_laurant.num_max_orders)
    print(new_laurant.denom_orders)
    print(new_laurant.num.to_str(pt.CStrVec(['x', 'y', 'z', 'w'])))
    val[flip_id] = 1/val[flip_id]
    print(new_laurant.eval(val))
    print()


if __name__ == '__main__':
    # main_test_partial_eval()
    # test_Laurant_basic()
    # test_laurant_reduction()
    test_laurant_flip()