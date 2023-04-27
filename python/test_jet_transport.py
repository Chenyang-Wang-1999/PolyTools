'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2023-04-21
Copyright © Department of Physics, Tsinghua University. All rights reserved
'''

import py_jet_transport as jet
import numpy as np

def test_ode_fun():
    ''' 
        test input and output of ode functions
        example:
            dx/dt = - x + 2*y + x*z
            dy/dt = - y - 2*x
            dz/dt = - z + 3(1 + deltaR) + x^2
    '''

    coeffs0 = jet.CScalarVec([
        -1, 2, 1
    ])
    orders0 = jet.CIndexVec([
        1,0,0,0, 
        0,1,0,0,
        1,0,1,0
    ])
    poly0 = jet.CPolyLinkedList(4)
    poly0.batch_add_elements(coeffs0, orders0)

    coeffs1 = jet.CScalarVec([
        -1, -2
    ])
    orders1 = jet.CIndexVec([
        0,1,0,0, 
        1,0,0,0
    ])
    poly1 = jet.CPolyLinkedList(4)
    poly1.batch_add_elements(coeffs1, orders1)

    coeffs2 = jet.CScalarVec([
        -1, 3, 3, 1
    ])
    orders2 = jet.CIndexVec([
        0,0,1,0,
        0,0,0,0,
        0,0,0,1,
        2,0,0,0
    ])
    poly2 = jet.CPolyLinkedList(4)
    poly2.batch_add_elements(coeffs2, orders2)

    jet_transport_test = jet.cJetTransport(3, 1, 5)
    # print(jet_transport_test.var_dim)
    # print(jet_transport_test.param_dim)
    # print(jet_transport_test.max_order)
    jet_transport_test.set_sol(0, poly0)
    jet_transport_test.set_sol(1, poly1)
    jet_transport_test.set_sol(2, poly2)

    for dim in range(3):
        print("ode fun %d:"%(dim))
        for order in range(4):
            print("order:",order)
            jet_transport_test.get_sol(dim, order).print_info()

def test_conversion_coeffs():

    sol_coeffs = jet.CScalarVec([
        0, 0, 0, 1, 0, -2, 0, 0, 0, 0, 0, 0, 1, 0, 0, 
        0, -1, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 1,
        3, 0, 0, 0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0
    ])
    jet_transport_test = jet.cJetTransport(3, 1, 5)
    # coeffs = jet_transport_test.get_dvars_as_scalars(3)
    # print(coeffs)
    for j in range(len(sol_coeffs)//3):
        print(jet_transport_test.get_series_indices(j))
    jet_transport_test.set_sol_by_scalars(sol_coeffs, 3)
    for j in range(3):
        print("var: ",j)
        for order in range(3):
            jet_transport_test.get_sol(j, order).print_info()
        print()

def test_conversion_coeffs2():
    coeffs0 = jet.CScalarVec([
        -1, 2, 1
    ])
    orders0 = jet.CIndexVec([
        1,0,0,0, 
        0,1,0,0,
        1,0,1,0
    ])
    poly0 = jet.CPolyLinkedList(4)
    poly0.batch_add_elements(coeffs0, orders0)

    coeffs1 = jet.CScalarVec([
        -1, -2
    ])
    orders1 = jet.CIndexVec([
        0,1,0,0, 
        1,0,0,0
    ])
    poly1 = jet.CPolyLinkedList(4)
    poly1.batch_add_elements(coeffs1, orders1)

    coeffs2 = jet.CScalarVec([
        -1, 3, 3, 1
    ])
    orders2 = jet.CIndexVec([
        0,0,1,0,
        0,0,0,0,
        0,0,0,1,
        2,0,0,0
    ])
    poly2 = jet.CPolyLinkedList(4)
    poly2.batch_add_elements(coeffs2, orders2)

    jet_transport_test = jet.cJetTransport(3, 1, 5)
    # print(jet_transport_test.var_dim)
    # print(jet_transport_test.param_dim)
    # print(jet_transport_test.max_order)
    jet_transport_test.set_sol(0, poly0)
    jet_transport_test.set_sol(1, poly1)
    jet_transport_test.set_sol(2, poly2)

    for dim in range(3):
        print("ode fun %d:"%(dim))
        for order in range(4):
            print("order:",order)
            jet_transport_test.get_sol(dim, order).print_info()
    
    coeffs = jet_transport_test.get_sol_as_scalars(4)
    print(np.reshape(coeffs, (3,-1)))
    for j in range(len(coeffs)//3):
        print(jet_transport_test.get_series_indices(j))

def test_calculate_dvars():
    coeffs0 = jet.CScalarVec([
        -1, 2, 1, 2
    ])
    orders0 = jet.CIndexVec([
        1,0,0,0, 
        0,1,0,0,
        1,0,1,0,
        0,0,0,0
    ])
    poly0 = jet.CPolyLinkedList(4)
    poly0.batch_add_elements(coeffs0, orders0)
    print(poly0.to_str(['x','y','z','lam']))
    print()

    coeffs1 = jet.CScalarVec([
        -1, -2
    ])
    orders1 = jet.CIndexVec([
        0,1,0,0, 
        1,0,0,0
    ])
    poly1 = jet.CPolyLinkedList(4)
    poly1.batch_add_elements(coeffs1, orders1)
    print(poly1.to_str(['x','y','z','lam']))
    print()

    coeffs2 = jet.CScalarVec([
        -1, 3, 3, 1
    ])
    orders2 = jet.CIndexVec([
        0,0,1,0,
        0,0,0,0,
        0,0,0,1,
        2,0,0,0
    ])
    poly2 = jet.CPolyLinkedList(4)
    poly2.batch_add_elements(coeffs2, orders2)
    print(poly2.to_str(['x','y','z','lam']))
    print()

    jet_transport_test = jet.cJetTransport(3, 1, 5)
    # print(jet_transport_test.var_dim)
    # print(jet_transport_test.param_dim)
    # print(jet_transport_test.max_order)
    jet_transport_test.set_sol(0, poly0)
    jet_transport_test.set_sol(1, poly1)
    jet_transport_test.set_sol(2, poly2)

    poly0 = jet.CPolyLinkedList(4)
    poly0.batch_add_elements(jet.CScalarVec([1,1]),jet.CIndexVec([1,0,0,0, 2,0,0,0]))
    poly1 = jet.CPolyLinkedList(4)
    poly1.batch_add_elements(jet.CScalarVec([1]),jet.CIndexVec([0,1,0,0]))
    poly2 = jet.CPolyLinkedList(4)
    poly2.batch_add_elements(jet.CScalarVec([1,1,1]),jet.CIndexVec([0,0,0,0, 0,0,1,0, 1,0,2,0]))
    jet_transport_test.set_ode_fun(0,poly0)
    jet_transport_test.set_ode_fun(1,poly1)
    jet_transport_test.set_ode_fun(2,poly2)

    jet_transport_test.calculate_dvars(3)
    coeffs = jet_transport_test.get_dvars_as_scalars(3)
    jet_transport_test.set_sol_by_scalars(coeffs, 3)
    for j in range(3):
        print(j)
        for order in range(3):
            # jet_transport_test.get_sol(j, order).print_info()
            print(jet_transport_test.get_sol(j, order).to_str(['x','y','z','lam']))
        print()
 
def test_Taylor_method():
    ''' 
        test Taylor expansion of the flow
        example:
            dx/dt = - x + 2*y + x*z
            dy/dt = - y - 2*x
            dz/dt = - z + 3(1 + deltaR) + x^2
    '''

    coeffs0 = jet.CScalarVec([
        -1, 2, 1
    ])
    orders0 = jet.CIndexVec([
        1,0,0,0, 
        0,1,0,0,
        1,0,1,0
    ])
    poly0 = jet.CPolyLinkedList(4)
    poly0.batch_add_elements(coeffs0, orders0)
    print(poly0.to_str(['x','y','z','lam']))

    coeffs1 = jet.CScalarVec([
        -1, -2
    ])
    orders1 = jet.CIndexVec([
        0,1,0,0, 
        1,0,0,0
    ])
    poly1 = jet.CPolyLinkedList(4)
    poly1.batch_add_elements(coeffs1, orders1)
    print(poly1.to_str(['x','y','z','lam']))

    coeffs2 = jet.CScalarVec([
        -1, 3, 3, 1
    ])
    orders2 = jet.CIndexVec([
        0,0,1,0,
        0,0,0,0,
        0,0,0,1,
        2,0,0,0
    ])
    poly2 = jet.CPolyLinkedList(4)
    poly2.batch_add_elements(coeffs2, orders2)
    print(poly2.to_str(['x','y','z','lam']))

    jet_transporter = jet.cJetTransport(3, 1, 6)
    jet_transporter.set_ode_fun(0,poly0)
    jet_transporter.set_ode_fun(1,poly1)
    jet_transporter.set_ode_fun(2,poly2)

    ini_poly0 = jet.CPolyLinkedList(4)
    ini_poly0.batch_add_elements([2, 1], [0,0,0,0,1,0,0,0])
    jet_transporter.set_sol(0, ini_poly0)
    ini_poly1 = jet.CPolyLinkedList(4)
    ini_poly1.batch_add_elements([1, 1], [0,0,0,0,0,1,0,0])
    jet_transporter.set_sol(1, ini_poly1)
    ini_poly2 = jet.CPolyLinkedList(4)
    ini_poly2.batch_add_elements([-1, 1], [0,0,0,0,0,0,1,0])
    jet_transporter.set_sol(2, ini_poly2)
    print()
    for var_id in range(3):
        sol_poly = jet.CPolyLinkedList(4)
        for curr_order in range(3):
            sol_poly.destructive_add_self(jet_transporter.get_sol(var_id, curr_order))
        print(sol_poly.to_str(['x','y','z','lam']))
    print()
    # jet_transporter.set_time_variable()
    # for var_id in range(3):
    #     sol_poly = jet.CPolyLinkedList(5)
    #     sol_poly2 = jet.CPolyLinkedList(4)
    #     for k in range(6):
    #         curr_poly = jet_transporter.get_sol_time(var_id, k)
    #         sol_poly.destructive_add_self(curr_poly)
    #         curr_poly = jet_transporter.get_sol(var_id, k)
    #         sol_poly2.destructive_add_self(curr_poly)
    #     sol_poly.print_info()
    #     sol_poly2.print_info()
    # for j in range(7):
    #     jet_transporter.derivative_iterate()
    #     for var_id in range(3):
    #         sol_poly = jet.CPolyLinkedList(5)
    #         for k in range(6):
    #             curr_poly = jet_transporter.get_sol_time(var_id, k)
    #             sol_poly.destructive_add_self(curr_poly)
    #         print(sol_poly.to_str(['x','y','z','lam','t']))
    #     print()

    jet_transporter.Poincare_projection(3, jet.CScalarVec([0,1,0]))

    for var_id in range(3):
        sol_poly = jet.CPolyLinkedList(5)
        for k in range(6):
            sol_poly.destructive_add_self(jet_transporter.get_sol_time(var_id, k))
        print(sol_poly.to_str(['x','y','z','lam','t']))
    print()
        
    for var_id in range(3):
        sol_poly = jet.CPolyLinkedList(4)
        for curr_order in range(3):
            sol_poly.destructive_add_self(jet_transporter.get_sol(var_id, curr_order))
        print(sol_poly.to_str(['x','y','z','lam']))


if __name__ == '__main__':
    # test_ode_fun()
    # test_conversion_coeffs2()
    # test_calculate_dvars()
    test_Taylor_method()
    