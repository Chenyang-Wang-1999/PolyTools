'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2022-12-26
Copyright © Department of Physics, Tsinghua University.  All rights reserved
'''

import py_invariant_manifold as pyIM
import numpy as np

def main_test_batch_add():
    N_terms = 20
    dim = 3
    max_order = 5
    coeffs = np.random.rand(N_terms)
    orders = np.random.randint(max_order, size=(N_terms*dim))
    p = pyIM.CPolyLinkedList(dim)
    p.batch_add_elements(coeffs, orders)
    p.print_info()

def main_test_save_and_write():
    N_terms = 20
    dim = 3
    max_order = 5
    coeffs = np.random.rand(N_terms)
    orders = np.random.randint(max_order, size=(N_terms*dim))
    p = pyIM.CPolyLinkedList(dim)
    p.batch_add_elements(coeffs, orders)
    print("p:")
    p.print_info()   
    p.save_to_file("./data/p.mat")
    q = pyIM.poly_load_from_file("./data/p.mat")
    input()
    print("q:")
    q.print_info()

def main_test_add_and_eval():
    N_terms = 20
    dim = 3
    max_order = 5
    coeffs = np.random.rand(N_terms)
    orders = np.random.randint(max_order, size=(N_terms*dim))
    p = pyIM.CPolyLinkedList(dim)
    p.batch_add_elements(coeffs, orders)

    coeffs = np.random.rand(N_terms)
    orders = np.random.randint(max_order, size=(N_terms*dim))
    q = pyIM.CPolyLinkedList(dim)
    q.batch_add_elements(coeffs, orders)

    r = p.subs(q)

    p.print_info()
    q.print_info()
    r.print_info()

    x = np.random.rand(dim, 100000)
    diff_order = pyIM.CIndexVec([0] * dim)
    p_val = p.batch_eval(diff_order, x)
    q_val = q.batch_eval(diff_order, x)
    r_val = r.batch_eval(diff_order, x)
    
    print(p_val)
    print(q_val)
    print(r_val)
    print(np.max(np.abs(r_val - p_val + q_val)))
    
def main_test_scale_var():
    N_terms = 20
    dim = 3
    max_order = 5
    scale_factor = 5.0

    coeffs = np.random.rand(N_terms)
    orders = np.random.randint(max_order, size=(N_terms*dim))
    p = pyIM.CPolyLinkedList(dim)
    p.batch_add_elements(coeffs, orders)   

    q = p.scale_var(pyIM.CScalarVec([scale_factor] * dim))

    x = np.random.rand(dim)
    p_x = p.eval(pyIM.CScalarVec(x))
    p_kx = p.eval(pyIM.CScalarVec(x*scale_factor))
    q_x = q.eval(pyIM.CScalarVec(x))
    print('p(x) = ', p_x)
    print('p(kx) = ', p_kx)
    print('q(x) = ', q_x)



if __name__ == '__main__':
    # main_test_batch_add()
    # main_test_save_and_write()
    # main_test_add_and_eval()
    main_test_scale_var()
