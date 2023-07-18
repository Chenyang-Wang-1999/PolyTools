'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2023-07-10
Copyright © Department of Physics, Tsinghua University. All rights reserved
'''

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


if __name__ == '__main__':
    main_test_partial_eval()