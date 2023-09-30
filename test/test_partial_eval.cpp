/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2023-09-01
 * Copyright © Department of Physics, Tsinghua University. All rights reserved
 */

#include "test_common.h"

int main()
{
    PolyTools::PolyLinkedList poly(3);
    size_t N_terms = 10;

    for(size_t j = 0; j < N_terms; j ++)
    {
        poly.add_term(random_term(3, 5));
    }

    std::cout << poly.to_str({"x", "y", "z"}) << '\n';

    PolyTools::PolyLinkedList poly_another(3);
    poly.scale_var({1, 1, PolyTools::Scalar(0, 1)}, poly_another);
    std::cout << poly_another.to_str({"x", "y", "z"}) << '\n';
    PolyTools::PolyLinkedList poly_2(2);
    PolyTools::ScalarVec var_val = {PolyTools::Scalar(0,1)};
    PolyTools::IndexVec eval_dim = {0}, dof_map = {1,2};
    poly.partial_eval(var_val, eval_dim, dof_map, &poly_2);
    std::cout << poly_2.to_str({"x", "y"}) << '\n';
}