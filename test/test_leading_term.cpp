/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2023-09-30
 * Copyright © Department of Physics, Tsinghua University. All rights reserved
 */

#include "test_common.h"

int test_leading_term()
{
    PolyTools::IndexType N_vars = 3, N_terms = 10;
    PolyTools::PolyLinkedList f(N_vars);
    random_poly(N_terms, 30, N_vars, f);
    std::cout << f.to_str({"x", "y", "z"}) << '\n';
    // f.update_leading_term();
    if(f.leading_term != NULL)
    {
        std::cout << f.leading_term->coeff << '\n';
        for(auto j = 0; j < N_vars ; j++)
        {
            std::cout << f.leading_term->var_order(j) <<',';
        }
        std::cout << '\n';
    }
    else
    {
        std::cout << "NULL\n";
    }
    return 0;
}

int main()
{
    srand(time(NULL));
    return test_leading_term();
}