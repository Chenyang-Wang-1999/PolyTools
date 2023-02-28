/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2022-12-20
 * Copyright © Department of Physics, Tsinghua University. All rights reserved
 */

#include "test_common.h"
#include <stdio.h>

int main()
{
    IndexType N_terms=10, order=5, dim=2;
    Homogen f(dim, order), g(dim, order), h(dim, order);
    random_homogen(N_terms - 5, order, dim, f);
    random_homogen(N_terms, order, dim, g);
    STDOUT << "f:\n";
    f.print_info();
    STDOUT << "g:\n";
    g.print_info();

    PolyMulSweeper sweeper(f, g);

    PolyTerm* curr_term = sweeper.first_term();
    STDOUT << "First term:\n";
    curr_term -> print_info();
    delete curr_term;
    while(! sweeper.is_finished())
    {
        sweeper.print_stack();
        sweeper.print_tick_matrix();
        curr_term = sweeper.next_term();
        if(curr_term == NULL)
        {
            STDOUT << "poped term is NULL\n";
        }
        else
        {
            STDOUT << "poped term:\n";
            curr_term -> print_info();
            delete curr_term;
        }
    }
    delete f;
    delete g;
    return 0;
}