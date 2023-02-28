/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2023-02-28
 * Copyright © Department of Physics, Tsinghua University. All rights reserved
 */

#include "test_common.h"

int main()
{
     // Set total number of terms, total order and dim
    IndexType N_terms=100, order=50, dim=2;
    Homogen f(dim, order), g(dim, order), h(dim, order);

    srand((unsigned int) time(NULL));

    // STDOUT << "f:\n";
    // f.print_info();
    // STDOUT << '\n';

    // test reinit
    for(unsigned int i = 0; i < 1000; i++)
    // while(true)
    {
        // random_homogen(N_terms, order, dim, f);
        // random_homogen(N_terms, order, dim, g);
        // f.add(g, h);
        for(unsigned int j = 0; j < 10; j++)
        {
            PolyTerm * new_term_f = new PolyTerm(3.0, {j, 50-j});
            PolyTerm * new_term_g = new PolyTerm(3.0, {j, 50-j});
            f.add_term(new_term_f);
            g.add_term(new_term_g);
        }
        f.reinit(dim, order);
        g.reinit(dim, order);
    }
}