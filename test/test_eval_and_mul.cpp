/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2022-12-21
 * Copyright © Department of Physics, Tsinghua University.  All rights reserved
 */
#include "test_common.h"
#include <stdio.h>

int main()
{
    IndexType N_terms=4, order=0, dim=2;
    Homogen f(dim, order), g(dim, order), h(dim, order);
    ScalarVec x_val = {Scalar(((double)rand())/RAND_MAX), Scalar(((double)rand())/RAND_MAX)};
    // ScalarVec x_val ={};

    srand(time(NULL));
    random_homogen(N_terms, order, dim, f);
    random_homogen(N_terms, order, dim, g);

    STDOUT << f.term_tree ->eval(x_val) <<'\n';

    STDOUT << "f:\n";
    f.print_info();
    STDOUT << "f value:" << f.eval(x_val) << "\n\n";
    
    STDOUT << "g:\n";
    g.print_info();
    STDOUT << "g value:" << g.eval(x_val) << "\n\n";

    f.add(g, h);
    STDOUT << "eval(f+g):" << h.eval(x_val) << "\n";
    STDOUT << "f_val + g_val:" << f.eval(x_val) + g.eval(x_val) << "\n\n";

    f.subs(g,h);
    STDOUT << "eval(f-g):" << h.eval(x_val) << "\n";
    STDOUT << "f_val - g_val:" << f.eval(x_val) - g.eval(x_val) << "\n\n";

    homogen_multiplication(f, g, h);
    STDOUT << "eval(f*g):" << h.eval(x_val) << "\n";
    STDOUT << "f_val * g_val:" << f.eval(x_val) * g.eval(x_val) << "\n\n";
}