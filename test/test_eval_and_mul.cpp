/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2022-12-21
 * Copyright © Department of Physics, Tsinghua University. All rights reserved
 */
#include "test_common.h"
#include <stdio.h>

int main()
{
    IndexType N_terms=1000, order=30, dim=2, max_diff_order=4;
    Homogen f(dim, order), g(dim, order), h(dim, order);
    ScalarVec x_val = {Scalar(((double)rand())/RAND_MAX), 
                        Scalar(((double)rand())/RAND_MAX),
                        Scalar(((double)rand())/RAND_MAX)};
    IndexVec diff_order = {rand()%(max_diff_order), rand()%(max_diff_order)};
    // ScalarVec x_val ={};

    srand(time(NULL));
    random_homogen(N_terms, order, dim, f);
    random_homogen(N_terms, order, dim, g);

    STDOUT << "x:" << x_val[0] << ", " << x_val[1] <<'\n';
    STDOUT << "diff order:" << diff_order[0] << ", " << diff_order[1] <<'\n';

    // f.print_info();
    STDOUT << "f eval new:" << f.eval_diff(diff_order, x_val) << '\n';
    
    Scalar fval = Scalar(0.0);
    PolyTerm * curr_term = f.term_tree;
    while(curr_term != NULL)
    {
        Monomial curr_monom = curr_term ->copy();
        for (IndexType var_id = 0; var_id < dim; var_id++)
        {
            for (IndexType j = 0; j < diff_order[var_id]; j++)
            {
                curr_monom = curr_monom.derivative(var_id);
            }
        }
        // fval += curr_term -> eval(x_val);
        fval += curr_monom.eval(x_val);
        curr_term = curr_term -> next;
    }
    STDOUT << "f eval old:" << fval << '\n';
    // STDOUT << f.term_tree ->eval(x_val) <<'\n';

    // STDOUT << "f:\n";
    // f.print_info();
    // STDOUT << "f value:" << f.eval(x_val) << "\n\n";
    
    // STDOUT << "g:\n";
    // g.print_info();
    // STDOUT << "g value:" << g.eval(x_val) << "\n\n";

    // f.add(g, h);
    // // h.print_info();
    // STDOUT << "eval(f+g):" << h.eval(x_val) << "\n";
    // STDOUT << "f_val + g_val:" << f.eval(x_val) + g.eval(x_val) << "\n\n";

    // f.subs(g,h);
    // // h.print_info();
    // STDOUT << "eval(f-g):" << h.eval(x_val) << "\n";
    // STDOUT << "f_val - g_val:" << f.eval(x_val) - g.eval(x_val) << "\n\n";

    // homogen_multiplication(f, g, h);
    // STDOUT << "eval(f*g):" << h.eval(x_val) << "\n";
    // STDOUT << "f_val * g_val:" << f.eval(x_val) * g.eval(x_val) << "\n\n";
}