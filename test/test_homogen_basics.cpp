/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2022-12-18
 * Copyright © Department of Physics, Tsinghua University.  All rights reserved
 */

#include "../src/polys.h"
#include <cstdlib>
#include <ctime>
#include <complex>

void random_homogen(IndexType N_terms, IndexType order, IndexType dim, Homogen & homogen_fun)
{
    IndexVec order_var(dim);
    homogen_fun.reinit(dim, order);
    // set random
    for(IndexType term_id=0; term_id<N_terms; term_id ++)
    {
        IndexType remaining_order = order;
        for(IndexType var_id = 0; var_id<dim-1; var_id ++)
        {
            IndexType curr_order = rand() % (remaining_order + 1);
            order_var[var_id] = curr_order;
            remaining_order -= curr_order;
        }
        order_var[dim-1] = remaining_order;
        Scalar coeff(rand()/((double) RAND_MAX), rand()/((double) RAND_MAX));
        homogen_fun.add_term(Monomial(coeff, order_var));
    }
 
}

int main()
{
    // Set total number of terms, total order and dim
    IndexType N_terms=10, order=5, dim=3;
    Homogen f(dim, order), g(dim, order), h(dim, order);

    srand((unsigned int) time(NULL));
    random_homogen(N_terms, order, dim, f);
    random_homogen(N_terms, order, dim, g);

    // STDOUT << "f:\n";
    // f.print_info();
    // STDOUT << '\n';

    // test reinit
    // STDOUT << "------Test Reinit-------\n";
    // f.reinit(dim, order);
    // f.print_info();
    // STDOUT << "\n";

    // test copy
    // STDOUT << "------Test Copy--------\n";
    // f.print_info();
    // g.print_info();
    // STDOUT << "------copy------\n";
    // f.copy(g);
    // g.print_info();
    // STDOUT << '\n';

    // negative:
    // STDOUT << "======Test neg=======\n";
    // STDOUT << "f:\n";
    // f.print_info();

    // f.neg_self();
    // STDOUT << "f after neg_self:\n";
    // f.print_info();

    // f.neg(g);
    // STDOUT << "f after neg:\n";
    // f.print_info();
    // STDOUT << "g after neg:\n";
    // g.print_info();
    // STDOUT << "------Test neg-------\n";

    // // scalar mult
    // STDOUT << "========Test Scalar Mult=========\n";
    // Scalar k(0.0, 1.0);
    // STDOUT << "f:\n";
    // f.print_info();

    // f.scalar_mul_self(k);
    // STDOUT << "f after scalar_mul_self\n";
    // f.print_info();

    // f.scalar_mul(k, g);
    // STDOUT << "f after scalar_mul(k, g)\n";
    // f.print_info();
    // STDOUT << "g after scalar_mul(k, g)\n";
    // g.print_info();

    // STDOUT << "--------Test Scalar Mult---------\n";

    // STDOUT << "========Test Add Self=========\n";
    // STDOUT << "f\n";
    // f.print_info();
    // STDOUT << "g\n";
    // g.print_info();

    // f.add_self(g);
    // STDOUT << "f after add_self(g)\n";
    // f.print_info();
    // STDOUT << "g after add_self(g)\n";
    // g.print_info();

    // f.destructive_add_self(g);
    // STDOUT << "f after destructive_add_self(g)\n";
    // f.print_info();
    // STDOUT << "g after destructive_add_self(g)\n";
    // g.print_info();
    // STDOUT << "--------Test Add Self---------\n";

    STDOUT << "========Test Add=========\n";
    STDOUT << "f\n";
    f.print_info();
    STDOUT << "g\n";
    g.print_info();

    f.add(g, h);
    STDOUT << "f after f.add(g,h)\n";
    f.print_info();
    STDOUT << "g after f.add(g,h)\n";
    g.print_info();
    STDOUT << "h after f.add(g,h)\n";
    h.print_info();

    f.destructive_add(g, h);
    STDOUT << "f after f.destructive_add(g,h)\n";
    f.print_info();
    STDOUT << "g after f.destructive_add(g,h)\n";
    g.print_info();
    STDOUT << "h after f.destructive_add(g,h)\n";
    h.print_info();

    STDOUT << "--------Test Add---------\n";

    // STDOUT << "========Test Subs Self=========\n";
    // STDOUT << "f\n";
    // f.print_info();
    // STDOUT << "g\n";
    // g.print_info();

    // f.subs_self(g);
    // STDOUT << "f after subs_self(g)\n";
    // f.print_info();
    // STDOUT << "g after subs_self(g)\n";
    // g.print_info();

    // f.destructive_subs_self(g);
    // STDOUT << "f after destructive_subs_self(g)\n";
    // f.print_info();
    // STDOUT << "g after destructive_subs_self(g)\n";
    // g.print_info();
    // STDOUT << "--------Test Subs Self---------\n";

    // STDOUT << "========Test Subs=========\n";
    // STDOUT << "f\n";
    // f.print_info();
    // STDOUT << "g\n";
    // g.print_info();

    // f.subs(g, h);
    // STDOUT << "f after f.subs(g,h)\n";
    // f.print_info();
    // STDOUT << "g after f.subs(g,h)\n";
    // g.print_info();
    // STDOUT << "h after f.subs(g,h)\n";
    // h.print_info();

    // f.destructive_subs(g, h);
    // STDOUT << "f after f.destructive_subs(g,h)\n";
    // f.print_info();
    // STDOUT << "g after f.destructive_subs(g,h)\n";
    // g.print_info();
    // STDOUT << "h after f.destructive_subs(g,h)\n";
    // h.print_info();

    // STDOUT << "--------Test Subs---------\n";

    // STDOUT << "=========Test Insert Zero==========\n";
    // f.neg(g);
    // f.print_info();
    // g.print_info();

    // PolyTerm * term_ptr = g.pop_first_term();
    // f.add_term(term_ptr);
    // f.print_info();
    // g.print_info();


    // term_ptr = g.pop_first_term();
    // term_ptr -> ~PolyTerm();
    // term_ptr = g.pop_first_term();
    // f.add_term(term_ptr);
    // f.print_info();
    // g.print_info();
    // STDOUT << "---------Test Insert Zero----------\n";
}