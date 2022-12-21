/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2022-12-20
 * Copyright © Department of Physics, Tsinghua University.  All rights reserved
 */

#include "../src/polys.h"
#include <cstdlib>
#include <ctime>
#include <complex>

Monomial random_term(IndexType dim, IndexType order)
{
    IndexVec order_var(dim);
    IndexType remaining_order = order;
    for(IndexType var_id = 0; var_id<dim-1; var_id ++)
    {
        IndexType curr_order = rand() % (remaining_order + 1);
        order_var[var_id] = curr_order;
        remaining_order -= curr_order;
    }
    order_var[dim-1] = remaining_order;
    Scalar coeff(rand()/((double) RAND_MAX), rand()/((double) RAND_MAX));

    return Monomial(coeff, order_var);
    
}

void random_homogen(IndexType N_terms, IndexType order, IndexType dim, Homogen & homogen_fun)
{
    IndexVec order_var(dim);
    homogen_fun.reinit(dim, order);
    // set random
    for(IndexType term_id=0; term_id<N_terms; term_id ++)
    {
        // IndexType remaining_order = order;
        // for(IndexType var_id = 0; var_id<dim-1; var_id ++)
        // {
        //     IndexType curr_order = rand() % (remaining_order + 1);
        //     order_var[var_id] = curr_order;
        //     remaining_order -= curr_order;
        // }
        // order_var[dim-1] = remaining_order;
        // Scalar coeff(rand()/((double) RAND_MAX), rand()/((double) RAND_MAX));
        homogen_fun.add_term(random_term(dim, order));
    }
 
}