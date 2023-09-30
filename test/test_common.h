/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2022-12-20
 * Copyright © Department of Physics, Tsinghua University. All rights reserved
 */

#include "../src/polys.hpp"
#include <cstdlib>
#include <ctime>
#include <complex>

PolyTools::Monomial random_term(PolyTools::IndexType dim, PolyTools::IndexType order)
{
    PolyTools::IndexVec order_var(dim);
    PolyTools::IndexType remaining_order = order;
    for(PolyTools::IndexType var_id = 0; var_id<dim-1; var_id ++)
    {
        PolyTools::IndexType curr_order = rand() % (remaining_order + 1);
        order_var[var_id] = curr_order;
        remaining_order -= curr_order;
    }
    order_var[dim-1] = remaining_order;
    PolyTools::Scalar coeff(rand()/((double) RAND_MAX), rand()/((double) RAND_MAX));

    return PolyTools::Monomial(coeff, order_var);
    
}

void random_homogen(PolyTools::IndexType N_terms, PolyTools::IndexType order, PolyTools::IndexType dim, PolyTools::Homogen & homogen_fun)
{
    PolyTools::IndexVec order_var(dim);
    homogen_fun.reinit(dim, order);
    // set random
    for(PolyTools::IndexType term_id=0; term_id<N_terms; term_id ++)
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