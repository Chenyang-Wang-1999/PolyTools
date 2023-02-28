/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2022-12-21
 * Copyright © Department of Physics, Tsinghua University. All rights reserved
 */

#include "test_common.h"

int main()
{
    IndexType N_terms=4, order=1, dim=2, n_seq = 3;
    IndexType K_max =20;
    IndexType k_target = 0;
    std::vector<Series*> series_ptr(n_seq);

    // initialize
    for(auto it = series_ptr.begin(); it != series_ptr.end(); it++)
    {
        (*it) = new Series(K_max, dim);
    }
    series_ptr[0] -> add_term(Monomial(Scalar(1.0),{1,0}));
    series_ptr[0] -> add_term(Monomial(Scalar(1.0),{0,1}));
    series_ptr[1] -> add_term(Monomial(Scalar(1.0),{1,0}));
    series_ptr[1] -> add_term(Monomial(Scalar(-1.0),{0,1}));
    series_ptr[2] -> add_term(Monomial(Scalar(1.0),{1,0}));
    series_ptr[2] -> add_term(Monomial(Scalar(2.0),{0,2}));

    // Monomial f(Scalar(1.0,0.0), {0,0,1});
    Series f(20, 3);
    f.add_term(Monomial(Scalar(3.0,2.0),{2,1,1}));
    f.add_term(Monomial(Scalar(-2.0,0.0),{1,0,1}));
    Homogen res(2, 4);

    for (k_target = 0; k_target < 6; k_target++)
    {
        STDOUT << "k target:" << k_target << '\n';
        series_comp(f, series_ptr, k_target, res);
        // term_comp(f, series_ptr, k_target, res);
        res.print_info();
    }

    // destroy
    for(auto it = series_ptr.begin(); it != series_ptr.end(); it++)
    {
        delete (*it);
    }
    
    return 0;   
}