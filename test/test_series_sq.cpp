/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2022-12-23
 * Copyright © Department of Physics, Tsinghua University.  All rights reserved
 */

#include "test_common.h"

int main()
{
    srand(time(NULL));

    IndexType dim = 4, Kmax = 100, N_terms = 1000;
    Series x(Kmax, dim);
    Homogen xsq(dim, 0), xsq_mul(dim,0);

    // randomly generate x
    for(IndexType j = 0; j < N_terms; j++)
    {
        x.add_term(random_term(dim, rand()%Kmax));
    }

    x.print_info();

    // calculate 
    for(IndexType k = 0; k<Kmax; k++)
    {
        STDOUT << "Current k:\n";
        series_mul(x, x, k, xsq_mul);
        // xsq_mul.print_info();
        series_renew_xqr_term(x, xsq, k);        
        // xsq.print_info();
        xsq.destructive_subs_self(xsq_mul);
        xsq.print_info();
    }

    return 0;
}