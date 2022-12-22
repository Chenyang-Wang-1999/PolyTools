/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2022-12-21
 * Copyright © Department of Physics, Tsinghua University.  All rights reserved
 */

#include "test_common.h"

int main()
{
    Series f(30, 3);
    SeriesVec g(3, 3, 30);

    f.add_term(Monomial(Scalar(3.0, 2.0), {2, 1, 1}));
    f.add_term(Monomial(Scalar(0.0, -2.0), {1, 0, 1}));
    f.add_term(Monomial(Scalar(3.0, 0.0), {0, 2, 2}));
    // f.add_term(Monomial(Scalar(1.0, 0.0), {0, 0, 1}));

    g.series_vec[0] -> add_term(Monomial(Scalar(1.0,0.0), {1,0,0}));
    g.series_vec[0] -> add_term(Monomial(Scalar(1.0,0.0), {0,1,0}));
    g.series_vec[1] -> add_term(Monomial(Scalar(1.0,0.0), {1,0,0}));
    g.series_vec[1] -> add_term(Monomial(Scalar(-1.0,0.0), {0,0,1}));
    g.series_vec[2] -> add_term(Monomial(Scalar(1.0,0.0), {0,1,0}));
    g.series_vec[2] -> add_term(Monomial(Scalar(2.0,0.0), {0,0,2}));

    for(IndexType k = 0; k < 10; k++)
    {
        Homogen h(3, k);
        STDOUT << "k:" << k << '\n';
        series_DW_dot_f(f, g.series_vec, k, h);
        h.print_info();
    }
}