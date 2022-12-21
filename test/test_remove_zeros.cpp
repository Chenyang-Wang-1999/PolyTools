/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2022-12-21
 * Copyright © Department of Physics, Tsinghua University.  All rights reserved
 */

#include "test_common.h"

int main()
{
    Homogen f(3, 2);

    f.add_term(Monomial(Scalar(1.0,0.0),{0,1,1}));
    f.add_term(Monomial(Scalar(1.0,1.0),{1,0,1}));
    f.add_term(Monomial(Scalar(1e-16,0.0),{0,0,2}));
    f.add_term(Monomial(Scalar(0.0,1e-15),{1,1,0}));
    f.add_term(Monomial(Scalar(2e-16,2e-16),{0,2,0}));

    f.print_info();
    f.remove_zeros();
    f.print_info();
    return 0;
}