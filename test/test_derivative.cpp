/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2022-12-21
 * Copyright © Department of Physics, Tsinghua University.  All rights reserved
 */

#include "test_common.h"

int main()
{
    srand(time(NULL));
    // Monomial f = random_term(3, 10);
    // f.print_info();
    // f.derivative(0).print_info();
    // f.derivative(1).print_info();
    // f.derivative(2).print_info();

    Homogen f(3, 10), g(3,9);
    // random_homogen(10, 0, 3, f);
    f.print_info();
    f.derivative(0, g);
    g.print_info();
    f.derivative(1, g);
    g.print_info();
    f.derivative(2, g);
    g.print_info();
}