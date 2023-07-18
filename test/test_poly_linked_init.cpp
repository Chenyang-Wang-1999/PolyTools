/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2023-03-27
 * Copyright © Department of Physics, Tsinghua University. All rights reserved
 */

#include "../src/polys.hpp"
#include <iostream>

int main()
{
    PolyLinkedList f(2);
    f.print_info();
    std::cout << "f.dim:" <<  f.dim << ", f.n_terms:" << f.n_terms <<'\n';
    return 0;
}