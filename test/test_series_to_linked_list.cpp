/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2022-12-26
 * Copyright © Department of Physics, Tsinghua University.  All rights reserved
 */

#include "test_common.h"


int main()
{
    srand(time(NULL));
    IndexType dim = 3, max_order = 500, N_terms = 1000000;
    Series poly_series(max_order, dim);
    PolyLinkedList poly_list(dim);

    for(IndexType j = 0; j < N_terms; j++)
    {
        poly_series.add_term(random_term(dim, rand()%max_order));
    }
    STDOUT << "Series:" << '\n';
    // poly_series.print_info();

    series_to_linklist(poly_series, poly_list);
    STDOUT << "Linked list" <<'\n';
    // poly_list.print_info();

}