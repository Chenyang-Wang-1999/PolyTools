/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2022-12-21
 * Copyright © Department of Physics, Tsinghua University.  All rights reserved
 */

#include "../src/poly_utils.h"
#include <vector>
#include <iostream>

int main()
{
    IndexType total_sum = 50;
    IndexVec sep_vec = {0,0,0}, index_vec(4);

    do
    {
        const_sum_sepper_to_data(sep_vec, index_vec, total_sum);
        for(auto it=index_vec.begin(); it != index_vec.end(); ++it)
        {
            STDOUT << *it << ',';
        }
        STDOUT << '\n';
    }while(! const_sum_next(sep_vec, total_sum));
    return 0;
}