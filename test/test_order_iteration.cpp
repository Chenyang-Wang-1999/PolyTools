/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2022-12-21
 * Copyright © Department of Physics, Tsinghua University. All rights reserved
 */

#include "../src/poly_utils.h"
#include <vector>
#include <iostream>

int main()
{
    IndexType total_sum = 0;
    IndexVec sep_vec = {0,0,0,0,0}, index_vec(4);
    IndexVec ball_min = {0,0,0,1}, ball_max = {30, 30, 30, 30};

    boards_left_most(sep_vec, 0, ball_min, ball_max, total_sum);

    do
    {
        const_sum_boards_to_data(sep_vec, index_vec);
        for(auto it=index_vec.begin(); it != index_vec.end(); ++it)
        {
            STDOUT << *it << ',';
        }
        STDOUT << '\n';
    }while(! const_sum_next(sep_vec, ball_min, ball_max));
    return 0;
}