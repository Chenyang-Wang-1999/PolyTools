/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2022-12-21
 * Copyright © Department of Physics, Tsinghua University.  All rights reserved
 */

#include "poly_utils.h"

// tools for constant sum indices
    /* algorithm: consider k dots seperated by k-1 boards
            |oo||o
        sepper = {0,2,2}, data = {0,2,0,1}
    */ 
void const_sum_sepper_to_data(std::vector<IndexType> &sepper, std::vector<IndexType> &data, IndexType total_sum)
{
    data.clear();
    data.resize(sepper.size()+1);
    data[0] = sepper[0];
    for(IndexType sepper_id = 1; sepper_id < sepper.size(); sepper_id ++)
    {
        data[sepper_id] = sepper[sepper_id] - sepper[sepper_id-1];
    }
    data[sepper.size()] = total_sum - sepper[sepper.size() -1];
}

bool const_sum_next(std::vector<IndexType> &sepper, IndexType total_sum)
{
    bool ptr_hanged = true;
    IndexType curr_id = sepper.size() - 1;
    do {
        sepper[curr_id] ++;
        if(sepper[curr_id] > total_sum)
        {
            if(curr_id > 0)
            {
                curr_id --;
            }
            else
            {
                return true;
            }
        }
        else
        {
            ptr_hanged = false;
        }
    }while (ptr_hanged);
    for(IndexType remaining_id = curr_id + 1; remaining_id < sepper.size(); remaining_id ++)
    {
        sepper[remaining_id] = sepper[curr_id];
    }

    return false;
}