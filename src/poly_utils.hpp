/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2022-12-21
 * Copyright © Department of Physics, Tsinghua University. All rights reserved
 */

#ifndef POLY_UTILS_H
#define POLY_UTILS_H

#include "basic_defs.h"

/* algorithm: consider k dots seperated by k-1 boards
        |oo||o
    boards = {0,0,2,2,3}, data = {0,2,0,1}
*/ 
void const_sum_boards_to_data(std::vector<IndexType> &boards, std::vector<IndexType> &data)
{
    // data.clear();
    // data.resize(boards.size()-1);
    if(data.size() != boards.size() - 1)
    {
        data.clear();
        data.resize(boards.size() - 1);
    }
    for(IndexType boards_id = 0; boards_id < boards.size()-1; boards_id ++)
    {
        data[boards_id] = boards[boards_id + 1] - boards[boards_id];
    }
}

bool boards_left_most(IndexVec & boards, IndexType n_freeze, IndexVec & ball_min, IndexVec & ball_max)
{
    assert(boards.size() == ball_min.size() + 1);
    // sweep from left to right
    for(IndexType board_id = n_freeze + 1; board_id < boards.size() - 1; ++board_id)
    {
        boards[board_id] = boards[board_id - 1] + ball_min[board_id - 1];
    }

    if((boards[boards.size() - 1]) < ball_min[boards.size() - 2] + (boards[boards.size() - 2]))
    {
        return false;
    }

    // sweep from right to left
    for(IndexType board_id = boards.size() - 2; board_id > n_freeze; board_id -- )
    {
        if(boards[board_id + 1] >= boards[board_id] + ball_max[board_id])
        {
            boards[board_id] = boards[board_id + 1] - ball_max[board_id] + 1;
        }
        else
        {
            return true;
        }
    }

    if((boards[n_freeze+1])>= ball_max[n_freeze] + boards[n_freeze])
    {
        return false;
    }
    else
    {
        return true;
    }
}

bool boards_left_most(IndexVec & boards, IndexType n_freeze, IndexVec & ball_min, IndexVec & ball_max, IndexType total_sum)
{
    // check whether ball_min >= ball_max
    for(IndexType ball_id=0; ball_id < ball_min.size(); ball_id ++)
    {
        assert(ball_min[ball_id] < ball_max[ball_id]);
    }

    // initialize
    if(boards.size() != ball_min.size() + 1)
    {
        if(n_freeze)
        {
            STDERR << "Warning: resize boards\n";
            n_freeze = 0;
        }
        boards.clear();
        boards.resize(ball_min.size() + 1);
    }

    // set left and right boards
    boards[0] = 0; boards[boards.size() - 1] = total_sum;
    return boards_left_most(boards, n_freeze, ball_min, ball_max);
}

IndexType curr_right_bound(IndexVec & boards, IndexVec & ball_min, IndexVec & ball_max, IndexType curr_id)
{
    IndexType right_bound_1 = boards[curr_id - 1] + ball_max[curr_id - 1];
    IndexType right_bound_2 = boards[curr_id + 1] - ball_min[curr_id] + 1;
    return (right_bound_1 > right_bound_2)? right_bound_2:right_bound_1;
}

bool const_sum_next(std::vector<IndexType> &boards, IndexVec & ball_min, IndexVec & ball_max)
{
    bool ptr_hanged = true;
    IndexType curr_id = boards.size() - 2;
    do {
        boards[curr_id] ++;
        if(boards[curr_id] >= (boards[curr_id - 1] + ball_max[curr_id - 1]))
        {
            curr_id --;
            if(curr_id == 0)
            {
                return true;
            }
        }
        else
        {
            bool reinit_successful = boards_left_most(boards, curr_id, ball_min, ball_max);
            if(reinit_successful)
            {
                ptr_hanged = false;
            }
            else
            {
                curr_id --;
                if(curr_id == 0)
                {
                    return true;
                }
            }
        }
    }while (ptr_hanged);

    return false;
}

#endif