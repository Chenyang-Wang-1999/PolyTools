/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2022-12-21
 * Copyright © Department of Physics, Tsinghua University.  All rights reserved
 */

#ifndef POLY_UTILS_H
#define POLY_UTILS_H

#include "basic_defs.h"
bool const_sum_next(std::vector<IndexType> &boards, IndexVec & ball_min, IndexVec & ball_max);
void const_sum_boards_to_data(std::vector<IndexType> &boards, std::vector<IndexType> &data);
bool boards_left_most(IndexVec & boards, IndexType n_freeze, IndexVec & ball_min, IndexVec & ball_max, IndexType total_sum);
bool boards_left_most(IndexVec & boards, IndexType n_freeze, IndexVec & ball_min, IndexVec & ball_max);
IndexType curr_right_bound(IndexVec & boards, IndexVec & ball_min, IndexVec & ball_max, IndexType curr_id);


#endif