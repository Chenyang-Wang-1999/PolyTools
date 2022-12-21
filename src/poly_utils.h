/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2022-12-21
 * Copyright © Department of Physics, Tsinghua University.  All rights reserved
 */

#ifndef POLY_UTILS_H
#define POLY_UTILS_H

#include "basic_defs.h"
bool const_sum_next(std::vector<IndexType> &sepper, IndexType total_sum);
void const_sum_sepper_to_data(std::vector<IndexType> &sepper, std::vector<IndexType> &data, IndexType total_sum);


#endif