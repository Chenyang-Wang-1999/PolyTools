/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2022-12-20
 * Copyright © Department of Physics, Tsinghua University.  All rights reserved
 */

#ifndef BASIC_DEFS_H
#define BASIC_DEFS_H

#include <complex>
#include <vector>
typedef std::complex<double> Scalar;
typedef unsigned int IndexType;
typedef std::vector<Scalar> ScalarVec;
typedef std::vector<IndexType> IndexVec; 
typedef enum{GT, EQ, LT, GEQ, LEQ} CompareResult;
#define STDOUT std::cout
#define EPS 1e-16
#define ABS_FUN abs

#endif