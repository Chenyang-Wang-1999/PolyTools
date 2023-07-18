/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2022-12-20
 * Copyright © Department of Physics, Tsinghua University. All rights reserved
 */

#ifndef BASIC_DEFS_H
#define BASIC_DEFS_H


#include <complex>
#include <vector>
#include <eigen3/Eigen/Core>
#include <iostream>

namespace PolyTools{

// Define variable type
#ifndef SCALAR_MODE
#define SCALAR_MODE 0
#endif


#if (SCALAR_MODE == 0)
    // complex coeff, complex var
    typedef std::complex<double> Scalar;
    typedef std::complex<double> VarScalar;
#elif (SCALAR_MODE == 1)
    // real coeff, complex var
    typedef double Scalar;
    typedef std::complex<double> VarScalar;
#else
    // real coeff, real var
    typedef double Scalar;
    typedef double VarScalar;
#endif

typedef unsigned int IndexType;
typedef std::vector<Scalar> ScalarVec;
typedef std::vector<VarScalar> VarScalarVec;
typedef std::vector<IndexType> IndexVec; 
typedef enum{GT, EQ, LT, GEQ, LEQ} CompareResult;
typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> EigenMatrixX;
typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> EigenVectorX;
typedef Eigen::Matrix<VarScalar, Eigen::Dynamic, Eigen::Dynamic> EigenVarMatrixX;
typedef Eigen::Matrix<VarScalar, Eigen::Dynamic, 1> EigenVarVectorX;
#define STDOUT std::cout
#define STDERR std::cerr
#define EPS 1e-16
#define EXP_BASE 1
#define ABS_FUN abs
}// namespace PolyTools

#endif