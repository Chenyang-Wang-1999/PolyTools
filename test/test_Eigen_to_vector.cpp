/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2022-12-22
 * Copyright © Department of Physics, Tsinghua University.  All rights reserved
 */

#include "../src/basic_defs.h"
#include <vector>
#include <eigen3/Eigen/Core>
#include <iostream>

int main()
{
    Eigen::MatrixXd mat = Eigen::MatrixXd::Random(5,7);
    std::cout << mat << '\n';

    for(unsigned int j = 0; j < mat.cols(); ++j)
    {
        std::vector<double> double_vec(mat.col(j).data(), mat.col(j).data() + mat.rows());
        std::cout << "col:" << j << '\n';
        for(unsigned int k = 0; k < mat.rows(); ++k)
        {
            std::cout << double_vec[k] << '\n';
        }
    }
    return 0;
}