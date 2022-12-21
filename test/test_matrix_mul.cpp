/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2022-12-21
 * Copyright © Department of Physics, Tsinghua University.  All rights reserved
 */

#include "test_common.h"
#include <eigen3/Eigen/Core>

int main()
{
    srand(time(NULL));
    Eigen::MatrixXcd A = Eigen::MatrixXcd::Random(3,2);
    Eigen::VectorXcd v(2);
    Eigen::VectorXcd x = Eigen::VectorXcd::Random(3);
    ScalarVec x_vec(3);
    for(IndexType i = 0; i<3; i++)
    {
        x_vec[i] = x(i);
    }

    // create series vector
    SeriesVec input_series(3, 2, 20), output_series(3,3,20);
    for(IndexType i = 0; i < 2; i++)
    {
        for(IndexType j = 0; j < 5; j++)
        {
            IndexType order = rand() % 10;
            Homogen * g_ptr = new Homogen(3, order);
            random_homogen(10, order, 3, *g_ptr);
            input_series.series_vec[i] -> destructive_add_homogen(*g_ptr);
            g_ptr->~Homogen();
        }
        v(i) = input_series.series_vec[i] -> eval(x_vec);
    }
    STDOUT << "v vec\n";
    STDOUT << v <<'\n';

    // calculate matrix mult
    scalar_matrix_mul(A, input_series, output_series);
    STDOUT << "A * v:\n";
    STDOUT << A*v << '\n';

    STDOUT << "series product:\n";
    for(IndexType i = 0; i<3; i++)
    {
        STDOUT << output_series.series_vec[i] -> eval(x_vec) <<'\n';
    }
    return 0;
}