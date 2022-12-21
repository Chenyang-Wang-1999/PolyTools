/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2022-12-21
 * Copyright © Department of Physics, Tsinghua University.  All rights reserved
 */

#include "test_common.h"
#include <vector>

int main()
{
    IndexType N_terms=4, order=5, dim=2, n_seq = 5;
    std::vector<Homogen *> homog_seq(n_seq);
    std::vector<Scalar> homog_val(n_seq);
    ScalarVec x_val = {Scalar(((double)rand())/RAND_MAX), Scalar(((double)rand())/RAND_MAX)};
    // ScalarVec x_val ={};

    srand(time(NULL));
    // Initialize
    for(IndexType homog_id = 0; homog_id < n_seq ; homog_id ++)
    {
        homog_seq[homog_id] = new Homogen(dim, order);
        random_homogen(N_terms, order,dim, *(homog_seq[homog_id]));
        homog_val[homog_id] = homog_seq[homog_id] -> eval(x_val);
    }

    Homogen result(dim, order * n_seq);
    homogen_mul_seq(homog_seq, order*n_seq, result);
    Scalar res_val(1.0);
    for(IndexType homog_id = 0; homog_id < n_seq; homog_id ++)
    {
        STDOUT << "homog "<<homog_id <<": " << homog_val[homog_id] <<'\n';
        res_val *= homog_val[homog_id];
    }

    STDOUT << "mul seq:\n";
    result.print_info();
    STDOUT << "product: " << res_val <<'\n';
    STDOUT << "val:" << result.eval(x_val) <<'\n';

    for(IndexType homog_id = 0; homog_id < n_seq ; homog_id ++)
    {
        homog_seq[homog_id] -> ~Homogen();
    }
    return 0;
}