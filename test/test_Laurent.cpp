/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2024-01-20
 * Copyright © Department of Physics, Tsinghua University. All rights reserved
 */

#include "../src/polys.hpp"
#include "test_common.h"

using namespace PolyTools;

int test_Laurent_derivative()
{
    srand(time(NULL));
    IndexType dim = 3, N_terms = 100, max_deg = 5;
    ScalarVec coeffs(N_terms);
    IndexVec degrees(dim * N_terms);

    // generate random Laurent
    PolyLinkedList f_num(dim);
    random_poly(N_terms, max_deg, dim, f_num);
    Monomial f_denom = random_term(dim, max_deg);
    Laurent f(dim);
    f.set_Laurent(f_num, f_denom);

    // calculate df, df(x)
    Laurent df(dim);
    IndexType diff_dim = 1;
    f.derivative(diff_dim, df);

    // generate x, dx
    VarScalarVec x(dim);
    Scalar dx((double)rand()/RAND_MAX, (double)rand()/RAND_MAX);
    dx *= 1e-3;
    for(auto x_it = x.begin(); x_it != x.end(); x_it++)
    {
        (*x_it) = Scalar(((double)rand())/RAND_MAX,((double)rand())/RAND_MAX) + 0.001;
    }


    // calculate f(x), f(x + dx)
    VarScalar fx = f.eval(x), dfx = df.eval(x);
    x[diff_dim] += dx;
    VarScalar fdx = f.eval(x);

    // compare f(x+dx) - f(x) with df(x) * dx
    std::cout << "x:";
    for(auto x_it : x)
    {
        std::cout << x_it << ',';
    }
    std::cout << '\n';

    std::cout << "f(x):" << fx << '\n';
    std::cout << "f(x+dx):" << fdx << '\n';
    std::cout << "f(x+dx) - f(x)" << fdx - fx << '\n';
    std::cout << "df*dx" << dfx * dx << '\n';
    std::cout << "error:" << fdx - fx - dfx*dx << '\n';
    return 0;
}

int test_Laurent_partial_eval()
{
    srand(time(NULL));
    IndexType dim = 3, N_terms = 100, max_deg = 5;
    ScalarVec coeffs(N_terms);
    IndexVec degrees(dim * N_terms);

    // generate random Laurent
    PolyLinkedList f_num(dim);
    random_poly(N_terms, max_deg, dim, f_num);
    Monomial f_denom = random_term(dim, max_deg);
    Laurent f(dim);
    f.set_Laurent(f_num, f_denom);

    // generate x, dx
    VarScalarVec x(dim);
    for(auto x_it = x.begin(); x_it != x.end(); x_it++)
    {
        (*x_it) = Scalar(((double)rand())/RAND_MAX,((double)rand())/RAND_MAX) + 0.001;
    }

    // generate f1
    Laurent f1(dim-1);
    ScalarVec vals = {x[1]};
    IndexVec dofs = {1}, dof_map = {0,2};
    f.partial_eval(vals, dofs, dof_map, f1);

    // evaluation
    std::cout << f.eval(x) << '\n';
    std::cout << f1.eval({x[0], x[2]}) << '\n';

    return 0;
}

int main()
{
    // return test_Laurent_derivative();
    return test_Laurent_partial_eval();
}