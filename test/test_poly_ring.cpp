/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2023-09-30
 * Copyright © Department of Physics, Tsinghua University. All rights reserved
 */

#include "../src/poly_ring.hpp"
#include "../src/polys.hpp"
#include "test_common.h"

int test_monomial_division()
{
    PolyTools::Monomial a(3.0, {1,2,3}), b(3.0, {2,1,3}), c(0.0, {0,0,0});
    bool result = PolyTools::monomial_division(b, a, c);
    std::cout << result << '\n';
    c.print_info();
    return 0;
}

int test_poly_division()
{
    PolyTools::PolyLinkedList f(2), g(2), h(2);
    f.batch_add_elements(
        {1, -1},
        {3, 1, 0, 0}
    );
    g.batch_add_elements(
        {1, -1},
        {1, 2, 0, 1}
    );
    h.batch_add_elements(
        {1, -1},
        {2, 1, 0, 1}
    );
    // f.update_leading_term();
    // g.update_leading_term();
    std::cout << f.to_str({"x", "y"}) << '\n';
    std::cout << g.to_str({"x", "y"}) << '\n';
    std::cout << h.to_str({"x", "y"}) << '\n';
    // PolyTools::remainder_division(f, g, 1e-12);

    PolyTools::PolyDeque poly_deque;
    poly_deque.push_back(f);
    poly_deque.push_back(g);
    poly_deque.push_back(h);
    // std::cout << f.to_str({"x", "y"}) << '\n';
    PolyTools::simplify_basis(poly_deque, 1e-12);
    std::cout << "Results:\n";
    for(auto poly_ptr:poly_deque.data)
    {
        std::cout << poly_ptr->to_str({"x", "y"}) << '\n';
    }
    return 0;
}

int test_S_poly()
{
    PolyTools::PolyLinkedList f(2), g(2), h(2);
    f.batch_add_elements(
        {1, -1},
        {3, 1, 0, 0}
    );
    g.batch_add_elements(
        {1, -1},
        {1, 2, 0, 1}
    );

    std::cout << "f:" << f.to_str({"x", "y"}) << '\n';
    std::cout << "g:" << g.to_str({"x", "y"}) << '\n';
    PolyTools::get_S_poly(f, g, h);
    std::cout << "h:" << h.to_str({"x", "y"}) << '\n';
    return 0;
}

int test_Groebner()
{
    PolyTools::PolyLinkedList f(2), g(2), h(2);
    f.batch_add_elements(
        {1, -1},
        {3, 1, 0, 0}
    );
    g.batch_add_elements(
        {1, -1},
        {1, 2, 0, 1}
    );

    std::cout << "f:" << f.to_str({"x", "y"}) << '\n';
    std::cout << "g:" << g.to_str({"x", "y"}) << '\n';

    PolyTools::PolyDeque Groebner_basis;
    Groebner_basis.push_back(f);
    Groebner_basis.push_back(g);

    PolyTools::buchberger(Groebner_basis, 1e-12);

    std::cout << "Groebner:\n";
    for(auto poly_ptr : Groebner_basis.data)
    {
        std::cout << poly_ptr -> to_str ({"x", "y"}) << '\n';
    }
    return 0;
}

int main()
{
    srand(time(NULL));
    // return test_monomial_division();
    // return test_poly_division();
    // return test_S_poly();
    return test_Groebner();
}