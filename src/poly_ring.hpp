/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2023-09-30
 * Copyright © Department of Physics, Tsinghua University. All rights reserved
 */

#include "polys.hpp"
#include <deque>

namespace PolyTools{

/* PolyDeque object */
class PolyDeque
{
public:
    std::deque<PolyLinkedList*> data;
    void push_back(PolyLinkedList & new_poly)
    {
        PolyLinkedList * new_poly_ptr = new PolyLinkedList(new_poly.dim);
        new_poly.copy_to(*new_poly_ptr);
        data.push_back(new_poly_ptr);
    }
    void push_front(PolyLinkedList & new_poly)
    {
        PolyLinkedList * new_poly_ptr = new PolyLinkedList(new_poly.dim);
        new_poly.copy_to(*new_poly_ptr);
        data.push_front(new_poly_ptr);
    }
    ~PolyDeque()
    {
        for(auto poly_ptr:data)
        {
            delete poly_ptr;
        }
    }
};

/* Divisions */

// division of monomials
bool monomial_division(Monomial & dividend, Monomial & divisor, Monomial & quotient)
{
    assert(dividend.dim == divisor.dim);
    assert(dividend.dim == quotient.dim);
    // assert(ABS_FUN(divisor.coeff) > EPS);

    // 1. check each dim 
    for(IndexType var_id = 0; var_id < dividend.dim ; var_id ++)
    {
        if((dividend.var_order(var_id)) >= (divisor.var_order(var_id)))
        {
            quotient.var_order(var_id) = dividend.var_order(var_id) - divisor.var_order(var_id);
        }
        else
        {
            return false;
        }
    }

    // 2. coefficients
    quotient.coeff = dividend.coeff / divisor.coeff;

    // 3. total order
    quotient.order = dividend.order - divisor.order;

    return true;
}

/* 
    remainder of a polynomial by another polynomial 
    assumption: the leading term of the divisor is updated
*/ 
bool remainder_division(PolyLinkedList & dividend, PolyLinkedList & divisor, double tol)
{
    // For each term in the dividend, check whether the term is divided by LT(divisor)
    // If dividend is divided by the divisor, subtract the divisor * quotient from the dividend, and move the term pointer to the begining
    PolyTerm * curr_ptr = dividend.term_tree;
    Monomial quo(0.0, IndexVec(dividend.dim));
    bool zero_quotient = true;
    if(! divisor.leading_term_flag)
    {
        divisor.update_leading_term();
    }
    while (curr_ptr != NULL)
    {
        if(monomial_division(*curr_ptr, *(divisor.leading_term), quo))
        {
            PolyLinkedList new_poly(divisor.dim);
            divisor.monomial_mul(quo, new_poly);
            dividend.destructive_subs_self(new_poly);
            dividend.remove_zeros_with_tol(tol, 1.0);
            curr_ptr = dividend.term_tree;

            zero_quotient = false;
        }
        else
        {
            curr_ptr = curr_ptr -> next;
        }
    }

    // // update the leading term of the dividend
    // dividend.update_leading_term();

    return zero_quotient;
}

/* 
    remainder of a polynomial by another group of polynomials 
    assumption: the leading term of all the divisor is updated
*/ 
bool remainder_division(PolyLinkedList & dividend, std::deque<PolyLinkedList*> & divisor, 
        double tol)
{
    bool nothing_changed = true;
    bool curr_turn_changed;
    do
    {
        curr_turn_changed = false;
        for(auto curr_iter = divisor.begin(); curr_iter != divisor.end(); curr_iter++)
        {
            bool div_result = remainder_division(dividend, *(*curr_iter), tol);
            if(! div_result)
            {
                nothing_changed = false;
                curr_turn_changed = true;
            }
        }
    }while (curr_turn_changed);
    

    return nothing_changed;
}

/* Groebner */
bool simplify_basis(PolyDeque & poly_basis, double tol)
{
    // divide the polynomial by others, until no poly is divisible by the others
    bool nothing_changed, basis_simplified = true;
    do
    {
        IndexType n_polys = poly_basis.data.size();
        nothing_changed = true;
        for(IndexType j = 0; j < n_polys; j++)
        {
            PolyLinkedList * curr_poly = poly_basis.data.front();            
            poly_basis.data.pop_front();
            bool div_result = remainder_division(*curr_poly, poly_basis.data, tol);
            if(!div_result) // something has changed
            {
                nothing_changed = false;
                basis_simplified = false;
                if(curr_poly->term_tree == NULL)
                {
                    // empty
                    delete curr_poly;
                }
                else
                {
                    curr_poly->update_leading_term();
                    poly_basis.data.push_back(curr_poly);
                }
            }
            else
            {
                poly_basis.data.push_back(curr_poly);
            }
        }
    } while (! nothing_changed);
    
    return basis_simplified;
}

void get_S_poly(PolyLinkedList & f, PolyLinkedList & g, PolyLinkedList & h)
{
    assert(f.dim == g.dim);
    h.reinit(f.dim);
    // 1. update leading term
    if(! f.leading_term_flag)
    {
        f.update_leading_term();
    }
    if(! g.leading_term_flag)
    {
        g.update_leading_term();
    }

    // 2. get the cross term
    IndexVec cross_term_order(f.dim);
    for(IndexType j = 0; j<f.dim; j++)
    {
        cross_term_order[j] = (f.leading_term->var_order(j) > g.leading_term->var_order(j))? 
                    f.leading_term->var_order(j): g.leading_term->var_order(j);
    }
    Monomial cross_term(1.0, cross_term_order);

    // 3. get the multiplier of f and g
    Monomial f_mul(0.0, IndexVec(f.dim)), g_mul(0.0, IndexVec(f.dim));
    monomial_division(cross_term, *(f.leading_term), f_mul);
    monomial_division(cross_term, *(g.leading_term), g_mul);
    g_mul.coeff *= -1;
    PolyLinkedList middle_poly(f.dim);
    f.monomial_mul(f_mul, h);
    g.monomial_mul(g_mul, middle_poly);
    h.destructive_add_self(middle_poly);
}

// buchberger algorithm to solve the Groebner basis
void buchberger(PolyDeque & poly_basis, double tol)
{
    std::deque<PolyLinkedList*> last_basis;
    if(poly_basis.data.size() < 2)
    {
        // no poly or one poly
        return;
    }
    bool basis_closed = true;
    do
    {
        basis_closed = true;

        simplify_basis(poly_basis, tol);

        last_basis = poly_basis.data;
        for(auto poly_ptr1 = last_basis.begin(); poly_ptr1 + 1 != last_basis.end(); poly_ptr1 ++)
        {
            for(auto poly_ptr2 = poly_ptr1 + 1; poly_ptr2 != last_basis.end(); poly_ptr2 ++)
            {
                PolyLinkedList * S_poly_ptr = new PolyLinkedList(poly_basis.data[0]->dim);

                // 1. calculate S_poly
                get_S_poly(*(*poly_ptr1), *(*poly_ptr2), *S_poly_ptr);

                // 2. simplify the S_poly under the basis of the new poly basis
                remainder_division(*S_poly_ptr, poly_basis.data, tol);

                // 3. if S_poly != 0, then basis_closed = false, append the S_poly into poly_basis. Else, delete the S_poly
                if(S_poly_ptr -> term_tree == NULL)
                {
                    delete S_poly_ptr;
                }
                else
                {
                    basis_closed = false;
                    poly_basis.data.push_back(S_poly_ptr);
                }
            }
        }
    }while(!basis_closed);

}


};// namespace PolyTool