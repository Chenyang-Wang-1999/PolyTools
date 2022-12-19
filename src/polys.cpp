/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2022-12-18
 * Copyright © Department of Physics, Tsinghua University.  All rights reserved
 */

#include "polys.h"
#include <math.h>

/* Monomial */

// monomial evaluation
Scalar Monomial::eval(ScalarVec x)
{
    // for each variable, calculate its order-th pow
    Scalar val = coeff;
    for(IndexType var_id = 0; var_id < dim; var_id++)
    {
        val *= pow(x[var_id], _order_var[var_id]);
    }
    return val;
}

// monomial product
Monomial Monomial::operator*(Monomial monomial_another)
{
    assert(dim == monomial_another.dim);
    // The new coefficient equals the product of the two monomial
    Scalar new_coeff = coeff * monomial_another.coeff;

    // The new order equals the sum of the two monomial
    IndexVec new_order_var(dim);
    for(IndexType var_id = 0; var_id < dim; var_id ++)
    {
        new_order_var[var_id] = _order_var[var_id] + monomial_another._order_var[var_id];
    }
    
    // create a new one
    return Monomial(new_coeff, new_order_var);
}

// compare different monomial
CompareResult Monomial::comp(Monomial monomial_another)
{
    assert(dim == monomial_another.dim);
    // 1. compare the total order. If not equal, return.
    if(order > monomial_another.order)
    {
        return GT;
    }
    else if(order < monomial_another.order)
    {
        return LT;
    }

    // 2. compare variable order by dictionary order.
    for(IndexType var_id =0; var_id < dim; var_id++)
    {
        if(_order_var[var_id] > monomial_another._order_var[var_id])
        {
            return GT;
        }
        else if(_order_var[var_id] < monomial_another._order_var[var_id])
        {
            return LT;
        }
    }

    // 3. all orders are equal
    return EQ;
}


// monomial power
Monomial Monomial::power(IndexType n)
{
    Scalar new_coeff = pow(coeff, n);
    IndexVec new_order(dim);
    for(IndexType var_id = 0; var_id < dim; var_id ++)
    {
        new_order[var_id] = _order_var[var_id] * n;
    }
    return Monomial(new_coeff, new_order);
}

void Monomial::print_info()
{
    STDOUT << "Variable dimension: " << dim <<'\n';
    STDOUT << "Total order: " << order << '\n';
    STDOUT << "Coefficient: " << coeff << '\n';
    STDOUT << "Power for each variable:\n";
    for(IndexType var_id = 0; var_id < dim; var_id ++) 
    {
        STDOUT << _order_var[var_id] << ", ";
    }
    STDOUT << '\n';
    STDOUT << "-------------------\n";
}

/* PolyTerm */
void polyterm_turn_neg(PolyTerm * term){
    term -> coeff = -(term->coeff);
}
void polyterm_scalar_mul(PolyTerm * term, Scalar * k){
    term -> coeff *= *k;
}
void polyterm_print_info(PolyTerm * term){
    term -> print_info();
}

/* Homogen */
Homogen::Homogen(std::vector<Monomial> & monomial_vec)
{
    // initialize
    n_terms = monomial_vec.size();
    dim = monomial_vec[0].dim;
    order = monomial_vec[0].order;
    term_tree = new PolyTerm(monomial_vec[0]);
    term_tree->next = NULL;
    PolyTerm * curr_ptr = term_tree;

    // sweep list
    for(IndexType curr_term = 1; curr_term < monomial_vec.size(); curr_term ++)
    {
        curr_ptr->next = new PolyTerm(monomial_vec[curr_term]);
        curr_ptr = curr_ptr->next;
        curr_ptr->next = NULL;
    }
}

void Homogen::copy(Homogen & new_homog)
{
    new_homog.reinit(dim, order);

    PolyTerm * new_term_ptr = NULL;
    PolyTerm * curr_term_ptr = term_tree;
    while(curr_term_ptr != NULL)
    {
        PolyTerm * new_term = new PolyTerm(curr_term_ptr->copy());
        new_term->next = NULL;
        if(new_homog.n_terms == 0)
        {
            new_homog.insert_at_head(new_term);
            new_term_ptr = new_homog.term_tree;
        }
        else
        {
            new_homog.add_term_after_ptr(new_term_ptr, new_term);
            new_term_ptr = new_term_ptr->next;
        }
        curr_term_ptr = curr_term_ptr->next;
    }
}

void Homogen::copy_call(void (*funcall)(PolyTerm *), Homogen & new_homog)
{
    new_homog.reinit(dim, order);

    PolyTerm * new_term_ptr = NULL;
    PolyTerm * curr_term_ptr = term_tree;
    while(curr_term_ptr != NULL)
    {
        PolyTerm* new_term = new PolyTerm(curr_term_ptr->copy());
        funcall(new_term);
        if(new_homog.n_terms == 0)
        {
            new_homog.insert_at_head(new_term);
            new_term_ptr = new_homog.term_tree;
        }
        else
        {
            new_homog.add_term_after_ptr(new_term_ptr, new_term);
            new_term_ptr = new_term_ptr->next;
        }
        curr_term_ptr = curr_term_ptr->next;
    }
}

// add term to correct position (increasing dictionary order)
void Homogen::add_term(PolyTerm * new_term)
{
    assert(new_term != NULL);
    CompareResult res = GT;
    // 0 if the tree is empty, add directly to the first term
    if(term_tree == NULL){
        insert_at_head(new_term);
        return;
    }

    // 1 check the first term
    res = new_term->comp(*term_tree);
    if(res == LT)
    {
        insert_at_head(new_term);
        return;
    }
    else if (res == EQ)
    {
        term_tree -> coeff += new_term ->coeff;
        if(ABS_FUN(term_tree->coeff) < EPS)
        {
            PolyTerm * junk_ptr =  pop_first_term();
            junk_ptr -> ~PolyTerm();
        }
        new_term -> ~PolyTerm();
        return;
    }

    // 2 for each node, check the new term and the next term
    //   until the term is less then the next term. Then add
    //   the term to the next position of the current pointer
    PolyTerm * curr_ptr = term_tree;
    while(curr_ptr->next != NULL)
    {
        res = new_term->comp(*(curr_ptr->next));
        if(res == LT)
        {
            add_term_after_ptr(curr_ptr, new_term);
            return;
        }
        else if (res == EQ)
        {
            curr_ptr->next->coeff += new_term -> coeff;
            new_term -> ~PolyTerm();
            if(ABS_FUN(curr_ptr->next->coeff) < EPS)
            {
                PolyTerm* junk_ptr =  pop_next_term(curr_ptr);
                junk_ptr -> ~PolyTerm();
            }
            return;
        }
        else
        {
            curr_ptr = curr_ptr->next;
        }
    }

    // 3 if the new term is larger than all terms, append the
    //   new term to the end of the link list
    add_term_after_ptr(curr_ptr, new_term);
}

// add another to this, but the terms of another will be destroyed
void Homogen::destructive_add_self(Homogen & another)
{
    assert(dim == another.dim);
    assert(order == another.order);

    PolyTerm * LHS_ptr=term_tree;
    bool is_first_term = true;// whether LHS_ptr points to the first term

    // check whether current term_tree is empty
    if(LHS_ptr == NULL)
    {
        term_tree = another.term_tree;
        n_terms = another.n_terms;
        another.term_tree = NULL;
        another.n_terms = 0;
        return;
    }

    // exhaust the other homogen
    while(another.term_tree != NULL)
    {
        PolyTerm * new_term_ptr;
        if(is_first_term)
        {
            CompareResult res = term_tree ->comp(*another.term_tree);
            switch (res)
            {
            case LT:
                is_first_term = false;
                LHS_ptr = term_tree;
                break;
            case EQ:
                new_term_ptr = another.pop_first_term();
                LHS_ptr->coeff += new_term_ptr->coeff;
                new_term_ptr -> ~PolyTerm();
                if(ABS_FUN(LHS_ptr->coeff) < EPS)
                {
                    PolyTerm * junk_ptr = pop_first_term();
                    junk_ptr -> ~PolyTerm();
                    LHS_ptr = term_tree;
                }
                break;
            case GT:
                new_term_ptr = another.pop_first_term();
                insert_at_head(new_term_ptr);
                LHS_ptr = term_tree;
                is_first_term = false;
                /* here the second term of "another" is larger than the poped
                   term, i.e., the new header of *this. So there's no need to 
                   check the first term of *this.
                   */
                break;
            }
        } // if(is_first_term)
        else
        {
            // consider the term after LHS_ptr
            if(LHS_ptr->next == NULL)
            {
                LHS_ptr->next = another.term_tree;
                n_terms += another.n_terms;
                another.n_terms = 0;
                another.term_tree = NULL;
                return;
            } //if(LHS_ptr -> next == NULL)
            else
            {

                CompareResult res = LHS_ptr->next ->comp(*another.term_tree);
                switch (res)
                {
                case LT:
                    LHS_ptr = LHS_ptr->next;
                    break;
                case GT:
                    new_term_ptr = another.pop_first_term();
                    add_term_after_ptr(LHS_ptr, new_term_ptr);
                    LHS_ptr = new_term_ptr;
                    break;
                case EQ:
                    new_term_ptr = another.pop_first_term();
                    LHS_ptr->next->coeff += new_term_ptr->coeff;
                    new_term_ptr->~PolyTerm();
                    if(ABS_FUN(LHS_ptr->next->coeff) < EPS)
                    {
                        PolyTerm * junk_ptr = pop_next_term(LHS_ptr);
                        junk_ptr -> ~PolyTerm();
                    }
                    break;
                }
            } //if(LHS_ptr -> next == NULL)else
        } // if(is_first_term)else

    } //while(another.term_tree != NULL)
}


// add up to the third object, and all the terms of self and the other operand
// will all be destroyed.
void Homogen::destructive_add(Homogen & another, Homogen & new_homog)
{
    assert(dim == another.dim);
    assert(order == another.order);
    // Compare the pointer at the list head, and append one by one to the new Homogen
    new_homog.reinit(dim, order);

    PolyTerm * new_term_ptr=NULL, *insert_pos_ptr=NULL;
    while((term_tree!=NULL) || (another.term_tree != NULL))
    {
        if(term_tree == NULL)
        {
            if(insert_pos_ptr == NULL)
            {
                new_homog.term_tree = another.term_tree;
            }
            else
            {
                insert_pos_ptr->next = another.term_tree;
            }
            new_homog.n_terms += another.n_terms;
            another.term_tree = NULL;
            another.n_terms = 0;
        }
        else if (another.term_tree == NULL)
        {
            if(insert_pos_ptr == NULL)
            {
                new_homog.term_tree = term_tree;
            }
            else
            {
                insert_pos_ptr->next = term_tree;
            }
            new_homog.n_terms += n_terms;
            term_tree = NULL;
            n_terms = 0;
        }
        else
        {
            // both LHS and RHS has content
            CompareResult res = term_tree -> comp(*(another.term_tree));
            switch (res)
            {
            case LT:
                new_term_ptr = pop_first_term();
                break;
            case GT:
                new_term_ptr = another.pop_first_term();
                break;
            case EQ:
                new_term_ptr = pop_first_term();
                PolyTerm * another_term_ptr = another.pop_first_term();
                new_term_ptr->coeff += another_term_ptr -> coeff;
                another_term_ptr->~PolyTerm();
                break;
            }

            // Insert to new 
            if(insert_pos_ptr==NULL)
            {
                if(ABS_FUN(new_term_ptr->coeff) > EPS)
                {
                    new_homog.insert_at_head(new_term_ptr);
                }
                else
                {
                    new_term_ptr -> ~PolyTerm();
                }
                // bug, n_terms is not updated. new_homog.term_tree = new_term_ptr;
            }
            else
            {
                if(ABS_FUN(new_term_ptr->coeff) > EPS)
                {
                    new_homog.add_term_after_ptr(insert_pos_ptr, new_term_ptr);
                }
                else
                {
                    new_term_ptr -> ~PolyTerm();
                }               
            }

            insert_pos_ptr = new_term_ptr;
        } // if(term_tree == NULL)
    } //while
    assert(n_terms == 0);
    assert(another.n_terms == 0);
}

// add to self without destroy the other one
void Homogen::add_self(Homogen & another)
{
    // copy another, and destructively add to self
    Homogen another_copy(dim, order);
    another.copy(another_copy);
    destructive_add_self(another_copy);
    another_copy.~Homogen();
}

// add to another Homogen
void Homogen::add(Homogen & another, Homogen & result_copy)
{
    // copy self
    copy(result_copy);
    result_copy.add_self(another);
}