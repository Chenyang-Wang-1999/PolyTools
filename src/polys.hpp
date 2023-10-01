/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2022-12-18
 * Copyright © Department of Physics, Tsinghua University. All rights reserved
 */

#ifndef POLYS_H
#define POLYS_H 

#include "basic_defs.h"

#include <complex>
#include <vector>
#include <list>
#include <string>
#include <assert.h>
#include "linklist_tools.hpp"
#include "poly_utils.hpp"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <eigen3/Eigen/Core>

namespace PolyTools{

/* Monomial: type for multivarite monomial */
class Monomial
{
public:
    Monomial(Scalar coeff, IndexVec order_var):
     dim(order_var.size()),coeff(coeff), _order_var(order_var){
        order = 0;
        for(IndexType var_id = 0; var_id < dim; var_id ++)
        {
            order += _order_var[var_id];
        }
    }

    // ~Monomial()
    // {
    //     _order_var.clear();
    // }
    IndexType order; // total order of the monomial
    IndexType dim;  // dimension of the variable 
    Scalar coeff;   // coefficient


    /* Mathematics */
    // evaluation
    VarScalar eval(VarScalarVec & x)
    {
        // for each variable, calculate its order-th pow
        VarScalar val = coeff;
        for(IndexType var_id = 0; var_id < dim; var_id++)
        {
            if(_order_var[var_id] != 0)
            {
                val *= pow(x[var_id], _order_var[var_id]);
            }
        }
        return val;
    }

    // monomial product
    Monomial operator*(Monomial monomial_another)
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

    void mul_self(Monomial monomial_another)
    {
        assert(dim == monomial_another.dim);
        // The new coefficient equals the product of the two monomial
        coeff *= monomial_another.coeff;

        // The new order equals the sum of the two monomial
        for(IndexType var_id = 0; var_id < dim; var_id ++)
        {
            _order_var[var_id] += monomial_another._order_var[var_id];
        }
        order += monomial_another.order;
    }

    // monomial power
    Monomial power(IndexType n)
    {
        Scalar new_coeff = pow(coeff, n);
        IndexVec new_order(dim);
        for(IndexType var_id = 0; var_id < dim; var_id ++)
        {
            new_order[var_id] = _order_var[var_id] * n;
        }
        return Monomial(new_coeff, new_order);
    }

    // compare different monomial
    CompareResult comp(Monomial monomial_another)
    {
        assert(dim == monomial_another.dim);
        // 1. compare variable order by dictionary order.
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

        // 2. all orders are equal
        return EQ;
    }

    bool operator>(Monomial monomial_another){
        CompareResult res = comp(monomial_another);
        return (res == GT);
    }
    bool operator<(Monomial monomial_another){
        CompareResult res = comp(monomial_another);
        return (res == LT);
    }
    bool operator==(Monomial monomial_another){
        CompareResult res = comp(monomial_another);
        return (res == EQ);
    }
    bool operator>=(Monomial monomial_another){
        CompareResult res = comp(monomial_another);
        return ((res == EQ) || (res == GT));
    }
    bool operator<=(Monomial monomial_another){
        CompareResult res = comp(monomial_another);
        return ((res == EQ) || (res == LT));
    }

    Monomial derivative(IndexType var_id)
    {
        IndexVec der_order(_order_var);
        // if order = 0, set coefficient = 0
        if(der_order[var_id] == 0)
        {
            return Monomial(Scalar(0.0), der_order);
        }
        else
        {
            // multiply coefficient by order
            Scalar new_coeff = coeff * Scalar(der_order[var_id]);
            // decrease the current order by 1
            der_order[var_id] --;
            return Monomial(new_coeff, der_order);
        }
    }

    void integrate_self(IndexType var_id)
    {
        _order_var[var_id] ++;
        coeff /= (_order_var[var_id]);
        order ++;
    }

    // Division
    bool divided_by(Monomial denom)
    {
        /* check whether *this is divided by denom */
        assert(dim == denom.dim);
        for(unsigned int var_id = 0; var_id < dim; var_id ++)
        {
            if(_order_var[var_id] < denom._order_var[var_id])
            {
                return false;
            }
        }
        return true;
    }
    bool get_division(Monomial denom, Monomial & result)
    {
        /* get division to result */
        assert(dim == denom.dim);
        for(unsigned int var_id = 0; var_id < dim; var_id ++)
        {
            if(_order_var[var_id] < denom._order_var[var_id])
            {
                return false;
            }
            else
            {
                result._order_var[var_id] = _order_var[var_id] - denom._order_var[var_id];
            }
        }
        result.coeff = coeff / denom.coeff;
        result.order = order - denom.order;
        return true;
    }


    /* Utilities */
    IndexType & var_order(IndexType var_id){return _order_var[var_id];}

    // copy
    Monomial copy(){
        IndexVec new_order_var(dim);
        for(IndexType var_id = 0; var_id < dim; var_id ++)
        {
            new_order_var[var_id] = _order_var[var_id];
        }
        return Monomial(coeff, new_order_var);}

    void print_info()
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

    // partial evaluation
    Monomial partial_eval(VarScalarVec & var_vals, IndexVec & eval_id_list, IndexVec & new_dof_map){
        // evaluate the values of a part of variables and restore the remaining variables
        assert(var_vals.size() == eval_id_list.size());
        assert((var_vals.size() + new_dof_map.size()) == this->dim);

        // 1. calculate new coefficient
        Scalar new_coeff = this->coeff;
        auto var_vals_it = var_vals.begin();
        auto eval_id_list_it = eval_id_list.begin();
        while(var_vals_it != var_vals.end())
        {
            if(this->_order_var[(*eval_id_list_it)])
            {
                new_coeff *= pow((*var_vals_it), this->_order_var[(*eval_id_list_it)]);
            }
            var_vals_it ++;
            eval_id_list_it++;
        }

        // 2. get new var orders
        IndexVec new_var_orders(new_dof_map.size()); 
        auto new_dof_map_it = new_dof_map.begin();
        auto new_var_orders_it = new_var_orders.begin();
        while(new_dof_map_it != new_dof_map.end())
        {
            (*new_var_orders_it) = _order_var[(*new_dof_map_it)];
            new_dof_map_it ++;
            new_var_orders_it ++;
        }

        return Monomial(new_coeff, new_var_orders);
    }
    

private:
    IndexVec _order_var;
};


/* Term for polynomial. Monomial + chain list */
class PolyTerm : public Monomial
{
public:
    PolyTerm(Monomial monomial):
        Monomial(monomial){}
    PolyTerm(Scalar coeff, IndexVec order_var):
        Monomial(coeff, order_var){}
    PolyTerm * next = NULL;
};

/* Functions for PolyTerm object manipulation */
void polyterm_turn_neg(PolyTerm * term){
    term -> coeff = -(term->coeff);
}

void polyterm_scalar_mul(PolyTerm * term, Scalar * k){
    term -> coeff *= *k;
}

void polyterm_print_info(PolyTerm * term){
    term -> print_info();
}

void polyterm_accumulate_eval(PolyTerm *term, VarScalarVec * x_and_res)
{
    VarScalar res = term->eval(*x_and_res);
    (*x_and_res)[term->dim] += res;
}

void polyterm_scale_term(PolyTerm * term, ScalarVec * scale_factor)
{
    for(IndexType var_id = 0; var_id < term->dim; var_id ++)
    {
        IndexType var_order = term -> var_order(var_id);        
        term->coeff *= pow((*scale_factor)[var_id], var_order);
    }
}


/* Tools for homogen multiplication */
typedef struct MatrixIndexType
{
    IndexType row_id;
    IndexType col_id;
} MatrixIndexType;


/* Polynomials implemented by linked list */
class PolyLinkedList
{
public: 
    IndexType n_terms, dim;
    bool increasing_order = true;
    PolyTerm * term_tree = NULL;
    PolyTerm * leading_term = NULL; // leading term
    PolyLinkedList(IndexType dim, bool increasing_order)
        :dim(dim), increasing_order(increasing_order){
        n_terms = 0;
        term_tree = NULL;
    };
    PolyLinkedList(IndexType dim): PolyLinkedList(dim, true)
    {}
    PolyLinkedList(std::vector<Monomial> & monomial_vec)
    {
        // initialize
        n_terms = monomial_vec.size();
        dim = monomial_vec[0].dim;
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

    ~PolyLinkedList()
    {
        destroy_tree<PolyTerm>(term_tree);
    }

    /* Utilities */
    void reinit(IndexType dim)
    {
        reinit(dim, true);
    }

    void reinit(IndexType dim, bool increasing_order)
    {
        this->dim = dim;
        this->increasing_order = increasing_order;
        this->n_terms = 0;
        destroy_tree<PolyTerm>(this->term_tree);
        this->term_tree = NULL;
        leading_term_flag = false;
    }

    // remove terms less then EPS
    void remove_zeros()
    {
        if(term_tree == NULL)
        {
            return;
        }
        
        // remove iterms at term_tree
        while(term_tree != NULL)
        {
            if(ABS_FUN(term_tree->coeff) * pow(EXP_BASE, term_tree->order) < EPS )
            {
                PolyTerm * junk_ptr = pop_first_term();
                delete junk_ptr;
            }
            else
            {
                break;
            }
        }
        if(term_tree == NULL)
        {
            return;
        }

        PolyTerm* curr_ptr = term_tree;
        while(curr_ptr->next != NULL)
        {
            if(ABS_FUN(curr_ptr->next->coeff * pow(EXP_BASE, curr_ptr->next->order))<EPS)
            {
                PolyTerm * junk_ptr = pop_next_term(curr_ptr);
                delete junk_ptr;
            }
            else
            {
                curr_ptr = curr_ptr->next;
            }
        }
    }

    void remove_zeros_with_tol(double tol, double exp_base)
    {
        if(term_tree == NULL)
        {
            return;
        }
        
        // remove iterms at term_tree
        while(term_tree != NULL)
        {
            if(ABS_FUN(term_tree->coeff) * pow(exp_base, term_tree->order) < tol )
            {
                PolyTerm * junk_ptr = pop_first_term();
                delete junk_ptr;
            }
            else
            {
                break;
            }
        }
        if(term_tree == NULL)
        {
            return;
        }

        PolyTerm* curr_ptr = term_tree;
        while(curr_ptr->next != NULL)
        {
            if(ABS_FUN(curr_ptr->next->coeff * pow(exp_base, curr_ptr->next->order))<tol)
            {
                PolyTerm * junk_ptr = pop_next_term(curr_ptr);
                delete junk_ptr;
            }
            else
            {
                curr_ptr = curr_ptr->next;
            }
        }
    }

    void copy_to(PolyLinkedList & new_homog)
    {
        new_homog.reinit(dim, increasing_order);

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

    void copy_and_call(void (*funcall)(PolyTerm *), PolyLinkedList & new_homog)
    {
        new_homog.reinit(dim, increasing_order);

        PolyTerm * new_term_ptr = NULL;
        PolyTerm * curr_term_ptr = term_tree;
        while(curr_term_ptr != NULL)
        {
            PolyTerm* new_term = new PolyTerm(curr_term_ptr->copy());
            funcall(new_term);
            if(new_homog.term_tree == NULL)
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

    
    template <typename DataType>
    void copy_and_call(void (*funcall)(PolyTerm*, DataType*), DataType* data, PolyLinkedList & new_homog)
    {
        new_homog.reinit(dim, increasing_order);

        PolyTerm * new_term_ptr = NULL;
        PolyTerm * curr_term_ptr = term_tree;
        while(curr_term_ptr != NULL)
        {
            PolyTerm * new_term = new PolyTerm(curr_term_ptr->copy());
            funcall(new_term, data);
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

    void print_info()
    {
        STDOUT << "=================\n";
        STDOUT << "Total terms: " << n_terms << '\n';
        STDOUT << "-----------------\n";
        traverse_from_node<PolyTerm>(term_tree, polyterm_print_info);
    }

    void get_max_order(IndexVec & max_orders)
    {
        max_orders.clear();
        max_orders.resize(dim);
        max_orders.assign(dim, 0);
        traverse_from_node<PolyTerm, IndexVec>(term_tree, [](PolyTerm * curr_term, IndexVec * max_order_ptr){
            for(IndexType curr_var = 0; curr_var < curr_term->dim; curr_var++)
            {
                if(curr_term->var_order(curr_var) > (*max_order_ptr)[curr_var])
                {
                    (*max_order_ptr)[curr_var] = curr_term->var_order(curr_var);
                }
            }
        }, &max_orders);
    }

    void get_min_order(const IndexVec & max_orders, IndexVec & min_orders)
    {
        min_orders.clear();
        min_orders.resize(dim);
        min_orders.assign(max_orders.begin(), max_orders.end());
        traverse_from_node<PolyTerm, IndexVec>(term_tree, [](PolyTerm * curr_term, IndexVec * min_order_ptr){
            for(IndexType curr_var = 0; curr_var < curr_term->dim; curr_var++)
            {
                if(curr_term->var_order(curr_var) < (*min_order_ptr)[curr_var])
                {
                    (*min_order_ptr)[curr_var] = curr_term->var_order(curr_var);
                }
            }
        }, &min_orders);
    }

    void get_max_and_min_order(IndexVec & max_orders, IndexVec & min_orders)
    {
        get_max_order(max_orders);
        get_min_order(max_orders, min_orders);
    }

    /* Link list manipulation*/

    void append_new_variable(PolyLinkedList & res)
    {
        traverse_from_node<PolyTerm, PolyLinkedList>(
            term_tree, [](PolyTerm * term, PolyLinkedList * new_poly)
            {
                IndexVec new_orders(term->dim+1);
                for(IndexType var_id = 0; var_id < term->dim; var_id++)
                {
                    new_orders[var_id] = term->var_order(var_id);
                }
                new_orders[term->dim] = 0;
                new_poly->add_term(Monomial(term->coeff, new_orders));
            }, &res
        );
        res.leading_term_flag = false;
    }

    /* add term to the current homogenerous polynomial, 
      the new term is referenced by pointer */     
    void add_term(PolyTerm * new_term)
    {
        assert(new_term != NULL);
        if(ABS_FUN(new_term->coeff) * pow(EXP_BASE, new_term->order) < EPS)
        {
            delete new_term;
            return;
        }

        CompareResult res = GT;
        // 0 if the tree is empty, add directly to the first term
        if(term_tree == NULL){
            insert_at_head(new_term);
            return;
        }

        // 1 check the first term
        res = new_term->comp(*term_tree);
        if(res == ((increasing_order)?LT:GT))
        {
            insert_at_head(new_term);
            return;
        }
        else if (res == EQ)
        {
            term_tree -> coeff += new_term ->coeff;
            if(ABS_FUN(term_tree->coeff) * pow(EXP_BASE, term_tree->order) < EPS)
            {
                PolyTerm * junk_ptr =  pop_first_term();
                delete junk_ptr;
            }
            delete new_term;
            return;
        }

        // 2 for each node, check the new term and the next term
        //   until the term is less then the next term. Then add
        //   the term to the next position of the current pointer
        PolyTerm * curr_ptr = term_tree;
        while(curr_ptr->next != NULL)
        {
            res = new_term->comp(*(curr_ptr->next));
            if(res == ((increasing_order)?LT:GT))
            {
                add_term_after_ptr(curr_ptr, new_term);
                return;
            }
            else if (res == EQ)
            {
                curr_ptr->next->coeff += new_term -> coeff;
                delete new_term;
                if(ABS_FUN(curr_ptr->next->coeff) * pow(EXP_BASE, curr_ptr->next->order) < EPS)
                {
                    PolyTerm* junk_ptr =  pop_next_term(curr_ptr);
                    delete junk_ptr;
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

    void add_term(Monomial new_term){
        PolyTerm * new_term_ptr = new PolyTerm(new_term);
        add_term(new_term_ptr);
    }
    void add_term(Scalar coeff, const IndexVec & order_var)
    {
        PolyTerm * new_term_ptr = new PolyTerm(coeff, order_var);
        add_term(new_term_ptr);
    }

    void add_term_after_ptr(PolyTerm * curr_term, PolyTerm * new_term)
    {
        assert(dim == new_term->dim);
        n_terms ++;
        insert_after_ptr<PolyTerm>(curr_term, new_term);
        leading_term_flag = false;
    }
    void insert_at_head(PolyTerm * new_term)
    {
        assert(dim == new_term->dim);
        n_terms ++;
        new_term->next = term_tree;
        term_tree = new_term;
        leading_term_flag = false;
    }

    PolyTerm* pop_first_term()
    {
        assert(term_tree != NULL);
        PolyTerm * curr_term = term_tree;
        term_tree = term_tree->next;
        n_terms--;
        curr_term -> next = NULL;
        leading_term_flag = false;
        return curr_term;
    }

    /* pop the next term of current term ptr */
    PolyTerm* pop_next_term(PolyTerm * curr_term_ptr)
    {
        PolyTerm* next_ptr = pop_next_node<PolyTerm>(curr_term_ptr);
        n_terms--;
        leading_term_flag = false;
        return next_ptr;
    }

    // update leading term
    bool leading_term_flag = false;
    void update_leading_term()
    {
        leading_term = term_tree;
        if(leading_term == NULL)
        {
            return;
        }
        while(leading_term->next != NULL)
        {
            leading_term = leading_term->next;
        }
        leading_term_flag = true;
    }


    void reverse()
    {
        // 1. Change direction
        increasing_order = ! increasing_order;

        // 2. push all the ptrs to a std::vector
        std::vector<PolyTerm*> all_terms;
        PolyTerm* curr_ptr = term_tree;
        while(curr_ptr != NULL)
        {
            all_terms.push_back(curr_ptr);
            curr_ptr = curr_ptr -> next;
        }

        // 3. pop all the ptrs and append to term_tree
        term_tree = *(all_terms.end());
        all_terms.pop_back();
        curr_ptr = term_tree;
        curr_ptr->next = NULL;
        while(all_terms.size())
        {
            curr_ptr -> next = *(all_terms.end());
            all_terms.pop_back();
            curr_ptr = curr_ptr->next;
            curr_ptr -> next = NULL;
        }
    }

    // to phcpy style string
    std::string to_str(const std::vector<std::string> & var_name_list)
    {
        PolyTerm * curr_ptr = term_tree;
        std::ostringstream poly_string;
        poly_string << std::setprecision(18);
        while(curr_ptr != NULL)
        {
            // 1. coefficient
            #if SCALAR_MODE == 0
                // complex coefficients
                poly_string << '+' << '(' << curr_ptr->coeff.real();
                if(curr_ptr->coeff.imag() >=0)
                {
                    poly_string << '+' << abs(curr_ptr->coeff.imag()) << "*i)";
                }
                else
                {
                    poly_string << '-' << abs(curr_ptr->coeff.imag()) << "*i)";
                }
            #else
                // real coefficients
                poly_string << '+' << '(' << curr_ptr->coeff << ')';
            #endif

            // 2. variables
            for(IndexType curr_var = 0; curr_var < dim; curr_var++)
            {
                if(curr_ptr->var_order(curr_var) > 1)
                {
                    poly_string << "*(" << var_name_list[curr_var] <<"**" << curr_ptr->var_order(curr_var) <<')';
                }
                else if(curr_ptr->var_order(curr_var) == 1)
                {
                    poly_string << "*" << var_name_list[curr_var];
                }
            }
            curr_ptr = curr_ptr->next;
        }
        poly_string << ';';
        return poly_string.str();
    }

    /* Mathematics */

    // add another to this, but the terms of another will be destroyed
    void destructive_add_self(PolyLinkedList & another)
    {
        assert(dim == another.dim);
        assert(increasing_order == another.increasing_order);
        leading_term_flag = false;
        another.leading_term_flag = false;

        PolyTerm * LHS_ptr=term_tree;
        bool is_first_term = true;// whether LHS_ptr points to the first term

        

        // exhaust the other homogen
        while(another.term_tree != NULL)
        {
            // check whether current term_tree is empty
            if(LHS_ptr == NULL)
            {
                term_tree = another.term_tree;
                n_terms = another.n_terms;
                another.term_tree = NULL;
                another.n_terms = 0;
                return;
            }

            PolyTerm * new_term_ptr;
            if(is_first_term)
            {
                CompareResult res = term_tree ->comp(*another.term_tree);
                if(res == ((increasing_order)?LT:GT))
                {
                    is_first_term = false;
                    LHS_ptr = term_tree;
                }
                else if (res == EQ)
                {
                    new_term_ptr = another.pop_first_term();
                    LHS_ptr->coeff += new_term_ptr->coeff;
                    delete new_term_ptr;
                    if(ABS_FUN(LHS_ptr->coeff) * pow(EXP_BASE, LHS_ptr->order) < EPS)
                    {
                        PolyTerm * junk_ptr = pop_first_term();
                        delete junk_ptr;
                        LHS_ptr = term_tree;
                    }
                }
                else
                {
                    new_term_ptr = another.pop_first_term();
                    insert_at_head(new_term_ptr);
                    LHS_ptr = term_tree;
                    is_first_term = false;
                    /* here the second term of "another" is larger than the poped
                    term, i.e., the new header of *this. So there's no need to 
                    check the first term of *this.
                    */
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
                    if(res == ((increasing_order)?LT:GT))
                    {
                        LHS_ptr = LHS_ptr->next;
                    }
                    else if (res == EQ)
                    {
                        new_term_ptr = another.pop_first_term();
                        LHS_ptr->next->coeff += new_term_ptr->coeff;
                        delete new_term_ptr;
                        if(ABS_FUN(LHS_ptr->next->coeff) * pow(EXP_BASE, LHS_ptr->next->order) < EPS)
                        {
                            PolyTerm * junk_ptr = pop_next_term(LHS_ptr);
                            delete junk_ptr;
                        }
                    }
                    else
                    {
                        new_term_ptr = another.pop_first_term();
                        add_term_after_ptr(LHS_ptr, new_term_ptr);
                        LHS_ptr = new_term_ptr;
                    }
                } //if(LHS_ptr -> next == NULL)else
            } // if(is_first_term)else

        } //while(another.term_tree != NULL)
    }


    // add up to the third object, and all the terms of self and the other operand
    // will all be destroyed.
    void destructive_add(PolyLinkedList & another, PolyLinkedList& new_homog)
    {
        assert(dim == another.dim);
        assert(increasing_order == another.increasing_order);
        // Compare the pointer at the list head, and append one by one to the new Homogen
        new_homog.reinit(dim, increasing_order);
        leading_term_flag = false;
        another.leading_term_flag = false;

        PolyTerm * new_term_ptr=NULL, *insert_pos_ptr=NULL;
        PolyTerm * another_term_ptr = NULL;
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
                if(res == ((increasing_order)?LT:GT))
                {
                    new_term_ptr = pop_first_term();
                }
                else if (res == EQ)
                {
                    new_term_ptr = pop_first_term();
                    another_term_ptr = another.pop_first_term();
                    new_term_ptr->coeff += another_term_ptr -> coeff;
                    delete another_term_ptr;
                }
                else
                {
                    new_term_ptr = another.pop_first_term();
                }

                // Insert to new 
                if(insert_pos_ptr==NULL)
                {
                    if(ABS_FUN(new_term_ptr->coeff) * pow(EXP_BASE, new_term_ptr->order) > EPS)
                    {
                        new_homog.insert_at_head(new_term_ptr);
                    }
                    else
                    {
                        delete new_term_ptr;
                    }
                    // bug, n_terms is not updated. new_homog.term_tree = new_term_ptr;
                }
                else
                {
                    if(ABS_FUN(new_term_ptr->coeff) * pow(EXP_BASE, new_term_ptr->order) > EPS)
                    {
                        new_homog.add_term_after_ptr(insert_pos_ptr, new_term_ptr);
                    }
                    else
                    {
                        delete new_term_ptr;
                    }               
                }

                insert_pos_ptr = new_term_ptr;
            } // if(term_tree == NULL)
        } //while
        assert(n_terms == 0);
        assert(another.n_terms == 0);
    }


    // add to self without destroy the other one
    void add_self(PolyLinkedList & another)
    {
        // copy another, and destructively add to self
        PolyLinkedList another_copy(dim, increasing_order);
        another.copy_to(another_copy);
        destructive_add_self(another_copy);
    }

    // add to the third poly
    void add(PolyLinkedList & another, PolyLinkedList & result_copy)
    {
        // copy self
        copy_to(result_copy);
        result_copy.add_self(another);
    }



    // scalar multiplication
    void scalar_mul_self(Scalar k){
        if(ABS_FUN(k) < EPS)
        {
            reinit(dim, increasing_order);
        }
        else{
            traverse_from_node<PolyTerm, Scalar>(term_tree, polyterm_scalar_mul, &k);
        }
    }
    void scalar_mul(Scalar k, PolyLinkedList & new_homog){
        if(ABS_FUN(k) < EPS)
        {
            new_homog.reinit(dim, increasing_order);
        }
        else{
            copy_and_call<Scalar>(polyterm_scalar_mul, &k, new_homog);
        }
    }

    // scale variable
    void scale_var_self(ScalarVec k){
        traverse_from_node<PolyTerm, ScalarVec>(term_tree, polyterm_scale_term, &k);
    }

    void scale_var(ScalarVec k, PolyLinkedList & new_poly)
    {
        copy_and_call<ScalarVec>(polyterm_scale_term, &k, new_poly);
    }

    // negative and substraction
    void neg_self(){
        traverse_from_node<PolyTerm>(term_tree, polyterm_turn_neg);
    }
    void neg(PolyLinkedList & new_homog){
        copy_and_call(polyterm_turn_neg, new_homog);
    }
    void destructive_subs_self(PolyLinkedList & another){
        another.neg_self();
        destructive_add_self(another);
    }
    void subs_self(PolyLinkedList & another){
        another.neg_self();
        add_self(another);
        another.neg_self();
    }
    void subs(PolyLinkedList & another, PolyLinkedList & new_homog){
        another.neg_self();
        add(another, new_homog);
        another.neg_self();
    }

    void destructive_subs(PolyLinkedList & another, PolyLinkedList & new_homog){
        another.neg_self();
        destructive_add(another, new_homog);
    }

    // derivative
    void derivative(IndexType var_id, PolyLinkedList & res)
    {
        res.reinit(dim, increasing_order);
        PolyTerm * new_term_ptr = NULL;
        PolyTerm * curr_term_ptr = term_tree;
        while(curr_term_ptr != NULL)
        {
            Monomial term_der = curr_term_ptr->derivative(var_id);
            if(ABS_FUN(term_der.coeff) * pow(EXP_BASE, term_der.order)>EPS)
            {
                PolyTerm* new_term = new PolyTerm(term_der);
                if(res.n_terms == 0)
                {
                    res.insert_at_head(new_term);
                    new_term_ptr = res.term_tree;
                }
                else
                {
                    res.add_term_after_ptr(new_term_ptr, new_term);
                    new_term_ptr = new_term_ptr->next;
                }
            }

            curr_term_ptr = curr_term_ptr->next;
        }

    }

    void integrate_self(IndexType var_id)
    {
        traverse_from_node<PolyTerm, IndexType>(term_tree, [](PolyTerm * term, IndexType * var_id){
            term -> integrate_self(*var_id);
        }, & var_id);
    }

    /* Evaluation */
    VarScalar eval(const VarScalarVec & x)
    {
        IndexVec diff_order(dim);
        for(auto it = diff_order.begin(); it != diff_order.end(); it++)
        {
            (*it) = 0;
        }
        return eval_diff(diff_order, x);
    }

    VarScalar eval_diff(const IndexVec & diff_order, const VarScalarVec & x)
    {
        return eval_diff_variable(diff_order, x, term_tree, NULL, 0);
    }

    void partial_eval(VarScalarVec & var_vals, IndexVec & eval_id_list, IndexVec & new_dof_map, PolyLinkedList * new_poly){
        new_poly->reinit(new_dof_map.size(), this->increasing_order);

        PolyTerm* curr_ptr = term_tree;
        while(curr_ptr != NULL)
        {
            new_poly->add_term(curr_ptr->partial_eval(var_vals, eval_id_list, new_dof_map));
            curr_ptr = curr_ptr->next;
        }
    }


    /* batch operation */

    // evaluation
    EigenVarVectorX batch_eval(const IndexVec & diff_order, const EigenVarMatrixX x)
    {
        assert(dim == x.rows());
        EigenVarVectorX result_vec(x.cols());
        for(IndexType col_id = 0; col_id < x.cols(); col_id ++)
        {
            VarScalarVec curr_col(x.col(col_id).data(), x.col(col_id).data() + dim);
            result_vec(col_id) = eval_diff(diff_order, curr_col);
        }
        return result_vec;
    }

    // set data
    void batch_add_elements(const ScalarVec & coeff, const IndexVec & orders)
    {
        assert(orders.size() == (dim * coeff.size()));
        IndexVec::const_iterator order_it = orders.begin();
        for(ScalarVec::const_iterator coeff_it = coeff.begin(); coeff_it != coeff.end(); ++coeff_it)
        {
            Scalar curr_coeff = *coeff_it;
            IndexVec curr_order_var(dim);
            for(IndexVec::iterator curr_order_it = curr_order_var.begin(); 
                curr_order_it != curr_order_var.end(); curr_order_it++)
            {
                (*curr_order_it) = (*order_it);
                order_it ++;
            }

            // add new terms
            add_term(curr_coeff, curr_order_var);
        }
    }

    // get data
    void batch_get_data(ScalarVec & coeff, IndexVec & orders)
    {
        coeff.clear(); orders.clear();
        coeff.resize(n_terms); orders.resize(n_terms * dim);
        IndexVec::iterator order_it = orders.begin();
        ScalarVec::iterator coeff_it = coeff.begin();

        // traverse all terms
        PolyTerm * curr_ptr = term_tree;
        while(curr_ptr != NULL)
        {
            (*coeff_it) = curr_ptr -> coeff;
            coeff_it ++;
            for(IndexType var_id = 0; var_id < dim; var_id++)
            {
                (*order_it) = curr_ptr -> var_order(var_id);
                order_it ++;
            }

            curr_ptr = curr_ptr->next;
        }

    }
    
    // initialize with data
    void init_with_data(const ScalarVec & coeff, const IndexVec & orders)
    {
        reinit(dim, increasing_order);
        assert(orders.size() == (dim * coeff.size()));
        IndexVec::const_iterator order_it = orders.begin();
        PolyTerm * curr_insert_ptr = term_tree;
        for(ScalarVec::const_iterator coeff_it = coeff.begin(); coeff_it != coeff.end(); ++coeff_it)
        {
            Scalar curr_coeff = *coeff_it;
            IndexVec curr_order_var(dim);
            for(IndexVec::iterator curr_order_it = curr_order_var.begin(); 
                curr_order_it != curr_order_var.end(); curr_order_it++)
            {
                (*curr_order_it) = (*order_it);
                order_it ++;
            }

            PolyTerm * curr_term = new PolyTerm(curr_coeff, curr_order_var);
            // add new terms
            if(term_tree == NULL)
            {
                insert_at_head(curr_term);
                curr_insert_ptr = term_tree;
            }
            else
            {
                add_term_after_ptr(curr_insert_ptr, curr_term);
                curr_insert_ptr = curr_insert_ptr->next;
            }
        }
    }

    // multiply with a monomial
    void monomial_mul_self(Monomial & new_monomial)
    {
        PolyTerm * curr_ptr = term_tree;
        while(curr_ptr != NULL)
        {
            curr_ptr -> mul_self(new_monomial);
            curr_ptr = curr_ptr -> next;
        }
    }

    void monomial_mul(Monomial & new_monomial, PolyLinkedList & new_poly)
    {
        new_poly.reinit(dim);
        copy_and_call<Monomial>([](PolyTerm * curr_ptr, Monomial * new_monomial_ptr){curr_ptr->mul_self(*new_monomial_ptr);}, 
                &new_monomial, new_poly);
    }

private:
    VarScalar eval_diff_variable(const IndexVec & diff_order, const VarScalarVec & x, 
                PolyTerm* start_ptr, PolyTerm* end_ptr, IndexType var_id)
    {
        PolyTerm* curr_ptr = start_ptr;
        // move to the first ptr that the order >= diff order
        while(curr_ptr != end_ptr)
        {
            if(curr_ptr -> var_order(var_id) >= diff_order[var_id])
            {
                break;
            }
            curr_ptr = curr_ptr -> next;
        }
        VarScalar val = VarScalar(0.0);

        // all order is less then diff order
        if(curr_ptr == end_ptr)
        {
            return val;
        }

        IndexType curr_order = (curr_ptr->var_order(var_id) - diff_order[var_id]);
        VarScalar curr_val = 1.0;
        if(curr_order > 0)
        {
            curr_val = pow(x[var_id], curr_order);
        }
        else
        {
            curr_val = 1.0;
        }

        // if var_id == dim - 1, collect the values from start to end.
        if(var_id == dim - 1)
        {
            while(curr_ptr != end_ptr)
            {
                // collect values
                VarScalar curr_coeff = curr_ptr->coeff;
                // derivative part: n * (n - 1) * ... * (n - m + 1)
                for(IndexType diff_coeff = curr_ptr-> var_order(var_id); 
                    diff_coeff > curr_order; diff_coeff--)
                    {
                        curr_coeff *= diff_coeff;
                    }
                
                val += curr_coeff * curr_val;

                curr_ptr = curr_ptr->next;
                if(curr_ptr != end_ptr)
                {
                    IndexType next_order = (curr_ptr -> var_order(var_id)) - diff_order[var_id];
                    assert(next_order > curr_order);
                    curr_val *= pow(x[var_id], (next_order-curr_order));
                    curr_order = next_order;
                }
            }
        }
        // else, calculate by a similar way except that the coefficients are calculated recursively by the other variables
        else
        {
            PolyTerm * next_start_term = curr_ptr;
            VarScalar coeff = 0.0;
            while(curr_ptr != end_ptr)
            {
                if((curr_ptr->var_order(var_id) - diff_order[var_id]) == curr_order)
                {
                    // move to next
                    curr_ptr = curr_ptr -> next;
                }
                else
                {
                    // calculate coefficients
                    coeff = eval_diff_variable(diff_order, x, next_start_term, curr_ptr, var_id + 1);
                    // calculate derivatives
                    for(IndexType diff_coeff = curr_order + diff_order[var_id]; diff_coeff > curr_order; diff_coeff--)
                    {
                        coeff *= diff_coeff;
                    }
                    val += coeff * curr_val;

                    // update curr_ variables
                    IndexType next_order = (curr_ptr->var_order(var_id) - diff_order[var_id]);
                    assert(next_order > curr_order);
                    curr_val *= pow(x[var_id], next_order - curr_order);
                    next_start_term = curr_ptr;
                    curr_order = next_order;

                    // move ptr
                    curr_ptr = curr_ptr -> next;
                }
            }        

            // The last solve from the last "next_start_term" to "end_ptr"
            coeff = eval_diff_variable(diff_order, x, next_start_term, end_ptr, var_id + 1);
            for(IndexType diff_coeff = curr_order + diff_order[var_id]; diff_coeff > curr_order; diff_coeff--)
            {
                coeff *= diff_coeff;
            }
            val += coeff * curr_val;
        }

        return val;
    }


};

class Laurant
{
public:
    IndexType dim;
    PolyLinkedList num;
    // denom: 1, max_negative_orders
    IndexVec num_max_orders;
    IndexVec denom_orders;

    Laurant(IndexType dim):dim(dim), num(dim){
        num_max_orders.assign(dim, 0);
        denom_orders.assign(dim, 0);
    }

    Laurant(Laurant& laurant_another):dim(laurant_another.dim), num(laurant_another.dim), 
        num_max_orders(laurant_another.num_max_orders), denom_orders(laurant_another.denom_orders)
    {
        laurant_another.num.copy_to(num);
    }

    void reinit()
    {
        num.reinit(dim, num.increasing_order);
        num_max_orders.assign(dim, 0);
        denom_orders.assign(dim, 0);
    }
    
    void reduction()
    {
        // conduct reduction of the Laurant polynomial
        IndexVec num_min_orders;
        
        num.get_min_order(num_max_orders, num_min_orders);

        // 1. find the common factor
        IndexVec orders_common_factor(dim);
        for(IndexType var_id = 0; var_id < dim; var_id++)
        {
            orders_common_factor[var_id] = (num_min_orders[var_id] < denom_orders[var_id])?
                         num_min_orders[var_id] : denom_orders[var_id];
            denom_orders[var_id] -= orders_common_factor[var_id];
            num_max_orders[var_id] -= orders_common_factor[var_id];
        }


        // 2. divided by the common factor
        traverse_from_node<PolyTerm, IndexVec>(
            num.term_tree,
            [](PolyTerm * curr_term, IndexVec* common_factor_ptr)
            {
                for(IndexType var_id = 0; var_id < curr_term->dim; var_id ++)
                {
                    curr_term->var_order(var_id) -= (*common_factor_ptr)[var_id];
                }
            },
            &orders_common_factor
        );
    }

    void set_Laurant(PolyLinkedList & num, Monomial & denom)
    {
        assert(num.dim == dim);
        assert(denom.dim == dim);
        num.scalar_mul(1.0/denom.coeff, this->num);
        num.get_max_order(num_max_orders);
        for(IndexType var_id = 0; var_id < dim; var_id ++)
        {
            denom_orders[var_id] = denom.var_order(var_id);
        }

        reduction();
    }

    void set_Laurant_by_terms(ScalarVec coeffs, std::vector<int> orders)
    {
        // assertion
        assert(coeffs.size() * dim == orders.size());
        // Find minimal orders
        num_max_orders.assign(dim, 0);
        denom_orders.assign(dim, 0);
        for(IndexType term_id = 0; term_id < coeffs.size(); term_id ++)
        {
            for(IndexType var_id = 0; var_id < dim; var_id ++)
            {
                if(orders[term_id * dim + var_id] < 0)
                {
                    if(abs(orders[term_id * dim + var_id]) > denom_orders[var_id])
                    {
                        denom_orders[var_id] = abs(orders[term_id * dim + var_id]);
                    }
                }
            }
        }

        // Add terms
        for(IndexType term_id = 0; term_id < coeffs.size(); term_id ++)
        {
            IndexVec curr_orders(dim);
            for(IndexType var_id = 0; var_id < dim; var_id ++)
            {
                curr_orders[var_id] = orders[term_id * dim + var_id] + denom_orders[var_id];
            }
            num.add_term(Monomial(coeffs[term_id], curr_orders));
        }
        num.get_max_order(num_max_orders);
    }

    VarScalar eval(VarScalarVec x_arr)
    {
        VarScalar val = num.eval(x_arr);
        for(IndexType var_id = 0; var_id < dim; var_id ++)
        {
            if(denom_orders[var_id] > 0)
            {
                val /= pow(x_arr[var_id], denom_orders[var_id]);
            }
        }
        return val;
    }

    /* 
        f(x1, x2, ..., xi, ..., xn) -> g(x) = f(x1, x2, ..., 1/xi, ..., xn) 
    */
    void flip_variable(IndexType var_id, Laurant & new_laurant)
    {
        PolyTerm * curr_ptr = num.term_tree;
        new_laurant.denom_orders.assign(denom_orders.begin(), denom_orders.end());
        bool larger_denom = (denom_orders[var_id] > num_max_orders[var_id]);
        if(larger_denom)
        {
            new_laurant.denom_orders[var_id] = 0;
        }
        else
        {
            new_laurant.denom_orders[var_id] = num_max_orders[var_id] - denom_orders[var_id];
        }
        while(curr_ptr != NULL)
        {
            // new order
            Monomial new_term(*curr_ptr);
            if(larger_denom)
            {
                new_term.var_order(var_id) = denom_orders[var_id] - curr_ptr->var_order(var_id);
            }
            else
            {
                new_term.var_order(var_id) = num_max_orders[var_id] - curr_ptr->var_order(var_id);
            }
            
            new_laurant.num.add_term(new_term);
            curr_ptr = curr_ptr -> next;
        }
        new_laurant.num.get_max_order(new_laurant.num_max_orders);
        new_laurant.reduction();
    }
};

/* Homogeneous polymial*/
class Homogen: public PolyLinkedList
{
public:
    IndexType order;
    Homogen(IndexType dim, IndexType order, bool increasing_order)
        :PolyLinkedList(dim, increasing_order), order(order){
    };
    Homogen(IndexType dim, IndexType order)
        :PolyLinkedList(dim), order(order){
    };
    Homogen(std::vector<Monomial> & monomial_vec): PolyLinkedList(monomial_vec){};

    /* Utilities */
    void reinit(IndexType dim, IndexType order, bool increasing_order)
    {
        this->order = order;
        PolyLinkedList::reinit(dim, increasing_order);
    }
    void reinit(IndexType dim, bool increasing_order)
    {
        PolyLinkedList::reinit(dim, increasing_order);
    }

    void copy_to(Homogen & new_homog){
        new_homog.reinit(dim, order, increasing_order);
        PolyLinkedList::copy_to(new_homog);
    };

    void copy_and_call(void (*funcall)(PolyTerm*), Homogen & new_homog){
        new_homog.reinit(dim, order, increasing_order);
        PolyLinkedList::copy_and_call(funcall, new_homog);
    };
    
    template <typename DataType>
    void copy_and_call(void (*funcall)(PolyTerm*, DataType*), DataType* data, Homogen & new_homog)
    {
        new_homog.reinit(dim, order);
        PolyLinkedList::copy_and_call<DataType>(funcall, data, new_homog);
    }

    void print_info()
    {
        STDOUT << "=================\n";
        STDOUT << "Order: " << order << '\n';
        STDOUT << "Total terms: " << n_terms << '\n';
        STDOUT << "-----------------\n";
        traverse_from_node<PolyTerm>(term_tree, polyterm_print_info);
    }

    /* Link list manipulation*/
    /* add term to the current homogenerous polynomial, 
      the new term is referenced by pointer */     
    void add_term(PolyTerm * new_term)
    {
        assert(new_term->order == order);
        PolyLinkedList::add_term(new_term);
    }
    void add_term(Monomial new_term){
        assert(new_term.order == order);
        PolyLinkedList::add_term(new_term);
    }

    void add_term_after_ptr(PolyTerm * curr_term, PolyTerm * new_term)
    {
        assert(order == new_term->order);
        PolyLinkedList::add_term_after_ptr(curr_term, new_term);
    }
    void insert_at_head(PolyTerm * new_term)
    {
        assert(order == new_term->order);
        PolyLinkedList::insert_at_head(new_term);
    }

    void cast_from_poly(PolyLinkedList & poly)
    {
        assert(poly.term_tree != NULL);
        cast_from_poly(poly, poly.term_tree->order);
    }
    void cast_from_poly(PolyLinkedList & poly, IndexType order)
    {
        this->reinit(this->dim, this->increasing_order);
        if(poly.term_tree == NULL)
        {
            this->order = order;
        }
        else
        {
            assert(order == (poly.term_tree->order));
            this->order = order;
        }
        this->dim = dim;
        this->term_tree = poly.term_tree;
        this->increasing_order = poly.increasing_order;
        this->n_terms = poly.n_terms;
        poly.term_tree = NULL;
    }

    /* Mathematics */

    // add another to this, but the terms of another will be destroyed
    void destructive_add_self(Homogen & another)
    {
        assert(another.order == order);
        PolyLinkedList::destructive_add_self(another);
    }

    // add up to the third object, and all the terms of self and the other operand
    // will all be destroyed.
    void destructive_add(Homogen & another, Homogen& new_homog)
    {
        assert(another.order == order);
        new_homog.reinit(dim, order, increasing_order);
        PolyLinkedList::destructive_add(another, new_homog);
    }

    // add to self without destroy the other one
    void add_self(Homogen & another)
    {
        assert(another.order == order);
        PolyLinkedList::add_self(another);
    }

    // add to the third Homogen
    void add(Homogen & another, Homogen & result_copy)
    {
        assert(another.order == order);
        result_copy.reinit(dim, order, increasing_order);
        PolyLinkedList::add(another, result_copy);
    }

    // scalar multiplication
    void scalar_mul_self(Scalar k){
        traverse_from_node<PolyTerm, Scalar>(term_tree, polyterm_scalar_mul, &k);
    }

    void scalar_mul(Scalar k, Homogen & new_homog){
        copy_and_call<Scalar>(polyterm_scalar_mul, &k, new_homog);
    }

    // negative and substraction
    void neg_self(){
        traverse_from_node<PolyTerm>(term_tree, polyterm_turn_neg);
    }
    void neg(Homogen & new_homog){
        copy_and_call(polyterm_turn_neg, new_homog);
    }
    void destructive_subs_self(Homogen & another){
        another.neg_self();
        destructive_add_self(another);
    }
    void subs_self(Homogen & another){
        another.neg_self();
        add_self(another);
        another.neg_self();
    }
    void subs(Homogen & another, Homogen & new_homog){
        another.neg_self();
        add(another, new_homog);
        another.neg_self();
    }

    void destructive_subs(Homogen & another, Homogen & new_homog){
        another.neg_self();
        destructive_add(another, new_homog);
    }

    // // derivative
    void derivative(IndexType var_id, Homogen & res)
    {
        if(order == 0)
        {
            res.reinit(dim,0, increasing_order);
            return;
        }
        else
        {
            res.reinit(dim,order-1, increasing_order);
            PolyLinkedList::derivative(var_id, res);
        }
    }
};

class PolyMulSweeper
{
public:
    bool increasing_order;

    PolyMulSweeper(const PolyLinkedList & f, const PolyLinkedList & g)
    {
        N_terms.row_id = f.n_terms;
        N_terms.col_id = g.n_terms;

        // clear the stack
        indices_stack.clear();

        // initialize matrices
        ptr_matrix.resize(N_terms.row_id);
        for(IndexType curr_row_id = 0; curr_row_id < N_terms.row_id; curr_row_id ++)
        {
            ptr_matrix[curr_row_id].resize(N_terms.col_id);
        }
        
        tick_matrix.resize(N_terms.row_id);
        for(IndexType curr_row_id = 0; curr_row_id < N_terms.row_id; curr_row_id ++)
        {
            tick_matrix[curr_row_id].resize(N_terms.col_id);
            for(IndexType curr_col_id = 0; curr_col_id < N_terms.col_id; curr_col_id ++)
            {
                tick_matrix[curr_row_id][curr_col_id] = false;
            }
        }

        // fill the ptr matrices by the product
        PolyTerm *f_ptr, *g_ptr;
        for(IndexType row = 0; row < N_terms.row_id; row++)
        {
            // renew ptr
            if(row == 0)
            {
                f_ptr = f.term_tree;
            }
            else
            {
                f_ptr = f_ptr -> next;
            }
            for(IndexType col = 0; col < N_terms.col_id; col++)
            {
                // renew ptr
                if(col == 0)
                {
                    g_ptr = g.term_tree;
                }
                else
                {
                    g_ptr = g_ptr -> next;
                }

                // get a new poly term by the product
                PolyTerm * new_term = new PolyTerm(f_ptr->operator*(*g_ptr));
                ptr_matrix[row][col] = new_term;
            }
        }
    }


    MatrixIndexType N_terms;

    // get the term (0,0), tick (0,0), (1,0), (0,1) and push the later two
    PolyTerm* first_term()
    {
        // return the term (0,0), then push back (0,1), (1,0) 
        tick_matrix[0][0] = true;

        MatrixIndexType new_index1, new_index2;
        new_index1.row_id = 1; new_index1.col_id = 0;
        new_index2.row_id = 0; new_index2.col_id = 1;
        push_new_term(new_index1);
        push_new_term(new_index2);
        return ptr_matrix[0][0];
    }

    // get next term. If next term is NULL, return NULL
    PolyTerm* next_term()
    {
        // find minimize iterators of current stack
        std::vector<std::list<MatrixIndexType>::iterator> minimal_terms;
        find_min_term(minimal_terms);

        // pick the first ptr, and push the lower and right ones
        MatrixIndexType curr_indices = *(minimal_terms[0]);
        indices_stack.erase(minimal_terms[0]);
        PolyTerm * first_ptr = ptr_matrix[
                curr_indices.row_id][
                curr_indices.col_id];
        MatrixIndexType new_mat_id_low, new_mat_id_right;
        new_mat_id_low.row_id = curr_indices.row_id + 1;
        new_mat_id_low.col_id = curr_indices.col_id;
        new_mat_id_right.row_id = curr_indices.row_id;
        new_mat_id_right.col_id = curr_indices.col_id + 1;
        push_new_term(new_mat_id_low);
        push_new_term(new_mat_id_right);

        // if only one element, pop the ptr and return directly
        if(minimal_terms.size() == 1)
        {
            return first_ptr;
        }
        else
        {
            // else, sweep the second to the last ones, pop the indices, 
            // add the coefficient to the first ptr, and destroy the terms
            for(auto min_term_it = minimal_terms.begin() + 1; min_term_it!=minimal_terms.end(); min_term_it ++)
            {
                // current indices ptr is (*min_term_it)
                // pop current indices
                curr_indices = *(*min_term_it);
                indices_stack.erase((*min_term_it));
                PolyTerm * curr_ptr = ptr_matrix[curr_indices.row_id][curr_indices.col_id];
                first_ptr->coeff += curr_ptr->coeff;
                delete curr_ptr;

                // push new indices
                new_mat_id_low.row_id = curr_indices.row_id + 1;
                new_mat_id_low.col_id = curr_indices.col_id;
                new_mat_id_right.row_id = curr_indices.row_id;
                new_mat_id_right.col_id = curr_indices.col_id + 1;
                push_new_term(new_mat_id_low);
                push_new_term(new_mat_id_right);
            }

            // check zero. If zero, destroy the ptr and return NULL
            if(ABS_FUN(first_ptr->coeff) * pow(EXP_BASE, first_ptr->order) < EPS)
            {
                delete first_ptr;
                return NULL;
            }
            else
            {
                return first_ptr;
            }

        }
    }


    /* check whether the sweep is finished */
    bool is_finished()
    {
        // if is not initialized, return false
        if(! tick_matrix[0][0])
        {
            return false;
        }

        // check by the length of stack
        if(indices_stack.size() > 0)
        {
            return false;
        }
        else
        {
            return true;
        }
    }

    /* print the terms in stack for debug */
    void print_stack()
    {
        STDOUT << "=======Current terms in stacks=======\n";
        for(std::list<MatrixIndexType>::iterator it = indices_stack.begin(); it != indices_stack.end(); it++)
        {
            PolyTerm* curr_ptr = ptr_matrix[it->row_id][it->col_id];
            STDOUT << "index:(" << it->row_id << ',' << it->col_id <<")\n"; 
            curr_ptr->print_info();
        }
        STDOUT << "=====================================\n";
    }


    /* print the tick matrix for debug */
    void print_tick_matrix()
    {
        STDOUT << "=======Current ticks=======\n";
        for(IndexType row = 0; row < N_terms.row_id; row ++)
        {
            for(IndexType col = 0; col < N_terms.col_id; col ++)
            {
                STDOUT << tick_matrix[row][col] <<',';
            }
            STDOUT << '\n';
        }
        STDOUT << "=====================================\n";
    }

private:
    void push_new_term(MatrixIndexType new_index)
    {
        // if the new index exceed the bound, do nothing
        if(new_index.row_id >= N_terms.row_id)
        {
            return;
        }
        if(new_index.col_id >= N_terms.col_id)
        {
            return;
        }
        // if the matrix is ticked, return
        if(tick_matrix[new_index.row_id][new_index.col_id])
        {
            return;
        }

        // push back the new index, and tick
        indices_stack.push_back(new_index);
        tick_matrix[new_index.row_id][new_index.col_id] = true;
    }

    // find min term, and add all the min term to the vector
    void find_min_term(std::vector<std::list<MatrixIndexType>::iterator> & min_id_vec)
    {
        // initialize the vector
        min_id_vec.clear(); 
        min_id_vec.push_back(indices_stack.begin());

        // sweep the whole stack
        std::list<MatrixIndexType>::iterator initial_it = indices_stack.begin();
        initial_it ++;
        for(std::list<MatrixIndexType>::iterator stack_it = initial_it; 
            stack_it != indices_stack.end(); stack_it ++)
            {
                PolyTerm * min_poly_term = ptr_matrix[min_id_vec[0]->row_id][min_id_vec[0]->col_id];
                PolyTerm * curr_poly_term = ptr_matrix[stack_it->row_id][stack_it->col_id];
                CompareResult res = min_poly_term -> comp(*curr_poly_term);
                if(res == EQ)
                {
                    min_id_vec.push_back(stack_it);
                }
                else if(res == GT)
                {
                    min_id_vec.clear();
                    min_id_vec.push_back(stack_it);
                }
            }
    }

    std::list<MatrixIndexType> indices_stack;
    std::vector<std::vector<PolyTerm*>> ptr_matrix;
    std::vector<std::vector<bool>> tick_matrix;
};

/* Polynomial arithmatics */
void poly_multiplication(const PolyLinkedList & f, const PolyLinkedList & g, PolyLinkedList & h)
{
    // assert(f.increasing_order == g.increasing_order);
    assert(f.dim == g.dim);
    assert(f.increasing_order);
    assert(g.increasing_order);

    h.reinit(f.dim, f.increasing_order);

    if(f.term_tree == NULL || g.term_tree == NULL)
    {
        return;
    }

    // create sweeper
    PolyMulSweeper sweeper(f, g);

    // add first term
    h.insert_at_head(sweeper.first_term());
    PolyTerm * curr_ptr = h.term_tree;

    // add remaining terms
    while(! sweeper.is_finished())
    {

        PolyTerm * next_term = sweeper.next_term();

        if(next_term != NULL)
        {
            h.add_term_after_ptr(curr_ptr, next_term);
            curr_ptr = curr_ptr->next;
        }
    }
}

void homogen_multiplication(const Homogen & f, const Homogen & g, Homogen & h)
{

    // total order  = f.order + g.order, dim = f.dim
    h.reinit(f.dim, f.order + g.order, f.increasing_order);
    poly_multiplication(f, g, h);
}

// get the multiplication of the homogenerous polynomials in f_seq
void poly_mul_seq(std::vector<PolyLinkedList*> & f_seq, PolyLinkedList & res)
{
    PolyLinkedList temp_poly0 = PolyLinkedList(f_seq[0]->dim, f_seq[0]->increasing_order);
    PolyLinkedList temp_poly1 = PolyLinkedList(f_seq[1]->dim, f_seq[1]->increasing_order);
    PolyLinkedList* temp_poly_ptr[2] = {&temp_poly0, &temp_poly1};
    f_seq[0] -> copy_to(temp_poly0);
    for (IndexType f_id = 1; f_id < f_seq.size(); f_id ++)
    {
        PolyLinkedList * curr_ptr = temp_poly_ptr[((f_id + 1)%2)];
        PolyLinkedList * res_ptr = temp_poly_ptr[(f_id % 2)];
        poly_multiplication(*curr_ptr, *(f_seq[f_id]), * res_ptr);
    }

    // return the result
    temp_poly_ptr[(f_seq.size()+1)%2]->copy_to(res);
}

void homogen_mul_seq(std::vector<Homogen*> & f_seq, Homogen & res)
{
    PolyLinkedList res_poly(f_seq[0]->dim, f_seq[0]->increasing_order);
    std::vector<PolyLinkedList*> f_seq_poly(f_seq.size());
    for(IndexType j = 0; j<f_seq.size();j++)
    {
        f_seq_poly[j] = (PolyLinkedList*) f_seq[j];
    }
    poly_mul_seq(f_seq_poly, res_poly);
    res.cast_from_poly(res_poly);
}

// series of homogenerous polynomial, with maximum order recorded by Kmax
class Series
{
public:
    IndexType Kmax;
    IndexType dim;
    IndexType curr_kmin;
    IndexType curr_kmax;
    bool increasing_order;
    Series(IndexType Kmax, IndexType dim, bool increasing_order)
    : Kmax(Kmax), dim(dim), increasing_order(increasing_order)
    {
        homogen_terms.resize(Kmax);
        for(IndexType k = 0; k < Kmax; k++)
        {
            homogen_terms[k] = new Homogen(dim, k, increasing_order);
        }
        curr_kmin = Kmax;
        curr_kmax = 0;
    }
    Series(IndexType Kmax, IndexType dim):
        Series(Kmax, dim, true){}
    ~Series(){
        for(IndexType k=0; k<Kmax; k++)
        {
            delete homogen_terms[k];
        }
        homogen_terms.clear();
    }

    // reinitialize to zeros
    void reinit()
    {
        for(IndexType k = 0; k < Kmax; k++)
        {
            homogen_terms[k] -> reinit(dim, k, homogen_terms[k]->increasing_order);
        }
        curr_kmin = Kmax;
        curr_kmax = 0;
    }

    // print information
    void print_info()
    {
        for(IndexType k=curr_kmin; k<curr_kmax; k++)
        {
            if(homogen_terms[k]->term_tree != NULL)
            {
                homogen_terms[k] -> print_info();
            }
        }
    }

    /* mathematics */
    void update_order_bound(IndexType new_order)
    {
        if(new_order < curr_kmin)
        {
            curr_kmin = new_order;
        }
        if(new_order >= curr_kmax)
        {
            curr_kmax = new_order + 1;
        }
    }

    // add term to current series
    void add_term(Monomial term)
    {
        if(term.order < Kmax)
        {
            homogen_terms[term.order] -> add_term(term);
            update_order_bound(term.order);
        }
    }

    void add_term(PolyTerm* term)
    {
        if(term->order < Kmax)
        {
            update_order_bound(term->order);
            homogen_terms[term->order] -> add_term(term);
        }
    }

    // add homogen to current series
    void add_homogen(Homogen & homog)
    {
        if(homog.order < Kmax)
        {
            homogen_terms[homog.order] -> add_self(homog); 
            update_order_bound(homog.order);
        }
    }

    // destructive add homogen to current series
    void destructive_add_homogen(Homogen & homog)
    {
        if(homog.order<Kmax)
        {
            homogen_terms[homog.order] -> destructive_add_self(homog);
            update_order_bound(homog.order);
        }
        else
        {
            homog.reinit(homog.dim, homog.order, homog.increasing_order);
        }
    }

    void destructive_add_poly(PolyLinkedList & poly)
    {
        while(poly.term_tree != NULL)
        {
            PolyTerm * curr_term = poly.pop_first_term();
            add_term(curr_term);
        }
    }

    void add_series(const Series & new_series)
    {
        for(IndexType k = new_series.curr_kmin; (k < Kmax && k < new_series.curr_kmax); k++)
        {
            homogen_terms[k]->add_self(*(new_series.homogen_terms[k]));
        }
        update_order_bound(new_series.curr_kmin);
        update_order_bound(new_series.curr_kmax);
    }

    // destructive add series
    void destructive_add_series(Series & new_series)
    {
        for(IndexType k = new_series.curr_kmin; (k< Kmax && k < new_series.curr_kmax); k++)
        {
            homogen_terms[k] -> destructive_add_self(*(new_series.homogen_terms[k]));
        }
        update_order_bound(new_series.curr_kmin);
        update_order_bound(new_series.curr_kmax);
        new_series.curr_kmax = 0; new_series.curr_kmin = new_series.Kmax;

    }

    void scalar_mul_self(Scalar k)
    {
        for(auto it=homogen_terms.begin(); it != homogen_terms.end(); it++)
        {
            (*it) -> scalar_mul_self(k);
        }
    }

    // derivative
    void derivative(IndexType var_id, Series & res)
    {
        res.reinit();
        for(IndexType k = 1; k < curr_kmax && k < (res.Kmax+1); k++)
        {
            homogen_terms[k]->derivative(var_id, *(res.homogen_terms[k-1]));
        }
        res.update_order_bound(curr_kmin - 1);
        res.update_order_bound(curr_kmax - 1);
    }

    VarScalar eval(const VarScalarVec & x)
    {
        VarScalar res = VarScalar(0.0);
        for(IndexType k = curr_kmin; k < curr_kmax; k++)
        {
            res += homogen_terms[k] -> eval(x);
        }
        return res;
    }
    void copy_to(Series & res)
    {
        res.reinit();
        for(IndexType k = curr_kmin; k<curr_kmax; k++)
        {
            res.add_homogen(*(homogen_terms[k]));
        }
    }
    void append_new_variable(Series & res)
    {
        res.reinit();
        for(IndexType k = curr_kmin; k<curr_kmax; k++)
        {
            Homogen new_homog(dim+1, k);
            homogen_terms[k] -> append_new_variable(new_homog);
            res.destructive_add_homogen(new_homog);
        }
    }

    std::vector<Homogen*> homogen_terms;
};

// series vector
class SeriesVec
{
public:
    IndexType var_dim, val_dim;
    IndexType Kmax;
    bool increasing_order = true;
    std::vector<Series*> series_vec;

    SeriesVec(IndexType var_dim, IndexType val_dim, IndexType Kmax, bool increasing_order):
        var_dim(var_dim), val_dim(val_dim), Kmax(Kmax), increasing_order(increasing_order)
        {
            series_vec.resize(val_dim);
            for(auto it = series_vec.begin(); it!= series_vec.end(); it++)
            {
                (*it) = new Series(Kmax, var_dim, increasing_order);
            }
        }
    SeriesVec(IndexType var_dim, IndexType val_dim, IndexType Kmax):
        SeriesVec(var_dim, val_dim, Kmax, true){}
    
    ~SeriesVec()
    {
        for(auto it = series_vec.begin(); it != series_vec.end(); it ++)
        {
            delete (*it);
        }
        series_vec.clear();
    }

    void reinit()
    {
        for(auto it = series_vec.begin(); it != series_vec.end(); it ++)
        {
            (*it) -> reinit();
        }
    }

    void print_info()
    {
        for(IndexType val_id = 0; val_id < val_dim; val_id ++)
        {
            STDOUT << "########################################\n";
            STDOUT << "Current value:" << val_id << '\n';
            series_vec[val_id] -> print_info();
            STDOUT << "########################################\n";
        }
    }

    void copy_to(SeriesVec & res)
    {
        res.reinit();
        for(IndexType j = 0; j < val_dim; j++)
        {
            series_vec[j] -> copy_to(*(res.series_vec[j]));
        }
    }

    void append_new_variable(SeriesVec & res)
    {
        res.reinit();
        for(IndexType j = 0; j < val_dim; j++)
        {
            series_vec[j] -> append_new_variable(*(res.series_vec[j]));
        }
    }

    void destructive_add_poly(IndexType val_id, PolyLinkedList & poly)
    {
        series_vec[val_id] -> destructive_add_poly(poly);
    }

    void add_term(IndexType val_id, Monomial term)
    {
        series_vec[val_id] -> add_term(term);
    }

    void get_poly(IndexType val_id, IndexType curr_order, PolyLinkedList & poly)
    {
        series_vec[val_id] -> homogen_terms[curr_order]->PolyLinkedList::copy_to(poly);
    }
};

// Homogen vector
class HomogenVec
{
public:
    IndexType var_dim, val_dim;
    IndexType order;
    bool increasing_order;
    std::vector<Homogen*> homog_vec;

    HomogenVec(IndexType var_dim, IndexType val_dim, IndexType order, bool increasing_order):
        var_dim(var_dim), val_dim(val_dim), order(order), increasing_order(increasing_order)
        {
            homog_vec.resize(val_dim);
            for(auto it = homog_vec.begin(); it!= homog_vec.end(); it++)
            {
                (*it) = new Homogen(var_dim, order, increasing_order);
            }
        }
    
    HomogenVec(IndexType var_dim, IndexType val_dim, IndexType order):
        HomogenVec(var_dim, val_dim, order, true){}
    ~HomogenVec()
    {
        for(auto it = homog_vec.begin(); it != homog_vec.end(); it ++)
        {
            delete (*it);
        }
        homog_vec.clear();
    }

    void reinit()
    {
        for(auto it = homog_vec.begin(); it != homog_vec.end(); it ++)
        {
            (*it) -> reinit(var_dim, order, increasing_order);
        }
    }    
    
    void reinit(IndexType new_order)
    {
        this->order = new_order;
        for(auto it = homog_vec.begin(); it != homog_vec.end(); it ++)
        {
            (*it) -> reinit(var_dim, order, increasing_order);
        }
    }

    void print_info()
    {
        for(IndexType val_id = 0; val_id < val_dim; val_id ++)
        {
            STDOUT << "########################################\n";
            STDOUT << "Current value:" << val_id << '\n';
            homog_vec[val_id] -> print_info();
            STDOUT << "########################################\n";
        }
    }
};



// arithmatics by order
void series_mul(Series & f, Series & g, IndexType k, Homogen& res)
{
    assert(f.dim == g.dim);
    assert(f.increasing_order == g.increasing_order);
    res.reinit(f.dim, k, f.increasing_order);

    // Sweep the all combinations with total order k
    Homogen current_mul(f.dim, k, f.increasing_order);
    for(IndexType k_f = f.curr_kmin; (k_f <= k && k_f < f.curr_kmax); k_f ++)
    {
        IndexType k_g = k - k_f;
        if(k_f >= f.Kmax)
        {
            break;
        }
        if(k_g >= g.curr_kmax)
        {
            continue;
        }
        if(k_g < g.curr_kmin)
        {
            break;
        }
        homogen_multiplication(*(f.homogen_terms[k_f]), *(g.homogen_terms[k_g]), current_mul);
        assert(current_mul.order == k);
        res.destructive_add_self(current_mul);
    }
}



// composition of a monomial f and a series vec, and return the index k term to the result
void term_comp(Monomial & f, std::vector<Series*> & series_vec, IndexType k, Homogen& res)
{
    assert(f.dim == series_vec.size());
    res.reinit(series_vec[0]->dim, k, series_vec[0]->increasing_order);

    if(f.order == 0)
    {
        if(k > 0)
        {
            return;
        }
        else
        {
            IndexVec new_term_orders(series_vec[0]->dim);
            for(auto it = new_term_orders.begin(); it!=new_term_orders.end(); it++)
            {
                (*it) = 0;
            }
            res.add_term(Monomial(f.coeff, new_term_orders));
            return;
        }
    }

    if(f.order == 1)
    {
        for(IndexType var_id = 0; var_id < f.dim; var_id ++)
        {
            if(f.var_order(var_id) == 1)
            {
                Homogen homog_temp(series_vec[0]->dim, k, series_vec[0]->increasing_order);
                series_vec[var_id] -> homogen_terms[k] -> scalar_mul(f.coeff, homog_temp);
                res.destructive_add_self(homog_temp);
                return;
            }
        }
    }

    // initialize the vectors for series ptr and homogen order
    std::vector<IndexType> homog_order(f.order);
    std::vector<Series*> series_ptr(f.order);
    std::vector<Series*>::iterator series_it = series_ptr.begin();
    for(IndexType new_term_var_id = 0; new_term_var_id < f.dim; new_term_var_id++)
    {
        for(IndexType j = 0; j < f.var_order(new_term_var_id); j++)
        {
            (*series_it) = series_vec[new_term_var_id];  
            series_it ++;
        }
    }

    // set k max and k min
    IndexVec k_min_list(f.order), k_max_list(f.order);
    for(IndexType series_id = 0; series_id < f.order; series_id ++)
    {
        if(series_ptr[series_id] -> curr_kmax <= series_ptr[series_id] -> curr_kmin)
        {
            return;
        }
        else
        {
            k_min_list[series_id] = series_ptr[series_id] -> curr_kmin;
            k_max_list[series_id] = series_ptr[series_id] -> curr_kmax;
        }
    }

    std::vector<Homogen*> homog_ptr(f.order);

    // sweep all combination with total sum k
    IndexVec sep_vec(f.order + 1);
    bool init_successful = boards_left_most(sep_vec, 0, k_min_list, k_max_list, k);
    if(! init_successful)
    {
        return;
    }

    Homogen homog_temp(series_vec[0]->dim, k, series_vec[0]->increasing_order);
    do
    {
        // sweep combinations
        const_sum_boards_to_data(sep_vec, homog_order);


        // check whether to dump the term
        bool dump = false;
        for(IndexType homog_id = 0; homog_id < f.order; homog_id ++)
        {
            if(homog_order[homog_id] >= series_ptr[homog_id] -> curr_kmax)            
            {
                dump = true;
                break;
            }
            homog_ptr[homog_id] = series_ptr[homog_id] -> homogen_terms[homog_order[homog_id]];

            if(homog_ptr[homog_id]->term_tree == NULL)
            {
                // empty term
                dump = true;
                break;
            }
        }


        // if not dump, sum and add to the result
        if(! dump)
        {
            homogen_mul_seq(homog_ptr, homog_temp);
            assert(homog_temp.order == k);

          
            res.destructive_add_self(homog_temp);
        }
    }while(! const_sum_next(sep_vec, k_min_list, k_max_list));

    res.scalar_mul_self(f.coeff);
}

void poly_comp(PolyLinkedList &f, std::vector<Series*> & series_vec, IndexType k, Homogen & res)
{
    assert(f.dim == series_vec.size());
    assert(f.increasing_order == series_vec[0]->increasing_order);
    res.reinit(series_vec[0]->dim, k, series_vec[0]->increasing_order);
    Homogen curr_result(series_vec[0]->dim, k, f.increasing_order);
    PolyTerm* curr_term = f.term_tree;
    while(curr_term != NULL)
    {
        term_comp(*curr_term, series_vec, k, curr_result);
        res.destructive_add_self(curr_result);
        curr_term = curr_term -> next;
    }
}

void poly_comp(PolyLinkedList &f, SeriesVec & series_vec, IndexType k, PolyLinkedList & res)
{
    res.reinit(series_vec.var_dim, true);
    Homogen res_homog(series_vec.var_dim, k);
    poly_comp(f, series_vec.series_vec, k, res_homog);
    res.destructive_add_self(res_homog);
}

void series_comp(Series &f, std::vector<Series*> & series_vec, IndexType k, Homogen& res)
{

    assert(f.dim == series_vec.size());
    assert(f.increasing_order == series_vec[0]->increasing_order);
    res.reinit(series_vec[0]->dim, k, series_vec[0]->increasing_order);
    if(f.curr_kmin >= f.curr_kmax)
    {
        return;
    }

    Homogen curr_result(series_vec[0]->dim, k, f.increasing_order);

    for(IndexType f_k = f.curr_kmin; f_k < f.curr_kmax; f_k ++)
    {
        Homogen & curr_homogen = *f.homogen_terms[f_k];
        PolyTerm * curr_term = curr_homogen.term_tree;
        while(curr_term != NULL)
        {
            term_comp(*curr_term, series_vec, k, curr_result);
            res.destructive_add_self(curr_result);
            curr_term = curr_term -> next;
        }
    }
}


// calculate nabla W dot f
void series_DW_dot_f(Series & W, std::vector<Series*> & f_vec, IndexType k, Homogen & res)
{
    assert(W.dim == f_vec.size());
    assert(W.dim == f_vec[0]->dim);
    assert(W.increasing_order == f_vec[0]->increasing_order);

    res.reinit(W.dim, k, W.increasing_order);

    // for each var_id, add them together
    for(IndexType var_id = 0; var_id < W.dim; var_id ++)
    {
        // W: k_W = 1, 2, ..., k+1, f: k_f = k + 1 - k_W
        for(IndexType k_W = 1; k_W <= k+1; k_W++)
        {
            IndexType k_f = k + 1 - k_W;
            Homogen dWdxi(W.dim, k_W - 1, W.increasing_order);
            Homogen temp_mul(W.dim, k, W.increasing_order);
            W.homogen_terms[k_W]->derivative(var_id, dWdxi);
            homogen_multiplication(dWdxi, *(f_vec[var_id]->homogen_terms[k_f]), temp_mul);
            res.destructive_add_self(temp_mul);
        }
    }
}


// multiply a matrix to a series vector
void scalar_matrix_mul(Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> scalar_matrix,
        const SeriesVec & input_series, SeriesVec & output_series)
{
    assert(scalar_matrix.cols() == input_series.val_dim);
    assert(scalar_matrix.rows() == output_series.val_dim);
    assert(input_series.var_dim == output_series.var_dim);

    output_series.reinit();
    Series temp_series(input_series.Kmax, input_series.var_dim);
    temp_series.reinit();
    // for each row
    for(IndexType row_id = 0; row_id < scalar_matrix.rows(); row_id ++)
    {
        for(IndexType col_id = 0; col_id < scalar_matrix.cols(); col_id++)
        {
            temp_series.add_series(*(input_series.series_vec[col_id]));
            temp_series.scalar_mul_self(scalar_matrix(row_id, col_id));
            (output_series.series_vec[row_id]) -> destructive_add_series(temp_series);
        }
    }
}

void scalar_matrix_mul(Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> scalar_matrix,
        const HomogenVec & input_series, HomogenVec & output_series)
{
    assert(scalar_matrix.cols() == input_series.val_dim);
    assert(scalar_matrix.rows() == output_series.val_dim);
    assert(input_series.var_dim == output_series.var_dim);

    output_series.reinit();
    Homogen temp_homogen(input_series.var_dim, input_series.order, input_series.increasing_order);
    // for each row
    for(IndexType row_id = 0; row_id < scalar_matrix.rows(); row_id ++)
    {
        for(IndexType col_id = 0; col_id < scalar_matrix.cols(); col_id++)
        {
            (input_series.homog_vec[col_id]) -> copy_to(temp_homogen);
            temp_homogen.scalar_mul_self(scalar_matrix(row_id, col_id));
            (output_series.homog_vec[row_id]) -> destructive_add_self (temp_homogen);
        }
    }       
}


void scalar_matrix_mul(Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> scalar_matrix,
        const std::vector<Homogen*> & input_series, HomogenVec & output_series)
{
    assert(scalar_matrix.cols() == input_series.size());
    assert(scalar_matrix.rows() == output_series.val_dim);
    assert(input_series[0]->dim == output_series.var_dim);

    output_series.reinit();
    Homogen temp_homogen(output_series.var_dim, input_series[0] -> order, input_series[0]->increasing_order);
    // for each row
    for(IndexType row_id = 0; row_id < scalar_matrix.rows(); row_id ++)
    {
        for(IndexType col_id = 0; col_id < scalar_matrix.cols(); col_id++)
        {
            (input_series[col_id]) -> copy_to(temp_homogen);
            temp_homogen.scalar_mul_self(scalar_matrix(row_id, col_id));
            (output_series.homog_vec[row_id]) -> destructive_add_self (temp_homogen);
        }
    }   
}


// calculate the integer log2 of x
IndexType int_log2(IndexType x)
{
    IndexType p = 0;
    while((((IndexType)1)<<p) < x)
    {
        p++;
    }
    return p;
}


/* suppose we have x_1, x_2, ..., x_(k-1), and x_1^2, x_2^2, x_3^2, ..., x_(k-1)^2
 calculate x_k^2
 [1] À. Haro, M. Canadell, J.-L. Figueras, A. Luque, and J. M. Mondelo, The Parameterization Method for Invariant Manifolds: From Rigorous Results to Effective Computations, Vol. 195 (Springer International Publishing, Cham, 2016).
    (k%2) == 1: 2 * sum(j=0 to (k-1)/2) x_j * x_{k-j}
    (k%2) == 0: 2 * sum(j=0 to (k-2)/2) x_j * x_{k-j}  +  (x_{k/2})^2
*/

void series_renew_xqr_term(const Series & x, Homogen & x_sqr_k, IndexType next_k)
{
    

    if(next_k < (x.curr_kmin * 2))
    {
        return;
    }

    x_sqr_k.reinit(x.dim, next_k, x.increasing_order);
    if(next_k % 2 == 1)
    {
        for(IndexType j = x.curr_kmin; 2*j + 1 <= next_k; j++)
        {
            if(j >= x.curr_kmax)
            {
                break;
            }
            if((next_k - j) >= x.curr_kmax)
            {
                continue;
            }

            if((next_k - j) < x.curr_kmin)
            {
                break;
            }
            Homogen curr_homog(x.dim, next_k, x.increasing_order);
            homogen_multiplication(*(x.homogen_terms[j]), *(x.homogen_terms[next_k - j]), curr_homog);
            x_sqr_k.destructive_add_self(curr_homog);
        }
        x_sqr_k.scalar_mul_self(Scalar(2.0));
    }
    else
    {
        for(IndexType j = x.curr_kmin; 2*j+2<= next_k; j++)
        {
            if(j >= x.curr_kmax)
            {
                break;
            }
            if((next_k - j) >= x.curr_kmax)
            {
                continue;
            }

            if((next_k - j) < x.curr_kmin)
            {
                break;
            }
            Homogen curr_homog(x.dim, next_k, x.increasing_order);
            homogen_multiplication(*(x.homogen_terms[j]), *(x.homogen_terms[next_k - j]), curr_homog);
            x_sqr_k.destructive_add_self(curr_homog);
        }
        x_sqr_k.scalar_mul_self(Scalar(2.0));

        // x(k/2)^2 term
        Homogen x_khalf_sqr(x.dim, next_k, x.increasing_order);
        homogen_multiplication(*(x.homogen_terms[next_k/2]), *(x.homogen_terms[next_k/2]), x_khalf_sqr);
        x_sqr_k.destructive_add_self(x_khalf_sqr);
    }
}

void series_renew_xqr_term(const Series & x, Series & x_sqr, IndexType next_k)
{
    // assert(x.homogen_terms[0]->term_tree == NULL);
    if(next_k >= x_sqr.Kmax)
    {
        return;
    }
    Homogen new_term(x_sqr.dim, next_k, x_sqr.increasing_order);
    series_renew_xqr_term(x, new_term, next_k);
    x_sqr.destructive_add_homogen(new_term);
}

/* power seq: terms for w^1, w^2, w^4, w^8, ..., where w^1 is a constant pointer
    pointing to a certain series
*/ 
class SeriesPowerSeq: public SeriesVec
{
public:
    IndexType curr_k;  // current k for all ws
    
    SeriesPowerSeq(IndexType max_order, IndexType Kmax, Series * series_w, bool increasing_order):
        SeriesVec(series_w->dim, int_log2(max_order), Kmax, increasing_order)
    {
        curr_k = 1;
        assert(series_w->homogen_terms[0]->term_tree == NULL);
        delete series_vec[0];
        series_vec[0] = series_w;
    }

    SeriesPowerSeq(IndexType max_order, IndexType Kmax, Series * series_w):
            SeriesVec(series_w->dim, int_log2(max_order), Kmax)
    {
        curr_k = 1;
        assert(series_w->homogen_terms[0]->term_tree == NULL);
        delete series_vec[0];
        series_vec[0] = series_w;
    }

    ~SeriesPowerSeq()
    {
        for(auto it = series_vec.begin() + 1; it != series_vec.end(); it ++)
        {
            delete (*it);
        }
        series_vec.clear();
    }

    void reinit()
    {
        for(auto it = series_vec.begin() + 1; it != series_vec.end(); it ++)
        {
            (*it) -> reinit();
        }
    }

    // renew term k 
    void renew_term_k(IndexType k)
    {
        if(k<= curr_k)
        {
            STDERR << "Warning (renew_term_k): k:"<<k<<" is smaller than current k:"<<curr_k<<". Nothing done.\n";
            return;
        }
        IndexType base_min_k = series_vec[0]->curr_kmin;
        
        for(IndexType pow_id = 1; pow_id < val_dim; pow_id++)
        {
            IndexType curr_min_k = (base_min_k << pow_id);
            if(curr_min_k > k)
            {
                break;
            }
            series_renew_xqr_term(*(series_vec[pow_id-1]), *(series_vec[pow_id]), k);

        }
        curr_k = k;
    }



    /* for w^{order} = w^{d1}*w^{d2}*..., append the pointers 
       &w^{d1}, &w^{d2},... to the series_ptr */
    void append_pow_ptr(IndexType order, std::vector<Series*> & series_ptr)
    {
        IndexType seq_id = 0; // order[seq_id] = 2^{seq_id}
        while(order)
        {
            if(order % 2)
            {
                if(seq_id >= val_dim)
                {
                    STDERR << "Error (append_pow_ptr) order exceed the max order of power sequence\n";
                    STDERR << "seq_id: " << seq_id << ", max pow: " << val_dim << '\n';
                    assert(0);
                }
                series_ptr.push_back(series_vec[seq_id]);
            }
            seq_id ++;
            order /= 2;
        }
    }
};

/* composition with the help of series power sequence */
void term_comp(Monomial & f, std::vector<SeriesPowerSeq*> & series_seq_vec, IndexType k, Homogen& res)
{
    assert(f.dim == series_seq_vec.size());
    res.reinit(series_seq_vec[0]->var_dim, k, series_seq_vec[0]->increasing_order);

    if(f.order == 0)
    {
        if(k > 0)
        {
            return;
        }
        else
        {
            IndexVec new_term_orders(series_seq_vec[0]->var_dim);
            for(auto it = new_term_orders.begin(); it!=new_term_orders.end(); it++)
            {
                (*it) = 0;
            }
            res.add_term(Monomial(f.coeff, new_term_orders));
            return;
        }
    }

    // generate series ptr by the power sequences
    std::vector<Series*> series_ptr;
    series_ptr.clear();
    for(IndexType var_id = 0; var_id < f.dim; var_id ++)
    {
        IndexType curr_order = f.var_order(var_id);
        series_seq_vec[var_id] -> append_pow_ptr(curr_order, series_ptr);
    }
    
    if(series_ptr.size() == 0)
    {
        return;
    }
    else if(series_ptr.size() == 1)
    {
        // only one term
        series_ptr[0] -> homogen_terms[k] -> copy_to(res);
        res.scalar_mul_self(f.coeff);
        return;
    }

    // initialize the vectors for series ptr and homogen order
    std::vector<IndexType> homog_order(series_ptr.size());

    // set k max and k min
    IndexVec k_min_list(series_ptr.size()), k_max_list(series_ptr.size());
    for(IndexType series_id = 0; series_id < series_ptr.size(); series_id ++)
    {
        if(series_ptr[series_id] -> curr_kmax <= series_ptr[series_id] -> curr_kmin)
        {
            return;
        }
        else
        {
            k_min_list[series_id] = series_ptr[series_id] -> curr_kmin;
            k_max_list[series_id] = series_ptr[series_id] -> curr_kmax;
        }
    }

    std::vector<Homogen*> homog_ptr(series_ptr.size());

    // sweep all combination with total sum k
    IndexVec sep_vec(series_ptr.size() + 1);
    bool init_successful = boards_left_most(sep_vec, 0, k_min_list, k_max_list, k);
    if(! init_successful)
    {
        return;
    }

    Homogen homog_temp(series_seq_vec[0]->var_dim, k, series_seq_vec[0]->increasing_order);
    do
    {
        // sweep combinations
        const_sum_boards_to_data(sep_vec, homog_order);


        // check whether to dump the term
        bool dump = false;
        for(IndexType homog_id = 0; homog_id < series_ptr.size(); homog_id ++)
        {
            if(homog_order[homog_id] >= series_ptr[homog_id] -> curr_kmax)            
            {
                dump = true;
                break;
            }
            homog_ptr[homog_id] = series_ptr[homog_id] -> homogen_terms[homog_order[homog_id]];

            if(homog_ptr[homog_id]->term_tree == NULL)
            {
                // empty term
                dump = true;
                break;
            }
        }


        // if not dump, sum and add to the result
        if(! dump)
        {
            homogen_mul_seq(homog_ptr, homog_temp);
            assert(homog_temp.order == k);

          
            res.destructive_add_self(homog_temp);
        }
    }while(! const_sum_next(sep_vec, k_min_list, k_max_list));

    res.scalar_mul_self(f.coeff);
}

/* series composition with the help of power sequences */
void series_comp(Series &f, std::vector<SeriesPowerSeq*> & series_seq_vec, IndexType k, Homogen& res)
{
    assert(f.dim == series_seq_vec.size());
    assert(f.increasing_order == series_seq_vec[0]->increasing_order);
    res.reinit(series_seq_vec[0]->var_dim, k, f.increasing_order);

    if(f.curr_kmin >= f.curr_kmax)
    {
        // f empty
        return;
    }

    Homogen curr_result(series_seq_vec[0]->var_dim, k, f.increasing_order);
    for(IndexType f_k = f.curr_kmin; f_k < f.curr_kmax; f_k ++)
    {
        Homogen & curr_homogen = *(f.homogen_terms[f_k]);
        PolyTerm * curr_term = curr_homogen.term_tree;
        while(curr_term != NULL)
        {
            term_comp(*curr_term, series_seq_vec, k, curr_result);
            res.destructive_add_self(curr_result);
            curr_term = curr_term -> next;
        }
    }
}



/* series and linked list */
void series_to_linklist(Series & poly_series, PolyLinkedList & poly_list)
{
    /* prepare a list of ptrs for each homogenerous terms, 
        use std::list to record the non-empty terms of them. */ 
    std::vector<PolyTerm*> homog_term_ptrs(poly_series.curr_kmax - poly_series.curr_kmin);
    std::list<IndexType> non_empty_ptrs;
    for(IndexType k_id = 0; k_id < (poly_series.curr_kmax - poly_series.curr_kmin); k_id++)
    {
        homog_term_ptrs[k_id] = poly_series.homogen_terms[k_id + poly_series.curr_kmin]->term_tree;
        if(homog_term_ptrs[k_id] != NULL)
        {
            non_empty_ptrs.push_back(k_id);
        }
    }

    /* sort term with the help of poly linked list object */
    poly_list.reinit(poly_series.dim, poly_series.homogen_terms[0]->increasing_order);
    PolyTerm * insert_ptr = poly_list.term_tree;
    while(non_empty_ptrs.size())
    {
        IndexType min_id = non_empty_ptrs.front();
        for(auto it = non_empty_ptrs.begin(); it != non_empty_ptrs.end(); it++)
        {
            PolyTerm* curr_ptr = homog_term_ptrs[(*it)];
            CompareResult res = curr_ptr->comp(*homog_term_ptrs[min_id]);
            if(res == ((poly_list.increasing_order)?LT:GT))
            {
                min_id = *it;
            }
        }
        


        // pop the first item
        PolyTerm * new_term = new PolyTerm(*homog_term_ptrs[min_id]);
        homog_term_ptrs[min_id] = homog_term_ptrs[min_id] -> next;

        // erase NULL ptrs
        std::vector<std::list<IndexType>::iterator> list_erase_its;
        list_erase_its.clear();
        for(auto it = non_empty_ptrs.begin(); it != non_empty_ptrs.end(); it++)
        {
            if(homog_term_ptrs[*it] == NULL)
            {
                list_erase_its.push_back(it);
            }
        }
        // erase empty elements
        for(auto it = list_erase_its.begin(); it != list_erase_its.end(); it++)
        {
            non_empty_ptrs.erase((*it));
        }

        if(poly_list.term_tree == NULL)
        {
            poly_list.insert_at_head(new_term);
            insert_ptr = poly_list.term_tree;
        }
        else
        {
            poly_list.add_term_after_ptr(insert_ptr, new_term);
            insert_ptr = insert_ptr->next;
        }
    }

}

} // namespace PolyTools
#endif