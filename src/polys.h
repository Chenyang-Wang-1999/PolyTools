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
#include <iostream>
#include <eigen3/Eigen/Core>


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
    VarScalar eval(VarScalarVec & x);

    // monomial product
    Monomial operator*(Monomial monomial_another);
    void mul_self(Monomial monomial_another);

    // monomial power
    Monomial power(IndexType n);

    // compare different monomial
    CompareResult comp(Monomial monomial_another);
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

    Monomial derivative(IndexType var_id);
    void integrate_self(IndexType var_id);

    // Division
    bool divided_by(Monomial denom);
    bool get_division(Monomial denom, Monomial & result);

    /* Utilities */
    IndexType var_order(IndexType var_id){return _order_var[var_id];}
    // copy
    Monomial copy(){
        IndexVec new_order_var(dim);
        for(IndexType var_id = 0; var_id < dim; var_id ++)
        {
            new_order_var[var_id] = _order_var[var_id];
        }
        return Monomial(coeff, new_order_var);}

    void print_info();
    

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
void polyterm_turn_neg(PolyTerm * term);
void polyterm_scalar_mul(PolyTerm * term, Scalar * k);
void polyterm_print_info(PolyTerm * term);
void polyterm_accumulate_eval(PolyTerm *term, VarScalarVec * x_and_res);
void polyterm_scale_term(PolyTerm * term, ScalarVec * scale_factor);

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
    PolyLinkedList(IndexType dim, bool increasing_order)
        :dim(dim), increasing_order(increasing_order){
        n_terms = 0;
        term_tree = NULL;
    };
    PolyLinkedList(IndexType dim): PolyLinkedList(dim, true)
    {}
    PolyLinkedList(std::vector<Monomial> & monomial_vec);
    ~PolyLinkedList()
    {
        destroy_tree<PolyTerm>(term_tree);
    }

    /* Utilities */
    void reinit(IndexType dim, bool increasing_order)
    {
        this->dim = dim;
        this->increasing_order = increasing_order;
        this->n_terms = 0;
        destroy_tree<PolyTerm>(this->term_tree);
        this->term_tree = NULL;
    }

    void remove_zeros();
    void remove_zeros_with_tol(double tol, double exp_base);
    void copy_to(PolyLinkedList &);
    void copy_and_call(void (*funcall)(PolyTerm*), PolyLinkedList &);
    
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
    }

    /* add term to the current homogenerous polynomial, 
      the new term is referenced by pointer */     
    void add_term(PolyTerm * new_term);
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
    }
    void insert_at_head(PolyTerm * new_term)
    {
        assert(dim == new_term->dim);
        n_terms ++;
        new_term->next = term_tree;
        term_tree = new_term;
    }

    PolyTerm* pop_first_term()
    {
        assert(term_tree != NULL);
        PolyTerm * curr_term = term_tree;
        term_tree = term_tree->next;
        n_terms--;
        curr_term -> next = NULL;
        return curr_term;
    }
    /* pop the next term of current term ptr */
    PolyTerm* pop_next_term(PolyTerm * curr_term_ptr)
    {
        PolyTerm* next_ptr = pop_next_node<PolyTerm>(curr_term_ptr);
        n_terms--;
        return next_ptr;
    }


    void reverse();

    /* Mathematics */

    // add another to this, but the terms of another will be destroyed
    void destructive_add_self(PolyLinkedList & another);

    // add up to the third object, and all the terms of self and the other operand
    // will all be destroyed.
    void destructive_add(PolyLinkedList & another, PolyLinkedList& new_homog);

    // add to self without destroy the other one
    void add_self(PolyLinkedList & another);

    // add to the third poly
    void add(PolyLinkedList & another, PolyLinkedList & result_copy);

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

    void derivative(IndexType var_id, PolyLinkedList & res);
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

private:
    VarScalar eval_diff_variable(const IndexVec & diff_order, const VarScalarVec & x, 
                PolyTerm* start_ptr, PolyTerm* end_ptr, IndexType var_id);

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

    PolyMulSweeper(const PolyLinkedList & f, const PolyLinkedList & g);

    MatrixIndexType N_terms;

    /* give the first term (0,0) of the ptr matrix */
    PolyTerm* first_term();

    /* give the next term of the ptr matrix*/
    PolyTerm* next_term();

    /* check whether the sweep is finished */
    bool is_finished();

    /* print the terms in stack for debug */
    void print_stack();

    /* print the tick matrix for debug */
    void print_tick_matrix();

private:
    void push_new_term(MatrixIndexType new_index);
    void find_min_term(std::vector<std::list<MatrixIndexType>::iterator> & min_id_vec);
    std::list<MatrixIndexType> indices_stack;
    std::vector<std::vector<PolyTerm*>> ptr_matrix;
    std::vector<std::vector<bool>> tick_matrix;
};

/* Homogen arithmatics */

// h = f * g
void poly_multiplication(const PolyLinkedList & f, const PolyLinkedList & g, PolyLinkedList & h);
void poly_mul_seq(std::vector<PolyLinkedList*> & f_seq, PolyLinkedList & res);
void homogen_multiplication(const Homogen & f, const Homogen & g, Homogen & h);
void homogen_mul_seq(std::vector<Homogen*> & f_seq, Homogen & res);

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
    void reinit();

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
    void add_term(Monomial term);
    void add_term(PolyTerm* term);
    void add_homogen(Homogen & homog);
    void destructive_add_homogen(Homogen & homog);
    void destructive_add_poly(PolyLinkedList & poly)
    {
        while(poly.term_tree != NULL)
        {
            PolyTerm * curr_term = poly.pop_first_term();
            add_term(curr_term);
        }
    }
    void add_series(const Series & new_series);
    void destructive_add_series(Series & new_series);
    void scalar_mul_self(Scalar k)
    {
        for(auto it=homogen_terms.begin(); it != homogen_terms.end(); it++)
        {
            (*it) -> scalar_mul_self(k);
        }
    }
    void derivative(IndexType var_id, Series & res);
    VarScalar eval(const VarScalarVec & x);
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

// arithmatics by order
void series_mul(Series & f, Series & g, IndexType k, Homogen& res);
void term_comp(Monomial & f, std::vector<Series*> & series_vec, IndexType k, Homogen& res);
void series_comp(Series &f, std::vector<Series*> & series_vec, IndexType k, Homogen& res);
void series_DW_dot_f(Series & W, std::vector<Series*> & f_vec, IndexType k, Homogen & res);

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


void scalar_matrix_mul(Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> scalar_matrix,
        const SeriesVec & input_series, SeriesVec & output_series);

void scalar_matrix_mul(Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> scalar_matrix,
        const HomogenVec & input_series, HomogenVec & output_series);


void scalar_matrix_mul(Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> scalar_matrix,
        const std::vector<Homogen*> & input_series, HomogenVec & output_series);

// calculate the integer log2 of x
IndexType int_log2(IndexType x);

/* suppose we have x_1, x_2, ..., x_(k-1), and x_1^2, x_2^2, x_3^2, ..., x_(k-1)^2
 calculate x_k^2*/
void series_renew_xqr_term(const Series & x, Series & x_sqr, IndexType next_k);
void series_renew_xqr_term(const Series & x, Homogen & x_sqr_k, IndexType next_k);

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

    void renew_term_k(IndexType k);

    /* for w^{order} = w^{d1}*w^{d2}*..., append the pointers 
       &w^{d1}, &w^{d2},... to the series_ptr */
    void append_pow_ptr(IndexType order, std::vector<Series*> & series_ptr);

};

/* composition with the help of series power sequence */
void term_comp(Monomial & f, std::vector<SeriesPowerSeq*> & series_seq_vec, IndexType k, Homogen& res);
void series_comp(Series &f, std::vector<SeriesPowerSeq*> & series_seq_vec, IndexType k, Homogen& res);
void poly_comp(PolyLinkedList &f, std::vector<Series*> & series_vec, IndexType k, Homogen & res);
void poly_comp(PolyLinkedList &f, SeriesVec & series_vec, IndexType k, PolyLinkedList & res);

/* series and linked list */
void series_to_linklist(Series & poly_series, PolyLinkedList & poly_list);

#endif