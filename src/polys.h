/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2022-12-18
 * Copyright © Department of Physics, Tsinghua University.  All rights reserved
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

    ~Monomial()
    {
        _order_var.clear();
    }
    IndexType order; // total order of the monomial
    IndexType dim;  // dimension of the variable 
    Scalar coeff;   // coefficient


    /* Mathematics */
    // evaluation
    Scalar eval(ScalarVec & x);

    // monomial product
    Monomial operator*(Monomial monomial_another);

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
void polyterm_accumulate_eval(PolyTerm *term, ScalarVec * x_and_res);

/* Tools for homogen multiplication */
typedef struct MatrixIndexType
{
    IndexType row_id;
    IndexType col_id;
} MatrixIndexType;



/* Homogeneous polymial*/
class Homogen
{
public:
    IndexType n_terms, dim, order;
    PolyTerm * term_tree = NULL;
    Homogen(IndexType dim, IndexType order):dim(dim), order(order){
        n_terms = 0;
        term_tree = NULL;
    };
    Homogen(std::vector<Monomial> & monomial_vec);
    ~Homogen()
    {
        destroy_tree<PolyTerm>(term_tree);
    }

    /* Utilities */
    void reinit(IndexType dim, IndexType order)
    {
        this->dim = dim;
        this->order = order;
        this->n_terms = 0;
        destroy_tree<PolyTerm>(this->term_tree);
        this->term_tree = NULL;
    }

    void remove_zeros();
    void copy(Homogen &);
    void copy_call(void (*funcall)(PolyTerm*), Homogen &);
    
    template <typename DataType>
    void copy_call(void (*funcall)(PolyTerm*, DataType*), DataType* data, Homogen & new_homog)
    {
        new_homog.reinit(dim, order);

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
        STDOUT << "Order: " << order << '\n';
        STDOUT << "Total terms: " << n_terms << '\n';
        STDOUT << "-----------------\n";
        traverse_from_node<PolyTerm>(term_tree, polyterm_print_info);
    }

    /* Link list manipulation*/
    /* add term to the current homogenerous polynomial, 
      the new term is referenced by pointer */     
    void add_term(PolyTerm * new_term);
    void add_term(Monomial new_term){
        PolyTerm * new_term_ptr = new PolyTerm(new_term);
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

    /* Mathematics */

    // add another to this, but the terms of another will be destroyed
    void destructive_add_self(Homogen & another);

    // add up to the third object, and all the terms of self and the other operand
    // will all be destroyed.
    void destructive_add(Homogen & another, Homogen& new_homog);

    // add to self without destroy the other one
    void add_self(Homogen & another);

    // add to the third Homogen
    void add(Homogen & another, Homogen & result_copy);

    // scalar multiplication
    void scalar_mul_self(Scalar k){
        traverse_from_node<PolyTerm, Scalar>(term_tree, polyterm_scalar_mul, &k);
    }
    void scalar_mul(Scalar k, Homogen & new_homog){
        copy_call<Scalar>(polyterm_scalar_mul, &k, new_homog);
    }

    // negative and substraction
    void neg_self(){
        traverse_from_node<PolyTerm>(term_tree, polyterm_turn_neg);
    }
    void neg(Homogen & new_homog){
        copy_call(polyterm_turn_neg, new_homog);
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

    // derivative
    void derivative(IndexType var_id, Homogen & res);

    /* Evaluation */
    Scalar eval(const ScalarVec & x);
};

class PolyMulSweeper
{
public:
    PolyMulSweeper(const Homogen & f, const Homogen & g);

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
void homogen_multiplication(const Homogen & f, const Homogen & g, Homogen & h);
void homogen_mul_seq(std::vector<Homogen*> & f_seq, IndexType total_order, Homogen & res);

// series of homogenerous polynomial, with maximum order recorded by Kmax
class Series
{
public:
    IndexType Kmax;
    IndexType dim;
    IndexType curr_kmin;
    IndexType curr_kmax;
    Series(IndexType Kmax, IndexType dim): Kmax(Kmax), dim(dim)
    {
        homogen_terms.resize(Kmax);
        for(IndexType k = 0; k < Kmax; k++)
        {
            homogen_terms[k] = new Homogen(dim, k);
        }
        curr_kmin = Kmax;
        curr_kmax = 0;
    }
    ~Series(){
        for(IndexType k=0; k<Kmax; k++)
        {
            homogen_terms[k] -> ~Homogen();
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
    void add_homogen(Homogen & homog);
    void destructive_add_homogen(Homogen & homog);
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
    Scalar eval(const ScalarVec & x);
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
    std::vector<Series*> series_vec;

    SeriesVec(IndexType var_dim, IndexType val_dim, IndexType Kmax):
        var_dim(var_dim), val_dim(val_dim), Kmax(Kmax)
        {
            series_vec.resize(val_dim);
            for(auto it = series_vec.begin(); it!= series_vec.end(); it++)
            {
                (*it) = new Series(Kmax, var_dim);
            }
        }
    
    ~SeriesVec()
    {
        for(auto it = series_vec.begin(); it != series_vec.end(); it ++)
        {
            (*it) -> ~Series();
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
};

// Homogen vector
class HomogenVec
{
public:
    IndexType var_dim, val_dim;
    IndexType order;
    std::vector<Homogen*> homog_vec;

    HomogenVec(IndexType var_dim, IndexType val_dim, IndexType order):
        var_dim(var_dim), val_dim(val_dim), order(order)
        {
            homog_vec.resize(val_dim);
            for(auto it = homog_vec.begin(); it!= homog_vec.end(); it++)
            {
                (*it) = new Homogen(var_dim, order);
            }
        }
    
    ~HomogenVec()
    {
        for(auto it = homog_vec.begin(); it != homog_vec.end(); it ++)
        {
            (*it) -> ~Homogen();
        }
        homog_vec.clear();
    }

    void reinit()
    {
        for(auto it = homog_vec.begin(); it != homog_vec.end(); it ++)
        {
            (*it) -> reinit(var_dim, order);
        }
    }    
    
    void reinit(IndexType new_order)
    {
        this->order = new_order;
        for(auto it = homog_vec.begin(); it != homog_vec.end(); it ++)
        {
            (*it) -> reinit(var_dim, order);
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

#endif