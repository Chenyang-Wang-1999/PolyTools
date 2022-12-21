/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2022-12-18
 * Copyright © Department of Physics, Tsinghua University.  All rights reserved
 */

#include "polys.h"
#include <math.h>

/* Monomial */

// monomial evaluation
Scalar Monomial::eval(ScalarVec & x)
{
    // for each variable, calculate its order-th pow
    Scalar val = coeff;
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

void polyterm_accumulate_eval(PolyTerm *term, ScalarVec * x_and_res)
{
    Scalar res = term->eval(*x_and_res);
#ifdef DEBUG
        STDOUT << res << '\n';
#endif
    (*x_and_res)[term->dim] += res;
}

/* PolyMulSweeper */
PolyMulSweeper::PolyMulSweeper(const Homogen & f, const Homogen & g)
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

// push new term and tick
void PolyMulSweeper::push_new_term(MatrixIndexType new_index)
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

// get the term (0,0), tick (0,0), (1,0), (0,1) and push the later two
PolyTerm * PolyMulSweeper::first_term()
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

// find min term, and add all the min term to the vector
void PolyMulSweeper::find_min_term(std::vector<std::list<MatrixIndexType>::iterator> & min_id_vec)
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

// get next term. If next term is NULL, return NULL
PolyTerm* PolyMulSweeper::next_term()
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
            curr_ptr -> ~PolyTerm();

            // push new indices
            new_mat_id_low.row_id = curr_indices.row_id + 1;
            new_mat_id_low.col_id = curr_indices.col_id;
            new_mat_id_right.row_id = curr_indices.row_id;
            new_mat_id_right.col_id = curr_indices.col_id + 1;
            push_new_term(new_mat_id_low);
            push_new_term(new_mat_id_right);
        }

        // check zero. If zero, destroy the ptr and return NULL
        if(ABS_FUN(first_ptr->coeff) < EPS)
        {
            first_ptr -> ~PolyTerm();
            return NULL;
        }
        else
        {
            return first_ptr;
        }

    }
}

// print all terms in stack
void PolyMulSweeper::print_stack()
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
void PolyMulSweeper::print_tick_matrix()
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

// check whether is finished
bool PolyMulSweeper::is_finished()
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

// remove terms less then EPS
void Homogen::remove_zeros()
{
    if(term_tree == NULL)
    {
        return;
    }
    
    // remove iterms at term_tree
    while(term_tree != NULL)
    {
        if(ABS_FUN(term_tree->coeff) < EPS )
        {
            PolyTerm * junk_ptr = pop_first_term();
            junk_ptr -> ~PolyTerm();
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
        if(ABS_FUN(curr_ptr->next->coeff)<EPS)
        {
            PolyTerm * junk_ptr = pop_next_term(curr_ptr);
            junk_ptr -> ~PolyTerm();
        }
        else
        {
            curr_ptr = curr_ptr->next;
        }
    }
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

// eval
Scalar Homogen::eval(const ScalarVec & x)
{
    // evaluation of the homogen
    ScalarVec x_and_res(dim+1);
    
    // initialize
    for(IndexType var_id = 0; var_id < dim; var_id++)
    {
        x_and_res[var_id] = x[var_id];
    }
    x_and_res[dim] = Scalar(0.0);

    // evaluate
    traverse_from_node<PolyTerm, ScalarVec>(term_tree, polyterm_accumulate_eval, &x_and_res);
    return x_and_res[dim];
}

void homogen_multiplication(const Homogen & f, const Homogen & g, Homogen & h)
{
    assert(f.dim == g.dim);

    // total order  = f.order + g.order, dim = f.dim
    h.reinit(f.dim, f.order + g.order);
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
        }
    }
}