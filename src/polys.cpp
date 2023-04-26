/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2022-12-18
 * Copyright © Department of Physics, Tsinghua University. All rights reserved
 */

#include "polys.h"
#include <math.h>
#include "poly_utils.h"

/* Monomial */

// monomial evaluation
VarScalar Monomial::eval(VarScalarVec & x)
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

void Monomial::mul_self(Monomial monomial_another)
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

// compare different monomial
CompareResult Monomial::comp(Monomial monomial_another)
{
    assert(dim == monomial_another.dim);
    // // 1. compare the total order. If not equal, return.
    // if(order > monomial_another.order)
    // {
    //     return GT;
    // }
    // else if(order < monomial_another.order)
    // {
    //     return LT;
    // }

    // directly compare by dictionary order

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

// derivative
Monomial Monomial::derivative(IndexType var_id)
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

void Monomial::integrate_self(IndexType var_id)
{
    _order_var[var_id] ++;
    coeff /= (_order_var[var_id]);
    order ++;
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

bool Monomial::divided_by(Monomial denom)
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

bool Monomial::get_division(Monomial denom, Monomial & result)
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

/* PolyMulSweeper */
PolyMulSweeper::PolyMulSweeper(const PolyLinkedList & f, const PolyLinkedList & g)
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


/* PolyLinkedList */
PolyLinkedList::PolyLinkedList(std::vector<Monomial> & monomial_vec)
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

void PolyLinkedList::copy_to(PolyLinkedList & new_homog)
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

void PolyLinkedList::copy_and_call(void (*funcall)(PolyTerm *), PolyLinkedList & new_homog)
{
    new_homog.reinit(dim, increasing_order);

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
void PolyLinkedList::add_term(PolyTerm * new_term)
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

// add another to this, but the terms of another will be destroyed
void PolyLinkedList::destructive_add_self(PolyLinkedList & another)
{
    assert(dim == another.dim);
    assert(increasing_order == another.increasing_order);

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

// remove terms less then EPS
void PolyLinkedList::remove_zeros()
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

void PolyLinkedList::remove_zeros_with_tol(double tol, double exp_base)
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

void PolyLinkedList::reverse()
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

// add up to the third object, and all the terms of self and the other operand
// will all be destroyed.
void PolyLinkedList::destructive_add(PolyLinkedList & another, PolyLinkedList & new_homog)
{
    assert(dim == another.dim);
    assert(increasing_order == another.increasing_order);
    // Compare the pointer at the list head, and append one by one to the new Homogen
    new_homog.reinit(dim, increasing_order);

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
void PolyLinkedList::add_self(PolyLinkedList & another)
{
    // copy another, and destructively add to self
    PolyLinkedList another_copy(dim, increasing_order);
    another.copy_to(another_copy);
    destructive_add_self(another_copy);
}

// add to another Homogen
void PolyLinkedList::add(PolyLinkedList & another, PolyLinkedList & result_copy)
{
    // copy self
    copy_to(result_copy);
    result_copy.add_self(another);
}




// // eval
// Scalar PolyLinkedList::eval(const ScalarVec & x)
// {
//     // evaluation of the homogen
//     ScalarVec x_and_res(dim+1);
    
//     // initialize
//     for(IndexType var_id = 0; var_id < dim; var_id++)
//     {
//         x_and_res[var_id] = x[var_id];
//     }
//     x_and_res[dim] = Scalar(0.0);

//     // evaluate
//     traverse_from_node<PolyTerm, ScalarVec>(term_tree, polyterm_accumulate_eval, &x_and_res);
//     return x_and_res[dim];
// }


VarScalar PolyLinkedList::eval_diff_variable(const IndexVec & diff_order, const VarScalarVec & x, 
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

// derivative
void PolyLinkedList::derivative(IndexType var_id, PolyLinkedList & res)
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



/* poly arithmatics */
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

/* Series */
VarScalar Series::eval(const VarScalarVec & x)
{
    VarScalar res = VarScalar(0.0);
    for(IndexType k = curr_kmin; k < curr_kmax; k++)
    {
        res += homogen_terms[k] -> eval(x);
    }
    return res;
}

// add term to current series
void Series::add_term(Monomial term)
{
    if(term.order < Kmax)
    {
        homogen_terms[term.order] -> add_term(term);
        update_order_bound(term.order);
    }
}

void Series::add_term(PolyTerm* term)
{
    if(term->order < Kmax)
    {
        homogen_terms[term->order] -> add_term(term);
        update_order_bound(term->order);
    }
}

// add homogen to current series
void Series::add_homogen(Homogen & homog)
{
    if(homog.order < Kmax)
    {
        homogen_terms[homog.order] -> add_self(homog); 
        update_order_bound(homog.order);
    }
}

// destructive add homogen to current series
void Series::destructive_add_homogen(Homogen & homog)
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

// add series
void Series::add_series(const Series & new_series)
{
    for(IndexType k = new_series.curr_kmin; (k < Kmax && k < new_series.curr_kmax); k++)
    {
        homogen_terms[k]->add_self(*(new_series.homogen_terms[k]));
    }
    update_order_bound(new_series.curr_kmin);
    update_order_bound(new_series.curr_kmax);
}

// destructive add series
void Series::destructive_add_series(Series & new_series)
{
    for(IndexType k = new_series.curr_kmin; (k< Kmax && k < new_series.curr_kmax); k++)
    {
        homogen_terms[k] -> destructive_add_self(*(new_series.homogen_terms[k]));
    }
    update_order_bound(new_series.curr_kmin);
    update_order_bound(new_series.curr_kmax);
    new_series.curr_kmax = 0; new_series.curr_kmin = new_series.Kmax;

}

// reinit
void Series::reinit()
{
    for(IndexType k = 0; k < Kmax; k++)
    {
        homogen_terms[k] -> reinit(dim, k, homogen_terms[k]->increasing_order);
    }
    curr_kmin = Kmax;
    curr_kmax = 0;
}

// derivative
void Series::derivative(IndexType var_id, Series & res)
{
    res.reinit();
    for(IndexType k = 1; k < curr_kmax && k < (res.Kmax+1); k++)
    {
        homogen_terms[k]->derivative(var_id, *(res.homogen_terms[k-1]));
    }
    res.update_order_bound(curr_kmin - 1);
    res.update_order_bound(curr_kmax - 1);
}

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

// renew term k 
void SeriesPowerSeq::renew_term_k(IndexType k)
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


void SeriesPowerSeq::append_pow_ptr(IndexType order, std::vector<Series*> & series_ptr)
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