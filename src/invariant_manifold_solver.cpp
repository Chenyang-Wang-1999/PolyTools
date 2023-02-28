/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2022-12-21
 * Copyright © Department of Physics, Tsinghua University. All rights reserved
 */

#include "invariant_manifold_solver.h"
#include "basic_defs.h"
#include <eigen3/Eigen/Eigenvalues>

void InvariantManifoldSolver::clear_all()
{
    // clear power sequences
    for(auto it = W_pow_seq.begin(); it != W_pow_seq.end(); it++)
    {
        delete (*it);
    }
    is_initialized = false;
    zero_T = false;
    resonance_tol = DEFAULT_RESONANCE_TOL;
    F.reinit(); 
    f.reinit();
    W.reinit();
}

/* evaluation */
template<typename PolyType>
EigenVarMatrixX eval_poly_vec(const std::vector<PolyType*> & poly_vec, EigenVarMatrixX x_arr, IndexType val_dim)
{
    assert(poly_vec[0] -> dim == x_arr.rows());
    EigenVarMatrixX result_matrix(val_dim, x_arr.cols());

    // evaluate for each row
    for(IndexType col_id = 0; col_id < x_arr.cols(); col_id ++)
    {
        VarScalarVec x_vec(x_arr.col(col_id).data(), x_arr.col(col_id).data() + x_arr.rows());
        for(IndexType row_id = 0; row_id < val_dim; row_id ++)
        {
            result_matrix(row_id, col_id) = poly_vec[row_id] -> eval(x_vec);
        }
    }
    return result_matrix;
}

template<typename PolyVecType>
IndexType poly_var_or_val_dim(const PolyVecType & poly_vec, bool is_var)
{
    if(is_var)
    {
        return poly_vec.var_dim;
    }
    else
    {
        return poly_vec.val_dim;
    }
}

/* Manipulation functions */
void InvariantManifoldSolver::init(EigenMatrixX P, EigenMatrixX T, EigenVectorX lam)
{
    this->P = P;
    this->P_inv = P.inverse();
    this-> T = T;
    this->lam = lam;
    is_initialized = true;
    zero_T = false;

    // fill first terms of W and f
    HomogenVec W1(manifold_dim, phys_dim, 1);
    Homogen curr_f1(manifold_dim, 1);
    scalar_matrix_mul(P.leftCols(manifold_dim), s_vec, W1);
    for(IndexType val_id = 0; val_id < phys_dim; val_id ++)
    {
        W.series_vec[val_id] -> destructive_add_homogen(*(W1.homog_vec[val_id]));
    }

    for(IndexType val_id = 0; val_id < manifold_dim; val_id ++)
    {
        s_vec.homog_vec[val_id] -> copy(curr_f1);
        curr_f1.scalar_mul_self(lam(val_id));
        f.series_vec[val_id] -> destructive_add_homogen(curr_f1);
    }

    // initialize power sequence
    init_power_seq_from_F();

    // set next_k
    next_k = 2;
}

// calculate jacobian by F
EigenMatrixX InvariantManifoldSolver::jacobian_by_F()
{
    // collect jacobian at 0
    EigenMatrixX J(phys_dim, phys_dim);
    Homogen der_homog(phys_dim, 0);
    Homogen * curr_homog = NULL;
    for(IndexType row_id = 0; row_id < phys_dim; row_id++)
    {
        for(IndexType col_id = 0; col_id < phys_dim; col_id++)
        {
            // derivative of the row_id-th series to the col_id-th variable
            curr_homog = F.series_vec[row_id] -> homogen_terms[1];
            curr_homog -> derivative(col_id, der_homog);
            assert(der_homog.order == 0);
            assert(der_homog.n_terms < 2);
            if(der_homog.n_terms)
            {
                J(row_id, col_id) = der_homog.term_tree->coeff;
            }
            else
            {
                J(row_id, col_id) = Scalar(0.0);
            }
        }
    }
    return J;
}

/*
// calculate the jacobian and eigen system by F
void InvariantManifoldSolver::init_by_F()
{
    // collect jacobian at 0
    EigenMatrixX J(phys_dim, phys_dim);
    Homogen der_homog(phys_dim, 0);
    Homogen * curr_homog = NULL;
    for(IndexType row_id = 0; row_id < phys_dim; row_id++)
    {
        for(IndexType col_id = 0; col_id < phys_dim; col_id++)
        {
            // derivative of the row_id-th series to the col_id-th variable
            curr_homog = F.series_vec[row_id] -> homogen_terms[1];
            curr_homog -> derivative(col_id, der_homog);
            assert(der_homog.order == 0);
            assert(der_homog.n_terms < 2);
            if(der_homog.n_terms)
            {
                J(row_id, col_id) = der_homog.term_tree->coeff;
            }
            else
            {
                J(row_id, col_id) = Scalar(0.0);
            }
        }
    }

    // calculate P and lambda
    Eigen::ComplexEigenSolver<EigenMatrixX> eigen_solver(J, true);
    lam = eigen_solver.eigenvalues();
    P = eigen_solver.eigenvectors();
}
*/

// add term with coeff and order_vec to the "val_id" of F
void InvariantManifoldSolver::add_term_F(IndexType val_id, IndexVec order_vec, Scalar coeff)
{
    if(is_initialized)
    {
        STDERR << "F must be edited before initialization. Abort\n";
        return;
    }
    F.series_vec[val_id] -> add_term(Monomial(coeff, order_vec));    
}

// get term order and coefficient
Scalar InvariantManifoldSolver::get_curr_term(IndexVec & order_vec)
{
    assert(acce_term != NULL);
    for(IndexType var_id = 0; var_id < acce_term->dim; var_id ++)
    {
        order_vec[var_id] = acce_term -> var_order(var_id);
    }
    return acce_term -> coeff;
}

// move all the indices to first nonzeros terms
bool InvariantManifoldSolver::move_to_nearest_nonzero_term()
{
    assert(acce_term == NULL);
    while(acce_term == NULL)
    {
        acce_k += 1;
        // all remaining homogen is empty
        if(acce_k >= Kceil)
        {
            return false;
        }

        // else, set the current pointer to the head of acce_k
        acce_term = acce_series_vec -> series_vec[acce_val_id] -> homogen_terms[acce_k] -> term_tree;
    }
    return true;
}

// set data accessor, find first non zero term. Return false if no term is found
bool InvariantManifoldSolver::set_data_accessor(IndexType which_val, IndexType which_poly, IndexType from_k)
{
    // set the labels
    acce_k = from_k; acce_val_id = which_val; acce_poly = which_poly;
    switch(which_poly)
    {
    case POLY_F:
        acce_series_vec = &F;
        break;
    case POLY_f:
        acce_series_vec = &f;
        break;
    case POLY_W:
        acce_series_vec = &W;
        break;
    default: 
        STDERR << "invalid poly index\n";
        assert(0);
    }

    // set the term pointer
    acce_term = acce_series_vec -> series_vec[acce_val_id] -> homogen_terms[acce_k] -> term_tree;
    if(acce_term != NULL)
    {
        return true;
    }
    else
    {
        return move_to_nearest_nonzero_term();
    }
}


// load data and move next
bool InvariantManifoldSolver::load_data_and_move()
{
    assert(acce_term != NULL);
    // term ptr move next
    acce_term = acce_term -> next;
    
    if(acce_term != NULL)
    {
        return true;
    }
    else
    {
        return move_to_nearest_nonzero_term();
    }
}

// get poly as linked list
void InvariantManifoldSolver::get_poly(IndexType which_poly, PolyLinkedList & data_poly, IndexType val_id)
{
    PolyLinkedList & s_vec_ref = *(s_vec.homog_vec[val_id]);
    PolyLinkedList & Ek_ref = *(Ek.homog_vec[val_id]);
    switch(which_poly)
    {
    case POLY_F:
        series_to_linklist(*(F.series_vec[val_id]), data_poly);        
        break;
    case POLY_W:
        series_to_linklist(*(W.series_vec[val_id]), data_poly);        
        break;
    case POLY_f:
        series_to_linklist(*(f.series_vec[val_id]), data_poly);        
        break;
    case POLY_s:
        s_vec_ref.copy(data_poly);
        break;
    case POLY_Ek:
        Ek_ref.copy(data_poly);
        break;
    default:
        STDERR << "invalid poly index\n";
        assert(0);
    }
}

// print poly info
void InvariantManifoldSolver::print_poly_info(IndexType which_poly)
{
    switch(which_poly)
    {
    case POLY_F:
        STDOUT << "Current vector: F\n";
        F.print_info();
        break;
    case POLY_W:
        STDOUT << "Current vector: W\n";
        W.print_info();
        break;
    case POLY_f:
        STDOUT << "Current vector: f\n";
        f.print_info();
        break;
    case POLY_s:
        STDOUT << "Current vector: s\n";
        s_vec.print_info();
        break;
    case POLY_Ek:
        STDOUT << "Current vector: Ek\n";
        Ek.print_info();
        break;
    default:
        STDERR << "invalid poly index\n";
        assert(0);
    }
}


IndexType InvariantManifoldSolver::get_poly_val_or_var_dim(IndexType which_poly, bool is_var)
{
    switch(which_poly)
    {
    case POLY_F:
        return poly_var_or_val_dim<SeriesVec>(F, is_var);
        break;
    case POLY_W:
        return poly_var_or_val_dim<SeriesVec>(W, is_var);
        break;

    case POLY_f:
        return poly_var_or_val_dim<SeriesVec>(f, is_var);
        break;
    case POLY_s:
        return poly_var_or_val_dim<HomogenVec>(s_vec, is_var);
        break;
    case POLY_Ek:
        return poly_var_or_val_dim<HomogenVec>(Ek, is_var);
        break;
    default:
        STDERR << "invalid poly index\n";
        assert(0);
    }
}

// initialize s vector
void InvariantManifoldSolver::init_s_vec()
{
    // s vector: the i-th value is s_i, that is, coeff = 1.0, _order_var[i] = 1, _order_var[j != i] = 0
    IndexVec order_var(manifold_dim);
    // set zero
    for(auto it = order_var.begin(); it != order_var.end(); it++)
    {
        (*it) = 0;
    }

    // for each var_id, add term to s_vec
    s_vec.reinit(1);
    for(IndexType var_id = 0; var_id < manifold_dim; var_id ++)
    {
        order_var[var_id] += 1;
        s_vec.homog_vec[var_id]->add_term(Monomial(Scalar(1.0), order_var));
        order_var[var_id] -= 1;
    }
}

void InvariantManifoldSolver::init_power_seq_from_F()
{
    // find maximum power for each physical variable
    IndexVec max_phys_var_order(phys_dim);
    for(auto it = max_phys_var_order.begin(); it != max_phys_var_order.end(); it++)
    {
        (*it) = 0;
    }

    for(IndexType phys_id=0; phys_id < phys_dim; phys_id++)
    {
        Series & F_curr = *F.series_vec[phys_id];
        for(IndexType f_k = F_curr.curr_kmin; f_k < F_curr.curr_kmax; f_k++)
        {
            PolyTerm * curr_term = F_curr.homogen_terms[f_k]->term_tree;
            while(curr_term != NULL)
            {
                // for each term, compare with max order
                for(IndexType var_id = 0; var_id < phys_dim; var_id ++)
                {
                    if(curr_term ->var_order(var_id) >= max_phys_var_order[var_id])
                    {
                        max_phys_var_order[var_id] = curr_term->var_order(var_id) + 1;
                    }
                }
                curr_term = curr_term -> next;
            }
        }
    }

    for(IndexType phys_id = 0; phys_id < phys_dim; phys_id++)
    {
        W_pow_seq[phys_id] = new SeriesPowerSeq(max_phys_var_order[phys_id], Kmax, W.series_vec[phys_id]);
    }
}

/* main algorithm */

// calculate Ek for "next_k"
void InvariantManifoldSolver::calculate_curr_Ek()
{
    
    // Ek can be calculated element wise. For each homogen term, fetch corresponding homogen terms from F and W
    for(IndexType val_id = 0; val_id < phys_dim; val_id ++)
    {
        Homogen FW(manifold_dim, next_k), DWf(manifold_dim, next_k);

        // calculate FW
        series_comp(*(F.series_vec[val_id]), W_pow_seq, next_k, FW);
        Ek.homog_vec[val_id] -> destructive_add_self(FW);


        series_DW_dot_f(*(W.series_vec[val_id]), f.series_vec, next_k, DWf);
        Ek.homog_vec[val_id] -> destructive_subs_self(DWf);
    }
}


// calculation for eta k
void InvariantManifoldSolver::calculate_eta_k()
{
    etak.reinit(next_k);
    scalar_matrix_mul(-P_inv, Ek, etak);
}

// calculate m\lambda_L = m_1 \lambda_1 + ... + m_d \lambda_d
Scalar InvariantManifoldSolver::get_term_mlam(PolyTerm * term)
{
    Scalar mlam = Scalar(0.0);
    for(IndexType var_id = 0; var_id < manifold_dim; var_id ++)
    {
        mlam += lam(var_id) * Scalar(term->var_order(var_id));
    }
    return mlam;
}

// update xi^N
void InvariantManifoldSolver::update_xiN()
{
    // element wise
    for(IndexType val_id = manifold_dim; val_id < phys_dim; val_id ++)
    {
        // traverse all terms for eta^j
        PolyTerm * curr_ptr = etak.homog_vec[val_id]->term_tree;
        while(curr_ptr != NULL)
        {
            Monomial curr_etak = curr_ptr->copy();
            curr_etak.coeff /= (lam(val_id) - get_term_mlam(curr_ptr));
            xik.homog_vec[val_id]->add_term(curr_etak);
            curr_ptr = curr_ptr -> next;
        }
    }
}

// update Txi^N, write result directly to eta_k
void InvariantManifoldSolver::update_TxiN()
{
    HomogenVec TxiN(manifold_dim, manifold_dim, next_k);
    std::vector<Homogen*> xiN(phys_dim - manifold_dim);
    for(IndexType var_id = manifold_dim; var_id < phys_dim; var_id++)
    {
        xiN[var_id - manifold_dim] = xik.homog_vec[var_id];
    }
    scalar_matrix_mul(T, xiN, TxiN);

    // substract this from etak
    for(IndexType var_id = 0; var_id < manifold_dim; var_id++)
    {
        etak.homog_vec[var_id] -> destructive_subs_self(*(TxiN.homog_vec[var_id]));
    }
}

// update tangent part (xi_L and f_k)
void InvariantManifoldSolver::update_tangent_part()
{
    // for each val_id, sweep all terms like normal part
    for(IndexType val_id = 0; val_id < manifold_dim; val_id++)
    {
        // traverse all terms for eta^j
        PolyTerm * curr_ptr = etak.homog_vec[val_id]->term_tree;
        while(curr_ptr != NULL)
        {
            Monomial curr_etak = curr_ptr->copy();
            Scalar cross_reson_coeff = lam(val_id) - get_term_mlam(curr_ptr);
            if(ABS_FUN(cross_reson_coeff) < resonance_tol)
            {
                // update f
                Monomial xim_term = curr_etak.copy();
                Monomial f_term = curr_etak.copy();
                xim_term.coeff = Scalar(1.0);
                f_term.coeff = cross_reson_coeff - curr_etak.coeff;
                f.series_vec[val_id] -> add_term(f_term);
                xik.homog_vec[val_id] -> add_term(xim_term);
            }
            else
            {
                // update xik
                curr_etak.coeff /= cross_reson_coeff;
                xik.homog_vec[val_id] -> add_term(curr_etak);
            }

            curr_ptr = curr_ptr -> next;
        }
    }
}

// update W
void InvariantManifoldSolver::update_W()
{
    HomogenVec Wk(manifold_dim, phys_dim, next_k);
    scalar_matrix_mul(P, xik, Wk);
    for(IndexType val_id = 0; val_id < phys_dim; val_id ++)
    {
        W.series_vec[val_id] -> destructive_add_homogen(*(Wk.homog_vec[val_id]));
        W_pow_seq[val_id] -> renew_term_k(next_k);
    }
}


// solve one step
int InvariantManifoldSolver::solve_step()
{
    if(is_initialized)
    {
        Ek.reinit(next_k);
        etak.reinit(next_k);
        xik.reinit(next_k);

        //DEBUG
        // STDOUT << "Begin calculate Ek\n";
        calculate_curr_Ek();
        // STDOUT << "Begin calculate etak\n";
        calculate_eta_k();
        // STDOUT << "Begin update xiN\n";
        update_xiN();
        if(!zero_T)
        {
            // STDOUT << "Begin update TxiN\n";
            update_TxiN();
        }
        // STDOUT << "Begin update tangent part\n";
        update_tangent_part();
        // STDOUT << "Begin update W\n";
        update_W();
        // STDOUT << "End\n";
        next_k ++;
    }
    else
    {
        STDERR << "Solver not initialized. Abort.\n";
        return 1;
    }
    return 0;
}

// solve to Kceil
int InvariantManifoldSolver::solve_to_Kceil(bool print_log)
{
    while(next_k < Kceil)
    {
        if(print_log)
        {
            STDOUT << "Solving order:" << next_k << '\n';
        }
        int err_code = solve_step();
        if(err_code)
        {
            return err_code;
        }
    }
    return 0;
}

// evaluation
EigenVarMatrixX InvariantManifoldSolver::eval(IndexType which_poly, EigenVarMatrixX x_arr)
{
    switch(which_poly)
    {
    case POLY_F:
        return eval_poly_vec<Series>(F.series_vec, x_arr, F.val_dim);
        break;
    case POLY_f:
        return eval_poly_vec<Series>(f.series_vec, x_arr, f.val_dim);
        break;
    case POLY_W:
        return eval_poly_vec<Series>(W.series_vec, x_arr, W.val_dim);
        break;
    case POLY_Ek:
        return eval_poly_vec<Homogen>(Ek.homog_vec, x_arr, Ek.val_dim);
        break;
    default: 
        STDERR << "invalid poly index\n";
        assert(0);
    }
}