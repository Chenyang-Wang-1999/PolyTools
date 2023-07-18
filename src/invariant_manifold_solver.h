/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2022-12-21
 * Copyright © Department of Physics, Tsinghua University. All rights reserved
 */

#ifndef INVARIANT_MANIFOLD_SOLVER_H
#define INVARIANT_MANIFOLD_SOLVER_H

#include "polys.hpp"
#include <eigen3/Eigen/Core>


#define POLY_F 0
#define POLY_f 1
#define POLY_W 2
#define POLY_s 3
#define POLY_Ek 4

#define DEFAULT_RESONANCE_TOL 0.05

class InvariantManifoldSolver
{
public: 
    InvariantManifoldSolver(IndexType phys_dim, IndexType manifold_dim, IndexType Kmax):
        Kmax(Kmax), phys_dim(phys_dim), manifold_dim(manifold_dim), F(phys_dim, phys_dim, Kmax),
        f(manifold_dim, manifold_dim, Kmax), W(manifold_dim, phys_dim, Kmax),
        Ek(manifold_dim, phys_dim, Kmax), etak(manifold_dim, phys_dim, Kmax),
        xik(manifold_dim, phys_dim, Kmax), s_vec(manifold_dim, manifold_dim, 1) {
            init_s_vec();
            Kceil = Kmax-1;
            W_pow_seq.resize(phys_dim);
        }

    IndexType Kmax;
    IndexType Kceil;
    IndexType phys_dim, manifold_dim;
    EigenMatrixX P, T, P_inv;
    EigenVectorX lam;
    SeriesVec F;    // vector field
    double resonance_tol = DEFAULT_RESONANCE_TOL; // tolerance for cross resonance

    /* initialization */

    // clear all data and reset initialization
    void clear_all();
    // initialize with T
    void init(EigenMatrixX P, EigenMatrixX T, EigenVectorX lam);
    // initialize without T
    void init_without_T(EigenMatrixX P, EigenVectorX lam)
    {
        init(P, EigenMatrixX::Zero(manifold_dim, phys_dim - manifold_dim), lam);
        zero_T = true;
    }

    ~InvariantManifoldSolver()
    {
        for(auto it = W_pow_seq.begin(); it != W_pow_seq.end(); it++)
        {
            delete (*it);
        }
    }

    // void init_by_F();
    EigenMatrixX jacobian_by_F();

    // solve one step
    int solve_step();

    // solve to Kceil
    int solve_to_Kceil(bool print_log = false);

    /* Utils for python access */
    // add term
    void add_term_F(IndexType val_id, IndexVec order_vec, Scalar coeff);

    // get term order and coefficient
    Scalar get_curr_term(IndexVec & order_vec);

    // set data accessor, find first non zero term. Return false if no term is found
    bool set_data_accessor(IndexType which_val, IndexType which_poly, IndexType from_k);

    // load data and move next
    bool load_data_and_move();

    // print poly info
    void print_poly_info(IndexType which_poly);

    // get poly value dim and variable dim
    IndexType get_poly_val_or_var_dim(IndexType which_poly, bool is_var);

    // evaluation
    EigenVarMatrixX eval(IndexType which_poly, EigenVarMatrixX x_arr);

    // get poly as linked list
    void get_poly(IndexType which_poly, PolyLinkedList & data_poly, IndexType val_id);

    IndexType get_k()
    {
        return next_k;
    }


    void calculate_err(IndexType err_k)
    {
        IndexType next_k_temp = next_k;
        next_k = err_k;
        Ek.reinit(err_k);
        calculate_curr_Ek();
        next_k = next_k_temp;
    }

private:
    std::vector<SeriesPowerSeq*> W_pow_seq;
    SeriesVec f, W;
    IndexType next_k;
    HomogenVec Ek;
    HomogenVec etak;
    HomogenVec xik;
    HomogenVec s_vec;

    /* flags */
    bool is_initialized = false;
    bool zero_T = false;

    // initialize s vector
    void init_s_vec();

    // initialize power sequence for current F
    void init_power_seq_from_F();

    // calculation for Ek
    void calculate_curr_Ek();

    // calculation for eta k
    void calculate_eta_k();

    Scalar get_term_mlam(PolyTerm * term);

    // update xi^N
    void update_xiN();

    // update Txi^N, write result directly to eta_k
    void update_TxiN();

    // update tangent part (xi_L and f_k)
    void update_tangent_part();

    // update W
    void update_W();

    /* utils for python access */
    IndexType acce_k, acce_val_id, acce_poly;
    PolyTerm * acce_term; // current term
    SeriesVec * acce_series_vec; // pointer to current vector

    // move all the indices to first nonzeros terms
    bool move_to_nearest_nonzero_term();
}; //InvariantManifoldSolver

#endif // INVARIANT_MANIFOLD_SOLVER_H