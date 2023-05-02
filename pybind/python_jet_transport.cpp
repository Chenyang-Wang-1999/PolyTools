/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2023-04-20
 * copy_toright © Department of Physics, Tsinghua University. All rights reserved
 */

#include "python_jet_transport.h"

void set_series_vector(SeriesVec & target_series_vec, IndexType dim, PolyLinkedList & poly)
{
    target_series_vec.series_vec[dim]->reinit();
    PolyTerm * term_ptr = poly.term_tree;
    while(term_ptr != NULL)
    {
        target_series_vec.series_vec[dim] -> add_term(*term_ptr);
        term_ptr = term_ptr -> next;
    }
}

void get_series_vector(const SeriesVec & source_series_vec, IndexType dim, IndexType order, PolyLinkedList & poly)
{
    source_series_vec.series_vec[dim] -> homogen_terms[order] -> PolyLinkedList::copy_to(poly);
}


void cJetTransport::set_sol_by_scalars(ScalarVec & coeffs, IndexType err_order)
{
    // traverse all coefficients, and batch add them to "sol"
    for(IndexType var_id = 0; var_id < var_dim; var_id++)
    {
        IndexType id_offset = var_id * order_start_indices[err_order];
        for(IndexType curr_order=0; curr_order < err_order; curr_order++)
        {
            IndexVec sol_indices;
            ScalarVec sol_coeffs;
            sol_indices.clear();
            sol_coeffs.clear();
            for(IndexType coeff_id = order_start_indices[curr_order];
             coeff_id < order_start_indices[curr_order + 1]; coeff_id++)
             {
                if(ABS_FUN(coeffs[id_offset + coeff_id]) * pow(EXP_BASE, curr_order)>EPS)
                {
                    sol_coeffs.push_back(coeffs[id_offset + coeff_id]);
                    for(auto it = all_series_orders[coeff_id].begin();
                        it != all_series_orders[coeff_id].end();
                        it++)
                    {
                        sol_indices.push_back((*it));
                    }
                }
             }
            sol.series_vec[var_id] -> homogen_terms[curr_order]->init_with_data(sol_coeffs, sol_indices);
            sol.series_vec[var_id] ->update_order_bound(curr_order);
        }
    }
}

void cJetTransport::get_sol_as_scalars(ScalarVec & coeffs, IndexType err_order)
{
    coeffs.clear();
    coeffs.resize(var_dim*order_start_indices[err_order]);
    // set zero
    for(auto it=coeffs.begin(); it!=coeffs.end();it++)
    {
        (*it) = 0.0;
    }

    for(IndexType var_id=0; var_id < var_dim; var_id++)
    {
        IndexType id_offset = var_id * order_start_indices[err_order];
        for(IndexType curr_order = 0; curr_order < err_order; curr_order++)
        {
            PolyTerm * term_ptr = sol.series_vec[var_id]->homogen_terms[curr_order]->term_tree;
            for(IndexType coeff_id = order_start_indices[curr_order];
             coeff_id < order_start_indices[curr_order + 1]; coeff_id++)
             {
                if(term_ptr == NULL)
                {
                    break;
                }

                CompareResult comp = term_ptr -> comp(Monomial(1, all_series_orders[coeff_id]));
                if(comp == EQ)
                {
                    coeffs[id_offset + coeff_id] = term_ptr -> coeff;
                    term_ptr = term_ptr->next;
                }
                else
                {
                    assert(comp == GT);
                }
             }
        }
    }
}

void cJetTransport::get_dvars_as_scalars(ScalarVec & coeffs, IndexType err_order)
{
    coeffs.clear();
    coeffs.resize(var_dim*order_start_indices[err_order]);
    // set zero
    for(auto it=coeffs.begin(); it!=coeffs.end();it++)
    {
        (*it) = 0.0;
    }

    for(IndexType var_id=0; var_id < var_dim; var_id++)
    {
        IndexType id_offset = var_id * order_start_indices[err_order];
        for(IndexType curr_order = 0; curr_order < err_order; curr_order++)
        {
            PolyTerm * term_ptr = dvars.series_vec[var_id]->homogen_terms[curr_order]->term_tree;
            for(IndexType coeff_id = order_start_indices[curr_order];
             coeff_id < order_start_indices[curr_order + 1]; coeff_id++)
             {
                if(term_ptr == NULL)
                {
                    break;
                }

                CompareResult comp = term_ptr -> comp(Monomial(1, all_series_orders[coeff_id]));
                if(comp == EQ)
                {
                    coeffs[id_offset + coeff_id] = term_ptr -> coeff;
                    term_ptr = term_ptr->next;
                }
                else
                {
                    assert(comp == GT);
                }
             }
        }
    }
}

void cJetTransport::calculate_dvars(IndexType err_order)
{
    // calculate dvars by composition
    std::vector<Series*> all_series(var_dim + param_dim);
    for(IndexType var_id = 0; var_id < var_dim; var_id++)
    {
        all_series[var_id] = sol.series_vec[var_id];
    }
    for(IndexType param_id = 0; param_id < param_dim; param_id++)
    {
        all_series[param_id + var_dim] = const_params.series_vec[param_id];
    }

    for(IndexType var_id = 0; var_id < var_dim; var_id++)
    {

        for(IndexType curr_order=0; curr_order<err_order; curr_order++)
        {
            series_comp(*(ode_fun.series_vec[var_id]), all_series, curr_order, *(dvars.series_vec[var_id]->homogen_terms[curr_order]));
            dvars.series_vec[var_id] -> update_order_bound(curr_order);
        }
    }

}

void cJetTransport::derivative_iterate()
{
    /*
        Do the iteration \Phi'_{n+1}(\tau, x_0 + \delta x; lambda_0 + \delta \lambda)
             = F(\Phi_n(\tau, x_0 + \delta x; \lambda_0 + \delta \lambda))
        If y_{n+1} == y_n, return true, else return false
    */

    // 1. prepare the series of (\Phi_n; \lambda)
    std::vector <Series*> all_series(var_dim + param_dim);
    for(IndexType var_id = 0; var_id < var_dim; var_id++)
    {
        all_series[var_id] = sol_time.series_vec[var_id];
    }
    for(IndexType param_id = 0; param_id < param_dim; param_id++)
    {
        all_series[var_dim + param_id] = new Series(max_order, var_dim + param_dim + 1);
        const_params.series_vec[param_id]->append_new_variable(*(all_series[var_dim + param_id]));
    }

    // 2. do composition and integrate
    dvars_time.reinit();
    for(IndexType var_id = 0; var_id < var_dim; var_id ++)
    {
        for(IndexType curr_order = 0; curr_order < max_order; curr_order++)
        {
            Homogen curr_homogen(var_dim + param_dim + 1, curr_order);
            series_comp(*(ode_fun.series_vec[var_id]), all_series, curr_order, 
                curr_homogen);
            while(curr_homogen.term_tree != NULL)
            {
                PolyTerm* curr_term = curr_homogen.pop_first_term();
                curr_term -> integrate_self(var_dim + param_dim);
                if(curr_term->order < max_order)
                {
                    dvars_time.series_vec[var_id] -> add_term(curr_term);
                }
            }
        }

        // add \tau = 0 term
        Series zero_time_pos(max_order, var_dim+param_dim+1);
        sol.series_vec[var_id]->append_new_variable(zero_time_pos);
        dvars_time.series_vec[var_id]->add_series(zero_time_pos);
    }
    
    // 3. move dvars_time to sol_time
    sol_time.reinit();
    dvars_time.copy_to(sol_time);

    // 4. free temp variables
    for(IndexType param_id = 0; param_id < param_dim; param_id++)
    {
        delete all_series[var_dim + param_id];
    }
}

void cJetTransport::Poincare_projection(IndexType err_order, ScalarVec & normal_vec)
{
    /*
        Calculate the projection of sol onto Poincare section
        defined by normal_vec to the err_order
    */
    assert(normal_vec.size() == var_dim);

    // 1. Preparation
    // get x1 and remove x1 terms in "sol"
    ScalarVec x1(var_dim);
    Scalar n_dot_F = 0.0;
    for(IndexType var_id = 0; var_id <var_dim; var_id ++)
    {
        if(sol.series_vec[var_id]->homogen_terms[0]->term_tree == NULL)
        {
            x1[var_id] = 0;
        }
        else
        {
            x1[var_id] = sol.series_vec[var_id]->homogen_terms[0]->term_tree->coeff;
            sol.series_vec[var_id]->homogen_terms[0]->reinit(var_dim + param_dim,0, true);
            sol.series_vec[var_id]->curr_kmin = 1;
        }
    }

    ScalarVec dx1(var_dim);
    for(IndexType var_id = 0; var_id < var_dim; var_id++)
    {
        ScalarVec temp_x0(var_dim + param_dim);
        for(IndexType k = 0; k < var_dim; k++)
        {
            temp_x0[k] = x1[k];
        }
        for(IndexType k = var_dim; k < var_dim + param_dim; k++)
        {
            temp_x0[k] = 0;
        }
        dx1[var_id] = ode_fun.series_vec[var_id] -> eval(temp_x0);
        n_dot_F += normal_vec[var_id] * dx1[var_id]; 
    }



    // calculate sol_time
    dvars.reinit();
    for(IndexType var_id = 0; var_id < var_dim; var_id ++)
    {
        dvars.series_vec[var_id] -> destructive_add_series(*(sol.series_vec[var_id]));
    }
    sol.reinit();

    IndexVec term_orders(var_dim + param_dim);
    for(auto it = term_orders.begin(); it != term_orders.end(); it++)
    {
        (*it) = 0;
    }
    for(IndexType var_id = 0; var_id < var_dim; var_id ++)
    {
        sol.series_vec[var_id] ->add_term(Monomial(x1[var_id], term_orders));
        term_orders[var_id] = 1;
        sol.series_vec[var_id] ->add_term(Monomial(1, term_orders));
        term_orders[var_id] = 0;
    }

    set_time_variable();
    for(IndexType curr_order = 0; curr_order < max_order; curr_order++)
    {
        derivative_iterate();
    }

    sol.reinit();
    for(IndexType var_id = 0; var_id < var_dim;var_id++)
    {
        sol.series_vec[var_id] -> destructive_add_series(*(dvars.series_vec[var_id]));
    }
    dvars.reinit();

    // set pointers of composition
    std::vector<Series*> series_ptr(var_dim + param_dim + 1);
    for(IndexType var_id = 0; var_id < var_dim; var_id ++ )
    {
        series_ptr[var_id] = sol.series_vec[var_id];
    }
    for(IndexType param_id = 0; param_id < param_dim; param_id ++)
    {
        series_ptr[var_dim + param_id] = const_params.series_vec[param_id];
    }
    series_ptr[var_dim + param_dim] = new Series(max_order, var_dim + param_dim);

    // 2. iteration
    for(IndexType curr_order = 1; curr_order < err_order; curr_order++)
    {
        // 2.1 calculate tau_k
        series_ptr[var_dim + param_dim] -> reinit();
        Homogen curr_tau(var_dim + param_dim, curr_order);
        for(IndexType var_id = 0; var_id < var_dim; var_id ++)
        {
            Homogen temp_homog(var_dim + param_dim, curr_order);
            sol.series_vec[var_id]->homogen_terms[curr_order]
                ->scalar_mul(normal_vec[var_id], temp_homog);
            curr_tau.destructive_add_self(temp_homog);
        }
        curr_tau.scalar_mul_self(- 1.0/n_dot_F);
        series_ptr[var_dim + param_dim] -> destructive_add_homogen(curr_tau);


        // 2.2 composition from sol to dvars
        dvars.reinit();
        for(IndexType var_id = 0; var_id < var_dim; var_id++)
        {
            for(IndexType comp_order = 1; comp_order < err_order; comp_order++)
            {
                Homogen temp_homog(var_dim + param_dim, comp_order);
                series_comp(*(sol_time.series_vec[var_id]), series_ptr, comp_order, 
                    temp_homog);
                dvars.series_vec[var_id] -> destructive_add_homogen(temp_homog);
            }
        }

        // 2.3 copy dvars to sol
        sol.reinit();
        for(IndexType var_id = 0; var_id < var_dim; var_id ++)
        {
            sol.series_vec[var_id]->destructive_add_series(*(dvars.series_vec[var_id]));
            for(IndexType comp_order = 1; comp_order < err_order; comp_order ++)
            {
                sol.series_vec[var_id] -> homogen_terms[comp_order]->remove_zeros();
            }

        }
    }

    // add constant term
    term_orders.clear();
    term_orders.resize(var_dim + param_dim);
    for(auto it = term_orders.begin(); it != term_orders.end(); it++)
    {
        (*it) = 0;
    }
    for(IndexType var_id =0; var_id < var_dim; var_id ++)
    {
        sol.series_vec[var_id] -> add_term(Monomial(x1[var_id], term_orders));
    }

    // free
    delete series_ptr[var_dim + param_dim];

}