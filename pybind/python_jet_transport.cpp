/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2023-04-20
 * Copyright © Department of Physics, Tsinghua University. All rights reserved
 */

#include "python_jet_transport.h"


// set ode fun by poly
void cJetTransport::set_ode_fun(IndexType dim, PolyLinkedList & poly)
{
    PolyTerm * term_ptr = poly.term_tree;
    while(term_ptr != NULL)
    {
        ode_fun.series_vec[dim]->add_term(*term_ptr);
        term_ptr = term_ptr -> next;
    }
}

void cJetTransport::set_sol(IndexType dim, PolyLinkedList & poly)
{
    PolyTerm * term_ptr = poly.term_tree;
    while(term_ptr != NULL)
    {
        sol.series_vec[dim]->add_term(*term_ptr);
        term_ptr = term_ptr -> next;
    }
}

// get ode fun to poly
void cJetTransport::get_ode_fun(IndexType dim, IndexType order, PolyLinkedList & poly)
{
    ode_fun.series_vec[dim]->homogen_terms[order]->PolyLinkedList::copy(poly);
}

void cJetTransport::get_sol(IndexType dim, IndexType order, PolyLinkedList & poly)
{
    sol.series_vec[dim]->homogen_terms[order]->PolyLinkedList::copy(poly);
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