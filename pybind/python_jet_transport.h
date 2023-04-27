/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2023-04-20
 * Copyright © Department of Physics, Tsinghua University. All rights reserved
 */

#ifndef PYTHON_JET_TRANSPORT_H
#define PYTHON_JET_TRANSPORT_H

#include "../src/polys.h"
#include "../src/poly_utils.h"
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl_bind.h>
// #include <pybind11/eigen.h>
#include <pybind11/iostream.h>

namespace py=pybind11;

void set_series_vector(SeriesVec & target_series_vec, IndexType dim, PolyLinkedList & poly);
void get_series_vector(const SeriesVec & source_series_vec, IndexType dim, IndexType order, PolyLinkedList & poly);

class cJetTransport
{
public:
    IndexType param_dim;
    IndexType var_dim;
    IndexType max_order;
    cJetTransport(IndexType var_dim, IndexType param_dim, IndexType max_order):
        param_dim(param_dim), var_dim(var_dim), max_order(max_order),
        // template_series(max_order, var_dim + param_dim), 
        const_params(var_dim+param_dim, param_dim, max_order),
        ode_fun(var_dim + param_dim, var_dim, max_order),
        sol(var_dim+param_dim, var_dim, max_order),
        sol_time(var_dim+param_dim + 1, var_dim, max_order),
        dvars(var_dim+param_dim, var_dim, max_order),
        dvars_time(var_dim+param_dim+1, var_dim, max_order),
        tau_function(max_order, var_dim + param_dim),
        order_start_indices(max_order+1)
        {
            // Initialize template series
            all_series_orders.clear();
            IndexType curr_start_id = 0;
            for(IndexType curr_order = 0; curr_order < max_order; curr_order ++)
            {
                order_start_indices[curr_order] = curr_start_id;
                IndexVec order_sep(var_dim + param_dim + 1), order_vec(var_dim + param_dim);
                IndexVec order_vec_min(var_dim + param_dim), order_vec_max(var_dim + param_dim);
                for(auto it = order_vec_min.begin(); it != order_vec_min.end(); it++)
                {
                    (*it) = 0;
                }
                for(auto it = order_vec_max.begin(); it != order_vec_max.end(); it++)
                {
                    (*it) = max_order;
                }
                
                boards_left_most(order_sep, 0, order_vec_min, order_vec_max, curr_order);
                do
                {
                    const_sum_boards_to_data(order_sep, order_vec);
                    // template_series.add_term(Monomial(1, order_vec));
                    all_series_orders.push_back(order_vec);
                    curr_start_id ++;
                }while(! const_sum_next(order_sep, order_vec_min, order_vec_max));
            }
            order_start_indices[max_order] = curr_start_id;
            assert(order_start_indices[max_order] == all_series_orders.size());

            // Initialize "const_params"
            IndexVec param_poly_indices(param_dim + var_dim);
            for(auto it = param_poly_indices.begin(); it!=param_poly_indices.end(); it++)
            {
                (*it) = 0;
            }
            for(IndexType param_id = 0; param_id < param_dim; param_id++ )
            {
                param_poly_indices[param_id + var_dim] = 1;
                const_params.series_vec[param_id]->add_term(Monomial(1, param_poly_indices));
                param_poly_indices[param_id + var_dim] = 0;
            }
        }
    
    // Accessor to solution and dvars
    void set_ode_fun(IndexType dim, PolyLinkedList & poly)
    {
        set_series_vector(ode_fun, dim, poly);
    }
    void get_ode_fun(IndexType dim, IndexType order, PolyLinkedList & poly)
    {
        get_series_vector(ode_fun, dim, order, poly);
    }
    void get_sol(IndexType dim, IndexType order, PolyLinkedList &poly)
    {
        get_series_vector(sol, dim, order, poly);
    }
    void get_dvars_time(IndexType dim, IndexType order, PolyLinkedList & poly)
    {
        get_series_vector(dvars_time, dim, order, poly);
    }
    void get_sol_time(IndexType dim, IndexType order, PolyLinkedList & poly)
    {
        get_series_vector(sol_time, dim, order, poly);
    }
    void set_sol(IndexType dim, PolyLinkedList & poly)
    {
        set_series_vector(sol, dim, poly);
    }
    void set_sol_by_scalars(ScalarVec & coeffs, IndexType err_order);
    void get_sol_as_scalars(ScalarVec & coeffs, IndexType err_order);
    void get_dvars_as_scalars(ScalarVec & coeffs, IndexType err_order);
    IndexVec get_series_indices(IndexType j){return all_series_orders[j];}
    
    // // calculation
    void calculate_dvars(IndexType err_order);

    void set_time_variable()
    {
        sol.append_new_variable(sol_time);
    }

    void derivative_iterate();

    void Poincare_projection(IndexType err_order, ScalarVec & normal_vec);

private:
    // Series template_series;
    std::vector<IndexVec> all_series_orders;
    SeriesVec const_params;
    SeriesVec ode_fun;
    SeriesVec sol;
    SeriesVec sol_time;
    SeriesVec dvars;
    SeriesVec dvars_time;
    Series tau_function;
    IndexVec order_start_indices;
};

PYBIND11_MAKE_OPAQUE(IndexVec);
PYBIND11_MAKE_OPAQUE(ScalarVec);
#if SCALAR_MODE == 1
    PYBIND11_MAKE_OPAQUE(VarScalarVec);
#endif 

#if SCALAR_MODE == 0
PYBIND11_MODULE(_jet_transport_cc, m)
#elif SCALAR_MODE == 1
PYBIND11_MODULE(_jet_transport_rc, m)
#else
PYBIND11_MODULE(_jet_transport_rr, m)
#endif
{
    py::scoped_ostream_redirect stream(
        STDOUT,
        py::module::import("sys").attr("stdout")
    );
    py::bind_vector<IndexVec>(m, "CIndexVec");
    py::bind_vector<ScalarVec>(m, "CScalarVec");
    #if SCALAR_MODE == 1
        py::bind_vector<VarScalarVec>(m, "CVarScalarVec");
    #endif
    py::class_<PolyLinkedList>(m, "_CPolyLinkedList")
        .def(py::init<IndexType>())
        .def("reinit", &PolyLinkedList::reinit)
        .def("remove_zeros", &PolyLinkedList::remove_zeros)
        .def("copy_to", &PolyLinkedList::copy_to)
        .def("print_info", &PolyLinkedList::print_info)
        .def("destructive_add_self", &PolyLinkedList::destructive_add_self)
        .def("destructive_add", &PolyLinkedList::destructive_add)
        .def("add_self", &PolyLinkedList::add_self)
        .def("add", &PolyLinkedList::add)
        .def("scalar_mul_self", &PolyLinkedList::scalar_mul_self)
        .def("scalar_mul", &PolyLinkedList::scalar_mul)
        .def("scale_var_self", &PolyLinkedList::scale_var_self)
        .def("scale_var", &PolyLinkedList::scale_var)
        .def("neg_self", &PolyLinkedList::neg_self)
        .def("neg", &PolyLinkedList::neg)
        .def("destructive_subs_self", &PolyLinkedList::destructive_subs_self)
        .def("destructive_subs", &PolyLinkedList::destructive_subs)
        .def("subs_self", &PolyLinkedList::subs_self)
        .def("subs", &PolyLinkedList::subs)
        .def("derivative", &PolyLinkedList::derivative)
        .def("eval", &PolyLinkedList::eval)
        .def("eval_diff", &PolyLinkedList::eval_diff)
        .def("batch_eval", &PolyLinkedList::batch_eval)
        .def("batch_add_elements", &PolyLinkedList::batch_add_elements)
        .def("batch_get_data", &PolyLinkedList::batch_get_data)
        .def("init_with_data", &PolyLinkedList::init_with_data)
        .def_readonly("dim", &PolyLinkedList::dim)
        .def_readonly("n_terms", &PolyLinkedList::n_terms);

    py::class_<Monomial>(m, "_Monomial")
        .def(py::init<Scalar, IndexVec>())
        .def("eval", &Monomial::eval)
        .def("__LT__", &Monomial::operator<)
        .def("__GT__", &Monomial::operator>)
        .def("__LE__", &Monomial::operator<=)
        .def("__GE__", &Monomial::operator>=)
        .def("__EQ__", &Monomial::operator==)
        .def("__mul__", &Monomial::operator*)
        .def("__pow__", &Monomial::power)
        .def("copy_to", &Monomial::copy)
        .def("get_order", &Monomial::var_order)
        .def("derivative", &Monomial::derivative)
        .def_readonly("dim", &Monomial::dim)
        .def_readwrite("coeff", &Monomial::coeff);
    
    py::class_<cJetTransport>(m, "_cJetTransport")
        .def(py::init<IndexType, IndexType, IndexType>())
        .def_readonly("param_dim", &cJetTransport::param_dim)
        .def_readonly("var_dim", &cJetTransport::var_dim)
        .def_readonly("max_order", &cJetTransport::max_order)
        .def("set_ode_fun", &cJetTransport::set_ode_fun)
        .def("get_ode_fun", &cJetTransport::get_ode_fun)
        .def("get_sol", &cJetTransport::get_sol)
        .def("set_sol", &cJetTransport::set_sol)
        .def("get_sol_time", &cJetTransport::get_sol_time)
        .def("get_dvars_time", &cJetTransport::get_dvars_time)
        .def("set_sol_by_scalars", &cJetTransport::set_sol_by_scalars)
        .def("get_dvars_as_scalars", &cJetTransport::get_dvars_as_scalars)
        .def("get_sol_as_scalars", &cJetTransport::get_sol_as_scalars)
        .def("get_series_indices", &cJetTransport::get_series_indices)
        .def("calculate_dvars", &cJetTransport::calculate_dvars)
        .def("set_time_variable", &cJetTransport::set_time_variable)
        .def("Poincare_projection", &cJetTransport::Poincare_projection)
        .def("derivative_iterate", &cJetTransport::derivative_iterate);
}

#endif