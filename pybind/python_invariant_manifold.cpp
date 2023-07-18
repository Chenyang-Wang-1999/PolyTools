/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2022-12-22
 * Copyright © Department of Physics, Tsinghua University. All rights reserved
 */

#include "../src/basic_defs.h"
#include "../src/invariant_manifold_solver.h"
#include "../src/polys.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl_bind.h>
#include <pybind11/eigen.h>

namespace py=pybind11;
PYBIND11_MAKE_OPAQUE(IndexVec);
PYBIND11_MAKE_OPAQUE(ScalarVec);

#if SCALAR_MODE == 1
    PYBIND11_MAKE_OPAQUE(VarScalarVec);
#endif 

#if SCALAR_MODE == 0
PYBIND11_MODULE(_invariant_manifold_cc, m)
#elif SCALAR_MODE == 1
PYBIND11_MODULE(_invariant_manifold_rc, m)
#else
PYBIND11_MODULE(_invariant_manifold_rr, m)
#endif
{
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

    py::class_<InvariantManifoldSolver>(m, "_CInvariantManifoldSolver")
        .def(py::init<IndexType, IndexType, IndexType>())
        .def_readwrite("phys_dim", &InvariantManifoldSolver::phys_dim)
        .def_readwrite("manifold_dim", &InvariantManifoldSolver::manifold_dim)
        .def_readwrite("Kmax", &InvariantManifoldSolver::Kmax)
        .def_readwrite("Kceil", &InvariantManifoldSolver::Kceil)
        .def_readwrite("resonance_tol", &InvariantManifoldSolver::resonance_tol)
        .def("init_without_T", &InvariantManifoldSolver::init_without_T)
        .def("jacobian_by_F", &InvariantManifoldSolver::jacobian_by_F)
        .def("add_term_F", &InvariantManifoldSolver::add_term_F)
        .def("get_curr_term", &InvariantManifoldSolver::get_curr_term)
        .def("set_data_accessor", &InvariantManifoldSolver::set_data_accessor)
        .def("load_data_and_move", &InvariantManifoldSolver::load_data_and_move)
        .def("print_poly_info", &InvariantManifoldSolver::print_poly_info)
        .def("get_poly_val_or_var_dim", &InvariantManifoldSolver::get_poly_val_or_var_dim)
        .def("solve_step", &InvariantManifoldSolver::solve_step)
        .def("solve_to_Kceil", &InvariantManifoldSolver::solve_to_Kceil)
        .def("eval", &InvariantManifoldSolver::eval)
        .def("get_k", &InvariantManifoldSolver::get_k)
        .def("calculate_err", &InvariantManifoldSolver::calculate_err)
        .def("clear_all", &InvariantManifoldSolver::clear_all)
        .def("get_poly", &InvariantManifoldSolver::get_poly)
        .def("init_with_T", &InvariantManifoldSolver::init);
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
}