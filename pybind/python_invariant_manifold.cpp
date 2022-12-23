/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2022-12-22
 * Copyright © Department of Physics, Tsinghua University.  All rights reserved
 */

#include "../src/basic_defs.h"
#include "../src/invariant_manifold_solver.h"
#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl_bind.h>
#include <pybind11/eigen.h>

namespace py=pybind11;
PYBIND11_MAKE_OPAQUE(IndexVec);
PYBIND11_MAKE_OPAQUE(ScalarVec);

PYBIND11_MODULE(_invariant_manifold, m)
{
    py::bind_vector<IndexVec>(m, "CIndexVec");
    py::bind_vector<ScalarVec>(m, "CScalarVec");
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
        .def("init_with_T", &InvariantManifoldSolver::init);
    // py::class_<Monomial>(m, "_Monomial")
    //     .def(py::init<Scalar, IndexVec>())
    //     .def("eval", &Monomial::eval)
    //     .def("__LT__", &Monomial::operator<)
    //     .def("__GT__", &Monomial::operator>)
    //     .def("__LE__", &Monomial::operator<=)
    //     .def("__GE__", &Monomial::operator>=)
    //     .def("__EQ__", &Monomial::operator==)
    //     .def("__mul__", &Monomial::operator*)
    //     .def("__pow__", &Monomial::power)
    //     .def("copy", &Monomial::copy)
    //     .def("get_order", &Monomial::var_order)
    //     .def_readonly("dim", &Monomial::dim)
    //     .def_property("_coeff", &Monomial::get_coeff, &Monomial::set_coeff);
}