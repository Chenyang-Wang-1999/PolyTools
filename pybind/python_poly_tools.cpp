/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2023-03-02
 * Copyright © Department of Physics, Tsinghua University.  All rights reserved
 */
#include "../src/polys.hpp"
#include "../src/poly_utils.hpp"
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl_bind.h>
#include <pybind11/eigen.h>
#include <pybind11/iostream.h>

namespace py=pybind11;
using namespace PolyTools;

PYBIND11_MAKE_OPAQUE(IndexVec);
PYBIND11_MAKE_OPAQUE(ScalarVec);
PYBIND11_MAKE_OPAQUE(std::vector<std::string>);
PYBIND11_MAKE_OPAQUE(std::vector<int>);
#if SCALAR_MODE == 1
    PYBIND11_MAKE_OPAQUE(VarScalarVec);
#endif 

void poly_comp_py(PolyLinkedList& f, SeriesVec & v, IndexType k, PolyLinkedList &res)
{
    poly_comp(f, v, k , res);
}

#if SCALAR_MODE == 0
PYBIND11_MODULE(_poly_tools_cc, m)
#elif SCALAR_MODE == 1
PYBIND11_MODULE(_poly_tools_rc, m)
#else
PYBIND11_MODULE(_poly_tools_rr, m)
#endif
{
    py::scoped_ostream_redirect stream(
        STDOUT,
        py::module::import("sys").attr("stdout")
    );
    py::bind_vector<std::vector<std::string>>(m, "CStrVec");
    py::bind_vector<IndexVec>(m, "CIndexVec");
    py::bind_vector<ScalarVec>(m, "CScalarVec");
    py::bind_vector<std::vector<int>>(m, "CLaurantIndexVec");
    py::class_<SeriesVec>(m,"_CSeriesVec")
        .def(py::init<IndexType, IndexType, IndexType>() )
        .def_readonly("var_dim", &SeriesVec::var_dim)
        .def_readonly("val_dim", &SeriesVec::val_dim)
        .def_readonly("Kmax", &SeriesVec::Kmax)
        .def("reinit", &SeriesVec::reinit)
        .def("copy_to", &SeriesVec::copy_to)
        .def("add_term", &SeriesVec::add_term)
        .def("get_poly", &SeriesVec::get_poly)
        .def("destructive_add_poly", &SeriesVec::destructive_add_poly);

    m.def("poly_comp", &poly_comp_py);
    m.def("poly_multiplication", &poly_multiplication);

    #if SCALAR_MODE == 1
        py::bind_vector<VarScalarVec>(m, "CVarScalarVec");
    #endif
    py::class_<PolyLinkedList>(m, "_CPolyLinkedList")
        .def(py::init<IndexType>())
        .def("reinit", &PolyLinkedList::reinit)
        .def("remove_zeros", &PolyLinkedList::remove_zeros)
        .def("copy_to", &PolyLinkedList::copy_to)
        .def("to_str", &PolyLinkedList::to_str)
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
        .def("partial_eval", &PolyLinkedList::partial_eval)
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
        .def("partial_eval", &Monomial::partial_eval)
        .def_readonly("dim", &Monomial::dim)
        .def_readwrite("coeff", &Monomial::coeff);
    
    py::class_<Laurant>(m, "_CLaurant")
        .def(py::init<IndexType>())
        .def_readonly("dim", &Laurant::dim)
        .def_readonly("num_max_orders", &Laurant::num_max_orders)
        .def_readonly("denom_orders", &Laurant::denom_orders)
        .def_readwrite("num", &Laurant::num)
        .def("reinit", &Laurant::reinit)
        .def("reduction", &Laurant::reduction)
        .def("set_Laurant", &Laurant::set_Laurant)
        .def("set_Laurant_by_terms", &Laurant::set_Laurant_by_terms)
        .def("eval", &Laurant::eval);
}