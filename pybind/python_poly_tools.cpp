/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2023-03-02
 * Copyright © Department of Physics, Tsinghua University. All rights reserved
 */
#include "../src/basic_defs.h"
#include "../src/polys.h"
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
PYBIND11_MODULE(_poly_tools_cc, m)
#elif SCALAR_MODE == 1
PYBIND11_MODULE(_poly_tools_rc, m)
#else
PYBIND11_MODULE(_poly_tools_rr, m)
#endif
{
    py::bind_vector<IndexVec>(m, "CIndexVec");
    py::bind_vector<ScalarVec>(m, "CScalarVec");
    #if SCALAR_MODE == 1
        py::bind_vector<VarScalarVec>(m, "CVarScalarVec");
    #endif
    py::class_<PolyLinkedList>(m, "_CPolyLinkedList")
        .def(py::init<IndexType, bool>())
        .def("reinit", &PolyLinkedList::reinit)
        .def("remove_zeros", &PolyLinkedList::remove_zeros)
        .def("copy", &PolyLinkedList::copy)
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
        .def("eval", &PolyLinkedList::eval)
        .def("eval_diff", &PolyLinkedList::eval_diff)
        .def("batch_eval", &PolyLinkedList::batch_eval)
        .def("batch_add_elements", &PolyLinkedList::batch_add_elements)
        .def("batch_get_data", &PolyLinkedList::batch_get_data)
        .def("init_with_data", &PolyLinkedList::init_with_data)
        .def_readonly("dim", &PolyLinkedList::dim)
        .def_readonly("increasing_order", &PolyLinkedList::increasing_order)
        .def_readonly("n_terms", &PolyLinkedList::n_terms);

}