/*
 * @author        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
 * @date          2022-12-18
 * Copyright © Department of Physics, Tsinghua University. All rights reserved
 */

#include "../src/polys.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl_bind.h>

namespace py=pybind11;
PYBIND11_MAKE_OPAQUE(IndexVec);
PYBIND11_MAKE_OPAQUE(ScalarVec);

PYBIND11_MODULE(_polymanip, m)
{
    py::bind_vector<IndexVec>(m, "PolymanipIndexVec");
    py::bind_vector<ScalarVec>(m, "PolymanipScalarVec");
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
        .def("copy", &Monomial::copy)
        .def("get_order", &Monomial::var_order)
        .def_readonly("dim", &Monomial::dim)
        .def_property("_coeff", &Monomial::get_coeff, &Monomial::set_coeff);
}