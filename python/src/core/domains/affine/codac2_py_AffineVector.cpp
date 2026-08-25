/**
 * \file codac2_py_AffineVector.cpp
 *
 * Python binding of codac2::AffineVector.
 *
 * AffineVector is an Eigen column vector whose scalar type is Affine.
 */

#include "codac2_py_Affine.h"

#include <sstream>

#include "codac2_Affine.h"

namespace py = pybind11;

namespace codac2 {

void export_AffineVector(py::module_& m)
{
  using V = AffineVector;

  py::class_<V> cls(m, "AffineVector");

  cls
    .def(py::init<>())
    .def(py::init<Index>(), py::arg("size"))
    .def("size", [](const V& x) { return static_cast<Index>(x.size()); })
    .def("rows", [](const V& x) { return static_cast<Index>(x.rows()); })
    .def("cols", [](const V& x) { return static_cast<Index>(x.cols()); })
    .def("__len__", [](const V& x) { return static_cast<Index>(x.size()); })
    .def("__getitem__", [](const V& x, Index i) -> const Affine& {
      if (i < 0 || i >= static_cast<Index>(x.size()))
        throw py::index_error();
      return x(i);
    }, py::return_value_policy::reference_internal)
    .def("__setitem__", [](V& x, Index i, const Affine& value) {
      if (i < 0 || i >= static_cast<Index>(x.size()))
        throw py::index_error();
      x(i) = value;
    })
    .def("__str__", [](const V& x) {
      std::ostringstream os; os << x; return os.str();
    })
    .def("__repr__", [](const V& x) {
      std::ostringstream os; os << x; return os.str();
    });
}

} // namespace codac2
