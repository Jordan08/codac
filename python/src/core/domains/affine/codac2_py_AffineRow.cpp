/**
 * \file codac2_py_AffineRow.cpp
 *
 * Python binding of codac2::AffineRow.
 */

#include "codac2_py_Affine.h"

#include <sstream>

#include "codac2_Affine.h"

namespace py = pybind11;

namespace codac2 {

void export_AffineRow(py::module_& m)
{
  using R = AffineRow;

  py::class_<R> cls(m, "AffineRow");

  cls
    .def(py::init<>())
    .def(py::init<Index>(), py::arg("size"))
    .def("size", [](const R& x) { return static_cast<Index>(x.size()); })
    .def("rows", [](const R& x) { return static_cast<Index>(x.rows()); })
    .def("cols", [](const R& x) { return static_cast<Index>(x.cols()); })
    .def("__len__", [](const R& x) { return static_cast<Index>(x.size()); })
    .def("__getitem__", [](const R& x, Index i) -> const Affine& {
      if (i < 0 || i >= static_cast<Index>(x.size()))
        throw py::index_error();
      return x(i);
    }, py::return_value_policy::reference_internal)
    .def("__setitem__", [](R& x, Index i, const Affine& value) {
      if (i < 0 || i >= static_cast<Index>(x.size()))
        throw py::index_error();
      x(i) = value;
    })
    .def("__str__", [](const R& x) {
      std::ostringstream os; os << x; return os.str();
    })
    .def("__repr__", [](const R& x) {
      std::ostringstream os; os << x; return os.str();
    });
}

} // namespace codac2
