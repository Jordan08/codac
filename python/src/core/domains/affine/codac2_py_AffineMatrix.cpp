/**
 * \file codac2_py_AffineMatrix.cpp
 *
 * Python binding of codac2::AffineMatrix.
 */

#include "codac2_py_Affine.h"

#include <sstream>

#include "codac2_Affine.h"

namespace py = pybind11;

namespace codac2 {

void export_AffineMatrix(py::module_& m)
{
  using M = AffineMatrix;

  py::class_<M> cls(m, "AffineMatrix");

  cls
    .def(py::init<>())
    .def(py::init<Index, Index>(), py::arg("rows"), py::arg("cols"))
    .def("rows", [](const M& x) { return static_cast<Index>(x.rows()); })
    .def("cols", [](const M& x) { return static_cast<Index>(x.cols()); })
    .def("size", [](const M& x) { return static_cast<Index>(x.size()); })
    .def("__getitem__", [](const M& x, std::pair<Index, Index> ij) -> const Affine& {
      if (ij.first < 0 || ij.second < 0 ||
          ij.first >= static_cast<Index>(x.rows()) ||
          ij.second >= static_cast<Index>(x.cols()))
        throw py::index_error();
      return x(ij.first, ij.second);
    }, py::return_value_policy::reference_internal)
    .def("__setitem__", [](M& x, std::pair<Index, Index> ij, const Affine& value) {
      if (ij.first < 0 || ij.second < 0 ||
          ij.first >= static_cast<Index>(x.rows()) ||
          ij.second >= static_cast<Index>(x.cols()))
        throw py::index_error();
      x(ij.first, ij.second) = value;
    })
    .def("__str__", [](const M& x) {
      std::ostringstream os; os << x; return os.str();
    })
    .def("__repr__", [](const M& x) {
      std::ostringstream os; os << x; return os.str();
    });
}

} // namespace codac2
