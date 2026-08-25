/**
 * \file codac2_py_AffineVariables.cpp
 *
 * Python binding of codac2::AffineVariables.
 *
 * AffineVariables is the public affine-symbol vector. The implementation
 * class AffineVarMain is intentionally not exported as a Python type.
 */

#include "codac2_py_Affine.h"

#include <sstream>

#include "codac2_Affine.h"

namespace py = pybind11;

namespace codac2 {

void export_AffineVariables(py::module_& m)
{
  using V = AffineVariables;

  py::class_<V> cls(m, "AffineVariables");

  cls
    .def(py::init<>())
    .def(py::init<Index>(), py::arg("size"))
    .def("size", [](const V& x) { return static_cast<Index>(x.size()); })
    .def("rows", [](const V& x) { return static_cast<Index>(x.rows()); })
    .def("cols", [](const V& x) { return static_cast<Index>(x.cols()); })
    .def("__len__", [](const V& x) { return static_cast<Index>(x.size()); })
    .def("__str__", [](const V& x) {
      std::ostringstream os; os << x; return os.str();
    })
    .def("__repr__", [](const V& x) {
      std::ostringstream os; os << x; return os.str();
    });
}

} // namespace codac2
