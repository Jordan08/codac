/** 
 *  Approx Python binding
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <limits>
#include <sstream>
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <codac2_Interval.h>
#include <codac2_Vector.h>
#include <codac2_IntervalVector.h>
#include <codac2_Matrix.h>
#include <codac2_IntervalMatrix.h>
#include <codac2_Affine.h>
#include <codac2_Approx.h>
#include "codac2_py_Approx_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py):
#include "codac2_py_doc.h"

using namespace std;
using namespace codac2;
namespace py = pybind11;
using namespace pybind11::literals;

template<typename T>
void _export_Approx(py::module& m, const string& class_name)
{
  py::class_<Approx<T>> exported_class(m, class_name.c_str(), APPROX_MAIN);
  exported_class
  
    .def(py::init<const T&,double>(),
      APPROX_T_APPROX_CONST_T_REF_DOUBLE,
      "x"_a, "eps"_a = std::numeric_limits<double>::epsilon()*10)

    .def("__eq__", [](const Approx<T>& x1, const T& x2) { return x1 == x2; }, 
      BOOL_OPERATOREQ_CONST_T_REF_CONST_APPROX_T_REF)

    .def("__repr__", [](const Approx<T>& x) {
          std::ostringstream stream;
          stream << x;
          return string(stream.str()); 
        },
      OSTREAM_REF_OPERATOROUT_OSTREAM_REF_CONST_APPROX_T_REF)
  ;
}

/*
 * Approx<AffineMain<T>> and Approx<AffineMainVector<T>> (codac2_Approx.h)
 * are hand-written template specializations, not instances of the generic
 * Approx<T> template above: their constructor takes the *expected*
 * Interval/IntervalVector directly (not an AffineMain<T>/AffineMainVector<T>),
 * since an affine form's own type is not what the caller wants to spell out
 * at the comparison site. _export_Approx<T>() cannot be reused for them
 * (it assumes a T-typed constructor), hence these two dedicated exporters.
 */

void _export_Approx_Affine(py::module& m)
{
  py::class_<Approx<Affine>> exported_class(m, "Approx_Affine", APPROX_MAIN);
  exported_class

    .def(py::init<const Interval&,double>(),
      DOC_TO_BE_DEFINED,
      "x"_a, "eps"_a = std::numeric_limits<double>::epsilon()*10)

    .def("__eq__", [](const Approx<Affine>& x1, const Affine& x2) { return x1 == x2; },
      DOC_TO_BE_DEFINED)

    .def("__repr__", [](const Approx<Affine>& x) {
          std::ostringstream stream;
          stream << x;
          return string(stream.str());
        },
      DOC_TO_BE_DEFINED)
  ;
}

void _export_Approx_AffineVector(py::module& m)
{
  py::class_<Approx<AffineVector>> exported_class(m, "Approx_AffineVector", APPROX_MAIN);
  exported_class

    .def(py::init<const IntervalVector&,double>(),
      DOC_TO_BE_DEFINED,
      "x"_a, "eps"_a = std::numeric_limits<double>::epsilon()*10)

    .def("__eq__", [](const Approx<AffineVector>& x1, const AffineVector& x2) { return x1 == x2; },
      DOC_TO_BE_DEFINED)

    .def("__repr__", [](const Approx<AffineVector>& x) {
          std::ostringstream stream;
          stream << x;
          return string(stream.str());
        },
      DOC_TO_BE_DEFINED)
  ;
}

void export_Approx(py::module& m)
{
  _export_Approx<double>(m, "Approx_double");
  _export_Approx<Interval>(m, "Approx_Interval");
  _export_Approx<Vector>(m, "Approx_Vector");
  _export_Approx<IntervalVector>(m, "Approx_IntervalVector");
  _export_Approx<Row>(m, "Approx_Row");
  _export_Approx<IntervalRow>(m, "Approx_IntervalRow");
  _export_Approx<Matrix>(m, "Approx_Matrix");
  _export_Approx<IntervalMatrix>(m, "Approx_IntervalMatrix");
  _export_Approx<Segment>(m, "Approx_Segment");
  _export_Approx<Polygon>(m, "Approx_Polygon");
  _export_Approx<std::pair<Interval,Interval>>(m, "Approx_pair_Interval");
  _export_Approx_Affine(m);
  _export_Approx_AffineVector(m);
}