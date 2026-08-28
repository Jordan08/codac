/**
 *  \file codac2_py_Affine_operations.cpp
 *
 *  Affine Python binding (arithmetic operators and functions)
 *
 * ----------------------------------------------------------------------------
 *  \date       2026
 *  \author     Jordan Ninin
 *  \copyright  Copyright 2026 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <sstream>
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <codac2_Affine.h>
#include <codac2_Interval.h>
#include "codac2_py_AffineMain_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py):
#include "codac2_py_AffineMain_operations_impl_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py):
#include "codac2_py_matlab.h"

using namespace std;
using namespace codac2;
namespace py = pybind11;
using namespace pybind11::literals;

void export_Affine_operations(py::module& m, py::class_<Affine>& py_Affine)
{
  // Members functions

  // Note: only the intersection and the hull of the interval enclosures are
  // provided (they return an Interval). The in-place operators &= and |= are
  // explicitly deleted for affine forms, and are therefore not exported.

  if constexpr(!FOR_MATLAB)
  {
    py_Affine

      .def(py::self & py::self,
        INTERVAL_OPERATORINTER_CONST_AFFINEMAIN_T_REF_CONST_AFFINEMAIN_T_REF,
        "x"_a)

      .def(py::self & Interval(),
        INTERVAL_OPERATORINTER_CONST_AFFINEMAIN_T_REF_CONST_INTERVAL_REF,
        "x"_a)

      .def(Interval() & py::self,
        INTERVAL_OPERATORINTER_CONST_INTERVAL_REF_CONST_AFFINEMAIN_T_REF,
        "x"_a)

      .def(py::self | py::self,
        INTERVAL_OPERATORUNION_CONST_AFFINEMAIN_T_REF_CONST_AFFINEMAIN_T_REF,
        "x"_a)

      .def(py::self | Interval(),
        INTERVAL_OPERATORUNION_CONST_AFFINEMAIN_T_REF_CONST_INTERVAL_REF,
        "x"_a)

      .def(Interval() | py::self,
        INTERVAL_OPERATORUNION_CONST_INTERVAL_REF_CONST_AFFINEMAIN_T_REF,
        "x"_a)
    ;
  }

  if constexpr(FOR_MATLAB)
  {
    // For MATLAB compatibility
    py_Affine

      .def("inter", (Interval(*)(const Affine&,const Affine&)) &codac2::operator&,
        INTERVAL_OPERATORINTER_CONST_AFFINEMAIN_T_REF_CONST_AFFINEMAIN_T_REF,
        "x"_a)

      .def("inter", (Interval(*)(const Affine&,const Interval&)) &codac2::operator&,
        INTERVAL_OPERATORINTER_CONST_AFFINEMAIN_T_REF_CONST_INTERVAL_REF,
        "x"_a)

      .def("union", (Interval(*)(const Affine&,const Affine&)) &codac2::operator|,
        INTERVAL_OPERATORUNION_CONST_AFFINEMAIN_T_REF_CONST_AFFINEMAIN_T_REF,
        "x"_a)

      .def("union", (Interval(*)(const Affine&,const Interval&)) &codac2::operator|,
        INTERVAL_OPERATORUNION_CONST_AFFINEMAIN_T_REF_CONST_INTERVAL_REF,
        "x"_a)
    ;
  }

  py_Affine

  .def(py::self + py::self,
    AFFINEMAIN_T_OPERATORPLUS_CONST_AFFINEMAIN_T_REF_CONST_AFFINEMAIN_T_REF,
    "x"_a)

  .def(py::self + double(),
    AFFINEMAIN_T_OPERATORPLUS_CONST_AFFINEMAIN_T_REF_DOUBLE,
    "d"_a)

  .def(double() + py::self,
    AFFINEMAIN_T_OPERATORPLUS_DOUBLE_CONST_AFFINEMAIN_T_REF,
    "d"_a)

  .def(py::self + Interval(),
    AFFINEMAIN_T_OPERATORPLUS_CONST_AFFINEMAIN_T_REF_CONST_INTERVAL_REF,
    "x"_a)

  .def(Interval() + py::self,
    AFFINEMAIN_T_OPERATORPLUS_CONST_INTERVAL_REF_CONST_AFFINEMAIN_T_REF,
    "x"_a)

  .def(py::self - py::self,
    AFFINEMAIN_T_OPERATORMINUS_CONST_AFFINEMAIN_T_REF_CONST_AFFINEMAIN_T_REF,
    "x"_a)

  .def(py::self - double(),
    AFFINEMAIN_T_OPERATORMINUS_CONST_AFFINEMAIN_T_REF_DOUBLE,
    "d"_a)

  .def(double() - py::self,
    AFFINEMAIN_T_OPERATORMINUS_DOUBLE_CONST_AFFINEMAIN_T_REF,
    "d"_a)

  .def(py::self - Interval(),
    AFFINEMAIN_T_OPERATORMINUS_CONST_AFFINEMAIN_T_REF_CONST_INTERVAL_REF,
    "x"_a)

  .def(Interval() - py::self,
    AFFINEMAIN_T_OPERATORMINUS_CONST_INTERVAL_REF_CONST_AFFINEMAIN_T_REF,
    "x"_a)

  .def(py::self * py::self,
    AFFINEMAIN_T_OPERATORMUL_CONST_AFFINEMAIN_T_REF_CONST_AFFINEMAIN_T_REF,
    "x"_a)

  .def(py::self * double(),
    AFFINEMAIN_T_OPERATORMUL_CONST_AFFINEMAIN_T_REF_DOUBLE,
    "d"_a)

  .def(double() * py::self,
    AFFINEMAIN_T_OPERATORMUL_DOUBLE_CONST_AFFINEMAIN_T_REF,
    "d"_a)

  .def(py::self * Interval(),
    AFFINEMAIN_T_OPERATORMUL_CONST_AFFINEMAIN_T_REF_CONST_INTERVAL_REF,
    "x"_a)

  .def(Interval() * py::self,
    AFFINEMAIN_T_OPERATORMUL_CONST_INTERVAL_REF_CONST_AFFINEMAIN_T_REF,
    "x"_a)

  .def(py::self / py::self,
    AFFINEMAIN_T_OPERATORDIV_CONST_AFFINEMAIN_T_REF_CONST_AFFINEMAIN_T_REF,
    "x"_a)

  .def(py::self / double(),
    AFFINEMAIN_T_OPERATORDIV_CONST_AFFINEMAIN_T_REF_DOUBLE,
    "d"_a)

  .def(double() / py::self,
    AFFINEMAIN_T_OPERATORDIV_DOUBLE_CONST_AFFINEMAIN_T_REF,
    "d"_a)

  .def(py::self / Interval(),
    AFFINEMAIN_T_OPERATORDIV_CONST_AFFINEMAIN_T_REF_CONST_INTERVAL_REF,
    "x"_a)

  // Note: operator/(const Interval&, const AffineMain<T>&) is declared in
  // codac2_AffineMain.h but cannot be instantiated (its implementation calls
  // the private method Ainv() without being declared as a friend, unlike
  // operator/(double, const AffineMain<T>&)). It is therefore not exported.

  ;

  // External functions outside the Affine class

  m

  .def("inv", (Affine(*)(const Affine&)) &codac2::inv,
    AFFINEMAIN_T_INV_CONST_AFFINEMAIN_T_REF,
    "x"_a)

  .def("sqr", (Affine(*)(const Affine&)) &codac2::sqr,
    AFFINEMAIN_T_SQR_CONST_AFFINEMAIN_T_REF,
    "x"_a)

  .def("sqrt", (Affine(*)(const Affine&)) &codac2::sqrt,
    AFFINEMAIN_T_SQRT_CONST_AFFINEMAIN_T_REF,
    "x"_a)

  .def("exp", (Affine(*)(const Affine&)) &codac2::exp,
    AFFINEMAIN_T_EXP_CONST_AFFINEMAIN_T_REF,
    "x"_a)

  .def("log", (Affine(*)(const Affine&)) &codac2::log,
    AFFINEMAIN_T_LOG_CONST_AFFINEMAIN_T_REF,
    "x"_a)

  .def("pow", (Affine(*)(const Affine&,int)) &codac2::pow,
    AFFINEMAIN_T_POW_CONST_AFFINEMAIN_T_REF_INT,
    "x"_a, "n"_a)

  .def("pow", (Affine(*)(const Affine&,double)) &codac2::pow,
    AFFINEMAIN_T_POW_CONST_AFFINEMAIN_T_REF_DOUBLE,
    "x"_a, "d"_a)

  .def("pow", (Affine(*)(const Affine&,const Interval&)) &codac2::pow,
    AFFINEMAIN_T_POW_CONST_AFFINEMAIN_T_REF_CONST_INTERVAL_REF,
    "x"_a, "y"_a)

  .def("pow", (Affine(*)(const Affine&,const Affine&)) &codac2::pow,
    AFFINEMAIN_T_POW_CONST_AFFINEMAIN_T_REF_CONST_AFFINEMAIN_T_REF,
    "x"_a, "y"_a)

  .def("pow", (Affine(*)(double,const Affine&)) &codac2::pow,
    AFFINEMAIN_T_POW_DOUBLE_CONST_AFFINEMAIN_T_REF,
    "d"_a, "x"_a)

  .def("pow", (Affine(*)(const Interval&,const Affine&)) &codac2::pow,
    AFFINEMAIN_T_POW_CONST_INTERVAL_REF_CONST_AFFINEMAIN_T_REF,
    "x"_a, "y"_a)

  .def("root", (Affine(*)(const Affine&,int)) &codac2::root,
    AFFINEMAIN_T_ROOT_CONST_AFFINEMAIN_T_REF_INT,
    "x"_a, "n"_a)

  .def("cos", (Affine(*)(const Affine&)) &codac2::cos,
    AFFINEMAIN_T_COS_CONST_AFFINEMAIN_T_REF,
    "x"_a)

  .def("sin", (Affine(*)(const Affine&)) &codac2::sin,
    AFFINEMAIN_T_SIN_CONST_AFFINEMAIN_T_REF,
    "x"_a)

  .def("tan", (Affine(*)(const Affine&)) &codac2::tan,
    AFFINEMAIN_T_TAN_CONST_AFFINEMAIN_T_REF,
    "x"_a)

  .def("acos", (Affine(*)(const Affine&)) &codac2::acos,
    AFFINEMAIN_T_ACOS_CONST_AFFINEMAIN_T_REF,
    "x"_a)

  .def("asin", (Affine(*)(const Affine&)) &codac2::asin,
    AFFINEMAIN_T_ASIN_CONST_AFFINEMAIN_T_REF,
    "x"_a)

  .def("atan", (Affine(*)(const Affine&)) &codac2::atan,
    AFFINEMAIN_T_ATAN_CONST_AFFINEMAIN_T_REF,
    "x"_a)

  .def("atan2", (Affine(*)(const Affine&,const Affine&)) &codac2::atan2,
    AFFINEMAIN_T_ATAN2_CONST_AFFINEMAIN_T_REF_CONST_AFFINEMAIN_T_REF,
    "y"_a, "x"_a)

  .def("cosh", (Affine(*)(const Affine&)) &codac2::cosh,
    AFFINEMAIN_T_COSH_CONST_AFFINEMAIN_T_REF,
    "x"_a)

  .def("sinh", (Affine(*)(const Affine&)) &codac2::sinh,
    AFFINEMAIN_T_SINH_CONST_AFFINEMAIN_T_REF,
    "x"_a)

  .def("tanh", (Affine(*)(const Affine&)) &codac2::tanh,
    AFFINEMAIN_T_TANH_CONST_AFFINEMAIN_T_REF,
    "x"_a)

  .def("acosh", (Affine(*)(const Affine&)) &codac2::acosh,
    AFFINEMAIN_T_ACOSH_CONST_AFFINEMAIN_T_REF,
    "x"_a)

  .def("asinh", (Affine(*)(const Affine&)) &codac2::asinh,
    AFFINEMAIN_T_ASINH_CONST_AFFINEMAIN_T_REF,
    "x"_a)

  .def("atanh", (Affine(*)(const Affine&)) &codac2::atanh,
    AFFINEMAIN_T_ATANH_CONST_AFFINEMAIN_T_REF,
    "x"_a)

  .def("abs", (Affine(*)(const Affine&)) &codac2::abs,
    AFFINEMAIN_T_ABS_CONST_AFFINEMAIN_T_REF,
    "x"_a)

  .def("sign", (Affine(*)(const Affine&)) &codac2::sign,
    AFFINEMAIN_T_SIGN_CONST_AFFINEMAIN_T_REF,
    "x"_a)

  .def("min", (Interval(*)(const Affine&,const Affine&)) &codac2::min,
    INTERVAL_MIN_CONST_AFFINEMAIN_T_REF_CONST_AFFINEMAIN_T_REF,
    "x"_a, "y"_a)

  .def("min", (Interval(*)(const Interval&,const Affine&)) &codac2::min,
    INTERVAL_MIN_CONST_INTERVAL_REF_CONST_AFFINEMAIN_T_REF,
    "x"_a, "y"_a)

  .def("min", (Interval(*)(const Affine&,const Interval&)) &codac2::min,
    INTERVAL_MIN_CONST_AFFINEMAIN_T_REF_CONST_INTERVAL_REF,
    "x"_a, "y"_a)

  .def("max", (Interval(*)(const Affine&,const Affine&)) &codac2::max,
    INTERVAL_MAX_CONST_AFFINEMAIN_T_REF_CONST_AFFINEMAIN_T_REF,
    "x"_a, "y"_a)

  .def("max", (Interval(*)(const Interval&,const Affine&)) &codac2::max,
    INTERVAL_MAX_CONST_INTERVAL_REF_CONST_AFFINEMAIN_T_REF,
    "x"_a, "y"_a)

  .def("max", (Interval(*)(const Affine&,const Interval&)) &codac2::max,
    INTERVAL_MAX_CONST_AFFINEMAIN_T_REF_CONST_INTERVAL_REF,
    "x"_a, "y"_a)

  .def("integer", (Interval(*)(const Affine&)) &codac2::integer,
    INTERVAL_INTEGER_CONST_AFFINEMAIN_T_REF,
    "x"_a)

  .def("floor", (Interval(*)(const Affine&)) &codac2::floor,
    INTERVAL_FLOOR_CONST_AFFINEMAIN_T_REF,
    "x"_a)

  .def("ceil", (Interval(*)(const Affine&)) &codac2::ceil,
    INTERVAL_CEIL_CONST_AFFINEMAIN_T_REF,
    "x"_a)

  // Note: the distance() functions declared in codac2_AffineMain.h cannot be
  // instantiated: they rely on a distance(Interval,Interval) function that is
  // not provided by codac2. They are therefore not exported here.

  .def("chi", (Affine(*)(const Affine&,const Affine&,const Affine&)) &codac2::chi,
    AFFINEMAIN_T_CHI_CONST_AFFINEMAIN_T_REF_CONST_AFFINEMAIN_T_REF_CONST_AFFINEMAIN_T_REF,
    "a"_a, "b"_a, "c"_a)

  .def("chi", (Affine(*)(const Interval&,const Affine&,const Affine&)) &codac2::chi,
    AFFINEMAIN_T_CHI_CONST_INTERVAL_REF_CONST_AFFINEMAIN_T_REF_CONST_AFFINEMAIN_T_REF,
    "a"_a, "b"_a, "c"_a)

  .def("chi", (Affine(*)(const Interval&,const Interval&,const Affine&)) &codac2::chi,
    AFFINEMAIN_T_CHI_CONST_INTERVAL_REF_CONST_INTERVAL_REF_CONST_AFFINEMAIN_T_REF,
    "a"_a, "b"_a, "c"_a)

  .def("chi", (Affine(*)(const Interval&,const Affine&,const Interval&)) &codac2::chi,
    AFFINEMAIN_T_CHI_CONST_INTERVAL_REF_CONST_AFFINEMAIN_T_REF_CONST_INTERVAL_REF,
    "a"_a, "b"_a, "c"_a)

  .def("chi", (Affine(*)(const Affine&,const Interval&,const Affine&)) &codac2::chi,
    AFFINEMAIN_T_CHI_CONST_AFFINEMAIN_T_REF_CONST_INTERVAL_REF_CONST_AFFINEMAIN_T_REF,
    "a"_a, "b"_a, "c"_a)

  .def("chi", (Affine(*)(const Affine&,const Affine&,const Interval&)) &codac2::chi,
    AFFINEMAIN_T_CHI_CONST_AFFINEMAIN_T_REF_CONST_AFFINEMAIN_T_REF_CONST_INTERVAL_REF,
    "a"_a, "b"_a, "c"_a)

  ;
}
