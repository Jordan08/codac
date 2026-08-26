/**
 *  \file codac2_py_AffineVariables.cpp
 *
 *  AffineVariables Python binding
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
#include <codac2_Vector.h>
#include <codac2_Row.h>
#include <codac2_Matrix.h>
#include <codac2_IntervalVector.h>
#include <codac2_IntervalRow.h>
#include <codac2_IntervalMatrix.h>
#include <codac2_Affine.h>

#include "codac2_py_doc.h"
#include "codac2_py_matlab.h"
#include "codac2_py_MatrixBase_addons_AffineMainVector_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py)
#include "codac2_py_AffineVar_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py)
#include "codac2_py_AffineVarVector_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py)

using namespace std;
using namespace codac2;
namespace py = pybind11;
using namespace pybind11::literals;

py::class_<AffineVariables> export_AffineVariables(py::module& m)
{
  py::class_<AffineVariables> exported_affinevariables_class(m, "AffineVariables",
    AFFINEVARMAINVECTOR_MAIN);

  exported_affinevariables_class

    .def(py::init(
        [](Index_type n)
        {
          matlab::test_integer(n);
          return std::make_unique<AffineVariables>((Index)n);
        }),
      AFFINEVARMAINVECTOR_T_AFFINEVARMAINVECTOR_INDEX,
      "n"_a)

    .def(py::init<const IntervalVector&>(),
      AFFINEVARMAINVECTOR_T_AFFINEVARMAINVECTOR_CONST_INTERVALVECTOR_REF,
      "x"_a)

    .def(py::init<const Vector&>(),
      AFFINEVARMAINVECTOR_T_AFFINEVARMAINVECTOR_CONST_VECTOR_REF,
      "x"_a)

    .def(py::init<const AffineVariables&>(),
      DOC_TO_BE_DEFINED,
      "x"_a)

    .def("__len__", [](const AffineVariables& x)
        {
          return x.size();
        },
      DOC_TO_BE_DEFINED)

    .def("size", [](const AffineVariables& x)
        {
          return x.size();
        },
      DOC_TO_BE_DEFINED)

    .def("rows", [](const AffineVariables& x)
        {
          return x.rows();
        },
      DOC_TO_BE_DEFINED)

    .def("cols", [](const AffineVariables& x)
        {
          return x.cols();
        },
      DOC_TO_BE_DEFINED)

    .def(
        #if FOR_MATLAB
          "__call__"
        #else
          "__getitem__"
        #endif
        ,
        [](const AffineVariables& x, Index_type i) -> Affine
        {
          matlab::test_integer(i);
          // An IndexError is required here: this class provides no __iter__,
          // and Python then iterates by calling __getitem__ until IndexError.
          if(matlab::input_index(i) < 0 || matlab::input_index(i) >= x.size())
            throw py::index_error();
          return Affine(x[matlab::input_index(i)]);
        },
      DOC_TO_BE_DEFINED)

    .def("get_item_0", [](const AffineVariables& x, Index_type i) -> Affine
        {
          matlab::test_integer(i);
          if(i < 0 || i >= x.size())
            throw py::index_error();
          return Affine(x[i]);
        },
      DOC_TO_BE_DEFINED)

    .def("noise_index", [](const AffineVariables& x, Index_type i)
        {
          matlab::test_integer(i);
          if(matlab::input_index(i) < 0 || matlab::input_index(i) >= x.size())
            throw py::index_error();
          return matlab::output_index(x[matlab::input_index(i)].noise_index());
        },
      INDEX_AFFINEVARMAIN_T_NOISE_INDEX_CONST,
      "i"_a)

    .def("resize", [](AffineVariables& x, Index_type n)
        {
          matlab::test_integer(n);
          x.resize((Index)n);
        },
      VOID_AFFINEVARMAINVECTOR_T_RESIZE_INDEX,
      "n"_a)

    .def("conservativeResize", [](AffineVariables& x, Index_type n)
        {
          matlab::test_integer(n);
          x.conservativeResize((Index)n);
        },
      VOID_AFFINEVARMAINVECTOR_T_CONSERVATIVERESIZE_INDEX,
      "n"_a)

    .def("init", [](AffineVariables& x, const Interval& i) -> AffineVariables&
        {
          return x.init(i);
        }, py::return_value_policy::reference_internal,
      AFFINEVARMAINVECTOR_T_REF_AFFINEVARMAINVECTOR_T_INIT_CONST_INTERVAL_REF,
      "x"_a)

    .def("itv", [](const AffineVariables& x)
        {
          return x.itv();
        },
      MATRIXBASE_ADDONS_AFFINEMAINVECTOR_AUTO_ITV_CONST)

    .def("is_empty", [](const AffineVariables& x)
        {
          return x.is_empty();
        },
      DOC_TO_BE_DEFINED)

    .def("__neg__", [](const AffineVariables& x) -> AffineVector
        {
          return -x;
        },
      AFFINEMAINVECTOR_T_AFFINEVARMAINVECTOR_T_OPERATORMINUS_CONST)

    // Products: an affine-variable vector is never modified in place, the
    // result of any product is a plain affine vector or matrix.

    .def("__mul__", [](const AffineVariables& x1, const Row& x2) -> AffineMatrix
        {
          return AffineVector(x1)*x2.cast<Affine>();
        },
      py::is_operator())

    .def("__mul__", [](const AffineVariables& x1, const IntervalRow& x2) -> AffineMatrix
        {
          return AffineVector(x1)*x2.cast<Affine>();
        },
      py::is_operator())

    .def("__mul__", [](const AffineVariables& x1, const AffineRow& x2) -> AffineMatrix
        {
          return AffineVector(x1)*x2;
        },
      py::is_operator())

    .def("__rmul__", [](const AffineVariables& x2, const Matrix& x1) -> AffineVector
        {
          return x1.cast<Affine>()*AffineVector(x2);
        },
      AUTO_OPERATORMUL_CONST_EIGEN_MATRIXBASE_OTHERDERIVED_REF_CONST_AFFINEVARMAINVECTOR_T_REF)

    .def("__rmul__", [](const AffineVariables& x2, const IntervalMatrix& x1) -> AffineVector
        {
          return x1.cast<Affine>()*AffineVector(x2);
        },
      AUTO_OPERATORMUL_CONST_EIGEN_MATRIXBASE_OTHERDERIVED_REF_CONST_AFFINEVARMAINVECTOR_T_REF)

    .def("__rmul__", [](const AffineVariables& x2, const AffineMatrix& x1) -> AffineVector
        {
          return x1*AffineVector(x2);
        },
      AUTO_OPERATORMUL_CONST_EIGEN_MATRIXBASE_OTHERDERIVED_REF_CONST_AFFINEVARMAINVECTOR_T_REF)

    .def("__repr__", [](const AffineVariables& x)
        {
          std::ostringstream s;
          s << x;
          return string(s.str());
        },
      OSTREAM_REF_OPERATOROUT_OSTREAM_REF_CONST_AFFINEVARMAINVECTOR_T_REF)
  ;

  return exported_affinevariables_class;
}
