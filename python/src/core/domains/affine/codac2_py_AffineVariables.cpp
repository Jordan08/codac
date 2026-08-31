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
#include "codac2_py_Matrix_addons_VectorBase_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py)
#include "codac2_py_MatrixBase_addons_AffineMainVector_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py)
#include "codac2_py_MatrixBase_addons_IntervalMatrixBase_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py)
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
        [](const AffineVariables& x, Index_type i) -> const Affine&
        {
          matlab::test_integer(i);
          // An IndexError is required here: this class provides no __iter__,
          // and Python then iterates by calling __getitem__ until IndexError.
          if(matlab::input_index(i) < 0 || matlab::input_index(i) >= x.size())
            throw py::index_error();
          // AffineVarMain<T> publicly inherits AffineMain<T> (Affine), so
          // this is a plain reference upcast (no slicing/copy): the
          // returned reference aliases the same component the container
          // owns, matching the reference semantics of __getitem__ on
          // AffineVector/AffineRow/AffineMatrix (see
          // codac2_py_AffineVector_templ.h / codac2_py_AffineMatrixBase.h).
          // Note that this also means the deleted AffineVarMain compound
          // assignments (+=,-=,*=,/=) and operator=(const AffineMain<T>&)
          // are reachable again through this reference: they are declared
          // non-virtual in AffineMain and only deleted on AffineVarMain
          // itself, so calling them through this AffineMain&-typed
          // reference bypasses that protection, same as it would in C++
          // through an AffineMain& alias of an AffineVarMain object.
          return x[matlab::input_index(i)];
        }, py::return_value_policy::reference_internal,
      MATRIX_ADDONS_VECTORBASE_CONST_SCALAR_REF_OPERATORCOMPO_INDEX_CONST)

    .def("get_item_0", [](const AffineVariables& x, Index_type i) -> const Affine&
        {
          matlab::test_integer(i);
          if(i < 0 || i >= x.size())
            throw py::index_error();
          return x[i];
        }, py::return_value_policy::reference_internal,
      MATRIX_ADDONS_VECTORBASE_CONST_SCALAR_REF_OPERATORCOMPO_INDEX_CONST)

    // The value type is deliberately Interval, not Affine: AffineVarMain
    // deletes operator=(const AffineMain<T>&) on purpose (see
    // codac2_AffineVar.h), so only assigning an Interval (or a double,
    // implicitly converted to Interval) reaches
    // AffineVarMain::operator=(const Interval&), which stores the radius
    // on this component's own dedicated noise symbol rather than in a
    // generic remainder error term.

    .def(
        #if FOR_MATLAB
          "set_item"
        #else
          "__setitem__"
        #endif
        , [](AffineVariables& x, Index_type i, const Interval& value)
        {
          matlab::test_integer(i);
          if(matlab::input_index(i) < 0 || matlab::input_index(i) >= x.size())
            throw py::index_error();
          x[matlab::input_index(i)] = value;
        },
      AFFINEVARMAIN_REF_AFFINEVARMAIN_T_OPERATORAFF_CONST_INTERVAL_REF)

    .def("set_item_0", [](AffineVariables& x, Index_type i, const Interval& value)
        {
          matlab::test_integer(i);
          if(i < 0 || i >= x.size())
            throw py::index_error();
          x[i] = value;
        },
      AFFINEVARMAIN_REF_AFFINEVARMAIN_T_OPERATORAFF_CONST_INTERVAL_REF)

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

    // lb()/ub()/mid() come from the generic IntervalMatrixBase Eigen addon:
    // AffineVarMain<T> is_interval_based (see codac2_AffineVar.h), so this
    // addon applies to AffineVarMainVector<T> exactly as it does to
    // AffineMainVector<T> (bound on AffineVector/AffineRow/AffineMatrix
    // through export_AffineMatrixBase(), see codac2_py_AffineMatrixBase.h).

    .def("lb", [](const AffineVariables& x)
        {
          return x.lb();
        },
      MATRIXBASE_ADDONS_INTERVALMATRIXBASE_AUTO_LB_CONST)

    .def("ub", [](const AffineVariables& x)
        {
          return x.ub();
        },
      MATRIXBASE_ADDONS_INTERVALMATRIXBASE_AUTO_UB_CONST)

    .def("mid", [](const AffineVariables& x)
        {
          return x.mid();
        },
      MATRIXBASE_ADDONS_INTERVALMATRIXBASE_AUTO_MID_CONST)

    .def("is_empty", [](const AffineVariables& x)
        {
          return x.is_empty();
        },
      MATRIXBASE_ADDONS_INTERVALMATRIXBASE_BOOL_IS_EMPTY_CONST)

    .def("__neg__", [](const AffineVariables& x) -> AffineVector
        {
          return -x;
        },
      AFFINEMAINVECTOR_T_AFFINEVARMAINVECTOR_T_OPERATORMINUS_CONST)

    // Products: an affine-variable vector is never modified in place, the
    // result of any product is a plain affine vector or matrix. These are
    // genuine matrix/vector products, so, matching Python/NumPy
    // convention, they are bound on "@" (__matmul__), not "*" (which is
    // reserved for scalar/elementwise multiplication).

    .def("__matmul__", [](const AffineVariables& x1, const Row& x2) -> AffineMatrix
        {
          return AffineVector(x1)*x2.cast<Affine>();
        },
      py::is_operator())

    .def("__matmul__", [](const AffineVariables& x1, const IntervalRow& x2) -> AffineMatrix
        {
          return AffineVector(x1)*x2.cast<Affine>();
        },
      py::is_operator())

    .def("__matmul__", [](const AffineVariables& x1, const AffineRow& x2) -> AffineMatrix
        {
          return AffineVector(x1)*x2;
        },
      py::is_operator())

    .def("__rmatmul__", [](const AffineVariables& x2, const Matrix& x1) -> AffineVector
        {
          return x1.cast<Affine>()*AffineVector(x2);
        },
      AUTO_OPERATORMUL_CONST_EIGEN_MATRIXBASE_OTHERDERIVED_REF_CONST_AFFINEVARMAINVECTOR_T_REF)

    .def("__rmatmul__", [](const AffineVariables& x2, const IntervalMatrix& x1) -> AffineVector
        {
          return x1.cast<Affine>()*AffineVector(x2);
        },
      AUTO_OPERATORMUL_CONST_EIGEN_MATRIXBASE_OTHERDERIVED_REF_CONST_AFFINEVARMAINVECTOR_T_REF)

    .def("__rmatmul__", [](const AffineVariables& x2, const AffineMatrix& x1) -> AffineVector
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
