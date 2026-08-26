/**
 *  \file codac2_py_AffineVector.cpp
 *
 *  AffineVector Python binding
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
#include "codac2_py_Matrix_addons_Base_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py)
#include "codac2_py_Matrix_addons_IntervalMatrixBase_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py)
#include "codac2_py_Matrix_addons_IntervalVector_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py)
#include "codac2_py_Matrix_addons_MatrixBase_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py)
#include "codac2_py_Matrix_addons_VectorBase_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py)
#include "codac2_py_MatrixBase_addons_Base_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py)
#include "codac2_py_MatrixBase_addons_IntervalMatrixBase_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py)
#include "codac2_py_MatrixBase_addons_VectorBase_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py)
#include "codac2_py_MatrixBase_addons_AffineMainVector_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py)
#include "codac2_py_AffineVector_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py)

#include "codac2_py_AffineVector_templ.h"

using namespace std;
using namespace codac2;
namespace py = pybind11;
using namespace pybind11::literals;

py::class_<AffineVector> export_AffineVector(py::module& m)
{
  py::class_<AffineVector> exported_affinevector_class(m, "AffineVector",
    USING_AFFINEMAINVECTOR_EQ_EIGEN_MATRIX_AFFINEMAIN_TMINUSONE1);
  export_AffineVector_<AffineVector,Vector,IntervalVector>(m, exported_affinevector_class);

  exported_affinevector_class

    // An affine-variable vector is converted into a plain affine vector

    .def(py::init(
        [](const AffineVariables& x)
        {
          return std::make_unique<AffineVector>(x);
        }),
      DOC_TO_BE_DEFINED,
      "x"_a)

    // Products with a row vector (the result is an affine matrix)

    .def("__mul__", [](const AffineVector& x1, const Row& x2) -> AffineMatrix
        {
          return x1*x2.cast<Affine>();
        },
      py::is_operator())

    .def("__mul__", [](const AffineVector& x1, const IntervalRow& x2) -> AffineMatrix
        {
          return x1*x2.cast<Affine>();
        },
      py::is_operator())

    .def("__mul__", [](const AffineVector& x1, const AffineRow& x2) -> AffineMatrix
        {
          return x1*x2;
        },
      py::is_operator())

    // Products with a matrix on the left-hand side

    .def("__rmul__", [](const AffineVector& x2, const Matrix& x1) -> AffineVector
        {
          return x1.cast<Affine>()*x2;
        },
      py::is_operator())

    .def("__rmul__", [](const AffineVector& x2, const IntervalMatrix& x1) -> AffineVector
        {
          return x1.cast<Affine>()*x2;
        },
      py::is_operator())

    .def("__repr__", [](const AffineVector& x)
        {
          std::ostringstream s;
          s << x;
          return string(s.str());
        },
      OSTREAM_REF_OPERATOROUT_OSTREAM_REF_CONST_AFFINEMAINVECTOR_T_REF)
  ;

  py::implicitly_convertible<AffineVariables,AffineVector>();

  return exported_affinevector_class;
}
