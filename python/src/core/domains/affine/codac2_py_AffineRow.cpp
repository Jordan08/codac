/**
 *  \file codac2_py_AffineRow.cpp
 *
 *  AffineRow Python binding
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
#include "codac2_py_AffineRow_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py)

#include "codac2_py_AffineVector_templ.h"

using namespace std;
using namespace codac2;
namespace py = pybind11;
using namespace pybind11::literals;

py::class_<AffineRow> export_AffineRow(py::module& m)
{
  py::class_<AffineRow> exported_affinerow_class(m, "AffineRow",
    USING_AFFINEMAINROW_EQ_EIGEN_MATRIX_AFFINEMAIN_T1MINUSONE);
  export_AffineVector_<AffineRow,Row,IntervalRow>(m, exported_affinerow_class);

  exported_affinerow_class

    // Products with a matrix (the result is an affine row vector). These
    // are genuine matrix products, so, matching Python/NumPy convention,
    // they are bound on "@" (__matmul__), not "*" (__mul__ is reserved for
    // scalar/elementwise multiplication, see codac2_py_AffineMatrixBase.h).

    .def("__matmul__", [](const AffineRow& x1, const Matrix& x2) -> AffineRow
        {
          return x1*x2.cast<Affine>();
        },
      py::is_operator())

    .def("__matmul__", [](const AffineRow& x1, const IntervalMatrix& x2) -> AffineRow
        {
          return x1*x2.cast<Affine>();
        },
      py::is_operator())

    .def("__matmul__", [](const AffineRow& x1, const AffineMatrix& x2) -> AffineRow
        {
          return x1*x2;
        },
      py::is_operator())

    // Dot products with a column vector (the result is an affine form):
    // a 1 x n row times an n x 1 column is still a matrix product, hence
    // "@" rather than "*".

    .def("__matmul__", [](const AffineRow& x1, const Vector& x2) -> Affine
        {
          assert_release(x1.size() == x2.size());
          return (x1*x2.cast<Affine>()).eval()(0,0);
        },
      py::is_operator())

    .def("__matmul__", [](const AffineRow& x1, const IntervalVector& x2) -> Affine
        {
          assert_release(x1.size() == x2.size());
          return (x1*x2.cast<Affine>()).eval()(0,0);
        },
      py::is_operator())

    .def("__matmul__", [](const AffineRow& x1, const AffineVector& x2) -> Affine
        {
          assert_release(x1.size() == x2.size());
          return (x1*x2).eval()(0,0);
        },
      py::is_operator())

    .def("__repr__", [](const AffineRow& x)
        {
          std::ostringstream s;
          s << x;
          return string(s.str());
        },
      OSTREAM_REF_OPERATOROUT_OSTREAM_REF_CONST_AFFINEMAINROW_T_REF)
  ;

  return exported_affinerow_class;
}
