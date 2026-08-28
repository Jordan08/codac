/**
 *  \file codac2_py_AffineMatrix.cpp
 *
 *  AffineMatrix Python binding
 *
 * ----------------------------------------------------------------------------
 *  \date       2026
 *  \author     Jordan Ninin
 *  \copyright  Copyright 2026 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <sstream>
#include <vector>
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
#include "codac2_py_Matrix_addons_MatrixBase_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py)
#include "codac2_py_MatrixBase_addons_Base_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py)
#include "codac2_py_MatrixBase_addons_IntervalMatrixBase_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py)
#include "codac2_py_MatrixBase_addons_AffineMainVector_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py)
#include "codac2_py_AffineMatrix_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py)

#include "codac2_py_AffineMatrixBase.h"

using namespace std;
using namespace codac2;
namespace py = pybind11;
using namespace pybind11::literals;

py::class_<AffineMatrix> export_AffineMatrix(py::module& m)
{
  py::class_<AffineMatrix> exported_affinematrix_class(m, "AffineMatrix",
    USING_AFFINEMAINMATRIX_EQ_TYPEDEF_EIGEN_MATRIX_AFFINEMAIN_TMINUSONEMINUSONE);
  export_AffineMatrixBase<AffineMatrix,Matrix,false>(m, exported_affinematrix_class);

  exported_affinematrix_class

    .def(py::init<>())

    .def(py::init(
        [](Index_type r, Index_type c)
        {
          matlab::test_integer(r,c);
          return std::make_unique<AffineMatrix>((int)r,(int)c);
        }),
      DOC_TO_BE_DEFINED,
      "r"_a, "c"_a)

    .def(py::init<const AffineMatrix&>(),
      "x"_a)

    .def(py::init(
        [](const Matrix& x)
        {
          auto am = std::make_unique<AffineMatrix>((int)x.rows(),(int)x.cols());
          for(Index i = 0 ; i < x.rows() ; i++)
            for(Index j = 0 ; j < x.cols() ; j++)
              (*am)(i,j) = x(i,j);
          return am;
        }),
      DOC_TO_BE_DEFINED,
      "x"_a)

    .def(py::init(
        [](const IntervalMatrix& x)
        {
          auto am = std::make_unique<AffineMatrix>((int)x.rows(),(int)x.cols());
          for(Index i = 0 ; i < x.rows() ; i++)
            for(Index j = 0 ; j < x.cols() ; j++)
              (*am)(i,j) = x(i,j);
          return am;
        }),
      DOC_TO_BE_DEFINED,
      "x"_a)

    .def(py::init(
        [](const AffineRow& x)
        {
          auto am = std::make_unique<AffineMatrix>(1,(int)x.size());
          am->row(0) = x;
          return am;
        }),
      DOC_TO_BE_DEFINED,
      "x"_a)

    .def(py::init(
        [](const AffineVector& x)
        {
          auto am = std::make_unique<AffineMatrix>((int)x.size(),1);
          am->col(0) = x;
          return am;
        }),
      DOC_TO_BE_DEFINED,
      "x"_a)

    .def(py::init( // this constructor must be the last one to be declared
        [](const std::vector<AffineVector>& v)
        {
          assert_release(!std::empty(v));
          auto am = std::make_unique<AffineMatrix>((int)v.size(),(int)v[0].size());
          for(size_t i = 0 ; i < v.size() ; i++)
          {
            assert_release(v[i].size() == am->cols() && "AffineVector objects of different size");
            am->row(i) = v[i].transpose();
          }
          return am;
        }),
      DOC_TO_BE_DEFINED,
      "v"_a)

    .def_static("empty", [](Index_type r, Index_type c)
        {
          matlab::test_integer(r,c);
          return AffineMatrix::constant((int)r,(int)c,Affine::empty());
        },
      MATRIX_ADDONS_INTERVALMATRIXBASE_STATIC_AUTO_EMPTY_INDEX_INDEX,
      "r"_a, "c"_a)

    // Products with vectors and matrices

    .def("__mul__", [](const AffineMatrix& x1, const Vector& x2) -> AffineVector
        {
          return x1*x2.cast<Affine>();
        },
      py::is_operator())

    .def("__mul__", [](const AffineMatrix& x1, const IntervalVector& x2) -> AffineVector
        {
          return x1*x2.cast<Affine>();
        },
      py::is_operator())

    .def("__mul__", [](const AffineMatrix& x1, const AffineVector& x2) -> AffineVector
        {
          return x1*x2;
        },
      py::is_operator())

    .def("__mul__", [](const AffineMatrix& x1, const Matrix& x2) -> AffineMatrix
        {
          return x1*x2.cast<Affine>();
        },
      py::is_operator())

    .def("__mul__", [](const AffineMatrix& x1, const IntervalMatrix& x2) -> AffineMatrix
        {
          return x1*x2.cast<Affine>();
        },
      py::is_operator())

    .def("__mul__", [](const AffineMatrix& x1, const AffineMatrix& x2) -> AffineMatrix
        {
          return x1*x2;
        },
      py::is_operator())

    .def("__rmul__", [](const AffineMatrix& x2, const Matrix& x1) -> AffineMatrix
        {
          return x1.cast<Affine>()*x2;
        },
      py::is_operator())

    .def("__rmul__", [](const AffineMatrix& x2, const IntervalMatrix& x1) -> AffineMatrix
        {
          return x1.cast<Affine>()*x2;
        },
      py::is_operator())

    .def("__repr__", [](const AffineMatrix& x)
        {
          std::ostringstream s;
          s << x;
          return string(s.str());
        },
      OSTREAM_REF_OPERATOROUT_OSTREAM_REF_CONST_AFFINEMAINMATRIX_T_REF)
  ;

  return exported_affinematrix_class;
}
