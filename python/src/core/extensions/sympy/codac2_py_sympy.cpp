/** 
 *  \file
 *  Codac binding (sympy)
 * ----------------------------------------------------------------------------
 *  \date       2026
 *  \author     Simon Rohou
 *  \copyright  Copyright 2026 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <pybind11/pybind11.h>
#include <codac2_AnalyticFunction.h>
#include <codac2_AnalyticType.h>
#include <codac2_sympy.h>
#include "codac2_py_sympy_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py):
#include "codac2_py_cast.h"

using namespace codac2;
namespace py = pybind11;


void export_sympy(py::module& m)
{
  m.def("sympy_diff_scalar", (AnalyticFunction<ScalarType> (*)(const AnalyticFunction<ScalarType>&))&codac2::sympy_diff,
    ANALYTICFUNCTION_SCALARTYPE_SYMPY_DIFF_CONST_ANALYTICFUNCTION_SCALARTYPE_REF,
    "f"_a);

  m.def("sympy_diff_vector", (AnalyticFunction<MatrixType> (*)(const AnalyticFunction<VectorType>&))&codac2::sympy_diff,
    ANALYTICFUNCTION_MATRIXTYPE_SYMPY_DIFF_CONST_ANALYTICFUNCTION_VECTORTYPE_REF,
    "f"_a);
}