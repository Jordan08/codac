/** 
 *  Codac binding (core)
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <codac2_template_tools.h>
#include <codac2_SlicedTube.h>
#include <codac2_CtcEval.h>
#include "codac2_py_Ctc.h"
#include "codac2_py_CtcEval_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py):

using namespace std;
using namespace codac2;
namespace py = pybind11;
using namespace pybind11::literals;

void export_CtcEval(py::module& m)
{
  py::class_<CtcEval> exported(m, "CtcEval", CTCEVAL_MAIN);
  exported

    .def(py::init<>(),
      CTCEVAL_CTCEVAL)

    .def("contract", (void (CtcEval::*)(Interval&,Interval&,SlicedTube<Interval>&,const SlicedTube<Interval>&) const)&CtcEval::contract,
      VOID_CTCEVAL_CONTRACT_INTERVAL_REF_T_REF_SLICEDTUBE_T_REF_CONST_SLICEDTUBE_T_REF_CONST,
      "t"_a, "y"_a, "x"_a, "v"_a)

    .def("contract", (void (CtcEval::*)(Interval&,IntervalVector&,SlicedTube<IntervalVector>&,const SlicedTube<IntervalVector>&) const)&CtcEval::contract,
      VOID_CTCEVAL_CONTRACT_INTERVAL_REF_T_REF_SLICEDTUBE_T_REF_CONST_SLICEDTUBE_T_REF_CONST,
      "t"_a, "y"_a, "x"_a, "v"_a)

  ;
}