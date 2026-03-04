/** 
 *  Codac binding (core)
 * ----------------------------------------------------------------------------
 *  \date       2026
 *  \author     Auguste Bourgois, Simon Rohou, Maël Godard
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <codac2_template_tools.h>
#include <codac2_SlicedTube.h>
#include <codac2_CtcLohner.h>
#include "codac2_py_Ctc.h"
#include "codac2_py_CtcLohner_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py):
#include "codac2_py_cast.h"

using namespace std;
using namespace codac2;
namespace py = pybind11;
using namespace pybind11::literals;

void export_CtcLohner(py::module& m)
{
  py::class_<CtcLohner> exported(m, "CtcLohner", CTCLOHNER_MAIN);
  exported

    .def(py::init(
      [](const py::object& f, int contractions, double eps)
          {
            if(!is_instance<AnalyticFunction<VectorType>>(f)) {
              assert_release("CtcLohner: invalid function type");
            }
            return std::make_unique<CtcLohner>(cast<AnalyticFunction<VectorType>>(f), contractions, eps);
          }
        ),
        CTCLOHNER_CTCLOHNER_CONST_ANALYTICFUNCTION_VECTORTYPE_REF_INT_DOUBLE,
        "f"_a, "contractions"_a = 5, "eps"_a = 0.1)

    .def("contract", [](const CtcLohner& c, py::object& tube, TimePropag t_propa) -> py::object&
        {
          if(!is_instance<SlicedTube<IntervalVector>>(tube)) {
            assert_release("contract: invalid tube type");
          }
          c.contract(cast<SlicedTube<IntervalVector>>(tube), t_propa);
          return tube;
        },
        VOID_CTCLOHNER_CONTRACT_SLICEDTUBE_INTERVALVECTOR_REF_TIMEPROPAG_CONST,
        "tube"_a, "t_propa"_a = TimePropag::FWD_BWD)
  ;
}