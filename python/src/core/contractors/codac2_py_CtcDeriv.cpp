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
#include <codac2_CtcDeriv.h>
#include "codac2_py_Ctc.h"
#include "codac2_py_CtcDeriv_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py):

using namespace std;
using namespace codac2;
namespace py = pybind11;
using namespace pybind11::literals;

void export_CtcDeriv(py::module& m)
{
  py::class_<CtcDeriv> exported(m, "CtcDeriv", CTCDERIV_MAIN);
  exported

    .def(py::init<const TimePropag&,bool>(),
      CTCDERIV_CTCDERIV_CONST_TIMEPROPAG_REF_BOOL,
      "time_propag"_a=TimePropag::FWD_BWD, "fast_mode"_a=true)

    .def("restrict_tdomain", &CtcDeriv::restrict_tdomain,
      VOID_CTCDERIV_RESTRICT_TDOMAIN_CONST_INTERVAL_REF,
      "tdomain"_a)

    // Contractions on Slice objects

    .def("contract", (void (CtcDeriv::*)(Slice<Interval>&,const Slice<Interval>&,const std::vector<Index>&) const)&CtcDeriv::contract,
      VOID_CTCDERIV_CONTRACT_SLICE_T_REF_CONST_SLICE_T_REF_CONST_VECTOR_INDEX_REF_CONST,
      "x"_a, "v"_a, "ctc_indices"_a=std::vector<Index>())

    .def("contract", [](const CtcDeriv& ctc, Slice<IntervalVector>& x, const Slice<IntervalVector>& v, const std::vector<Index_type>& ctc_indices)
        {
          ctc.contract(x, v, matlab::convert_indices(ctc_indices));
        },
      VOID_CTCDERIV_CONTRACT_SLICE_T_REF_CONST_SLICE_T_REF_CONST_VECTOR_INDEX_REF_CONST,
      "x"_a, "v"_a, "ctc_indices"_a=std::vector<Index>())

    // Contractions on SlicedTube objects
    
    .def("contract", (void (CtcDeriv::*)(SlicedTube<Interval>&,const SlicedTube<Interval>&,const std::vector<Index>&) const)&CtcDeriv::contract,
      VOID_CTCDERIV_CONTRACT_SLICEDTUBE_T_REF_CONST_SLICEDTUBE_T_REF_CONST_VECTOR_INDEX_REF_CONST,
      "x"_a, "v"_a, "ctc_indices"_a=std::vector<Index>())

    .def("contract", [](const CtcDeriv& ctc, SlicedTube<IntervalVector>& x, SlicedTube<IntervalVector>& v, const std::vector<Index_type>& ctc_indices)
        {
          return ctc.contract(x, v, matlab::convert_indices(ctc_indices));
        },
      VOID_CTCDERIV_CONTRACT_SLICEDTUBE_T_REF_CONST_SLICEDTUBE_T_REF_CONST_VECTOR_INDEX_REF_CONST,
      "x"_a, "v"_a, "ctc_indices"_a=std::vector<Index>())

  ;
}