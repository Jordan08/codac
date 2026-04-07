/** * Codac binding (core)
 * ----------------------------------------------------------------------------
 *  \date       2026
 *  \author     Quentin Brateau
 *  \copyright  Copyright 2026 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <codac2_CtcVisible.h>
#include "codac2_py_Ctc.h"
#include "codac2_py_CtcVisible_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py):

using namespace std;
using namespace codac2;
namespace py = pybind11;
using namespace pybind11::literals;

void export_CtcVisible(py::module& m, py::class_<CtcBase<IntervalVector>, pyCtcIntervalVector>& pyctc)
{
  py::class_<CtcVisible> vis(m, "CtcVisible", pyctc, CTCVISIBLE_MAIN);
  vis
    .def(py::init<const IntervalVector&, const Segment&>(), 
      CTCVISIBLE_CTCVISIBLE_CONST_INTERVALVECTOR_REF_CONST_SEGMENT_REF, "a"_a, "s"_a)
    .def(py::init<const IntervalVector&, const std::vector<Segment>&>(), 
      CTCVISIBLE_CTCVISIBLE_CONST_INTERVALVECTOR_REF_CONST_VECTOR_SEGMENT_REF, "a"_a, "l"_a)
    .def(py::init<const IntervalVector&, const Polygon&>(), 
      CTCVISIBLE_CTCVISIBLE_CONST_INTERVALVECTOR_REF_CONST_POLYGON_REF, "a"_a, "p"_a)
    .def(CONTRACT_BOX_METHOD(CtcVisible, VOID_CTCVISIBLE_CONTRACT_INTERVALVECTOR_REF_CONST));
}

void export_CtcNoVisible(py::module& m, py::class_<CtcBase<IntervalVector>, pyCtcIntervalVector>& pyctc)
{
  py::class_<CtcNoVisible> nvis(m, "CtcNoVisible", pyctc, CTCNOVISIBLE_MAIN);
  nvis
    .def(py::init<const IntervalVector&, const Segment&>(), 
      CTCNOVISIBLE_CTCNOVISIBLE_CONST_INTERVALVECTOR_REF_CONST_SEGMENT_REF, "a"_a, "s"_a)
    .def(py::init<const IntervalVector&, const std::vector<Segment>&>(), 
      CTCNOVISIBLE_CTCNOVISIBLE_CONST_INTERVALVECTOR_REF_CONST_VECTOR_SEGMENT_REF, "a"_a, "l"_a)
    .def(py::init<const IntervalVector&, const Polygon&>(), 
      CTCNOVISIBLE_CTCNOVISIBLE_CONST_INTERVALVECTOR_REF_CONST_POLYGON_REF, "a"_a, "p"_a)
    .def(CONTRACT_BOX_METHOD(CtcNoVisible, VOID_CTCNOVISIBLE_CONTRACT_INTERVALVECTOR_REF_CONST));
}
