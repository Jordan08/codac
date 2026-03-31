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
#include "codac2_py_CtcWrapper_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py):

#define CTCVISIBLE_MAIN "Contractor for visibility from a point relative to a segment."
#define CTCVISIBLE_INIT "Initialize CtcVisible with an observation point and an obstacle segment."
#define CTCVISIBLE_CONTRACT "Contract the box x based on visibility criteria."

#define CTCNOVISIBLE_MAIN "Contractor for the hidden zone (shadow) behind a segment."
#define CTCNOVISIBLE_INIT "Initialize CtcNoVisible with an observation point and an obstacle segment."
#define CTCNOVISIBLE_CONTRACT "Contract the box x to the hidden area."

using namespace std;
using namespace codac2;
namespace py = pybind11;
using namespace pybind11::literals;

void export_CtcVisibility(py::module& m, py::class_<CtcBase<IntervalVector>, pyCtcIntervalVector>& pyctc)
{
  // --- CtcVisible ---
  py::class_<CtcVisible> vis(m, "CtcVisible", pyctc, CTCVISIBLE_MAIN);
  vis
    .def(py::init<const Vector&, const Segment&>(), 
      CTCVISIBLE_INIT, "a"_a, "s"_a)
    .def(CONTRACT_BOX_METHOD(CtcVisible, CTCVISIBLE_CONTRACT));

  // --- CtcNoVisible ---
  py::class_<CtcNoVisible> nvis(m, "CtcNoVisible", pyctc, CTCNOVISIBLE_MAIN);
  nvis
    .def(py::init<const Vector&, const Segment&>(), 
      CTCNOVISIBLE_INIT, "a"_a, "s"_a)
    .def(CONTRACT_BOX_METHOD(CtcNoVisible, CTCNOVISIBLE_CONTRACT));
}