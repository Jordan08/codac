/** * Codac binding (core)
 * ----------------------------------------------------------------------------
 *  \date       2026
 *  \author     Quentin Brateau
 *  \copyright  Copyright 2026 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <codac2_SepVisible.h>
#include "codac2_py_Sep.h"
#include "codac2_py_SepUnion_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py):

#define SEPVISIBILITY_MAIN "Separator for visibility characterization (Inside=Hidden, Outside=Visible)."
#define SEPVISIBILITY_INIT "Initialize SepVisible with an observation point 'a' and an obstacle segment 's'."
#define SEPVISIBILITY_SEPARATE "Separate the box into hidden and visible parts, returning a tuple (x_in, x_out)."

using namespace std;
using namespace codac2;
namespace py = pybind11;
using namespace pybind11::literals;

void export_SepVisible(py::module& m, py::class_<SepBase, pySep>& pysep)
{
  py::class_<SepVisible> exported(m, "SepVisible", pysep, SEPVISIBILITY_MAIN);
  exported
    .def(py::init<const Vector&, const Segment&>(),
      "a"_a, "s"_a, 
      SEPVISIBILITY_INIT)
    
    .def("separate", &SepVisible::separate, 
      "x"_a, 
      SEPVISIBILITY_SEPARATE)
  ;
}