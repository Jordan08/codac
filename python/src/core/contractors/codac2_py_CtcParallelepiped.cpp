/** 
 *  Codac binding (core)
 * ----------------------------------------------------------------------------
 *  \date       2026
 *  \author     Maël Godard
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <codac2_template_tools.h>
#include <codac2_CtcParallelepiped.h>
#include "codac2_py_Ctc.h"
#include "codac2_py_CtcParallelepiped_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py):

using namespace std;
using namespace codac2;
namespace py = pybind11;
using namespace pybind11::literals;

void export_CtcParallelepiped(py::module& m, py::class_<CtcBase<IntervalVector>,pyCtcIntervalVector>& pyctc)
{
  py::class_<CtcParallelepiped> exported(m, "CtcParallelepiped", pyctc, CTCPARALLELEPIPED_MAIN);
  exported

    .def(py::init<const Parallelepiped&>(),
      CTCPARALLELEPIPED_CTCPARALLELEPIPED_CONST_PARALLELEPIPED_REF,
      "p"_a)

  ;

  CONTRACT_METHODS(exported, CtcParallelepiped,
    VOID_CTCPARALLELEPIPED_CONTRACT_INTERVALVECTOR_REF_CONST)
}