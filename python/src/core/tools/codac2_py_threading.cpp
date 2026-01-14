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
#include <codac2_threading.h>
#include "codac2_py_threading_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py):

using namespace std;
using namespace codac2;
namespace py = pybind11;
using namespace pybind11::literals;

void export_threading(py::module& m)
{
  m.def("max_threads", &codac2::max_threads,
        INT_MAX_THREADS);
  
  m.def("set_threads_used", &codac2::set_threads_used,
        VOID_SET_THREADS_USED_INT,
        "n"_a);

  m.def("get_threads_used", &codac2::get_threads_used,
        INT_GET_THREADS_USED);
}