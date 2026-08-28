/**
 *  \file codac2_py_Affine.cpp
 *
 *  Affine Python binding
 *
 * ----------------------------------------------------------------------------
 *  \date       2026
 *  \author     Jordan Ninin
 *  \copyright  Copyright 2026 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <sstream>
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <codac2_Affine.h>
#include <codac2_Interval.h>
#include "codac2_py_Affine_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py):
#include "codac2_py_AffineMain_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py):
#include "codac2_py_doc.h"
#include "codac2_py_matlab.h"

using namespace std;
using namespace codac2;
namespace py = pybind11;
using namespace pybind11::literals;

py::class_<Affine> export_Affine(py::module& m)
{
  py::class_<Affine> exported_affine_class(m, "Affine", AFFINEMAIN_MAIN);

  // Linearization mode of the affine forms (one mode per thread)

  py::enum_<Affine::Affine_Mode>(exported_affine_class, "Affine_Mode", DOC_TO_BE_DEFINED)

    .value("AF_Lin_Chebyshev", Affine::AF_Lin_Chebyshev)
    .value("AF_Lin_MinRange", Affine::AF_Lin_MinRange)
    .export_values()
  ;

  exported_affine_class

    .def(py::init<>(),
      AFFINEMAIN_T_AFFINEMAIN)

    .def(py::init<double>(),
      AFFINEMAIN_T_AFFINEMAIN_DOUBLE,
      "value"_a)

    .def(py::init<const Interval&>(),
      AFFINEMAIN_T_AFFINEMAIN_CONST_INTERVAL_REF,
      "itv"_a)

    .def(py::init<const Affine&>(),
      AFFINEMAIN_T_AFFINEMAIN_CONST_AFFINEMAIN_T_REF,
      "x"_a)

    .def_static("change_mode", &Affine::change_mode,
      STATIC_VOID_AFFINEMAIN_T_CHANGE_MODE_AFFINE_MODE,
      "tt"_a = Affine::AF_Lin_Chebyshev)

    .def_static("get_mode", &Affine::get_mode,
      STATIC_AFFINEMAIN_T_AFFINE_MODE_AFFINEMAIN_T_GET_MODE)

    .def(py::self == py::self,
      BOOL_AFFINEMAIN_T_OPERATOREQ_CONST_AFFINEMAIN_T_REF_CONST,
      "x"_a)

    .def(py::self == Interval(),
      BOOL_AFFINEMAIN_T_OPERATOREQ_CONST_INTERVAL_REF_CONST,
      "x"_a)

    .def(py::self == double(),
      BOOL_AFFINEMAIN_T_OPERATOREQ_CONST_DOUBLE_CONST,
      "x"_a)

    .def(py::self != py::self,
      BOOL_AFFINEMAIN_T_OPERATORNEQ_CONST_AFFINEMAIN_T_REF_CONST,
      "x"_a)

    .def(py::self != Interval(),
      BOOL_AFFINEMAIN_T_OPERATORNEQ_CONST_INTERVAL_REF_CONST,
      "x"_a)

    .def(py::self != double(),
      BOOL_AFFINEMAIN_T_OPERATORNEQ_CONST_DOUBLE_CONST,
      "x"_a)

    .def(py::self < py::self,
      BOOLINTERVAL_AFFINEMAIN_T_OPERATOR__CONST_AFFINEMAIN_T_REF_CONST,
      "x"_a)

    .def(py::self > py::self,
      BOOLINTERVAL_AFFINEMAIN_T_OPERATOR_CONST_AFFINEMAIN_T_REF_CONST,
      "x"_a)

    .def(py::self < Interval(),
      BOOLINTERVAL_AFFINEMAIN_T_OPERATOR__CONST_INTERVAL_REF_CONST,
      "x"_a)

    .def(py::self > Interval(),
      BOOLINTERVAL_AFFINEMAIN_T_OPERATOR_CONST_INTERVAL_REF_CONST,
      "x"_a)

    .def("lb", &Affine::lb,
      DOUBLE_AFFINEMAIN_T_LB_CONST)

    .def("ub", &Affine::ub,
      DOUBLE_AFFINEMAIN_T_UB_CONST)

    .def("mid", &Affine::mid,
      DOUBLE_AFFINEMAIN_T_MID_CONST)

    .def("mag", &Affine::mag,
      DOUBLE_AFFINEMAIN_T_MAG_CONST)

    .def("mig", &Affine::mig,
      DOUBLE_AFFINEMAIN_T_MIG_CONST)

    .def("smag", &Affine::smag,
      DOUBLE_AFFINEMAIN_T_SMAG_CONST)

    .def("smig", &Affine::smig,
      DOUBLE_AFFINEMAIN_T_SMIG_CONST)

    .def("rad", &Affine::rad,
      DOUBLE_AFFINEMAIN_T_RAD_CONST)

    .def("diam", &Affine::diam,
      DOUBLE_AFFINEMAIN_T_DIAM_CONST)

    .def("volume", &Affine::volume,
      DOUBLE_AFFINEMAIN_T_VOLUME_CONST)

    .def("size", [](const Affine& x)
        {
          return matlab::output_index(x.size());
        },
      INDEX_AFFINEMAIN_T_SIZE_CONST)

    .def("itv", &Affine::itv,
      CONST_INTERVAL_AFFINEMAIN_T_ITV_CONST)

    .def("noise_count", [](const Affine& x)
        {
          return matlab::output_index(x.noise_count());
        },
      INDEX_AFFINEMAIN_T_NOISE_COUNT_CONST)

    .def("noise", [](const Affine& x, Index_type i)
        {
          matlab::test_integer(i);
          return x.noise(matlab::input_index(i));
        },
      DOUBLE_AFFINEMAIN_T_NOISE_INDEX_CONST,
      "i"_a)

    .def("err", &Affine::err,
      DOUBLE_AFFINEMAIN_T_ERR_CONST)

    .def("set_empty", &Affine::set_empty,
      VOID_AFFINEMAIN_T_SET_EMPTY)

    .def("init", (Affine&(Affine::*)())&Affine::init,
      AFFINEMAIN_T_REF_AFFINEMAIN_T_INIT)

    .def("init", (Affine&(Affine::*)(const Interval&))&Affine::init,
      AFFINEMAIN_T_REF_AFFINEMAIN_T_INIT_CONST_INTERVAL_REF,
      "x"_a)

    .def("init_from_list", &Affine::init_from_list,
      AFFINEMAIN_T_REF_AFFINEMAIN_T_INIT_FROM_LIST_CONST_LIST_DOUBLE_REF,
      "l"_a)

    .def("is_active", &Affine::is_active,
      BOOL_AFFINEMAIN_T_IS_ACTIVE_CONST)

    .def("is_empty", &Affine::is_empty,
      BOOL_AFFINEMAIN_T_IS_EMPTY_CONST)

    .def("is_degenerated", &Affine::is_degenerated,
      BOOL_AFFINEMAIN_T_IS_DEGENERATED_CONST)

    .def("is_unbounded", &Affine::is_unbounded,
      BOOL_AFFINEMAIN_T_IS_UNBOUNDED_CONST)

    .def("compact", (void(Affine::*)(double))&Affine::compact,
      VOID_AFFINEMAIN_T_COMPACT_DOUBLE,
      "tol"_a)

    .def("compact", (void(Affine::*)())&Affine::compact,
      VOID_AFFINEMAIN_T_COMPACT)

    .def("contains", &Affine::contains,
      BOOL_AFFINEMAIN_T_CONTAINS_CONST_DOUBLE_REF_CONST,
      "d"_a)

    .def("__contains__", &Affine::contains,
      BOOL_AFFINEMAIN_T_CONTAINS_CONST_DOUBLE_REF_CONST,
      "d"_a)

    .def("interior_contains", &Affine::interior_contains,
      BOOL_AFFINEMAIN_T_INTERIOR_CONTAINS_CONST_DOUBLE_REF_CONST,
      "d"_a)

    .def("intersects", (bool(Affine::*)(const Interval&)const)&Affine::intersects,
      BOOL_AFFINEMAIN_T_INTERSECTS_CONST_INTERVAL_REF_CONST,
      "x"_a)

    .def("intersects", (bool(Affine::*)(const Affine&)const)&Affine::intersects,
      BOOL_AFFINEMAIN_T_INTERSECTS_CONST_AFFINEMAIN_T_REF_CONST,
      "x"_a)

    .def("is_disjoint", (bool(Affine::*)(const Interval&)const)&Affine::is_disjoint,
      BOOL_AFFINEMAIN_T_IS_DISJOINT_CONST_INTERVAL_REF_CONST,
      "x"_a)

    .def("is_disjoint", (bool(Affine::*)(const Affine&)const)&Affine::is_disjoint,
      BOOL_AFFINEMAIN_T_IS_DISJOINT_CONST_AFFINEMAIN_T_REF_CONST,
      "x"_a)

    .def("overlaps", (bool(Affine::*)(const Interval&)const)&Affine::overlaps,
      BOOL_AFFINEMAIN_T_OVERLAPS_CONST_INTERVAL_REF_CONST,
      "x"_a)

    .def("overlaps", (bool(Affine::*)(const Affine&)const)&Affine::overlaps,
      BOOL_AFFINEMAIN_T_OVERLAPS_CONST_AFFINEMAIN_T_REF_CONST,
      "x"_a)

    .def("is_subset", (bool(Affine::*)(const Interval&)const)&Affine::is_subset,
      BOOL_AFFINEMAIN_T_IS_SUBSET_CONST_INTERVAL_REF_CONST,
      "x"_a)

    .def("is_subset", (bool(Affine::*)(const Affine&)const)&Affine::is_subset,
      BOOL_AFFINEMAIN_T_IS_SUBSET_CONST_AFFINEMAIN_T_REF_CONST,
      "x"_a)

    .def("is_strict_subset", (bool(Affine::*)(const Interval&)const)&Affine::is_strict_subset,
      BOOL_AFFINEMAIN_T_IS_STRICT_SUBSET_CONST_INTERVAL_REF_CONST,
      "x"_a)

    .def("is_strict_subset", (bool(Affine::*)(const Affine&)const)&Affine::is_strict_subset,
      BOOL_AFFINEMAIN_T_IS_STRICT_SUBSET_CONST_AFFINEMAIN_T_REF_CONST,
      "x"_a)

    .def("is_interior_subset", (bool(Affine::*)(const Interval&)const)&Affine::is_interior_subset,
      BOOL_AFFINEMAIN_T_IS_INTERIOR_SUBSET_CONST_INTERVAL_REF_CONST,
      "x"_a)

    .def("is_interior_subset", (bool(Affine::*)(const Affine&)const)&Affine::is_interior_subset,
      BOOL_AFFINEMAIN_T_IS_INTERIOR_SUBSET_CONST_AFFINEMAIN_T_REF_CONST,
      "x"_a)

    // Note: is_relative_interior_subset() is declared in codac2_AffineMain.h but
    // cannot be instantiated (codac2::Interval provides no such method), it is
    // therefore not exported here.

    .def("is_strict_interior_subset", (bool(Affine::*)(const Interval&)const)&Affine::is_strict_interior_subset,
      BOOL_AFFINEMAIN_T_IS_STRICT_INTERIOR_SUBSET_CONST_INTERVAL_REF_CONST,
      "x"_a)

    .def("is_strict_interior_subset", (bool(Affine::*)(const Affine&)const)&Affine::is_strict_interior_subset,
      BOOL_AFFINEMAIN_T_IS_STRICT_INTERIOR_SUBSET_CONST_AFFINEMAIN_T_REF_CONST,
      "x"_a)

    .def("is_superset", (bool(Affine::*)(const Interval&)const)&Affine::is_superset,
      BOOL_AFFINEMAIN_T_IS_SUPERSET_CONST_INTERVAL_REF_CONST,
      "x"_a)

    .def("is_superset", (bool(Affine::*)(const Affine&)const)&Affine::is_superset,
      BOOL_AFFINEMAIN_T_IS_SUPERSET_CONST_AFFINEMAIN_T_REF_CONST,
      "x"_a)

    .def("is_strict_superset", (bool(Affine::*)(const Interval&)const)&Affine::is_strict_superset,
      BOOL_AFFINEMAIN_T_IS_STRICT_SUPERSET_CONST_INTERVAL_REF_CONST,
      "x"_a)

    .def("is_strict_superset", (bool(Affine::*)(const Affine&)const)&Affine::is_strict_superset,
      BOOL_AFFINEMAIN_T_IS_STRICT_SUPERSET_CONST_AFFINEMAIN_T_REF_CONST,
      "x"_a)

    .def("inflate", &Affine::inflate,
      AFFINEMAIN_T_REF_AFFINEMAIN_T_INFLATE_DOUBLE,
      "radd"_a)

    // Note: operator&= and operator|= are explicitly deleted for affine forms
    // (see codac2_AffineMain.h): the intersection or the hull of two affine
    // forms cannot be represented without discarding affine information.
    // They are therefore not exported, unlike for the Interval class.

    .def(-py::self,
      AFFINEMAIN_T_AFFINEMAIN_T_OPERATORMINUS_CONST)

    .def(py::self += double(),
      AFFINEMAIN_T_REF_AFFINEMAIN_T_OPERATORPLUSEQ_DOUBLE,
      "d"_a)

    .def(py::self += Interval(),
      AFFINEMAIN_T_REF_AFFINEMAIN_T_OPERATORPLUSEQ_CONST_INTERVAL_REF,
      "x"_a)

    .def(py::self += py::self,
      AFFINEMAIN_T_REF_AFFINEMAIN_T_OPERATORPLUSEQ_CONST_AFFINEMAIN_T_REF,
      "x"_a)

    .def(py::self -= double(),
      AFFINEMAIN_T_REF_AFFINEMAIN_T_OPERATORMINUSEQ_DOUBLE,
      "d"_a)

    .def(py::self -= Interval(),
      AFFINEMAIN_T_REF_AFFINEMAIN_T_OPERATORMINUSEQ_CONST_INTERVAL_REF,
      "x"_a)

    .def(py::self -= py::self,
      AFFINEMAIN_T_REF_AFFINEMAIN_T_OPERATORMINUSEQ_CONST_AFFINEMAIN_T_REF,
      "x"_a)

    .def(py::self *= double(),
      AFFINEMAIN_T_REF_AFFINEMAIN_T_OPERATORMULEQ_DOUBLE,
      "d"_a)

    .def(py::self *= Interval(),
      AFFINEMAIN_T_REF_AFFINEMAIN_T_OPERATORMULEQ_CONST_INTERVAL_REF,
      "x"_a)

    .def(py::self *= py::self,
      AFFINEMAIN_T_REF_AFFINEMAIN_T_OPERATORMULEQ_CONST_AFFINEMAIN_T_REF,
      "x"_a)

    .def(py::self /= double(),
      AFFINEMAIN_T_REF_AFFINEMAIN_T_OPERATORDIVEQ_DOUBLE,
      "d"_a)

    .def(py::self /= Interval(),
      AFFINEMAIN_T_REF_AFFINEMAIN_T_OPERATORDIVEQ_CONST_INTERVAL_REF,
      "x"_a)

    .def(py::self /= py::self,
      AFFINEMAIN_T_REF_AFFINEMAIN_T_OPERATORDIVEQ_CONST_AFFINEMAIN_T_REF,
      "x"_a)

    .def("__pow__", (Affine(*)(const Affine&,int)) &codac2::pow,
      AFFINEMAIN_T_POW_CONST_AFFINEMAIN_T_REF_INT,
      "n"_a)

    .def("__pow__", (Affine(*)(const Affine&,double)) &codac2::pow,
      AFFINEMAIN_T_POW_CONST_AFFINEMAIN_T_REF_DOUBLE,
      "n"_a)

    .def("__pow__", (Affine(*)(const Affine&,const Interval&)) &codac2::pow,
      AFFINEMAIN_T_POW_CONST_AFFINEMAIN_T_REF_CONST_INTERVAL_REF,
      "n"_a)

    .def("__pow__", (Affine(*)(const Affine&,const Affine&)) &codac2::pow,
      AFFINEMAIN_T_POW_CONST_AFFINEMAIN_T_REF_CONST_AFFINEMAIN_T_REF,
      "n"_a)

    .def_static("empty", &Affine::empty,
      STATIC_AFFINEMAIN_T_AFFINEMAIN_T_EMPTY)

    .def("__repr__", [](const Affine& x) {
          std::ostringstream stream;
          stream << x;
          return string(stream.str());
        },
      OSTREAM_REF_OPERATOROUT_OSTREAM_REF_CONST_AFFINEMAIN_T_REF)
  ;

  // Automatic cast of reals and intervals into affine forms
  py::implicitly_convertible<double,Affine>();
  py::implicitly_convertible<int,Affine>();
  py::implicitly_convertible<Interval,Affine>();

  return exported_affine_class;
}
