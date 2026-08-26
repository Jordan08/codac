/**
 *  \file codac2_py_AffineVector_templ.h
 *
 *  AffineVector/AffineRow binding (shared part)
 *
 * ----------------------------------------------------------------------------
 *  \date       2026
 *  \author     Jordan Ninin
 *  \copyright  Copyright 2026 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include <sstream>
#include <stdexcept>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include <codac2_Affine.h>

#include "codac2_py_AffineMatrixBase.h"

using namespace std;
using namespace codac2;
namespace py = pybind11;
using namespace pybind11::literals;

template<typename AV,typename V,typename IV>
void export_AffineVector_(py::module& m, py::class_<AV>& pyclass)
{
  export_AffineMatrixBase<AV,V,true>(m, pyclass);

  pyclass

    .def(py::init<>())

    .def(py::init(
        [](Index_type n)
        {
          matlab::test_integer(n);
          return std::make_unique<AV>((int)n);
        }),
      DOC_TO_BE_DEFINED,
      "n"_a)

    .def(py::init<const AV&>(),
      "x"_a)

    .def(py::init(
        [](const V& x)
        {
          auto v = std::make_unique<AV>((int)x.size());
          for(Index i = 0 ; i < x.size() ; i++)
            (*v)[i] = x[i];
          return v;
        }),
      DOC_TO_BE_DEFINED,
      "x"_a)

    .def(py::init(
        [](const IV& x)
        {
          auto v = std::make_unique<AV>((int)x.size());
          for(Index i = 0 ; i < x.size() ; i++)
            (*v)[i] = x[i];
          return v;
        }),
      DOC_TO_BE_DEFINED,
      "x"_a)

    .def(py::init(
        [](const std::vector<std::vector<double>>& v)
        {
          auto av = std::make_unique<AV>((int)v.size());
          for(size_t i = 0 ; i < v.size() ; i++)
          {
            if(v[i].size() == 1)
              (*av)[i] = Affine(Interval(v[i][0]));
            else if(v[i].size() == 2)
              (*av)[i] = Affine(Interval(v[i][0],v[i][1]));
            else
              throw invalid_argument("Affine form is not made of one or two values.");
          }
          return av;
        }),
      DOC_TO_BE_DEFINED,
      "v"_a)

    .def(py::init( // this constructor must be the last one to be declared
        [](const std::vector<Affine>& v)
        {
          auto av = std::make_unique<AV>((int)v.size());
          for(size_t i = 0 ; i < v.size() ; i++)
            (*av)[i] = v[i];
          return av;
        }),
      DOC_TO_BE_DEFINED,
      "v"_a)

    .def(
        #if FOR_MATLAB
          "__call__"
        #else
          "__getitem__"
        #endif
        ,
        [](const AV& x, Index_type i) -> const Affine&
        {
          matlab::test_integer(i);
          return x[matlab::input_index(i)];
        }, py::return_value_policy::reference_internal,
      MATRIX_ADDONS_VECTORBASE_CONST_SCALAR_REF_OPERATORCOMPO_INDEX_CONST)

    .def("get_item_0",
        [](const AV& x, Index_type i) -> const Affine&
        {
          matlab::test_integer(i);
          return x[i];
        }, py::return_value_policy::reference_internal,
      MATRIX_ADDONS_VECTORBASE_CONST_SCALAR_REF_OPERATORCOMPO_INDEX_CONST)

    .def(
        #if FOR_MATLAB
          "set_item"
        #else
          "__setitem__"
        #endif
        , [](AV& x, Index_type i, const Affine& a)
        {
          matlab::test_integer(i);
          x[matlab::input_index(i)] = a;
        },
      MATRIX_ADDONS_VECTORBASE_SCALAR_REF_OPERATORCOMPO_INDEX)

    .def("set_item_0", [](AV& x, Index_type i, const Affine& a)
        {
          matlab::test_integer(i);
          x[i] = a;
        },
      MATRIX_ADDONS_VECTORBASE_SCALAR_REF_OPERATORCOMPO_INDEX)

    .def("subvector", [](const AV& x, Index_type start_id, Index_type end_id) -> AV
        {
          matlab::test_integer(start_id, end_id);
          return x.subvector(matlab::input_index(start_id), matlab::input_index(end_id));
        },
      MATRIXBASE_ADDONS_VECTORBASE_AUTO_SUBVECTOR_INDEX_INDEX_CONST,
      "start_id"_a, "end_id"_a)

    .def("subvector_0", [](const AV& x, Index_type start_id, Index_type end_id) -> AV
        {
          matlab::test_integer(start_id, end_id);
          return x.subvector(start_id, end_id);
        },
      MATRIXBASE_ADDONS_VECTORBASE_AUTO_SUBVECTOR_INDEX_INDEX_CONST,
      "start_id"_a, "end_id"_a)

    .def("diag_matrix", [](const AV& x)
        {
          return x.diag_matrix().eval();
        },
      DOC_TO_BE_DEFINED)

    .def("resize", [](AV& x, Index_type n)
        {
          matlab::test_integer(n);
          x.resize(n);
        },
      DOC_TO_BE_DEFINED,
      "n"_a)

    .def("conservativeResize", [](AV& x, Index_type n)
        {
          matlab::test_integer(n);
          x.conservativeResize(n);
        },
        "Resize the vector while preserving the existing coefficients.",
      "n"_a)

    .def("put", [](AV& x, Index_type start_id, const AV& x1)
        {
          matlab::test_integer(start_id);
          x.put(matlab::input_index(start_id), x1);
        },
      MATRIX_ADDONS_VECTORBASE_VOID_PUT_INDEX_CONST_MATRIXBASE_OTHERDERIVED_REF,
      "start_id"_a, "x"_a)

    .def_static("empty", [](Index_type n)
        {
          matlab::test_integer(n);
          return AV::constant((int)n, Affine::empty());
        },
      MATRIX_ADDONS_INTERVALVECTOR_STATIC_AUTO_EMPTY_INDEX,
      "n"_a)

    .def_static("zero", [](Index_type n)
        {
          matlab::test_integer(n);
          return AV::zero(n);
        },
      MATRIX_ADDONS_VECTORBASE_STATIC_MATRIX_SCALARRC_ZERO_INDEX,
      "n"_a)

    .def_static("ones", [](Index_type n)
        {
          matlab::test_integer(n);
          return AV::ones(n);
        },
      MATRIX_ADDONS_VECTORBASE_STATIC_MATRIX_SCALARRC_ONES_INDEX,
      "n"_a)

    .def_static("constant", [](Index_type n, const Affine& x)
        {
          matlab::test_integer(n);
          return AV::constant(n,x);
        },
      MATRIX_ADDONS_VECTORBASE_STATIC_MATRIX_SCALARRC_CONSTANT_INDEX_CONST_SCALAR_REF,
      "n"_a, "x"_a)

    .def("__iter__", [](const AV& x)
        {
          return py::make_iterator(x.begin(), x.end());
        },
      py::keep_alive<0,1>() /* keep object alive while iterator exists */)
  ;

  py::implicitly_convertible<py::list,AV>();
  py::implicitly_convertible<V,AV>();
  py::implicitly_convertible<IV,AV>();
}
