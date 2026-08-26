/**
 *  \file codac2_py_AffineMatrixBase.h
 *
 *  AffineMatrixBase binding
 *
 * ----------------------------------------------------------------------------
 *  \date       2026
 *  \author     Jordan Ninin
 *  \copyright  Copyright 2026 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include <sstream>
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include <codac2_Affine.h>
#include <codac2_Interval.h>
#include <codac2_IntervalVector.h>
#include <codac2_IntervalMatrix.h>

#include "codac2_py_doc.h"
#include "codac2_py_matlab.h"

using namespace std;
using namespace codac2;
namespace py = pybind11;
using namespace pybind11::literals;

/*
 * Deliberately does not call export_MatrixBase(): the generic matrix binding
 * exports min_coeff(), max_coeff(), norm() and squared_norm(), which assume a
 * totally ordered scalar type. The comparison operators of AffineMain return
 * a BoolInterval, so these functions are not available for affine domains.
 * The affine-compatible part of the generic binding is reproduced here.
 */
template<typename S,typename V,bool VECTOR_INHERITANCE>
void export_AffineMatrixBase([[maybe_unused]] py::module& m, py::class_<S>& pyclass)
{
  // Interval counterpart of S (same shape, Interval scalar)
  using IS = Eigen::Matrix<Interval,(int)S::RowsAtCompileTime,(int)S::ColsAtCompileTime>;

  pyclass

    .def("__eq__", [](const S& x1, const S& x2)
        {
          return x1 == x2;
        },
      DOC_TO_BE_DEFINED)

    .def("__ne__", [](const S& x1, const S& x2)
        {
          return !(x1 == x2);
        },
      DOC_TO_BE_DEFINED)

    .def("__len__", [](const S& x)
        {
          return x.size();
        },
      DOC_TO_BE_DEFINED)

    .def("size", [](const S& x)
        {
          return x.size();
        },
      DOC_TO_BE_DEFINED)

    .def("rows", [](const S& x)
        {
          return x.rows();
        },
      DOC_TO_BE_DEFINED)

    .def("cols", [](const S& x)
        {
          return x.cols();
        },
      DOC_TO_BE_DEFINED)

    .def("init", [](S& x, const Affine& a)
        {
          x.init(a);
        },
      MATRIX_ADDONS_BASE_AUTO_REF_INIT_CONST_SCALAR_REF,
      "x"_a)

    .def("transpose", [](const S& x) -> Eigen::Matrix<Affine,-1,-1>
        {
          return x.transpose().eval();
        },
      DOC_TO_BE_DEFINED)

    .def("itv", [](const S& x)
        {
          return x.itv();
        },
      MATRIXBASE_ADDONS_AFFINEMAINVECTOR_AUTO_ITV_CONST)
  ;

  if constexpr(!VECTOR_INHERITANCE)
  {
    pyclass

      .def("is_squared", [](const S& x)
          {
            return x.is_squared();
          },
        MATRIXBASE_ADDONS_BASE_BOOL_IS_SQUARED_CONST)

      .def(
          #if FOR_MATLAB
            "__call__"
          #else
            "__getitem__"
          #endif
          , [](const S& x, const py::tuple& ij) -> const Affine&
          {
            if constexpr(FOR_MATLAB)
              assert_release(py::isinstance<py::int_>(ij[0]) && py::isinstance<py::int_>(ij[1]));

            int i = ij[0].cast<int>();
            int j = ij[1].cast<int>();

            return x(matlab::input_index(i), matlab::input_index(j));
          }, py::return_value_policy::reference_internal,
        MATRIX_ADDONS_BASE_CONST_SCALAR_REF_OPERATORCALL_INDEX_INDEX_CONST)

      .def(
          #if FOR_MATLAB
            "set_item"
          #else
            "__setitem__"
          #endif
          , [](S& x, const std::vector<Index_type>& ij, const Affine& a)
          {
            auto ij_conv = matlab::convert_indices(ij);
            x(ij_conv[0],ij_conv[1]) = a;
          },
        MATRIX_ADDONS_BASE_SCALAR_REF_OPERATORCALL_INDEX_INDEX)

      .def("__call__", [](S& x, Index_type i, Index_type j) -> Affine&
          {
            matlab::test_integer(i,j);
            return x(matlab::input_index(i),matlab::input_index(j));
          }, py::return_value_policy::reference_internal,
        MATRIX_ADDONS_BASE_SCALAR_REF_OPERATORCALL_INDEX_INDEX)

      .def("block", [](S& x, Index_type i, Index_type j, Index_type p, Index_type q) -> Eigen::Matrix<Affine,-1,-1>
          {
            matlab::test_integer(i,j);
            matlab::test_integer(p,q);
            return x.block(matlab::input_index(i),matlab::input_index(j),matlab::input_index(p),matlab::input_index(q));
          },
        py::keep_alive<0,1>(),
        DOC_TO_BE_DEFINED)

      .def("col", [](S& x, Index_type i) -> Eigen::Matrix<Affine,-1,1>
          {
            matlab::test_integer(i);
            return x.col(matlab::input_index(i)).eval();
          },
        DOC_TO_BE_DEFINED)

      .def("row", [](S& x, Index_type i) -> Eigen::Matrix<Affine,1,-1>
          {
            matlab::test_integer(i);
            return x.row(matlab::input_index(i)).eval();
          },
        DOC_TO_BE_DEFINED)

      .def("set_block", [](S& x, Index_type i, Index_type j, Index_type p, Index_type q, const S& y)
          {
            matlab::test_integer(i,j);
            matlab::test_integer(p,q);
            x.block(matlab::input_index(i),matlab::input_index(j),matlab::input_index(p),matlab::input_index(q)) = y;
          },
        py::keep_alive<0,1>(),
        DOC_TO_BE_DEFINED)

      .def("set_col", [](S& x, Index_type i, const Eigen::Matrix<Affine,-1,1>& y)
          {
            matlab::test_integer(i);
            x.col(matlab::input_index(i)) = y;
          },
        DOC_TO_BE_DEFINED)

      .def("set_row", [](S& x, Index_type i, const Eigen::Matrix<Affine,1,-1>& y)
          {
            matlab::test_integer(i);
            x.row(matlab::input_index(i)) = y;
          },
        DOC_TO_BE_DEFINED)

      .def("resize", [](S& x, Index_type nb_rows, Index_type nb_cols)
          {
            matlab::test_integer(nb_rows, nb_cols);
            x.resize(nb_rows, nb_cols);
          },
        DOC_TO_BE_DEFINED,
        "nb_rows"_a, "nb_cols"_a)

      .def("conservativeResize", [](S& x, Index_type nb_rows, Index_type nb_cols)
          {
            matlab::test_integer(nb_rows, nb_cols);
            x.conservativeResize(nb_rows, nb_cols);
          },
        "Resize the matrix while preserving the existing coefficients.",
        "nb_rows"_a, "nb_cols"_a)

      .def_static("zero", [](Index_type r, Index_type c)
          {
            matlab::test_integer(r,c);
            return S::zero(r,c);
          },
        MATRIX_ADDONS_MATRIXBASE_STATIC_MATRIX_SCALARRC_ZERO_INDEX_INDEX,
        "r"_a, "c"_a)

      .def_static("ones", [](Index_type r, Index_type c)
          {
            matlab::test_integer(r,c);
            return S::ones(r,c);
          },
        MATRIX_ADDONS_MATRIXBASE_STATIC_MATRIX_SCALARRC_ONES_INDEX_INDEX,
        "r"_a, "c"_a)

      .def_static("constant", [](Index_type r, Index_type c, const Affine& x)
          {
            matlab::test_integer(r,c);
            return S::constant(r,c,x);
          },
        MATRIX_ADDONS_MATRIXBASE_STATIC_MATRIX_SCALARRC_CONSTANT_INDEX_INDEX_CONST_SCALAR_REF,
        "r"_a, "c"_a, "x"_a)

      .def_static("eye", [](Index_type r, Index_type c)
          {
            matlab::test_integer(r,c);
            return S::eye(r,c);
          },
        MATRIX_ADDONS_MATRIXBASE_STATIC_MATRIX_SCALARRC_EYE_INDEX_INDEX,
        "r"_a, "c"_a)
    ;
  }

  // Methods inherited from the interval-domain part of Eigen::MatrixBase
  // (the affine scalar type satisfies the IsIntervalDomain concept)

  pyclass

    .def("volume", [](const S& x) { return x.volume(); },
      MATRIXBASE_ADDONS_INTERVALMATRIXBASE_DOUBLE_VOLUME_CONST)

    .def("is_empty", [](const S& x) { return x.is_empty(); },
      MATRIXBASE_ADDONS_INTERVALMATRIXBASE_BOOL_IS_EMPTY_CONST)

    .def("set_empty", [](S& x) { x.set_empty(); },
      MATRIX_ADDONS_INTERVALMATRIXBASE_VOID_SET_EMPTY)

    .def("lb", [](const S& x) { return x.lb(); },
      MATRIXBASE_ADDONS_INTERVALMATRIXBASE_AUTO_LB_CONST)

    .def("ub", [](const S& x) { return x.ub(); },
      MATRIXBASE_ADDONS_INTERVALMATRIXBASE_AUTO_UB_CONST)

    .def("mid", [](const S& x) { return x.mid(); },
      MATRIXBASE_ADDONS_INTERVALMATRIXBASE_AUTO_MID_CONST)

    .def("mag", [](const S& x) { return x.mag(); },
      MATRIXBASE_ADDONS_INTERVALMATRIXBASE_AUTO_MAG_CONST)

    .def("mig", [](const S& x) { return x.mig(); },
      MATRIXBASE_ADDONS_INTERVALMATRIXBASE_AUTO_MIG_CONST)

    .def("smag", [](const S& x) { return x.smag(); },
      MATRIXBASE_ADDONS_INTERVALMATRIXBASE_AUTO_SMAG_CONST)

    .def("smig", [](const S& x) { return x.smig(); },
      MATRIXBASE_ADDONS_INTERVALMATRIXBASE_AUTO_SMIG_CONST)

    .def("rad", [](const S& x) { return x.rad(); },
      MATRIXBASE_ADDONS_INTERVALMATRIXBASE_AUTO_RAD_CONST)

    .def("diam", [](const S& x) { return x.diam(); },
      MATRIXBASE_ADDONS_INTERVALMATRIXBASE_AUTO_DIAM_CONST)

    .def("min_rad", [](const S& x, const std::vector<Index_type>& among_indices)
        {
          return x.min_rad(matlab::convert_indices(among_indices));
        },
      MATRIXBASE_ADDONS_INTERVALMATRIXBASE_DOUBLE_MIN_RAD_CONST_VECTOR_INDEX_REF_CONST,
      "among_indices"_a=std::vector<Index_type>())

    .def("max_rad", [](const S& x, const std::vector<Index_type>& among_indices)
        {
          return x.max_rad(matlab::convert_indices(among_indices));
        },
      MATRIXBASE_ADDONS_INTERVALMATRIXBASE_DOUBLE_MAX_RAD_CONST_VECTOR_INDEX_REF_CONST,
      "among_indices"_a=std::vector<Index_type>())

    .def("min_diam", [](const S& x, const std::vector<Index_type>& among_indices)
        {
          return x.min_diam(matlab::convert_indices(among_indices));
        },
      MATRIXBASE_ADDONS_INTERVALMATRIXBASE_DOUBLE_MIN_DIAM_CONST_VECTOR_INDEX_REF_CONST,
      "among_indices"_a=std::vector<Index_type>())

    .def("max_diam", [](const S& x, const std::vector<Index_type>& among_indices)
        {
          return x.max_diam(matlab::convert_indices(among_indices));
        },
      MATRIXBASE_ADDONS_INTERVALMATRIXBASE_DOUBLE_MAX_DIAM_CONST_VECTOR_INDEX_REF_CONST,
      "among_indices"_a=std::vector<Index_type>())

    .def("min_diam_index", [](const S& x, const std::vector<Index_type>& among_indices)
        {
          return matlab::output_index(x.min_diam_index(matlab::convert_indices(among_indices)));
        },
      MATRIXBASE_ADDONS_INTERVALMATRIXBASE_INDEX_MIN_DIAM_INDEX_CONST_VECTOR_INDEX_REF_CONST,
      "among_indices"_a=std::vector<Index_type>())

    .def("max_diam_index", [](const S& x, const std::vector<Index_type>& among_indices)
        {
          return matlab::output_index(x.max_diam_index(matlab::convert_indices(among_indices)));
        },
      MATRIXBASE_ADDONS_INTERVALMATRIXBASE_INDEX_MAX_DIAM_INDEX_CONST_VECTOR_INDEX_REF_CONST,
      "among_indices"_a=std::vector<Index_type>())

    .def("extr_diam_index", [](const S& x, bool min, const std::vector<Index_type>& among_indices)
        {
          return matlab::output_index(x.extr_diam_index(min,matlab::convert_indices(among_indices)));
        },
      MATRIXBASE_ADDONS_INTERVALMATRIXBASE_INDEX_EXTR_DIAM_INDEX_BOOL_CONST_VECTOR_INDEX_REF_CONST,
      "min"_a, "among_indices"_a=std::vector<Index_type>())

    .def("__contains__", [](const S& x, const V& y) { return x.contains(y); },
      MATRIXBASE_ADDONS_INTERVALMATRIXBASE_BOOL_CONTAINS_CONST_MATRIXBASE_OTHERDERIVED_REF_CONST)

    .def("contains", [](const S& x, const V& y) { return x.contains(y); },
      MATRIXBASE_ADDONS_INTERVALMATRIXBASE_BOOL_CONTAINS_CONST_MATRIXBASE_OTHERDERIVED_REF_CONST)

    .def("interior_contains", [](const S& x, const V& y) { return x.interior_contains(y); },
      MATRIXBASE_ADDONS_INTERVALMATRIXBASE_BOOL_INTERIOR_CONTAINS_CONST_MATRIXBASE_OTHERDERIVED_REF_CONST)

    .def("is_unbounded", [](const S& x) { return x.is_unbounded(); },
      MATRIXBASE_ADDONS_INTERVALMATRIXBASE_BOOL_IS_UNBOUNDED_CONST)

    .def("is_degenerated", [](const S& x) { return x.is_degenerated(); },
      MATRIXBASE_ADDONS_INTERVALMATRIXBASE_BOOL_IS_DEGENERATED_CONST)

    .def("is_flat", [](const S& x) { return x.is_flat(); },
      MATRIXBASE_ADDONS_INTERVALMATRIXBASE_BOOL_IS_FLAT_CONST)

    .def("intersects", [](const S& x, const IS& y) { return x.intersects(y); },
      MATRIXBASE_ADDONS_INTERVALMATRIXBASE_BOOL_INTERSECTS_CONST_MATRIXBASE_OTHERDERIVED_REF_CONST)

    .def("is_disjoint", [](const S& x, const IS& y) { return x.is_disjoint(y); },
      MATRIXBASE_ADDONS_INTERVALMATRIXBASE_BOOL_IS_DISJOINT_CONST_MATRIXBASE_OTHERDERIVED_REF_CONST)

    .def("overlaps", [](const S& x, const IS& y) { return x.overlaps(y); },
      MATRIXBASE_ADDONS_INTERVALMATRIXBASE_BOOL_OVERLAPS_CONST_MATRIXBASE_OTHERDERIVED_REF_CONST)

    .def("is_subset", [](const S& x, const IS& y) { return x.is_subset(y); },
      MATRIXBASE_ADDONS_INTERVALMATRIXBASE_BOOL_IS_SUBSET_CONST_MATRIXBASE_OTHERDERIVED_REF_CONST)

    .def("is_strict_subset", [](const S& x, const IS& y) { return x.is_strict_subset(y); },
      MATRIXBASE_ADDONS_INTERVALMATRIXBASE_BOOL_IS_STRICT_SUBSET_CONST_MATRIXBASE_OTHERDERIVED_REF_CONST)

    .def("is_interior_subset", [](const S& x, const IS& y) { return x.is_interior_subset(y); },
      MATRIXBASE_ADDONS_INTERVALMATRIXBASE_BOOL_IS_INTERIOR_SUBSET_CONST_MATRIXBASE_OTHERDERIVED_REF_CONST)

    .def("is_strict_interior_subset", [](const S& x, const IS& y) { return x.is_strict_interior_subset(y); },
      MATRIXBASE_ADDONS_INTERVALMATRIXBASE_BOOL_IS_STRICT_INTERIOR_SUBSET_CONST_MATRIXBASE_OTHERDERIVED_REF_CONST)

    .def("is_superset", [](const S& x, const IS& y) { return x.is_superset(y); },
      MATRIXBASE_ADDONS_INTERVALMATRIXBASE_BOOL_IS_SUPERSET_CONST_MATRIXBASE_OTHERDERIVED_REF_CONST)

    .def("is_strict_superset", [](const S& x, const IS& y) { return x.is_strict_superset(y); },
      MATRIXBASE_ADDONS_INTERVALMATRIXBASE_BOOL_IS_STRICT_SUPERSET_CONST_MATRIXBASE_OTHERDERIVED_REF_CONST)

    .def("inflate", [](S& x, double r) { return x.inflate(r); },
      MATRIX_ADDONS_INTERVALMATRIXBASE_AUTO_REF_INFLATE_DOUBLE,
      "r"_a)

    .def("inflate", [](S& x, const V& r) { return x.inflate(r); },
      MATRIX_ADDONS_INTERVALMATRIXBASE_AUTO_REF_INFLATE_CONST_MATRIXBASE_OTHERDERIVED_REF,
      "r"_a)

  // Note: bisect(), bisect_largest(), is_bisectable(), has_integer_bounds(),
  // rand() and the set operators &=, |=, & and | are not exported: they are
  // not provided by the affine scalar type (see codac2_AffineMain.h).

  // Arithmetic operators

    .def("__add__", [](const S& x1, const S& x2) -> S { return (x1+x2).eval(); },
      py::is_operator())

    .def("__add__", [](const S& x1, const V& x2) -> S { return (x1+x2).eval(); },
      py::is_operator())

    .def("__add__", [](const S& x1, const IS& x2) -> S { return (x1+x2).eval(); },
      py::is_operator())

    .def("__radd__", [](const S& x1, const V& x2) -> S { return (x2+x1).eval(); },
      py::is_operator())

    .def("__radd__", [](const S& x1, const IS& x2) -> S { return (x2+x1).eval(); },
      py::is_operator())

    .def("__sub__", [](const S& x1, const S& x2) -> S { return (x1-x2).eval(); },
      py::is_operator())

    .def("__sub__", [](const S& x1, const V& x2) -> S { return (x1-x2).eval(); },
      py::is_operator())

    .def("__sub__", [](const S& x1, const IS& x2) -> S { return (x1-x2).eval(); },
      py::is_operator())

    .def("__rsub__", [](const S& x1, const V& x2) -> S { return (x2-x1).eval(); },
      py::is_operator())

    .def("__rsub__", [](const S& x1, const IS& x2) -> S { return (x2-x1).eval(); },
      py::is_operator())

    .def("__neg__", [](const S& x) -> S { return (-x).eval(); },
      py::is_operator())

    .def("__mul__", [](const S& x1, double x2) -> S { return (x1*x2).eval(); },
      py::is_operator())

    .def("__mul__", [](const S& x1, const Interval& x2) -> S { return (x1*x2).eval(); },
      py::is_operator())

    .def("__mul__", [](const S& x1, const Affine& x2) -> S { return (x1*x2).eval(); },
      py::is_operator())

    .def("__rmul__", [](const S& x1, double x2) -> S { return (x2*x1).eval(); },
      py::is_operator())

    .def("__rmul__", [](const S& x1, const Interval& x2) -> S { return (x2*x1).eval(); },
      py::is_operator())

    .def("__rmul__", [](const S& x1, const Affine& x2) -> S { return (x2*x1).eval(); },
      py::is_operator())

    .def("__truediv__", [](const S& x1, double x2) -> S { return (x1/x2).eval(); },
      py::is_operator())

    .def("__truediv__", [](const S& x1, const Interval& x2) -> S { return (x1/x2).eval(); },
      py::is_operator())

    .def("__truediv__", [](const S& x1, const Affine& x2) -> S { return (x1/x2).eval(); },
      py::is_operator())

    .def("__iadd__", [](S& x1, const S& x2) -> S& { x1 += x2; return x1; },
      py::is_operator())

    .def("__isub__", [](S& x1, const S& x2) -> S& { x1 -= x2; return x1; },
      py::is_operator())

    .def("__imul__", [](S& x1, double x2) -> S& { x1 *= x2; return x1; },
      py::is_operator())

    .def("__imul__", [](S& x1, const Interval& x2) -> S& { x1 *= x2; return x1; },
      py::is_operator())

    .def("__itruediv__", [](S& x1, double x2) -> S& { x1 /= x2; return x1; },
      py::is_operator())

    .def("__itruediv__", [](S& x1, const Interval& x2) -> S& { x1 /= x2; return x1; },
      py::is_operator())
  ;

  py::implicitly_convertible<V,S>();
}
