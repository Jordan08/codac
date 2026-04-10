/** 
 *  Codac binding (core)
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <sstream>
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <codac2_Interval.h>
#include <codac2_SlicedTube.h>
#include <codac2_SlicedTube_operations.h>
#include <codac2_TubeBase.h>
#include <codac2_SampledTraj.h>
#include "codac2_py_matlab.h"
#include "codac2_py_SlicedTubeBase_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py):
#include "codac2_py_SlicedTube_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py):
#include "codac2_py_SlicedTube_operations_docs.h" // Generated file from Doxygen XML (doxygen2docstring.py):

using namespace std;
using namespace codac2;
namespace py = pybind11;
using namespace pybind11::literals;


template<typename T>
py::class_<SlicedTube<T>,TubeBase> _export_SlicedTube(py::module& m, const std::string& name)
{
  py::class_<SlicedTube<T>,TubeBase> exported_slicedtube_class(m, name.c_str(), SLICEDTUBEBASE_MAIN);
  exported_slicedtube_class

    .def(py::init<const std::shared_ptr<TDomain>&,const T&>(),
      SLICEDTUBE_T_SLICEDTUBE_CONST_SHARED_PTR_TDOMAIN_REF_CONST_T_REF,
      "tdomain"_a, "codomain"_a)

    .def(py::init<const std::shared_ptr<TDomain>&,const AnalyticFunction<typename ExprType<T>::Type>&>(),
      SLICEDTUBE_T_SLICEDTUBE_CONST_SHARED_PTR_TDOMAIN_REF_CONST_ANALYTICFUNCTION_TYPENAME_EXPRTYPE_T_TYPE_REF,
      "tdomain"_a, "f"_a)

    .def(py::init<const std::shared_ptr<TDomain>&,const SampledTraj<typename ExprType<T>::Type::Scalar>&>(),
      SLICEDTUBE_T_SLICEDTUBE_CONST_SHARED_PTR_TDOMAIN_REF_CONST_SAMPLEDTRAJ_V_REF,
      "tdomain"_a, "x"_a)

    .def(py::init<const SlicedTube<T>&>(),
      SLICEDTUBE_T_SLICEDTUBE_CONST_SLICEDTUBE_T_REF,
      "x"_a)

    .def("__iter__", [](const SlicedTube<T>& x)
        {
          return py::make_iterator(x.begin(), x.end());
        },
        py::keep_alive<0,1>()) // essential: keep object alive while iterator exists

    // From SlicedTubeBase class:
    .def("nb_slices", &SlicedTube<T>::nb_slices,
      SIZET_SLICEDTUBEBASE_NB_SLICES_CONST)

    .def("size", &SlicedTube<T>::size,
      INDEX_SLICEDTUBE_T_SIZE_CONST)
    
    .def("volume", &SlicedTube<T>::volume,
      DOUBLE_SLICEDTUBE_T_VOLUME_CONST)
    
    .def("first_slice", [](SlicedTube<T>& x) -> Slice<T>&
        {
          return *x.first_slice();
        },
      py::return_value_policy::reference,
      SHARED_PTR_SLICE_T_SLICEDTUBE_T_FIRST_SLICE)
    
    .def("last_slice", [](SlicedTube<T>& x) -> Slice<T>&
        {
          return *x.last_slice();
        },
      py::return_value_policy::reference,
      SHARED_PTR_SLICE_T_SLICEDTUBE_T_LAST_SLICE)
    
    .def("slice", [](SlicedTube<T>& x, std::shared_ptr<TSlice> it) -> Slice<T>&
        {
          return *x.slice(it);
        },
      py::return_value_policy::reference,
      SHARED_PTR_SLICE_T_SLICEDTUBE_T_SLICE_SHARED_PTR_TSLICE,
      "it"_a)
    
    .def("is_empty", &SlicedTube<T>::is_empty,
      BOOL_SLICEDTUBE_T_IS_EMPTY_CONST)
    
    .def("is_unbounded", &SlicedTube<T>::is_unbounded,
      BOOL_SLICEDTUBE_T_IS_UNBOUNDED_CONST)
    
    .def("codomain", &SlicedTube<T>::codomain,
      T_SLICEDTUBE_T_CODOMAIN_CONST)
    
    .def("__call__", [](const SlicedTube<T>& x, const Interval& t)
        {
          return x(t);
        },
      T_SLICEDTUBE_T_OPERATORCALL_CONST_INTERVAL_REF_CONST,
      "t"_a)

    .def("enclosed_bounds", &SlicedTube<T>::enclosed_bounds,
      PAIR_TT_SLICEDTUBE_T_ENCLOSED_BOUNDS_CONST_INTERVAL_REF_CONST,
      "t"_a)

    .def("invert", (Interval (SlicedTube<T>::*)(const T&,const Interval&) const) &SlicedTube<T>::invert,
      INTERVAL_SLICEDTUBE_T_INVERT_CONST_T_REF_CONST_INTERVAL_REF_CONST,
      "y"_a, "t"_a=Interval())

    .def("invert", [](const SlicedTube<T>& x, const T& y, py::list& v_t, const Interval& t)
        {
          vector<Interval> vector_t;
          x.invert(y, vector_t, t);
          v_t.clear();
          for(const auto& ti : vector_t)
            v_t.append(ti);
        },
      INTERVAL_SLICEDTUBE_T_INVERT_CONST_T_REF_CONST_SLICEDTUBE_T_REF_CONST_INTERVAL_REF_CONST,
      "y"_a, "v_t"_a, "t"_a=Interval())
    
    .def("set", (void (SlicedTube<T>::*)(const T&)) &SlicedTube<T>::set,
      VOID_SLICEDTUBE_T_SET_CONST_T_REF,
      "codomain"_a)
    
    .def("set", (void (SlicedTube<T>::*)(const T&,double)) &SlicedTube<T>::set,
      VOID_SLICEDTUBE_T_SET_CONST_T_REF_DOUBLE,
      "codomain"_a, "t"_a)
    
    .def("set", (void (SlicedTube<T>::*)(const T&,const Interval&)) &SlicedTube<T>::set,
      VOID_SLICEDTUBE_T_SET_CONST_T_REF_CONST_INTERVAL_REF,
      "codomain"_a, "t"_a)
    
    .def("set_ith_slice", &SlicedTube<T>::set_ith_slice,
      VOID_SLICEDTUBE_T_SET_ITH_SLICE_CONST_T_REF_INDEX,
      "codomain"_a, "i"_a)
    
    .def("inflate", (const SlicedTube<T>& (SlicedTube<T>::*)(const double&)) &SlicedTube<T>::inflate,
      CONST_SLICEDTUBE_T_REF_SLICEDTUBE_T_INFLATE_CONST_V_REF,
      "rad"_a)
    
    .def("inflate", (const SlicedTube<T>& (SlicedTube<T>::*)(const SampledTraj<double>&)) &SlicedTube<T>::inflate,
      CONST_SLICEDTUBE_T_REF_SLICEDTUBE_T_INFLATE_CONST_SAMPLEDTRAJ_V_REF,
      "rad"_a)
  ;

  if constexpr(std::is_same_v<T,IntervalVector>)
  {
    exported_slicedtube_class

    .def("inflate", (const SlicedTube<T>& (SlicedTube<T>::*)(const Vector&)) &SlicedTube<T>::inflate,
      CONST_SLICEDTUBE_T_REF_SLICEDTUBE_T_INFLATE_CONST_V_REF,
      "rad"_a)
    
    .def("inflate", (const SlicedTube<T>& (SlicedTube<T>::*)(const SampledTraj<Vector>&)) &SlicedTube<T>::inflate,
      CONST_SLICEDTUBE_T_REF_SLICEDTUBE_T_INFLATE_CONST_SAMPLEDTRAJ_V_REF,
      "rad"_a)
    ;
  }

  exported_slicedtube_class

    .def(py::self == py::self,
      BOOL_SLICEDTUBE_T_OPERATOREQ_CONST_SLICEDTUBE_REF_CONST,
      "x"_a)

    .def(py::self != py::self,
      "x"_a)

    .def(py::self &= py::self,
      SLICEDTUBE_SLICEDTUBE_T_OPERATORINTEREQ_CONST_SLICEDTUBE_REF,
      "x"_a)

    // For MATLAB compatibility
    .def("self_inter", &SlicedTube<T>::operator&=,
      SLICEDTUBE_SLICEDTUBE_T_OPERATORINTEREQ_CONST_SLICEDTUBE_REF,
      "x"_a)

    .def("as_function", &SlicedTube<T>::as_function,
      ANALYTICFUNCTION_TYPENAME_EXPRTYPE_T_TYPE_SLICEDTUBE_T_AS_FUNCTION_CONST)

    .def("all_reals_value", &SlicedTube<T>::all_reals_value,
      T_SLICEDTUBE_T_ALL_REALS_VALUE_CONST)

    .def("empty_value", &SlicedTube<T>::empty_value,
      T_SLICEDTUBE_T_EMPTY_VALUE_CONST)

    .def("__repr__", [](const SlicedTube<T>& x) {
          std::ostringstream stream;
          stream << x;
          return string(stream.str()); 
        },
      OSTREAM_REF_OPERATOROUT_OSTREAM_REF_CONST_SLICEDTUBE_T_REF)
  ;

  if constexpr(std::is_same_v<T,Interval> || std::is_same_v<T,IntervalVector>)
  {
    exported_slicedtube_class
    
      .def("__call__", [](const SlicedTube<T>& x, const Interval& t, const SlicedTube<T>& v)
          {
            return x(t,v);
          },
        T_SLICEDTUBE_T_OPERATORCALL_CONST_INTERVAL_REF_CONST_SLICEDTUBE_T_REF_CONST,
        "t"_a, "v"_a)

      .def("integral", (T (SlicedTube<T>::*)(const Interval&) const) &SlicedTube<T>::integral,
        T_SLICEDTUBE_T_INTEGRAL_CONST_INTERVAL_REF_CONST,
        "t"_a)
      
      .def("integral", (T (SlicedTube<T>::*)(const Interval&,const Interval&) const) &SlicedTube<T>::integral,
        T_SLICEDTUBE_T_INTEGRAL_CONST_INTERVAL_REF_CONST_INTERVAL_REF_CONST,
        "t1"_a, "t2"_a)
      
      .def("partial_integral", (std::pair<T,T> (SlicedTube<T>::*)(const Interval&) const) &SlicedTube<T>::partial_integral,
        PAIR_TT_SLICEDTUBE_T_PARTIAL_INTEGRAL_CONST_INTERVAL_REF_CONST,
        "t"_a)
      
      .def("partial_integral", (std::pair<T,T> (SlicedTube<T>::*)(const Interval&,const Interval&) const) &SlicedTube<T>::partial_integral,
        PAIR_TT_SLICEDTUBE_T_PARTIAL_INTEGRAL_CONST_INTERVAL_REF_CONST_INTERVAL_REF_CONST,
        "t1"_a, "t2"_a)

      .def("primitive", (SlicedTube<T> (SlicedTube<T>::*)() const) &SlicedTube<T>::primitive,
        SLICEDTUBE_T_SLICEDTUBE_T_PRIMITIVE_CONST)

      .def("primitive", (SlicedTube<T> (SlicedTube<T>::*)(const T&) const) &SlicedTube<T>::primitive,
        SLICEDTUBE_T_SLICEDTUBE_T_PRIMITIVE_CONST_T_REF_CONST,
        "x0"_a)

      .def("invert", (Interval (SlicedTube<T>::*)(const T&,const Interval&) const) &SlicedTube<T>::invert,
        INTERVAL_SLICEDTUBE_T_INVERT_CONST_T_REF_CONST_INTERVAL_REF_CONST,
        "y"_a, "t"_a=Interval())
      
      .def("invert", (void (SlicedTube<T>::*)(const T&,std::vector<Interval>&,const Interval&) const) &SlicedTube<T>::invert,
        VOID_SLICEDTUBE_T_INVERT_CONST_T_REF_VECTOR_INTERVAL_REF_CONST_INTERVAL_REF_CONST,
        "y"_a, "v_t"_a, "t"_a=Interval())
      
      .def("invert", (Interval (SlicedTube<T>::*)(const T&,const SlicedTube<T>&,const Interval&) const) &SlicedTube<T>::invert,
        INTERVAL_SLICEDTUBE_T_INVERT_CONST_T_REF_CONST_SLICEDTUBE_T_REF_CONST_INTERVAL_REF_CONST,
        "y"_a, "v"_a, "t"_a=Interval())
      
      .def("invert", (void (SlicedTube<T>::*)(const T&,std::vector<Interval>&,const SlicedTube<T>&,const Interval&) const) &SlicedTube<T>::invert,
        VOID_SLICEDTUBE_T_INVERT_CONST_T_REF_VECTOR_INTERVAL_REF_CONST_SLICEDTUBE_T_REF_CONST_INTERVAL_REF_CONST,
        "y"_a, "v_t"_a, "v"_a, "t"_a=Interval())

    ;
  }

  if constexpr(std::is_same_v<T,IntervalVector>)
  {
    exported_slicedtube_class

      .def(
          #if FOR_MATLAB
            "get_item"
          #else
            "__getitem__"
          #endif
          ,
          [](const SlicedTube<IntervalVector>& x, Index_type i) -> SlicedTube<Interval>
          {
            matlab::test_integer(i);
            return x[matlab::input_index(i)];
          },
        SLICEDTUBE_INTERVAL_SLICEDTUBE_T_OPERATORCOMPO_INDEX_CONST,
        "i"_a)
        
      .def("get_item_0",
          [](const SlicedTube<IntervalVector>& x, Index_type i) -> SlicedTube<Interval>
          {
            matlab::test_integer(i);
            return x[i];
          },
        SLICEDTUBE_INTERVAL_SLICEDTUBE_T_OPERATORCOMPO_INDEX_CONST,
        "i"_a)

      .def("subvector",
          [](const SlicedTube<IntervalVector>& x, Index_type i, Index_type j) -> SlicedTube<IntervalVector>
          {
            matlab::test_integer(i,j);
            return x.subvector(matlab::input_index(i),matlab::input_index(j));
          },
        SLICEDTUBE_INTERVALVECTOR_SLICEDTUBE_T_SUBVECTOR_INDEX_INDEX_CONST,
        "i"_a, "j"_a)

      .def("subvector_0",
          [](const SlicedTube<IntervalVector>& x, Index_type i, Index_type j) -> SlicedTube<IntervalVector>
          {
            matlab::test_integer(i,j);
            return x.subvector(i,j);
          },
        SLICEDTUBE_INTERVALVECTOR_SLICEDTUBE_T_SUBVECTOR_INDEX_INDEX_CONST,
        "i"_a, "j"_a)
    ;
  }

  if constexpr(std::is_same_v<T,Interval>)
  {
    exported_slicedtube_class
      .def("mid", &SlicedTube<Interval>::mid<double>,
        AUTO_SLICEDTUBE_T_MID_CONST);
  }

  if constexpr(std::is_same_v<T,IntervalVector>)
  {
    exported_slicedtube_class
      .def("mid", &SlicedTube<IntervalVector>::mid<Vector>,
        AUTO_SLICEDTUBE_T_MID_CONST);
  }

  if constexpr(std::is_same_v<T,IntervalMatrix>)
  {
    exported_slicedtube_class
      .def("mid", &SlicedTube<IntervalMatrix>::mid<Matrix>,
        AUTO_SLICEDTUBE_T_MID_CONST);
  }

  return exported_slicedtube_class;
}

template<typename T>
void add_tube_operators(py::class_<SlicedTube<T>,TubeBase>& pyclass)
{
  pyclass

    .def("__add__", [](py::object& x1) { return cast<SlicedTube<T>>(x1); },
      CONST_SLICEDTUBE_T_REF_OPERATORPLUS_CONST_SLICEDTUBE_T_REF,
      py::is_operator())

    .def("__add__", [](py::object& x1, py::object& x2) { return cast<SlicedTube<T>>(x1)+cast<SlicedTube<T>>(x2); },
      SLICEDTUBE_T_OPERATORPLUS_CONST_SLICEDTUBE_T_REF_CONST_SLICEDTUBE_T_REF,
      py::is_operator())

    .def("__add__", [](py::object& x1, const T& x2) { return cast<SlicedTube<T>>(x1)+x2; },
      SLICEDTUBE_T_OPERATORPLUS_CONST_SLICEDTUBE_T_REF_CONST_Q_REF,
      py::is_operator())

    .def("__radd__", [](py::object& x2, const T& x1) { return x1+cast<SlicedTube<T>>(x2); },
      SLICEDTUBE_T_OPERATORPLUS_CONST_Q_REF_CONST_SLICEDTUBE_T_REF,
      py::is_operator())

    .def("__sub__", [](py::object& x1) { return -cast<SlicedTube<T>>(x1); },
      SLICEDTUBE_T_OPERATORMINUS_CONST_SLICEDTUBE_T_REF,
      py::is_operator())

    .def("__sub__", [](py::object& x1, py::object& x2) { return cast<SlicedTube<T>>(x1)-cast<SlicedTube<T>>(x2); },
      SLICEDTUBE_T_OPERATORMINUS_CONST_SLICEDTUBE_T_REF_CONST_SLICEDTUBE_T_REF,
      py::is_operator())

    .def("__sub__", [](py::object& x1, const T& x2) { return cast<SlicedTube<T>>(x1)-x2; },
      SLICEDTUBE_T_OPERATORMINUS_CONST_SLICEDTUBE_T_REF_CONST_Q_REF,
      py::is_operator())

    .def("__rsub__", [](py::object& x2, const T& x1) { return x1-cast<SlicedTube<T>>(x2); },
      SLICEDTUBE_T_OPERATORMINUS_CONST_Q_REF_CONST_SLICEDTUBE_T_REF,
      py::is_operator())

    .def("__mul__", [](py::object& x1, const Interval& x2) { return x2*cast<SlicedTube<T>>(x1); },
      SLICEDTUBE_T_OPERATORMUL_CONST_SLICEDTUBE_T_REF_CONST_INTERVAL_REF,
      py::is_operator())

    .def("__rmul__", [](py::object& x2, const Interval& x1) { return x1*cast<SlicedTube<T>>(x2); },
      SLICEDTUBE_T_OPERATORMUL_CONST_INTERVAL_REF_CONST_SLICEDTUBE_T_REF,
      py::is_operator())

    .def("__mul__", [](const Interval& x1, py::object& x2) { return x1*cast<SlicedTube<T>>(x2); },
      SLICEDTUBE_T_OPERATORMUL_CONST_INTERVAL_REF_CONST_SLICEDTUBE_T_REF,
      py::is_operator())

    .def("__truediv__", [](py::object& x1, const Interval& x2) { return cast<SlicedTube<T>>(x1)/x2; },
      SLICEDTUBE_T_OPERATORDIV_CONST_SLICEDTUBE_T_REF_CONST_Q_REF,
      py::is_operator())
  ;
}


void export_SlicedTube(py::module& m)
{
  auto py_SlicedTube_Interval = _export_SlicedTube<Interval>(m, "SlicedTube_Interval");
  auto py_SlicedTube_IntervalVector = _export_SlicedTube<IntervalVector>(m, "SlicedTube_IntervalVector");
  auto py_SlicedTube_IntervalMatrix = _export_SlicedTube<IntervalMatrix>(m, "SlicedTube_IntervalMatrix");

  add_tube_operators<Interval>(py_SlicedTube_Interval);

  py_SlicedTube_Interval

    .def(py::self * py::self,
      SLICEDTUBE_T_OPERATORMUL_CONST_SLICEDTUBE_T_REF_CONST_SLICEDTUBE_T_REF,
      py::is_operator())

    .def(py::self / py::self,
      SLICEDTUBE_T_OPERATORDIV_CONST_SLICEDTUBE_T_REF_CONST_SLICEDTUBE_T_REF,
      py::is_operator())

    .def(double() / py::self,
      SLICEDTUBE_T_OPERATORDIV_CONST_Q_REF_CONST_SLICEDTUBE_T_REF,
      py::is_operator())
  ;

  add_tube_operators<IntervalMatrix>(py_SlicedTube_IntervalMatrix);
  add_tube_operators<IntervalVector>(py_SlicedTube_IntervalVector);
  py_SlicedTube_IntervalMatrix

    .def("__mul__", [](const SlicedTube<IntervalMatrix>& x1, const SlicedTube<IntervalVector>& x2) { return x1*x2; },
      SLICEDTUBE_INTERVALVECTOR_OPERATORMUL_CONST_SLICEDTUBE_INTERVALMATRIX_REF_CONST_SLICEDTUBE_INTERVALVECTOR_REF,
      py::is_operator())

    .def("__mul__", [](const SlicedTube<IntervalMatrix>& x1, const SlicedTube<IntervalMatrix>& x2) { return x1*x2; },
      SLICEDTUBE_T_OPERATORMUL_CONST_SLICEDTUBE_T_REF_CONST_SLICEDTUBE_T_REF,
      py::is_operator())

  ;

  m

    .def("sqr", (SlicedTube<Interval> (*)(const SlicedTube<Interval>&)) &codac2::sqr,
      SLICEDTUBE_INTERVAL_SQR_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, py::is_operator())

    .def("sqrt", (SlicedTube<Interval> (*)(const SlicedTube<Interval>&)) &codac2::sqrt,
      SLICEDTUBE_INTERVAL_SQRT_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, py::is_operator())

    .def("pow", (SlicedTube<Interval> (*)(const SlicedTube<Interval>&,const Interval&)) &codac2::pow,
      SLICEDTUBE_INTERVAL_POW_CONST_SLICEDTUBE_INTERVAL_REF_CONST_INTERVAL_REF,
      "x1"_a, "x2"_a, py::is_operator())

    //.def("root", (SlicedTube<Interval> (*)(const SlicedTube<Interval>&,int)) &codac2::root,
    //  SAMPLEDTRAJ_DOUBLE_ROOT_CONST_SAMPLEDTRAJ_DOUBLE_REF_INT,
    //  "x1"_a, "x2"_a, py::is_operator())

    .def("exp", (SlicedTube<Interval> (*)(const SlicedTube<Interval>&)) &codac2::exp,
      SLICEDTUBE_INTERVAL_EXP_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, py::is_operator())

    .def("log", (SlicedTube<Interval> (*)(const SlicedTube<Interval>&)) &codac2::log,
      SLICEDTUBE_INTERVAL_LOG_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, py::is_operator())

    .def("cos", (SlicedTube<Interval> (*)(const SlicedTube<Interval>&)) &codac2::cos,
      SLICEDTUBE_INTERVAL_COS_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, py::is_operator())

    .def("sin", (SlicedTube<Interval> (*)(const SlicedTube<Interval>&)) &codac2::sin,
      SLICEDTUBE_INTERVAL_SIN_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, py::is_operator())

    .def("tan", (SlicedTube<Interval> (*)(const SlicedTube<Interval>&)) &codac2::tan,
      SLICEDTUBE_INTERVAL_TAN_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, py::is_operator())

    .def("acos", (SlicedTube<Interval> (*)(const SlicedTube<Interval>&)) &codac2::acos,
      SLICEDTUBE_INTERVAL_ACOS_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, py::is_operator())

    .def("asin", (SlicedTube<Interval> (*)(const SlicedTube<Interval>&)) &codac2::asin,
      SLICEDTUBE_INTERVAL_ASIN_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, py::is_operator())

    .def("atan", (SlicedTube<Interval> (*)(const SlicedTube<Interval>&)) &codac2::atan,
      SLICEDTUBE_INTERVAL_ATAN_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, py::is_operator())

    .def("atan2", (SlicedTube<Interval> (*)(const SlicedTube<Interval>&,const SlicedTube<Interval>&)) &codac2::atan2,
      SLICEDTUBE_INTERVAL_ATAN2_CONST_SLICEDTUBE_INTERVAL_REF_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, "x2"_a, py::is_operator())

    .def("atan2", (SlicedTube<Interval> (*)(const SlicedTube<Interval>&,const Interval&)) &codac2::atan2,
      SLICEDTUBE_INTERVAL_ATAN2_CONST_SLICEDTUBE_INTERVAL_REF_CONST_INTERVAL_REF,
      "x1"_a, "x2"_a, py::is_operator())

    .def("atan2", (SlicedTube<Interval> (*)(const Interval&,const SlicedTube<Interval>&)) &codac2::atan2,
      SLICEDTUBE_INTERVAL_ATAN2_CONST_INTERVAL_REF_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, "x2"_a, py::is_operator())

    .def("cosh", (SlicedTube<Interval> (*)(const SlicedTube<Interval>&)) &codac2::cosh,
      SLICEDTUBE_INTERVAL_COSH_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, py::is_operator())

    .def("sinh", (SlicedTube<Interval> (*)(const SlicedTube<Interval>&)) &codac2::sinh,
      SLICEDTUBE_INTERVAL_SINH_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, py::is_operator())

    .def("tanh", (SlicedTube<Interval> (*)(const SlicedTube<Interval>&)) &codac2::tanh,
      SLICEDTUBE_INTERVAL_TANH_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, py::is_operator())

    .def("acosh", (SlicedTube<Interval> (*)(const SlicedTube<Interval>&)) &codac2::acosh,
      SLICEDTUBE_INTERVAL_ACOSH_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, py::is_operator())

    .def("asinh", (SlicedTube<Interval> (*)(const SlicedTube<Interval>&)) &codac2::asinh,
      SLICEDTUBE_INTERVAL_ASINH_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, py::is_operator())

    .def("atanh", (SlicedTube<Interval> (*)(const SlicedTube<Interval>&)) &codac2::atanh,
      SLICEDTUBE_INTERVAL_ATANH_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, py::is_operator())

    .def("abs", (SlicedTube<Interval> (*)(const SlicedTube<Interval>&)) &codac2::abs,
      SLICEDTUBE_INTERVAL_ABS_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, py::is_operator())

    .def("min", (SlicedTube<Interval> (*)(const SlicedTube<Interval>&,const SlicedTube<Interval>&)) &codac2::min,
      SLICEDTUBE_INTERVAL_MIN_CONST_SLICEDTUBE_INTERVAL_REF_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, "x2"_a, py::is_operator())

    .def("min", (SlicedTube<Interval> (*)(const SlicedTube<Interval>&,const Interval&)) &codac2::min,
      SLICEDTUBE_INTERVAL_MIN_CONST_SLICEDTUBE_INTERVAL_REF_CONST_INTERVAL_REF,
      "x1"_a, "x2"_a, py::is_operator())

    .def("min", (SlicedTube<Interval> (*)(const Interval&,const SlicedTube<Interval>&)) &codac2::min,
      SLICEDTUBE_INTERVAL_MIN_CONST_INTERVAL_REF_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, "x2"_a, py::is_operator())

    .def("max", (SlicedTube<Interval> (*)(const SlicedTube<Interval>&,const SlicedTube<Interval>&)) &codac2::max,
      SLICEDTUBE_INTERVAL_MAX_CONST_SLICEDTUBE_INTERVAL_REF_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, "x2"_a, py::is_operator())

    .def("max", (SlicedTube<Interval> (*)(const SlicedTube<Interval>&,const Interval&)) &codac2::max,
      SLICEDTUBE_INTERVAL_MAX_CONST_SLICEDTUBE_INTERVAL_REF_CONST_INTERVAL_REF,
      "x1"_a, "x2"_a, py::is_operator())

    .def("max", (SlicedTube<Interval> (*)(const Interval&,const SlicedTube<Interval>&)) &codac2::max,
      SLICEDTUBE_INTERVAL_MAX_CONST_INTERVAL_REF_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, "x2"_a, py::is_operator())

    .def("sign", (SlicedTube<Interval> (*)(const SlicedTube<Interval>&)) &codac2::sign,
      SLICEDTUBE_INTERVAL_SIGN_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, py::is_operator())

    .def("integer", (SlicedTube<Interval> (*)(const SlicedTube<Interval>&)) &codac2::integer,
      SLICEDTUBE_INTERVAL_INTEGER_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, py::is_operator())

    .def("floor", (SlicedTube<Interval> (*)(const SlicedTube<Interval>&)) &codac2::floor,
      SLICEDTUBE_INTERVAL_FLOOR_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, py::is_operator())

    .def("ceil", (SlicedTube<Interval> (*)(const SlicedTube<Interval>&)) &codac2::ceil,
      SLICEDTUBE_INTERVAL_CEIL_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, py::is_operator())

  ;
}



#if 0

// Keeping the following just in case...
// Equivalent to the previous block of code, but using py::object for tubes

void export_SlicedTube(py::module& m)
{
  auto py_SlicedTube_Interval = _export_SlicedTube<Interval>(m, "SlicedTube_Interval");
  auto py_SlicedTube_IntervalVector = _export_SlicedTube<IntervalVector>(m, "SlicedTube_IntervalVector");
  auto py_SlicedTube_IntervalMatrix = _export_SlicedTube<IntervalMatrix>(m, "SlicedTube_IntervalMatrix");

  add_tube_operators<Interval>(py_SlicedTube_Interval);

  py_SlicedTube_Interval

    .def(py::self * py::self,
      SLICEDTUBE_T_OPERATORMUL_CONST_SLICEDTUBE_T_REF_CONST_SLICEDTUBE_T_REF,
      py::is_operator())

    .def(py::self / py::self,
      SLICEDTUBE_T_OPERATORDIV_CONST_SLICEDTUBE_T_REF_CONST_SLICEDTUBE_T_REF,
      py::is_operator())

    .def(double() / py::self,
      SLICEDTUBE_T_OPERATORDIV_CONST_Q_REF_CONST_SLICEDTUBE_T_REF,
      py::is_operator())
  ;

  add_tube_operators<IntervalMatrix>(py_SlicedTube_IntervalMatrix);
  add_tube_operators<IntervalVector>(py_SlicedTube_IntervalVector);
  py_SlicedTube_IntervalMatrix

    .def("__mul__", [](py::object& x1, py::object& x2) { return cast<SlicedTube<IntervalMatrix>>(x1)*cast<SlicedTube<IntervalVector>>(x2); },
      SLICEDTUBE_INTERVALVECTOR_OPERATORMUL_CONST_SLICEDTUBE_INTERVALMATRIX_REF_CONST_SLICEDTUBE_INTERVALVECTOR_REF,
      py::is_operator())

    .def("__mul__", [](py::object& x1, py::object& x2) { return cast<SlicedTube<IntervalMatrix>>(x1)*cast<SlicedTube<IntervalMatrix>>(x2); },
      SLICEDTUBE_T_OPERATORMUL_CONST_SLICEDTUBE_T_REF_CONST_SLICEDTUBE_T_REF,
      py::is_operator())

  ;

  m

    .def("sqr", [](py::object& x1) { return codac2::sqr(cast<SlicedTube<Interval>>(x1)); },
      SLICEDTUBE_INTERVAL_SQR_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, py::is_operator())

    .def("sqrt", [](py::object& x1) { return codac2::sqrt(cast<SlicedTube<Interval>>(x1)); },
      SLICEDTUBE_INTERVAL_SQRT_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, py::is_operator())

    .def("pow", [](py::object& x1, const Interval& x2) { return codac2::pow(cast<SlicedTube<Interval>>(x1),x2); },
      SLICEDTUBE_INTERVAL_POW_CONST_SLICEDTUBE_INTERVAL_REF_CONST_INTERVAL_REF,
      "x1"_a, "x2"_a, py::is_operator())

    //.def("root", (SlicedTube<Interval> (*)(const SlicedTube<Interval>&,int)) &codac2::root,
    //  SAMPLEDTRAJ_DOUBLE_ROOT_CONST_SAMPLEDTRAJ_DOUBLE_REF_INT,
    //  "x1"_a, "x2"_a, py::is_operator())

    .def("exp", [](py::object& x1) { return codac2::exp(cast<SlicedTube<Interval>>(x1)); },
      SLICEDTUBE_INTERVAL_EXP_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, py::is_operator())

    .def("log", [](py::object& x1) { return codac2::log(cast<SlicedTube<Interval>>(x1)); },
      SLICEDTUBE_INTERVAL_LOG_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, py::is_operator())

    .def("cos", [](py::object& x1) { return codac2::cos(cast<SlicedTube<Interval>>(x1)); },
      SLICEDTUBE_INTERVAL_COS_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, py::is_operator())

    .def("sin", [](py::object& x1) { return codac2::sin(cast<SlicedTube<Interval>>(x1)); },
      SLICEDTUBE_INTERVAL_SIN_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, py::is_operator())

    .def("tan", [](py::object& x1) { return codac2::tan(cast<SlicedTube<Interval>>(x1)); },
      SLICEDTUBE_INTERVAL_TAN_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, py::is_operator())

    .def("acos", [](py::object& x1) { return codac2::acos(cast<SlicedTube<Interval>>(x1)); },
      SLICEDTUBE_INTERVAL_ACOS_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, py::is_operator())

    .def("asin", [](py::object& x1) { return codac2::asin(cast<SlicedTube<Interval>>(x1)); },
      SLICEDTUBE_INTERVAL_ASIN_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, py::is_operator())

    .def("atan", [](py::object& x1) { return codac2::atan(cast<SlicedTube<Interval>>(x1)); },
      SLICEDTUBE_INTERVAL_ATAN_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, py::is_operator())

    .def("atan2", [](py::object& x1, py::object& x2) { return codac2::atan2(cast<SlicedTube<Interval>>(x1),cast<SlicedTube<Interval>>(x2)); },
      SLICEDTUBE_INTERVAL_ATAN2_CONST_SLICEDTUBE_INTERVAL_REF_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, "x2"_a, py::is_operator())

    .def("atan2", [](py::object& x1, const Interval& x2) { return codac2::atan2(cast<SlicedTube<Interval>>(x1),x2); },
      SLICEDTUBE_INTERVAL_ATAN2_CONST_SLICEDTUBE_INTERVAL_REF_CONST_INTERVAL_REF,
      "x1"_a, "x2"_a, py::is_operator())

    .def("atan2", [](const Interval& x1, py::object& x2) { return codac2::atan2(x1,cast<SlicedTube<Interval>>(x2)); },
      SLICEDTUBE_INTERVAL_ATAN2_CONST_INTERVAL_REF_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, "x2"_a, py::is_operator())

    .def("cosh", [](py::object& x1) { return codac2::cosh(cast<SlicedTube<Interval>>(x1)); },
      SLICEDTUBE_INTERVAL_COSH_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, py::is_operator())

    .def("sinh", [](py::object& x1) { return codac2::sinh(cast<SlicedTube<Interval>>(x1)); },
      SLICEDTUBE_INTERVAL_SINH_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, py::is_operator())

    .def("tanh", [](py::object& x1) { return codac2::tanh(cast<SlicedTube<Interval>>(x1)); },
      SLICEDTUBE_INTERVAL_TANH_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, py::is_operator())

    .def("acosh", [](py::object& x1) { return codac2::acosh(cast<SlicedTube<Interval>>(x1)); },
      SLICEDTUBE_INTERVAL_ACOSH_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, py::is_operator())

    .def("asinh", [](py::object& x1) { return codac2::asinh(cast<SlicedTube<Interval>>(x1)); },
      SLICEDTUBE_INTERVAL_ASINH_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, py::is_operator())

    .def("atanh", [](py::object& x1) { return codac2::atanh(cast<SlicedTube<Interval>>(x1)); },
      SLICEDTUBE_INTERVAL_ATANH_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, py::is_operator())

    .def("abs", [](py::object& x1) { return codac2::abs(cast<SlicedTube<Interval>>(x1)); },
      SLICEDTUBE_INTERVAL_ABS_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, py::is_operator())

    .def("min", [](py::object& x1, py::object& x2) { return codac2::min(cast<SlicedTube<Interval>>(x1),cast<SlicedTube<Interval>>(x2)); },
      SLICEDTUBE_INTERVAL_MIN_CONST_SLICEDTUBE_INTERVAL_REF_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, "x2"_a, py::is_operator())

    .def("min", [](py::object& x1, const Interval& x2) { return codac2::min(cast<SlicedTube<Interval>>(x1),x2); },
      SLICEDTUBE_INTERVAL_MIN_CONST_SLICEDTUBE_INTERVAL_REF_CONST_INTERVAL_REF,
      "x1"_a, "x2"_a, py::is_operator())

    .def("min", [](const Interval& x1, py::object& x2) { return codac2::min(x1,cast<SlicedTube<Interval>>(x2)); },
      SLICEDTUBE_INTERVAL_MIN_CONST_INTERVAL_REF_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, "x2"_a, py::is_operator())

    .def("max", [](py::object& x1, py::object& x2) { return codac2::max(cast<SlicedTube<Interval>>(x1),cast<SlicedTube<Interval>>(x2)); },
      SLICEDTUBE_INTERVAL_MAX_CONST_SLICEDTUBE_INTERVAL_REF_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, "x2"_a, py::is_operator())

    .def("max", [](py::object& x1, const Interval& x2) { return codac2::max(cast<SlicedTube<Interval>>(x1),x2); },
      SLICEDTUBE_INTERVAL_MAX_CONST_SLICEDTUBE_INTERVAL_REF_CONST_INTERVAL_REF,
      "x1"_a, "x2"_a, py::is_operator())

    .def("max", [](const Interval& x1, py::object& x2) { return codac2::max(x1,cast<SlicedTube<Interval>>(x2)); },
      SLICEDTUBE_INTERVAL_MAX_CONST_INTERVAL_REF_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, "x2"_a, py::is_operator())

    .def("sign", [](py::object& x1) { return codac2::sign(cast<SlicedTube<Interval>>(x1)); },
      SLICEDTUBE_INTERVAL_SIGN_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, py::is_operator())

    .def("integer", [](py::object& x1) { return codac2::integer(cast<SlicedTube<Interval>>(x1)); },
      SLICEDTUBE_INTERVAL_INTEGER_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, py::is_operator())

    .def("floor", [](py::object& x1) { return codac2::floor(cast<SlicedTube<Interval>>(x1)); },
      SLICEDTUBE_INTERVAL_FLOOR_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, py::is_operator())

    .def("ceil", [](py::object& x1) { return codac2::ceil(cast<SlicedTube<Interval>>(x1)); },
      SLICEDTUBE_INTERVAL_CEIL_CONST_SLICEDTUBE_INTERVAL_REF,
      "x1"_a, py::is_operator())

  ;
}
#endif