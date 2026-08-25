/**
 * \file codac2_py_Affine.cpp
 *
 * Python binding of codac2::Affine.
 */

#include "codac2_py_Affine.h"

#include <sstream>

#include <pybind11/operators.h>

#include "codac2_Affine.h"

namespace py = pybind11;

namespace codac2 {

void export_Affine(py::module_& m)
{
  using A = Affine;

  py::enum_<AffineStatus>(m, "AffineStatus")
    .value("LowerUnbounded", AffineStatus::LowerUnbounded)
    .value("UpperUnbounded", AffineStatus::UpperUnbounded)
    .value("AllReals", AffineStatus::AllReals)
    .value("Empty", AffineStatus::Empty)
    .value("Degenerate", AffineStatus::Degenerate)
    .value("Active", AffineStatus::Active);

  py::enum_<A::Affine_Mode>(m, "AffineMode")
    .value("AF_Lin_Chebyshev", A::AF_Lin_Chebyshev)
    .value("AF_Lin_MinRange", A::AF_Lin_MinRange);

  py::class_<A> cls(m, "Affine");

  cls
    .def(py::init<>())
    .def(py::init<double>(), py::arg("value"))
    .def(py::init<const Interval&>(), py::arg("interval"))
    .def(py::init<const A&>(), py::arg("other"))
    .def_static("empty", &A::empty)

    .def_static("change_mode", &A::change_mode,
                py::arg("mode") = A::AF_Lin_Chebyshev)
    .def_static("get_mode", &A::get_mode)

    .def("__neg__", [](const A& x) { return -x; })

    .def("__add__", [](const A& x, const A& y) { return x + y; })
    .def("__add__", [](const A& x, double y) { return x + y; })
    .def("__add__", [](const A& x, const Interval& y) { return x + y; })
    .def("__radd__", [](const A& x, double y) { return y + x; })
    .def("__radd__", [](const A& x, const Interval& y) { return y + x; })

    .def("__sub__", [](const A& x, const A& y) { return x - y; })
    .def("__sub__", [](const A& x, double y) { return x - y; })
    .def("__sub__", [](const A& x, const Interval& y) { return x - y; })
    .def("__rsub__", [](const A& x, double y) { return y - x; })
    .def("__rsub__", [](const A& x, const Interval& y) { return y - x; })

    .def("__mul__", [](const A& x, const A& y) { return x * y; })
    .def("__mul__", [](const A& x, double y) { return x * y; })
    .def("__mul__", [](const A& x, const Interval& y) { return x * y; })
    .def("__rmul__", [](const A& x, double y) { return y * x; })
    .def("__rmul__", [](const A& x, const Interval& y) { return y * x; })

    .def("__truediv__", [](const A& x, const A& y) { return x / y; })
    .def("__truediv__", [](const A& x, double y) { return x / y; })
    .def("__truediv__", [](const A& x, const Interval& y) { return x / y; })
    .def("__rtruediv__", [](const A& x, double y) { return y / x; })
    .def("__rtruediv__", [](const A& x, const Interval& y) { return y / x; })

    .def("__iadd__", [](A& x, const A& y) -> A& { x += y; return x; })
    .def("__iadd__", [](A& x, double y) -> A& { x += y; return x; })
    .def("__iadd__", [](A& x, const Interval& y) -> A& { x += y; return x; })
    .def("__isub__", [](A& x, const A& y) -> A& { x -= y; return x; })
    .def("__isub__", [](A& x, double y) -> A& { x -= y; return x; })
    .def("__isub__", [](A& x, const Interval& y) -> A& { x -= y; return x; })
    .def("__imul__", [](A& x, const A& y) -> A& { x *= y; return x; })
    .def("__imul__", [](A& x, double y) -> A& { x *= y; return x; })
    .def("__imul__", [](A& x, const Interval& y) -> A& { x *= y; return x; })
    .def("__itruediv__", [](A& x, const A& y) -> A& { x /= y; return x; })
    .def("__itruediv__", [](A& x, double y) -> A& { x /= y; return x; })
    .def("__itruediv__", [](A& x, const Interval& y) -> A& { x /= y; return x; })

    .def("__eq__", [](const A& x, const A& y) { return x == y; })
    .def("__eq__", [](const A& x, const Interval& y) { return x == y; })
    .def("__eq__", [](const A& x, double y) { return x == y; })
    .def("__ne__", [](const A& x, const A& y) { return x != y; })
    .def("__ne__", [](const A& x, const Interval& y) { return x != y; })
    .def("__ne__", [](const A& x, double y) { return x != y; })

    .def("__lt__", [](const A& x, const A& y) { return x < y; })
    .def("__lt__", [](const A& x, const Interval& y) { return x < y; })
    .def("__gt__", [](const A& x, const A& y) { return x > y; })
    .def("__gt__", [](const A& x, const Interval& y) { return x > y; })

    .def("lb", &A::lb)
    .def("ub", &A::ub)
    .def("mid", &A::mid)
    .def("rad", &A::rad)
    .def("diam", &A::diam)
    .def("mig", &A::mig)
    .def("mag", &A::mag)
    .def("smag", &A::smag)
    .def("smig", &A::smig)
    .def("volume", &A::volume)
    .def("noise_count", &A::noise_count)
    .def("noise", &A::noise, py::arg("index"))
    .def("err", &A::err)
    .def("itv", &A::itv)
    .def("inflate", &A::inflate, py::arg("radius"))
    .def("compact", py::overload_cast<>(&A::compact))
    .def("compact", py::overload_cast<double>(&A::compact), py::arg("tol"))
    .def("set_empty", &A::set_empty)
    .def("is_active", &A::is_active)
    .def("is_empty", &A::is_empty)
    .def("is_degenerated", &A::is_degenerated)
    .def("is_unbounded", &A::is_unbounded)
    .def("is_subset", py::overload_cast<const Interval&>(&A::is_subset, py::const_))
    .def("is_subset", py::overload_cast<const A&>(&A::is_subset, py::const_))
    .def("is_strict_subset", py::overload_cast<const Interval&>(&A::is_strict_subset, py::const_))
    .def("is_strict_subset", py::overload_cast<const A&>(&A::is_strict_subset, py::const_))
    .def("is_interior_subset", py::overload_cast<const Interval&>(&A::is_interior_subset, py::const_))
    .def("is_interior_subset", py::overload_cast<const A&>(&A::is_interior_subset, py::const_))
    .def("is_relative_interior_subset", py::overload_cast<const Interval&>(&A::is_relative_interior_subset, py::const_))
    .def("is_relative_interior_subset", py::overload_cast<const A&>(&A::is_relative_interior_subset, py::const_))
    .def("is_strict_interior_subset", py::overload_cast<const Interval&>(&A::is_strict_interior_subset, py::const_))
    .def("is_strict_interior_subset", py::overload_cast<const A&>(&A::is_strict_interior_subset, py::const_))
    .def("is_superset", py::overload_cast<const Interval&>(&A::is_superset, py::const_))
    .def("is_superset", py::overload_cast<const A&>(&A::is_superset, py::const_))
    .def("is_strict_superset", py::overload_cast<const Interval&>(&A::is_strict_superset, py::const_))
    .def("is_strict_superset", py::overload_cast<const A&>(&A::is_strict_superset, py::const_))
    .def("contains", &A::contains)
    .def("interior_contains", &A::interior_contains)
    .def("intersects", py::overload_cast<const Interval&>(&A::intersects, py::const_))
    .def("intersects", py::overload_cast<const A&>(&A::intersects, py::const_))
    .def("overlaps", py::overload_cast<const Interval&>(&A::overlaps, py::const_))
    .def("overlaps", py::overload_cast<const A&>(&A::overlaps, py::const_))
    .def("is_disjoint", py::overload_cast<const Interval&>(&A::is_disjoint, py::const_))
    .def("is_disjoint", py::overload_cast<const A&>(&A::is_disjoint, py::const_))

    .def("__str__", [](const A& x) {
      std::ostringstream os;
      os << x;
      return os.str();
    })
    .def("__repr__", [](const A& x) {
      std::ostringstream os;
      os << x;
      return os.str();
    });

  m.def("distance", [](const A& x, const A& y) { return distance<AF_Default>(x, y); });
  m.def("distance", [](const Interval& x, const A& y) { return distance<AF_Default>(x, y); });
  m.def("distance", [](const A& x, const Interval& y) { return distance<AF_Default>(x, y); });

#define CODAC_AFFINE_UNARY(NAME) \
  m.def(#NAME, [](const A& x) { return NAME<AF_Default>(x); })

  CODAC_AFFINE_UNARY(inv);
  CODAC_AFFINE_UNARY(sqr);
  CODAC_AFFINE_UNARY(sqrt);
  CODAC_AFFINE_UNARY(exp);
  CODAC_AFFINE_UNARY(log);
  CODAC_AFFINE_UNARY(root);
  CODAC_AFFINE_UNARY(sin);
  CODAC_AFFINE_UNARY(cos);
  CODAC_AFFINE_UNARY(tan);
  CODAC_AFFINE_UNARY(asin);
  CODAC_AFFINE_UNARY(acos);
  CODAC_AFFINE_UNARY(atan);
  CODAC_AFFINE_UNARY(sinh);
  CODAC_AFFINE_UNARY(cosh);
  CODAC_AFFINE_UNARY(tanh);
  CODAC_AFFINE_UNARY(asinh);
  CODAC_AFFINE_UNARY(acosh);
  CODAC_AFFINE_UNARY(atanh);
  CODAC_AFFINE_UNARY(abs);
  CODAC_AFFINE_UNARY(sign);
#undef CODAC_AFFINE_UNARY

  m.def("root", [](const A& x, int n) { return root<AF_Default>(x, n); },
        py::arg("x"), py::arg("n"));
  m.def("pow", [](const A& x, int n) { return pow<AF_Default>(x, n); });
  m.def("pow", [](const A& x, double n) { return pow<AF_Default>(x, n); });
  m.def("pow", [](const A& x, const Interval& y) { return pow<AF_Default>(x, y); });
  m.def("pow", [](const A& x, const A& y) { return pow<AF_Default>(x, y); });
  m.def("pow", [](double x, const A& y) { return pow<AF_Default>(x, y); });
  m.def("pow", [](const Interval& x, const A& y) { return pow<AF_Default>(x, y); });
  m.def("atan2", [](const A& y, const A& x) { return atan2<AF_Default>(y, x); });
  m.def("chi", [](const A& a, const A& b, const A& c) {
    return chi<AF_Default>(a, b, c);
  });
}

} // namespace codac2
