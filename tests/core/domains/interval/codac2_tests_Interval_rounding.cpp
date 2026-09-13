/**
 *  Codac tests
 *
 *  Rounding of the Interval operations, which Codac delegates to GAOL. The
 *  bounds are checked exactly, independently of GAOL and of the floating-point
 *  flags: they have to enclose the exact result of each operation, as tightly
 *  as possible for the arithmetic operations, and they have to go on doing so
 *  after any operation, GAOL computing with the rounding direction set upward
 *  once for all.
 * ----------------------------------------------------------------------------
 *  \date       2026
 *  \author     Jordan Ninin
 *  \copyright  Copyright 2026 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <catch2/catch_test_macros.hpp>
#include <codac2_Interval.h>
#include <codac2_Interval_operations.h>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

using namespace std;
using namespace codac2;

namespace
{
  // The checks below compute with rounding to nearest, which this object sets
  // for its lifetime, restoring the rounding direction it found afterwards.
  class RoundingToNearest
  {
    public:

      RoundingToNearest()
        : _saved(std::fegetround())
      {
        std::fesetround(FE_TONEAREST);
      }

      ~RoundingToNearest()
      {
        std::fesetround(_saved);
      }

    private:

      const int _saved;
  };

  // The exact result of an operation on two doubles, as the sum hi+lo of two
  // doubles, hi being that result rounded to nearest. It is computed by the
  // error-free transformations of Knuth (sum) and Dekker (product), which are
  // exact with rounding to nearest, far from overflow and underflow.
  //
  // The operands are read back from memory (volatile) once rounding to nearest
  // is set. Visual C++ for 32-bit x86 with /arch:AVX or /arch:AVX2 otherwise
  // reuses a+b or a*b as the interval operation just checked computed it, with
  // rounding upward, in spite of /fp:strict: about half of the random sums and
  // products then came out rounded the wrong way. Codac leaves /arch:AVX2 out
  // for that target (see the top-level CMakeLists.txt), and the checks stay
  // right for a build that would add it.
  struct Exact
  {
    double hi, lo;
  };

  Exact exact_sum(double a0, double b0)
  {
    RoundingToNearest nearest;
    volatile double va = a0, vb = b0;
    const double a = va, b = vb;
    const double s = a + b;
    const double bb = s - a;
    return { s, (a - (s - bb)) + (b - bb) };
  }

  Exact exact_product(double a0, double b0)
  {
    RoundingToNearest nearest;
    volatile double va = a0, vb = b0;
    const double a = va, b = vb;
    const double p = a * b;
    // 2^27+1 splits a double into two halves of at most 27 bits
    const double ta = 134217729. * a;
    const double ah = ta - (ta - a);
    const double al = a - ah;
    const double tb = 134217729. * b;
    const double bh = tb - (tb - b);
    const double bl = b - bh;
    return { p, ((ah * bh - p) + ah * bl + al * bh) + al * bl };
  }

  // x <= v, compared exactly
  bool is_below(double x, const Exact& v)
  {
    return x < v.hi || (x == v.hi && v.lo >= 0.);
  }

  // x >= v, compared exactly
  bool is_above(double x, const Exact& v)
  {
    return x > v.hi || (x == v.hi && v.lo <= 0.);
  }

  // x <= a/b, compared exactly through the product x*b
  bool is_below_quotient(double x, double a, double b)
  {
    const Exact xb = exact_product(x, b);
    return b > 0. ? is_above(a, xb) : is_below(a, xb);
  }

  // x >= a/b, compared exactly through the product x*b
  bool is_above_quotient(double x, double a, double b)
  {
    const Exact xb = exact_product(x, b);
    return b > 0. ? is_below(a, xb) : is_above(a, xb);
  }

  double next_double(double x)
  {
    RoundingToNearest nearest;
    return std::nextafter(x, numeric_limits<double>::infinity());
  }

  double previous_double(double x)
  {
    RoundingToNearest nearest;
    return std::nextafter(x, -numeric_limits<double>::infinity());
  }

  // Whether [x] is the tightest interval of doubles enclosing v
  bool is_tightest_enclosure(const Interval& x, const Exact& v)
  {
    return is_below(x.lb(), v) && !is_below(next_double(x.lb()), v)
      && is_above(x.ub(), v) && !is_above(previous_double(x.ub()), v);
  }

  // Whether [x] is the tightest interval of doubles enclosing a/b
  bool is_tightest_quotient(const Interval& x, double a, double b)
  {
    return is_below_quotient(x.lb(), a, b) && !is_below_quotient(next_double(x.lb()), a, b)
      && is_above_quotient(x.ub(), a, b) && !is_above_quotient(previous_double(x.ub()), a, b);
  }

  // A deterministic generator of doubles (xorshift64): a random sign and
  // mantissa, and an exponent between emin and emax, which keeps the exact
  // computations above far from overflow and underflow.
  class RandomDoubles
  {
    public:

      double operator()(int emin, int emax)
      {
        const uint64_t mantissa = next_bits() & 0x000FFFFFFFFFFFFFull;
        const uint64_t exponent = static_cast<uint64_t>(1023 + emin)
          + next_bits() % static_cast<uint64_t>(emax - emin + 1);
        const uint64_t sign = next_bits() & 1u;
        const uint64_t bits = (sign << 63) | (exponent << 52) | mantissa;
        double x;
        std::memcpy(&x, &bits, sizeof x);
        return x;
      }

    private:

      uint64_t next_bits()
      {
        _state ^= _state << 13;
        _state ^= _state >> 7;
        _state ^= _state << 17;
        return _state;
      }

      uint64_t _state = 0x9E3779B97F4A7C15ull;
  };

  string hex(double x)
  {
    ostringstream s;
    s << hexfloat << x;
    return s.str();
  }

  string hex(const Interval& x)
  {
    return "[" + hex(x.lb()) + ", " + hex(x.ub()) + "]";
  }

  // Counts the failures of a check repeated on many values, and describes the
  // first one
  struct Failures
  {
    int count = 0;
    string first;

    void add(bool ok, const string& description)
    {
      if(!ok && count++ == 0)
        first = description;
    }
  };

  const int nb_random_values = 5000;
}

TEST_CASE("Interval rounding - arithmetic operations")
{
  RandomDoubles random;
  Failures sum, difference, product, product_by_double, quotient, square, square_root, interval_product;

  for(int i = 0; i < nb_random_values; i++)
  {
    const double a = random(-30,30), b = random(-30,30);
    const string values = "a=" + hex(a) + " b=" + hex(b);

    const Interval s = Interval(a) + Interval(b);
    sum.add(is_tightest_enclosure(s, exact_sum(a,b)), values + ": a+b=" + hex(s));

    const Interval d = Interval(a) - Interval(b);
    difference.add(is_tightest_enclosure(d, exact_sum(a,-b)), values + ": a-b=" + hex(d));

    const Interval p = Interval(a) * Interval(b);
    product.add(is_tightest_enclosure(p, exact_product(a,b)), values + ": a*b=" + hex(p));

    const Interval pd = Interval(a) * b;
    product_by_double.add(is_tightest_enclosure(pd, exact_product(a,b)), values + ": a*b=" + hex(pd));

    const Interval q = Interval(a) / Interval(b);
    quotient.add(is_tightest_quotient(q, a, b), values + ": a/b=" + hex(q));

    const Interval sq = sqr(Interval(a));
    square.add(is_tightest_enclosure(sq, exact_product(a,a)), values + ": a^2=" + hex(sq));

    // GAOL encloses square roots, not always as tightly as possible:
    // lb^2 <= |a| <= ub^2 only
    const double m = std::fabs(a);
    const Interval r = sqrt(Interval(m));
    square_root.add(r.lb() >= 0. && is_above(m, exact_product(r.lb(),r.lb())) && is_below(m, exact_product(r.ub(),r.ub())),
      "|a|=" + hex(m) + ": sqrt(|a|)=" + hex(r));

    // Products of intervals, whose bounds have all the signs GAOL distinguishes:
    // the tightest enclosure of the four products of bounds
    const double c = random(-30,30), e = random(-30,30);
    const Interval x(a < c ? a : c, a < c ? c : a), y(b < e ? b : e, b < e ? e : b);
    const Interval xy = x * y;
    const vector<Exact> corners = {
      exact_product(x.lb(),y.lb()), exact_product(x.lb(),y.ub()),
      exact_product(x.ub(),y.lb()), exact_product(x.ub(),y.ub())
    };
    bool lb_below_all = true, next_lb_above_one = false, ub_above_all = true, previous_ub_below_one = false;
    for(const auto& v : corners)
    {
      lb_below_all = lb_below_all && is_below(xy.lb(), v);
      next_lb_above_one = next_lb_above_one || !is_below(next_double(xy.lb()), v);
      ub_above_all = ub_above_all && is_above(xy.ub(), v);
      previous_ub_below_one = previous_ub_below_one || !is_above(previous_double(xy.ub()), v);
    }
    interval_product.add(lb_below_all && next_lb_above_one && ub_above_all && previous_ub_below_one,
      "x=" + hex(x) + " y=" + hex(y) + ": x*y=" + hex(xy));
  }

  { INFO("a+b, first failure: " << sum.first); CHECK(sum.count == 0); }
  { INFO("a-b, first failure: " << difference.first); CHECK(difference.count == 0); }
  { INFO("a*b, first failure: " << product.first); CHECK(product.count == 0); }
  { INFO("[a]*b, first failure: " << product_by_double.first); CHECK(product_by_double.count == 0); }
  { INFO("a/b, first failure: " << quotient.first); CHECK(quotient.count == 0); }
  { INFO("sqr(a), first failure: " << square.first); CHECK(square.count == 0); }
  { INFO("sqrt(|a|), first failure: " << square_root.first); CHECK(square_root.count == 0); }
  { INFO("[x]*[y], first failure: " << interval_product.first); CHECK(interval_product.count == 0); }
}

TEST_CASE("Interval rounding - elementary functions")
{
  // The doubles immediately below and above the exact value of each function at
  // x, computed with 400 bits of precision (mpmath): the enclosure computed by
  // Codac has to contain both.
  struct Value
  {
    const char* function;
    Interval (*f)(const Interval&);
    double x, below, above;
  };

  const auto f_exp = [](const Interval& x) { return exp(x); };
  const auto f_log = [](const Interval& x) { return log(x); };
  const auto f_sin = [](const Interval& x) { return sin(x); };
  const auto f_cos = [](const Interval& x) { return cos(x); };
  const auto f_tan = [](const Interval& x) { return tan(x); };
  const auto f_atan = [](const Interval& x) { return atan(x); };
  const auto f_sqrt = [](const Interval& x) { return sqrt(x); };

  const vector<Value> values = {
    { "exp", f_exp,  0x1.0000000000000p+0,  0x1.5bf0a8b145769p+1,  0x1.5bf0a8b14576ap+1 },
    { "exp", f_exp,  0x1.0000000000000p-1,  0x1.a61298e1e069bp+0,  0x1.a61298e1e069cp+0 },
    { "exp", f_exp, -0x1.0000000000000p+0,  0x1.78b56362cef37p-2,  0x1.78b56362cef38p-2 },
    { "exp", f_exp,  0x1.4000000000000p+3,  0x1.5829dcf95055fp+14, 0x1.5829dcf950560p+14 },
    { "exp", f_exp,  0x1.0624dd2f1a9fcp-10, 0x1.0041919b7ee33p+0,  0x1.0041919b7ee34p+0 },
    { "exp", f_exp, -0x1.4000000000000p+4,  0x1.1b48655f37266p-29, 0x1.1b48655f37267p-29 },
    { "log", f_log,  0x1.0000000000000p+1,  0x1.62e42fefa39efp-1,  0x1.62e42fefa39f0p-1 },
    { "log", f_log,  0x1.4000000000000p+3,  0x1.26bb1bbb55515p+1,  0x1.26bb1bbb55516p+1 },
    { "log", f_log,  0x1.999999999999ap-4, -0x1.26bb1bbb55516p+1, -0x1.26bb1bbb55515p+1 },
    { "log", f_log,  0x1.8000000000000p+0,  0x1.9f323ecbf984bp-2,  0x1.9f323ecbf984cp-2 },
    { "log", f_log,  0x1.4f8b588e368f1p-17,-0x1.7069e2aa2aa5bp+3, -0x1.7069e2aa2aa5ap+3 },
    { "log", f_log,  0x1.2a05f20000000p+33, 0x1.7069e2aa2aa5ap+4,  0x1.7069e2aa2aa5bp+4 },
    { "sin", f_sin,  0x1.0000000000000p+0,  0x1.aed548f090ceep-1,  0x1.aed548f090cefp-1 },
    { "sin", f_sin,  0x1.0000000000000p-1,  0x1.eaee8744b05efp-2,  0x1.eaee8744b05f0p-2 },
    { "sin", f_sin,  0x1.8000000000000p+1,  0x1.210386db6d55bp-3,  0x1.210386db6d55cp-3 },
    { "sin", f_sin, -0x1.0000000000000p+1, -0x1.d18f6ead1b446p-1, -0x1.d18f6ead1b445p-1 },
    { "sin", f_sin,  0x1.9000000000000p+6, -0x1.03425b78c4db9p-1, -0x1.03425b78c4db8p-1 },
    { "sin", f_sin,  0x1.0624dd2f1a9fcp-10, 0x1.0624da5218a62p-10, 0x1.0624da5218a63p-10 },
    { "cos", f_cos,  0x1.0000000000000p+0,  0x1.14a280fb5068bp-1,  0x1.14a280fb5068cp-1 },
    { "cos", f_cos,  0x1.0000000000000p-1,  0x1.c1528065b7d4fp-1,  0x1.c1528065b7d50p-1 },
    { "cos", f_cos,  0x1.8000000000000p+1, -0x1.fae04be85e5d3p-1, -0x1.fae04be85e5d2p-1 },
    { "cos", f_cos, -0x1.0000000000000p+1, -0x1.aa22657537205p-2, -0x1.aa22657537204p-2 },
    { "cos", f_cos,  0x1.9000000000000p+6,  0x1.b981dbf665fdfp-1,  0x1.b981dbf665fe0p-1 },
    { "cos", f_cos,  0x1.8000000000000p+0,  0x1.21bd54fc5f9a7p-4,  0x1.21bd54fc5f9a8p-4 },
    { "tan", f_tan,  0x1.0000000000000p+0,  0x1.8eb245cbee3a5p+0,  0x1.8eb245cbee3a6p+0 },
    { "tan", f_tan,  0x1.0000000000000p-1,  0x1.17b4f5bf3474ap-1,  0x1.17b4f5bf3474bp-1 },
    { "tan", f_tan, -0x1.3333333333333p+0, -0x1.493c43acb164dp+1, -0x1.493c43acb164cp+1 },
    { "tan", f_tan,  0x1.8000000000000p+1, -0x1.23ef71254b870p-3, -0x1.23ef71254b86fp-3 },
    { "atan", f_atan,  0x1.0000000000000p+0,  0x1.921fb54442d18p-1,  0x1.921fb54442d19p-1 },
    { "atan", f_atan,  0x1.0000000000000p-1,  0x1.dac670561bb4fp-2,  0x1.dac670561bb50p-2 },
    { "atan", f_atan, -0x1.0000000000000p+1, -0x1.1b6e192ebbe45p+0, -0x1.1b6e192ebbe44p+0 },
    { "atan", f_atan,  0x1.f400000000000p+9,  0x1.91de2c0e658bcp+0,  0x1.91de2c0e658bdp+0 },
    { "sqrt", f_sqrt,  0x1.0000000000000p+1,  0x1.6a09e667f3bccp+0,  0x1.6a09e667f3bcdp+0 },
    { "sqrt", f_sqrt,  0x1.8000000000000p+1,  0x1.bb67ae8584caap+0,  0x1.bb67ae8584cabp+0 },
    { "sqrt", f_sqrt,  0x1.999999999999ap-4,  0x1.43d136248490fp-2,  0x1.43d1362484910p-2 },
    { "sqrt", f_sqrt,  0x1.ad7f29abcaf48p-24, 0x1.4b96be9c2da2bp-12, 0x1.4b96be9c2da2cp-12 },
  };

  for(const auto& v : values)
  {
    const Interval y = v.f(Interval(v.x));
    INFO(v.function << "(" << hex(v.x) << ") = " << hex(y));
    CHECK(y.lb() <= v.below);
    CHECK(y.ub() >= v.above);
  }

  // pi, whose neighbours scaled by a power of two are those of pi/2 and 2pi
  const double pi_below = 0x1.921fb54442d18p+1, pi_above = 0x1.921fb54442d19p+1;
  CHECK((Interval::pi().lb() <= pi_below && Interval::pi().ub() >= pi_above));
  CHECK((Interval::half_pi().lb() <= pi_below/2. && Interval::half_pi().ub() >= pi_above/2.));
  CHECK((Interval::two_pi().lb() <= pi_below*2. && Interval::two_pi().ub() >= pi_above*2.));
}

TEST_CASE("Interval rounding - after any operation")
{
  // Codac and GAOL compute with the rounding direction set upward, once for
  // all (see the initialisation of oo in codac2_math.h): every operation has to
  // leave it so, or the operations after it return wrong bounds.
  CHECK(std::fegetround() == FE_UPWARD);

  struct Operation
  {
    const char* name;
    void (*run)(const Interval& x, const Interval& y);
  };

  // x is within [0,1] and y above 1, which is the domain of every function below
  const Interval x(0.1,0.3), y(1.5,2.5);

  const vector<Operation> operations = {
    { "Interval(a,b)", [](const Interval&, const Interval&) { (void)Interval(0.1,1./3.); } },
    { "Interval::mid()", [](const Interval& x, const Interval&) { (void)x.mid(); } },
    { "Interval::rad()", [](const Interval& x, const Interval&) { (void)x.rad(); } },
    { "Interval::diam()", [](const Interval& x, const Interval&) { (void)x.diam(); } },
    { "Interval::mag()", [](const Interval& x, const Interval&) { (void)x.mag(); } },
    { "Interval::mig()", [](const Interval& x, const Interval&) { (void)x.mig(); } },
    { "Interval::rand()", [](const Interval& x, const Interval&) { (void)x.rand(); } },
    { "Interval::bisect()", [](const Interval& x, const Interval&) { (void)x.bisect(); } },
    { "Interval::inflate()", [](const Interval& x, const Interval&) { Interval z(x); z.inflate(0.1); } },
    { "operator<<", [](const Interval& x, const Interval&) { ostringstream s; s << x; } },
    { "operator+=", [](const Interval& x, const Interval& y) { Interval z(x); z += y; } },
    { "operator-=", [](const Interval& x, const Interval& y) { Interval z(x); z -= y; } },
    { "operator*=", [](const Interval& x, const Interval& y) { Interval z(x); z *= y; } },
    { "operator/=", [](const Interval& x, const Interval& y) { Interval z(x); z /= y; } },
    { "operator|=", [](const Interval& x, const Interval& y) { Interval z(x); z |= y; } },
    { "operator&=", [](const Interval& x, const Interval& y) { Interval z(x); z &= y; } },
    { "sqr", [](const Interval& x, const Interval&) { (void)sqr(x); } },
    { "sqrt", [](const Interval& x, const Interval&) { (void)sqrt(x); } },
    { "pow(x,int)", [](const Interval&, const Interval& y) { (void)pow(y,3); } },
    { "pow(x,double)", [](const Interval&, const Interval& y) { (void)pow(y,2.5); } },
    { "pow(x,Interval)", [](const Interval& x, const Interval& y) { (void)pow(y,x); } },
    { "root", [](const Interval&, const Interval& y) { (void)root(y,3); } },
    { "exp", [](const Interval& x, const Interval&) { (void)exp(x); } },
    { "log", [](const Interval& x, const Interval&) { (void)log(x); } },
    { "cos", [](const Interval& x, const Interval&) { (void)cos(x); } },
    { "sin", [](const Interval& x, const Interval&) { (void)sin(x); } },
    { "tan", [](const Interval& x, const Interval&) { (void)tan(x); } },
    { "acos", [](const Interval& x, const Interval&) { (void)acos(x); } },
    { "asin", [](const Interval& x, const Interval&) { (void)asin(x); } },
    { "atan", [](const Interval& x, const Interval&) { (void)atan(x); } },
    { "atan2", [](const Interval& x, const Interval& y) { (void)atan2(y,x); } },
    { "cosh", [](const Interval& x, const Interval&) { (void)cosh(x); } },
    { "sinh", [](const Interval& x, const Interval&) { (void)sinh(x); } },
    { "tanh", [](const Interval& x, const Interval&) { (void)tanh(x); } },
    { "acosh", [](const Interval&, const Interval& y) { (void)acosh(y); } },
    { "asinh", [](const Interval& x, const Interval&) { (void)asinh(x); } },
    { "atanh", [](const Interval& x, const Interval&) { (void)atanh(x); } },
    { "abs", [](const Interval& x, const Interval&) { (void)abs(x); } },
    { "min", [](const Interval& x, const Interval& y) { (void)min(x,y); } },
    { "max", [](const Interval& x, const Interval& y) { (void)max(x,y); } },
    { "sign", [](const Interval& x, const Interval&) { (void)sign(x); } },
    { "integer", [](const Interval&, const Interval& y) { (void)integer(y); } },
    { "floor", [](const Interval&, const Interval& y) { (void)floor(y); } },
    { "ceil", [](const Interval&, const Interval& y) { (void)ceil(y); } },
    { "chi", [](const Interval& x, const Interval& y) { (void)chi(x,y,x); } },
    { "Interval::pi()", [](const Interval&, const Interval&) { (void)Interval::pi(); } },
  };

  RandomDoubles random;
  for(const auto& op : operations)
  {
    op.run(x,y);
    INFO("after " << op.name);
    CHECK(std::fegetround() == FE_UPWARD);

    // A product and a sum whose exact results are not doubles, as the random
    // mantissas make almost certain: their bounds are the tightest ones only if
    // the rounding direction is still upward
    const double a = random(-30,30), b = random(-30,30);
    const Interval p = Interval(a) * Interval(b);
    INFO("a=" << hex(a) << " b=" << hex(b) << ": a*b=" << hex(p));
    CHECK(is_tightest_enclosure(p, exact_product(a,b)));
    const Interval s = Interval(a) + Interval(b);
    INFO("a+b=" << hex(s));
    CHECK(is_tightest_enclosure(s, exact_sum(a,b)));
  }
}
