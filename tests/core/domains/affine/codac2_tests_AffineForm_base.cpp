/* ============================================================================
 * I B E X - Test of Affine operations
 * ============================================================================
 * Copyright   : ENSTA Bretagne (FRANCE)
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file LICENSE
 *
 * Author(s)   : Jordan NININ
 * Created     : Juin 20, 2021
 * ----------------------------------------------------------------------------
 *
 * \brief Base tests for AffineMain<T>: construction, copy/assignment
 *        semantics, static_assert traits, comparisons and conversions
 *        to/from Interval, the Eigen comma-initializer, the
 *        Approx<AffineMain<T>> matcher, and the low-level AF_fAF2
 *        representation (coefficient storage, compact(), twoProd,
 *        ASan/UBSan robustness).
 * ---------------------------------------------------------------------------- */

#include <catch2/catch_test_macros.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <type_traits>
#include <vector>
#include "codac2_Affine.h"
#include "codac2_AffineMatrix.h"
#include "codac2_AffineRow.h"
#include "codac2_AffineVector.h"
#include "codac2_IntervalMatrix.h"
#include "codac2_IntervalRow.h"
#include "codac2_Approx.h"
#include "codac2_Matrix.h"

using namespace codac2;

using AA = AF_Default;
using AffineT = AffineMain<AA>;
using AffineTMatrix = AffineMainMatrix<AA>;
using AffineTVector = AffineMainVector<AA>;
using AffineTVarVector = AffineVarMainVector<AA>;
using AffineTRow = AffineMainRow<AA>;

using FixedAffineTMatrix = Eigen::Matrix<AffineT, 2, 3>;
using FixedAffineTVector = Eigen::Matrix<AffineT, 3, 1>;
using FixedAffineTRow    = Eigen::Matrix<AffineT, 1, 3>;

const double ERROR = std::numeric_limits<double>::epsilon()*10;

template<class T>
void CHECK_change_mode_MinRange()
{
  AffineMain<T>::change_mode(AffineMain<T>::AF_Lin_MinRange);
  CHECK(AffineMain<T>::AF_Lin_MinRange==AffineMain<T>::get_mode());

  AffineMain<T>::change_mode(AffineMain<T>::AF_Lin_Chebyshev);
  CHECK(AffineMain<T>::AF_Lin_Chebyshev==AffineMain<T>::get_mode());

}

template<class T>
void CHECK_affine_inclu(const AffineMain<T>& y_actual, const Interval& y_expected)
{
  if (y_expected.is_empty()) { CHECK(y_actual.is_empty()); return; }

  CHECK_FALSE(y_actual.is_empty());
  CHECK(y_expected.lb()>=y_actual.itv().lb());
  CHECK(y_expected.ub()<=y_actual.itv().ub());
}

namespace {

void check_point(const AffineT& actual, double expected)
{
    CAPTURE(actual.itv(), expected);
    CHECK(actual == Approx<AffineT>(Interval(expected)));
}

void check_interval_enclosure(const AffineT& actual, const Interval& expected)
{
    CAPTURE(actual.itv(), expected);
    CHECK(actual.itv().is_superset(expected));
}

} // namespace

TEST_CASE("Approx<AffineMain<T>>")
{

  // empty
  {
    AffineTVarVector ax(1);
    ax[0] = Interval::empty();
    CHECK(ax[0] == Approx<AffineT>(Interval::empty()));
  }

  // unbounded
  {
    AffineTVarVector ax(1);
    CHECK(ax[0] == Approx<AffineT>(Interval()));
  }
  {

	  AffineTVarVector ax(IntervalVector({{0,1},{0,1},{0,1}}));
	  CHECK(ax[0] == Approx<AffineT>(Interval(0,1)));
	  CHECK(ax[1] == Approx<AffineT>(Interval(0,1)));
	  CHECK(ax[2] == Approx<AffineT>(Interval(0,1)));
  }

  // simple linear cases: x in [0,1], 1 noise var
  {
    AffineTVarVector ax(IntervalVector({{0,1}}));
    CHECK(ax[0] == Approx<AffineT>(Interval(0,1)));
    CHECK((ax[0]+1.0) == Approx<AffineT>(Interval(1,2)));
    CHECK((-ax[0]) == Approx<AffineT>(Interval(-1,0)));
    CHECK((2.0*ax[0]) == Approx<AffineT>(Interval(0,2)));
  }

  // two affine variables summed together (2 noise variables)
  {
	  AffineTVarVector ax(IntervalVector({{0,1},{0,1}}));
	  AffineT y = ax[0] + ax[1];
	  CHECK(y == Approx<AffineT>(Interval(0,2)));
  }

  // must not match a wrong expected interval
  {
    AffineTVarVector bx(2);
    bx[0] =Interval(0,1);
    CHECK_FALSE(bx[0] == Approx<AffineT>(Interval(5,6)));
    bx[0] =Interval(5,6);
    CHECK(bx[0] == Approx<AffineT>(Interval(5,6)));
  }
}




TEST_CASE("AffineForm operations")
{

  CHECK_change_mode_MinRange<AA>();
  {
	  static_assert(
	    std::is_constructible_v<
	      AffineT,
	      double>
	  );

	  static_assert(
	    std::is_convertible_v<
	      double,
	      AffineT>
	  );

	  static_assert(
	    std::is_copy_constructible_v<
	      AffineVarMain<AA>>
	  );

	  static_assert(
	    std::is_copy_assignable_v<
	      AffineVarMain<AA>>
	  );

	  static_assert(
	    !std::is_assignable_v<
	      AffineVarMain<AA>&,
	      AffineT>);
  }
  /* test: operator-() */
  {
    AffineTVarVector ax(1);
    ax[0] =  Interval(0,1);
    CHECK(-(ax[0]) == Approx<AffineT>(Interval(-1,0), ERROR));
  }
  {
    AffineTVarVector ax(1);
    CHECK(-(ax[0]) == Approx<AffineT>(Interval(), ERROR));
  }
  {
    AffineTVarVector ax(1);
    ax[0] = Interval(-oo,0);
    CHECK(-(ax[0]) == Approx<AffineT>(Interval(0,oo), ERROR));
  }
  {
    AffineTVarVector ax(4);
    ax[0] = Interval(-oo,1);
    CHECK(-(ax[0]) == Approx<AffineT>(Interval(-1,oo), ERROR));
    CHECK(ax[2] == Approx<AffineT>(Interval(), ERROR));
  }
  {
    AffineT x(3.);
    AffineT y(x);

    CHECK(y.itv() == Interval(3.));

    AffineT z;
    z = x;

    CHECK(z.itv() == Interval(3.));
    CHECK((-z).itv() == Interval(-3.));
  }
  {

    AffineTVarVector x(IntervalVector({ {1.,2.}, {3.,4.}, {5.,6.} }) );

    AffineT source = x[0];
    source.set_empty();

    AffineT destination = x[1];
    destination = source;

    CHECK(destination.is_empty());
    CHECK(destination.noise_count() == source.noise_count());
  }
  {
    CHECK(std::is_nothrow_move_constructible_v<AF_fAF2> );
    CHECK(std::is_nothrow_move_assignable_v<AF_fAF2> );
    CHECK_FALSE(std::is_copy_constructible_v<AF_fAF2>);
    CHECK_FALSE(std::is_copy_assignable_v<AF_fAF2>);
  }
  {
    AffineTVarVector x(2);
    x[0] = 5.0;
    CHECK(x[0].itv() == Interval(5.));
  }
  {
    AffineTVarVector variables(IntervalVector({{1., 2.}}) );

    AffineT x = variables[0];
    AffineT residual = floor(x) - x;

    CHECK(residual.itv() != Interval(0.));
    CHECK(Interval(-1., 0.).is_subset(residual.itv()));
  }
}


// ============================================================================
// Comma-initializer with Interval literals on AffineT-scalar containers
// ============================================================================
//
// These are the exact call sites that used to fail with "invalid operands
// to binary expression" before AffineTMain(const Interval&) was added:
// Eigen's comma-initializer needs a converting constructor for each
// inserted value, operator=(const Interval&) alone isn't enough.

TEST_CASE("Eigen AffineT comma-initializer with Interval literals", "[Eigen][comma-initializer][AffineT]")
{
    AffineTVector vector(3);
    vector << Interval(1.0), Interval(-2.0), Interval(4.0);
    check_point(vector[0], 1.0);
    check_point(vector[1], -2.0);
    check_point(vector[2], 4.0);

    AffineTRow row(3);
    row << Interval(1.0), Interval(-2.0), Interval(4.0);
    check_point(row[0], 1.0);
    check_point(row[1], -2.0);
    check_point(row[2], 4.0);

    AffineTMatrix matrix(2, 3);
    matrix << Interval(1.0), Interval(-2.0), Interval(3.0),
              Interval(4.0), Interval(5.0), Interval(-6.0);
    check_point(matrix(0, 0), 1.0);
    check_point(matrix(0, 1), -2.0);
    check_point(matrix(1, 2), -6.0);

    FixedAffineTVector fvector;
    fvector << Interval(1.0), Interval(-2.0), Interval(4.0);
    check_point(fvector[0], 1.0);
    check_point(fvector[1], -2.0);
    check_point(fvector[2], 4.0);

    FixedAffineTRow frow;
    frow << Interval(1.0), Interval(-2.0), Interval(4.0);
    check_point(frow[0], 1.0);
    check_point(frow[1], -2.0);
    check_point(frow[2], 4.0);

    FixedAffineTMatrix fmatrix;
    fmatrix << Interval(1.0), Interval(-2.0), Interval(3.0),
               Interval(4.0), Interval(5.0), Interval(-6.0);
    check_point(fmatrix(0, 0), 1.0);
    check_point(fmatrix(0, 1), -2.0);
    check_point(fmatrix(1, 2), -6.0);
}


// ============================================================================
// 5. AffineTMain <-> Interval conversions, and operator==/!=(double)
// ============================================================================



TEST_CASE("Eigen AffineT equality with double", "[AffineT][conversion]")
{
    const AffineT a(4.0);

    CHECK(a == 4.0);
    CHECK_FALSE(a == 5.0);
    CHECK(a != 5.0);
    CHECK_FALSE(a != 4.0);
}



TEST_CASE("Eigen AffineT implicit conversion to and from Interval", "[AffineT][conversion]")
{
    const AffineT a(4.0);

    const Interval as_interval = a; // implicit operator Interval()
    CHECK(as_interval == Interval(4.0));

    const AffineT from_interval = Interval(2.0, 3.0); // implicit AffineTMain(const Interval&)
    check_interval_enclosure(from_interval, Interval(2.0, 3.0));

    // The two conversions being bidirectional and implicit is exactly the
    // kind of thing that can silently become ambiguous: pin down that
    // plain arithmetic between an AffineT and an Interval still resolves
    // cleanly instead of failing to compile.
    const Interval itv(1.0, 2.0);
    const AffineT sum = a + itv;
    check_interval_enclosure(sum, Interval(4.0) + itv);
}


// ============================================================================
// 6. AffineTVarVector::init(const Interval&): each component keeps its own
//    dedicated noise symbol after a broadcast initialization
// ============================================================================
//
// Regression test for the ex_affineform.cpp idiom `x.init(Interval(...))`.
// The generic Eigen::Matrix::init(const Scalar&) cannot be used here
// (Interval -> AffineVarMain<T> is intentionally not constructible, see
// test 3 above); AffineTVarVector needs its own init(const Interval&)
// overload assigning component-wise while preserving each component's
// dedicated noise symbol.



/*
 * Low-level AF_fAF2 robustness tests.
 *
 * The first two tests are especially useful with AddressSanitizer and
 * UndefinedBehaviorSanitizer. Without a sanitizer, an out-of-bounds access
 * or a null dereference may be silent, allocator-dependent, or may crash the
 * test process before Catch2 can report an assertion.
 */

TEST_CASE("AF_fAF2 inactive assignment keeps coefficient storage consistent")
{
  // destination initially owns storage for one affine variable
  AffineTVarVector small(IntervalVector({{0.0, 1.0}}));
  AffineT destination = small[0];

  // source is inactive, has a larger logical dimension, and normally owns no
  // coefficient array. Copy-assignment must not leave destination._n equal to
  // 16 while retaining its former two-element allocation.
  AffineTVarVector large_inactive(16);
  REQUIRE_FALSE(large_inactive[0].is_active());

  destination = large_inactive[0];
  CHECK(destination.noise_count() == 16);
  CHECK(destination.is_unbounded());

  // Re-activating the destination writes all coefficients from 0 through
  // noise_count(). ASan must not report a heap-buffer-overflow here.
  destination = Interval(-2.0, 3.0);

  CHECK(destination.noise_count() == 16);
  CHECK(destination.is_active());
  CHECK_affine_inclu<AA>(destination, Interval(-2.0, 3.0));

  for (Index i = 0; i < destination.noise_count(); ++i) {
    CAPTURE(i, destination.noise(i));
    CHECK(destination.noise(i) == 0.0);
  }
}


TEST_CASE("AF_fAF2 compact is safe for every inactive status")
{
  AffineTVarVector variables(12);
  AffineT value = variables[0];

  REQUIRE_FALSE(value.is_active());
  CHECK_NOTHROW(value.compact());

  value = Interval::empty();
  REQUIRE(value.is_empty());
  CHECK_NOTHROW(value.compact());

  value = Interval();
  REQUIRE(value.is_unbounded());
  CHECK_NOTHROW(value.compact());

  value = Interval(-oo, 2.0);
  REQUIRE(value.is_unbounded());
  CHECK_NOTHROW(value.compact());

  value = Interval(-2.0, oo);
  REQUIRE(value.is_unbounded());
  CHECK_NOTHROW(value.compact());
}


TEST_CASE("AF_fAF2 dimension changes preserve coefficients and inclusion")
{
  // Public arithmetic currently exposes dimension growth, not contraction.
  // This test exercises the observable grow path used by operator+= and
  // operator*=. A future public contraction path must transfer every removed
  // coefficient magnitude into the remainder instead of dropping it.
  AffineTVarVector small(IntervalVector({{1.0, 2.0}}));
  AffineTVarVector large(
    IntervalVector({
      {-1.0, 1.0}, {-2.0, 2.0}, {-3.0, 3.0}, {-4.0, 4.0},
      {-5.0, 5.0}, {-6.0, 6.0}, {-7.0, 7.0}, {-8.0, 8.0}
    })
  );

  AffineT sum = small[0];
  sum += large[7];
  CHECK(sum.noise_count() == 8);
  CHECK_affine_inclu<AA>(sum, small[0].itv() + large[7].itv());

  AffineT product = small[0];
  product *= large[7];
  CHECK(product.noise_count() == 8);
  CHECK_affine_inclu<AA>(product, small[0].itv() * large[7].itv());
}


TEST_CASE("AF_fAF2 multiplication handles opposite finite scales")
{
  const std::vector<std::pair<Interval, Interval>> inputs = {
    {Interval(1.e299, 1.e300), Interval(1.e-300, 1.e-299)},
    {Interval(-1.e300, -1.e299), Interval(1.e-300, 1.e-299)},
    {Interval(1.e200, 2.e200), Interval(1.e-200, 2.e-200)},
    {Interval(-1.e200, 2.e200), Interval(-2.e-200, 1.e-200)}
  };

  for (const auto& input : inputs) {
    CAPTURE(input.first, input.second);
    AffineTVarVector variables(
      IntervalVector({input.first, input.second})
    );

    const Interval reference = input.first * input.second;
    const AffineT result = variables[0] * variables[1];

    CHECK_FALSE(result.is_empty());
    CHECK_affine_inclu<AA>(result, reference);

    if (result.is_active()) {
      CHECK(std::isfinite(result.err()));
      CHECK(result.err() >= 0.0);
    }
  }
}


TEST_CASE("AF_fAF2 scalar products and sums handle extreme cancellation")
{
  const double large = 1.e300;
  const double small = 1.e-300;

  AffineTVarVector variables(
    IntervalVector({{large, codac2::next_float(large)}})
  );

  const AffineT scaled = variables[0] * small;
  CHECK_FALSE(scaled.is_empty());
  CHECK_affine_inclu<AA>(scaled, variables[0].itv() * small);

  const AffineT cancelled = variables[0] + (-large);
  CHECK_FALSE(cancelled.is_empty());
  CHECK_affine_inclu<AA>(cancelled, variables[0].itv() - large);

  if (scaled.is_active()) {
    CHECK(std::isfinite(scaled.err()));
    CHECK(scaled.err() >= 0.0);
  }
  if (cancelled.is_active()) {
    CHECK(std::isfinite(cancelled.err()));
    CHECK(cancelled.err() >= 0.0);
  }
}


TEST_CASE("AF_fAF2 active remainder stays finite and non-negative")
{
  AffineTVarVector variables(
    IntervalVector({{-2.0, 3.0}, {0.5, 2.0}, {-1.0, 1.0}})
  );

  const std::vector<AffineT> results = {
    variables[0] + variables[1],
    variables[0] * variables[1],
    sqr(variables[0]),
    exp(variables[2]),
    sin(variables[0]),
    asinh(variables[0]),
    inv(variables[1])
  };

  for (const AffineT& result : results) {
    CAPTURE(result.itv());
    if (result.is_active()) {
      CHECK(std::isfinite(result.err()));
      CHECK(result.err() >= 0.0);
    }
  }
}


TEST_CASE("AF_fAF2 repeated multiplication and status transitions remain stable")
{
  // This is a stress/regression test for temporary storage and repeated
  // active/inactive transitions. LeakSanitizer is needed to diagnose an
  // actual leak; the assertions check observable stability and inclusion.
  for (int iteration = 0; iteration < 256; ++iteration) {
    CAPTURE(iteration);
    AffineTVarVector variables(
      IntervalVector({{-1.0, 2.0}, {2.0, 3.0}, {-0.5, 0.5}})
    );

    AffineT value = variables[0];
    value *= variables[1];
    value += variables[2];
    CHECK_affine_inclu<AA>(
      value,
      variables[0].itv() * variables[1].itv() + variables[2].itv()
    );

    value = Interval();
    CHECK(value.is_unbounded());
    value = Interval(-1.0, 1.0);
    CHECK_affine_inclu<AA>(value, Interval(-1.0, 1.0));
    value.compact();
  }
}


TEST_CASE("AF_fAF2 interval extraction matches all status transitions")
{
  AffineT value;

  value = Interval::empty();
  CHECK(value.itv().is_empty());

  value = Interval();
  CHECK(value.itv() == Interval());

  value = Interval(-oo, 4.0);
  CHECK(value.itv() == Interval(-oo, 4.0));

  value = Interval(-3.0, oo);
  CHECK(value.itv() == Interval(-3.0, oo));

  value = Interval(2.0);
  CHECK(value.itv() == Interval(2.0));

  value = Interval(-2.0, 5.0);
  CHECK_affine_inclu<AA>(value, Interval(-2.0, 5.0));
}



TEST_CASE(  "AF_fAF2 square preserves small coefficients amplified by a large remainder")
{
  struct TestData {
    Interval variable;
    Interval uncertainty;
  };

  const std::vector<TestData> cases = {
    {
      Interval(-1.e-20, 1.e-20),
      Interval(-1.e20, 1.e20)
    },
    {
      Interval(-1.e-30, 1.e-30),
      Interval(-1.e30, 1.e30)
    },
    {
      Interval(-1.e-50, 1.e-50),
      Interval(-1.e50, 1.e50)
    }
  };

  for (const TestData& data : cases) {
    CAPTURE(data.variable, data.uncertainty);

    AffineTVarVector variables(
      IntervalVector({data.variable})
    );

    AffineT value = variables[0];
    value += data.uncertainty;

    const Interval before = value.itv();
    const AffineT result = sqr(value);
    const Interval reference = sqr(before);

    CHECK_FALSE(result.is_empty());
    CHECK_affine_inclu<AA>(result, reference);

    if (result.is_active()) {
      CHECK(std::isfinite(result.err()));
      CHECK(result.err() >= 0.0);
    }
  }
}

TEST_CASE(  "AF_fAF2 normalized square preserves mixed-scale contributions")
{
  const std::vector<double> scales = {
    1.e20,
    1.e30,
    1.e50,
    1.e100
  };

  for (const double scale : scales) {
    CAPTURE(scale);

    const double inverse_scale = 1.0/scale;

    AffineTVarVector variables(
      IntervalVector({
        {-inverse_scale, inverse_scale}
      })
    );

    AffineT value = variables[0];
    value += Interval(-scale, scale);

    const Interval reference =
      sqr(value.itv()) / (scale*scale);

    const AffineT normalized =
      sqr(value) / (scale*scale);

    CHECK_FALSE(normalized.is_empty());
    CHECK_affine_inclu<AA>(
      normalized,
      reference
    );

    if (normalized.is_active()) {
      CHECK(std::isfinite(normalized.err()));
      CHECK(normalized.err() >= 0.0);
    }
  }
}


TEST_CASE(  "AF_fAF2 scalar multiplication transfers discarded coefficients to remainder")
{
  AffineTVarVector variables(
    IntervalVector({
      {-1.e-20, 1.e-20}
    })
  );

  const AffineT result =
    variables[0] * 0.5;

  CHECK_affine_inclu<AA>(
    result,
    variables[0].itv() * 0.5
  );
}

TEST_CASE(  "AF_fAF2 addition transfers cancelled coefficients to remainder")
{
  AffineTVarVector variables(
    IntervalVector({
      {-1.e-20, 1.e-20}
    })
  );

  const AffineT x = variables[0];

  // Introduces a nearly complete cancellation while retaining
  // a representable small residual coefficient.
  const AffineT result =
    x + (-x * codac2::prev_float(1.0));

  const Interval reference =
    x.itv() +
    (-x.itv() * codac2::prev_float(1.0));

  CHECK_affine_inclu<AA>(
    result,
    reference
  );
}

TEST_CASE(  "AF_fAF2 compact transfers removed coefficients to remainder")
{
  AffineTVarVector variables(
    IntervalVector({
      {-1.e-8, 1.e-8}
    })
  );

  AffineT value = variables[0];
  const Interval before = value.itv();

  value.compact(1.e-6);

  CHECK_affine_inclu<AA>(value, before);

  if (value.is_active()) {
    CHECK(value.err() >= 1.e-8);
  }
}

TEST_CASE(  "AF_fAF2 software twoProd handles opposite scales")
{
if (!codac_fma_runtime()) {

  SKIP(
    "This test requires FMA option available on your computer and enabled at runtime"
  );
} else {

  const std::vector< std::pair<Interval, Interval> > cases = {
    {
      Interval(1.e299, 1.e300),
      Interval(1.e-300, 1.e-299)
    },
    {
      Interval(-1.e300, -1.e299),
      Interval(1.e-300, 1.e-299)
    },
    {
      Interval(1.e300),
      Interval(1.e-300)
    },
    {
      Interval(-1.e300),
      Interval(1.e-300)
    },
    {
      Interval(1.e200, 2.e200),
      Interval(1.e-200, 2.e-200)
    }
  };

  for (const auto& data : cases) {
    CAPTURE(data.first, data.second);

    AffineTVarVector variables(
      IntervalVector({
        data.first,
        data.second
      })
    );

    const Interval reference =
      data.first * data.second;

    const AffineT result =
      variables[0] * variables[1];

    CHECK_FALSE(result.is_empty());

    CHECK_affine_inclu<AA>(
      result,
      reference
    );

    if (result.is_active()) {
      CHECK(std::isfinite(result.err()));
      CHECK(result.err() >= 0.0);
    }
  }

}
}




TEST_CASE(
  "AF_fAF2 square handles a small coefficient with a large remainder"
)
{
  const double small =
    0.5 * 0x1p-55; // inférieur à AF_EC

  const std::vector<double> remainders = {
    1.0,
    2.0,
    10.0,
    1.e4,
    1.e8,
    1.e16
  };

  for (const double remainder : remainders) {
    CAPTURE(small, remainder);

    AffineTVarVector variables(
      IntervalVector({
        {-small, small}
      })
    );

    AffineT value = variables[0];

    // Le petit coefficient reste dans _val[1].
    // L'incertitude ajoutée est placée dans _err.
    value += Interval(-remainder, remainder);

    REQUIRE(value.is_active());
    REQUIRE(value.noise_count() == 1);

    CHECK(std::fabs(value.noise(0)) > 0.0);
    CHECK(std::fabs(value.noise(0)) < 0x1p-55);
    CHECK(value.err() >= remainder);

    const Interval input_enclosure = value.itv();
    const Interval expected = sqr(input_enclosure);

    const AffineT result = sqr(value);

    CAPTURE(
      input_enclosure,
      expected,
      result.itv()
    );

    CHECK_FALSE(result.is_empty());
    CHECK_affine_inclu<AA>(result, expected);

    if (result.is_active()) {
      CHECK(std::isfinite(result.err()));
      CHECK(result.err() >= 0.0);
    }
  }
}

TEST_CASE(  "AF_fAF2 square handles several small coefficients with a large remainder")
{
  const double small = 0.25 * 0x1p-55;

  AffineTVarVector variables(
    IntervalVector({
      {-small, small},
      {-small, small},
      {-small, small},
      {-small, small}
    })
  );

  AffineT value =
    variables[0] +
    variables[1] +
    variables[2] +
    variables[3];

  value += Interval(-1.e8, 1.e8);

  REQUIRE(value.is_active());

  const Interval before = value.itv();
  const Interval expected = sqr(before);
  const AffineT result = sqr(value);

  CAPTURE(
    before,
    expected,
    result.itv()
  );

  CHECK_FALSE(result.is_empty());
  CHECK_affine_inclu<AA>(result, expected);

  if (result.is_active()) {
    CHECK(std::isfinite(result.err()));
    CHECK(result.err() >= 0.0);
  }
}


TEST_CASE(
  "AF_fAF2 normalized square remains inclusive for mixed scales"
)
{
  struct TestData {
    double small;
    double remainder;
  };

  const std::vector<TestData> cases = {
    {0.5 * 0x1p-55, 2.0},
    {0.5 * 0x1p-55, 10.0},
    {0.25 * 0x1p-55, 1.e4},
    {0.125 * 0x1p-55, 1.e8}
  };

  for (const TestData& data : cases) {
    CAPTURE(data.small, data.remainder);

    AffineTVarVector variables(
      IntervalVector({
        {-data.small, data.small}
      })
    );

    AffineT value = variables[0];
    value += Interval(
      -data.remainder,
      data.remainder
    );

    const double normalization =
      data.remainder * data.remainder;

    REQUIRE(std::isfinite(normalization));
    REQUIRE(normalization > 0.0);

    const Interval expected =
      sqr(value.itv()) / normalization;

    const AffineT result =
      sqr(value) / normalization;

    CAPTURE(
      value.itv(),
      expected,
      result.itv()
    );

    CHECK_FALSE(result.is_empty());
    CHECK_affine_inclu<AA>(result, expected);

    if (result.is_active()) {
      CHECK(std::isfinite(result.err()));
      CHECK(result.err() >= 0.0);
    }
  }
}


TEST_CASE("AffineMain empty() factory matches an Interval::empty() assignment")
{
  const AffineT from_factory = AffineT::empty();
  AffineT from_assignment = Interval::empty();

  CHECK(from_factory.is_empty());
  CHECK(from_assignment.is_empty());
  CHECK(from_factory.itv().is_empty());
  CHECK(from_factory.noise_count() == from_assignment.noise_count());
  CHECK(from_factory == Approx<AffineT>(Interval::empty()));
}


TEST_CASE("AffineForm strict comparison operators return the documented BoolInterval")
{
  AffineTVarVector ax(2);
  ax[0] = Interval(0, 1);
  ax[1] = Interval(2, 3);

  // ub(ax[0]) < lb(ax[1]): certainly true, both ways and for every
  // Interval/AffineMain combination of the operands.
  CHECK((ax[0] < ax[1]) == BoolInterval::TRUE);
  CHECK((ax[1] > ax[0]) == BoolInterval::TRUE);
  CHECK((ax[0] < Interval(2, 3)) == BoolInterval::TRUE);
  CHECK((Interval(2, 3) > ax[0]) == BoolInterval::TRUE);

  // The reverse comparisons are certainly false.
  CHECK((ax[1] < ax[0]) == BoolInterval::FALSE);
  CHECK((ax[0] > ax[1]) == BoolInterval::FALSE);
  CHECK((Interval(2, 3) < ax[0]) == BoolInterval::FALSE);

  // Overlapping ranges: the sign is undetermined.
  AffineTVarVector bx(2);
  bx[0] = Interval(0, 2);
  bx[1] = Interval(1, 3);
  CHECK((bx[0] < bx[1]) == BoolInterval::UNKNOWN);
  CHECK((bx[0] > bx[1]) == BoolInterval::UNKNOWN);

  // An empty operand makes the comparison empty.
  AffineTVarVector ex(1);
  ex[0] = Interval::empty();
  CHECK((ex[0] < ax[1]) == BoolInterval::EMPTY);
  CHECK((ex[0] > ax[1]) == BoolInterval::EMPTY);
  CHECK((ex[0] < Interval(2, 3)) == BoolInterval::EMPTY);
  CHECK((Interval(2, 3) < ex[0]) == BoolInterval::EMPTY);
}


TEST_CASE("AffineForm set-relation predicates delegate to the interval enclosure")
{
  AffineTVarVector ax(4);
  ax[0] = Interval(1, 2);   // subset candidate
  ax[1] = Interval(0, 3);   // superset candidate
  ax[2] = Interval(5, 6);   // disjoint from ax[0]/ax[1]
  ax[3] = Interval(2, 4);   // overlaps ax[1] without being a subset

  // is_subset / is_superset, against both an Interval and an AffineForm.
  CHECK(ax[0].is_subset(Interval(0, 3)));
  CHECK(ax[0].is_subset(ax[1]));
  CHECK_FALSE(ax[1].is_subset(ax[0]));
  CHECK(ax[1].is_superset(ax[0]));
  CHECK(ax[1].is_superset(Interval(1, 2)));
  CHECK_FALSE(ax[0].is_superset(ax[1]));

  // Strict variants fail on touching boundaries, succeed strictly inside.
  CHECK(ax[0].is_strict_subset(ax[1]));
  CHECK_FALSE(ax[1].is_strict_subset(ax[1]));
  CHECK(ax[1].is_strict_superset(ax[0]));
  CHECK_FALSE(ax[1].is_strict_superset(ax[1]));

  // NOTE: is_relative_interior_subset() was removed from AffineMain<T>
  // (it forwarded to an Interval method that does not exist and could
  // never be instantiated); there is nothing left to test here.

  // Interior variant: ax[0]=[1,2] lies in the interior of [0,3], but its
  // own lower bound touches the boundary of [1,3].
  CHECK(ax[0].is_interior_subset(Interval(0, 3)));
  CHECK_FALSE(ax[0].is_interior_subset(Interval(1, 3)));
  CHECK(ax[0].is_strict_interior_subset(Interval(0, 3)));

  // contains / interior_contains.
  CHECK(ax[0].contains(1.0));
  CHECK(ax[0].contains(2.0));
  CHECK_FALSE(ax[0].interior_contains(1.0));
  CHECK(ax[0].interior_contains(1.5));

  // intersects / overlaps / is_disjoint, on overlapping and disjoint pairs.
  CHECK(ax[0].intersects(ax[1]));
  CHECK(ax[0].overlaps(ax[1]));
  CHECK_FALSE(ax[0].is_disjoint(ax[1]));

  CHECK_FALSE(ax[0].intersects(ax[2]));
  CHECK_FALSE(ax[0].overlaps(ax[2]));
  CHECK(ax[0].is_disjoint(ax[2]));

  // A boundary-touching pair intersects (non-empty intersection) but does
  // not overlap (zero-volume intersection), and is not disjoint.
  AffineTVarVector touching(2);
  touching[0] = Interval(0, 2);
  touching[1] = Interval(2, 4);
  CHECK(touching[0].intersects(touching[1]));
  CHECK_FALSE(touching[0].overlaps(touching[1]));
  CHECK_FALSE(touching[0].is_disjoint(touching[1]));

  // Every predicate must agree exactly with the equivalent Interval query
  // on the enclosure, since AffineMain simply forwards to itv().
  CHECK(ax[3].is_subset(ax[1]) == ax[3].itv().is_subset(ax[1].itv()));
  CHECK(ax[3].overlaps(ax[1]) == ax[3].itv().overlaps(ax[1].itv()));
}


TEST_CASE("AffineForm operator<< streams the interval enclosure and coefficients")
{
  {
    AffineTVarVector ax(1);
    ax[0] = Interval(1.0, 2.0);
    std::ostringstream stream;
    stream << ax[0];

    std::ostringstream expected_prefix;
    expected_prefix << ax[0].itv();

    CAPTURE(stream.str());
    CHECK(stream.str().substr(0, expected_prefix.str().size()) == expected_prefix.str());
    CHECK(stream.str().find(" : ") != std::string::npos);
    CHECK(stream.str().find("eps_0") != std::string::npos);
  }
  {
    const AffineT empty_form = AffineT::empty();
    std::ostringstream stream;
    stream << empty_form;
    CHECK(stream.str().find("not enabled") != std::string::npos);
  }
}


TEST_CASE("AffineMain size(), init() and init_from_list() reinitialize in place")
{
  // size() is always 1, matching Interval::size()'s template-uniformity role.
  CHECK(AffineT().size() == 1);
  CHECK(AffineT::empty().size() == 1);
  CHECK(AffineT(Interval(1.0, 2.0)).size() == 1);

  AffineTVarVector ax(1);
  ax[0] = Interval(1.0, 2.0);
  AffineT x = ax[0];
  CHECK(x.size() == 1);

  // init() resets to [-oo,oo].
  x.init();
  CHECK(x.itv() == Interval());
  CHECK(x.is_unbounded());

  // init(const Interval&) assigns like operator=(const Interval&): for a
  // plain AffineMain (with no dedicated noise symbol of its own) the whole
  // uncertainty is carried by the remainder error term, and every noise
  // coefficient is left at zero.
  x.init(Interval(3.0, 7.0));
  CHECK(x.itv() == Interval(3.0, 7.0));
  CHECK(x.mid() == 5.0);
  CHECK(x.err() == 2.0);
  for (Index i = 0; i < x.noise_count(); ++i) {
    CAPTURE(i);
    CHECK(x.noise(i) == 0.0);
  }

  // A degenerate interval carries no uncertainty at all.
  AffineT y;
  y.init(Interval(4.0));
  CHECK(y.itv() == Interval(4.0));
  CHECK(y.err() == 0.0);
  CHECK(y.is_degenerated());

  // init_from_list mirrors Interval::init_from_list: one value gives a
  // degenerate point, two values give the bounds of the resulting hull.
  AffineT single;
  single.init_from_list({4.0});
  CHECK(single.itv() == Interval(4.0));
  CHECK(single.is_degenerated());

  AffineT pair;
  pair.init_from_list({1.0, 3.0});
  CHECK(pair.itv() == Interval(1.0, 3.0));
  CHECK(pair.mid() == 2.0);
  CHECK(pair.err() == 1.0);
}


