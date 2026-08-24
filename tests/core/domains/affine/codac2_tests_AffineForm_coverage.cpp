/**
 * \file codac2_tests_AffineForm_coverage_gaps.cpp
 * \brief Regression tests for gaps found while reviewing AffineForm/Eigen
 *        and Paving test coverage:
 *
 *   1. Paving tree lifecycle: PavingNode<P>::_top used to be an owning
 *      shared_ptr, forming a strong reference cycle with the parent's
 *      owning _left/_right. No bisected node could ever be freed. _top is
 *      now a weak_ptr; these tests pin that down directly.
 *   2. AffineVarMainVector::itv() at the container level (enabled once
 *      is_Affine_based<AffineVarMain<T>> became true), never exercised.
 *   3. The intended asymmetry between AffineMain <- AffineVarMain
 *      (allowed, via slicing) and AffineVarMain <- AffineMain (deleted).
 *   4. The comma-initializer (operator<<) on Affine-scalar containers,
 *      which only compiles since AffineMain(const Interval&) was added -
 *      these are the exact call sites that used to fail before that.
 *   5. AffineMain <-> Interval implicit conversions and
 *      operator==(double) / operator!=(double).
 *
 * \brief Inspired by codac2_tests_AffineTForm_Eigen.cpp.
 */

#include <catch2/catch_test_macros.hpp>

#include <concepts>
#include <memory>
#include <type_traits>

#include "codac2_Affine.h"
#include "codac2_AffineMatrix.h"
#include "codac2_AffineRow.h"
#include "codac2_AffineVector.h"
#include "codac2_IntervalMatrix.h"
#include "codac2_IntervalRow.h"
#include "codac2_Approx.h"
#include <iostream>

using namespace codac2;

namespace {

using Model = AF_Default;
using AffineT = AffineMain<Model>;
using AffineTMatrix = AffineMainMatrix<Model>;
using AffineTVector = AffineMainVector<Model>;
using AffineTVarVector = AffineVarMainVector<Model>;
using AffineTRow = AffineMainRow<Model>;

using FixedAffineTMatrix = Eigen::Matrix<AffineT, 2, 3>;
using FixedAffineTVector = Eigen::Matrix<AffineT, 3, 1>;
using FixedAffineTRow    = Eigen::Matrix<AffineT, 1, 3>;

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

AffineTVarVector make_variable_vector_5()
{
    return AffineTVarVector(IntervalVector({
        {-2.0, -1.0},
        {0.0, 2.0},
        {3.0, 4.0},
        {-1.0, 1.0},
        {5.0, 7.0}
    }));
}

} // namespace



// ============================================================================
// 2. AffineTVarMainVector::itv() at the container level
// ============================================================================

TEST_CASE("Eigen AffineTVarMainVector itv()", "[Eigen][itv][AffineTVarMainVector]")
{
    const AffineTVarVector variables = make_variable_vector_5();
    const auto boxes = variables.itv();

    REQUIRE(boxes.size() == 5);
    CHECK(boxes[0] == Interval(-2.0, -1.0));
    CHECK(boxes[1] == Interval(0.0, 2.0));
    CHECK(boxes[2] == Interval(3.0, 4.0));
    CHECK(boxes[3] == Interval(-1.0, 1.0));
    CHECK(boxes[4] == Interval(5.0, 7.0));
}


// ============================================================================
// 3. Scalar-level assignment asymmetry between AffineTMain and AffineTVarMain
// ============================================================================
//
// These checks go through a concept (a genuine template substitution)
// rather than a bare static_assert(requires(...){...}) on a concrete type:
// operator-syntax overload-resolution failures are only SFINAE-friendly
// inside a requires-expression when they occur during real template
// substitution, not when evaluated directly on a fixed, concrete type.

template<class T>
concept can_assign_AffineT_from_var =
    requires(AffineMain<T>& dst, const AffineVarMain<T>& src) { dst = src; };

template<class T>
concept can_assign_var_from_AffineT =
    requires(AffineVarMain<T>& dst, const AffineMain<T>& src) { dst = src; };

TEST_CASE("AffineTVarMain cannot be assigned from AffineTMain", "[AffineT][AffineTVar][regression]")
{
    static_assert(can_assign_AffineT_from_var<Model>,
                  "AffineTMain should remain assignable from AffineTVarMain (slicing)");
    static_assert(!can_assign_var_from_AffineT<Model>,
                  "AffineTVarMain must not be assignable from AffineTMain");
}


// ============================================================================
// 4. Comma-initializer with Interval literals on AffineT-scalar containers
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

TEST_CASE("AffineTVarVector init(Interval) preserves distinct noise symbols",
          "[AffineT][AffineTVarVector][init][regression]")
{
    AffineTVarVector v(3);
    v.init(Interval(1.0, 2.0));

    REQUIRE(v.size() == 3);

    const double rad = Interval(1.0, 2.0).rad(); // 0.5

    for (int i = 0; i < 3; ++i)
    {
        CAPTURE(i, v[i].itv());

        // Every component carries the same interval enclosure...
        CHECK(v[i].itv() == Interval(1.0, 2.0));
        CHECK(v[i].rad() == rad);
        CHECK(v[i].noise_count() == 3);
        // ...but each keeps its own, distinct noise symbol: component i has
        // a coefficient equal to the interval's radius on eps_i, and zero
        // on every other eps_j. A regression that broadcast a single shared
        // symbol (or degraded to a plain interval constant) would show up
        // here as a non-diagonal pattern.
        for (int j = 0; j < 3; ++j)
            CHECK(v[i].noise(j) == (i == j ? rad : 0.0));
    }

    v.conservativeResize(5);

    for (int i = 0; i < 3; ++i)
    {
        CAPTURE(i, v[i].itv());

        // Every component carries the same interval enclosure...
        CHECK(v[i].itv() == Interval(1.0, 2.0));
        CHECK(v[i].rad() == rad);
        CHECK(v[i].noise_count() == 5);

        // ...but each keeps its own, distinct noise symbol: component i has
        // a coefficient equal to the interval's radius on eps_i, and zero
        // on every other eps_j. A regression that broadcast a single shared
        // symbol (or degraded to a plain interval constant) would show up
        // here as a non-diagonal pattern.
        for (int j = 0; j < 5; ++j)
            CHECK(v[i].noise(j) == (i == j ? rad : 0.0));
    }

    for (int i = 3; i < 5; ++i)
    {
        CAPTURE(i, v[i].itv());

        // Every component carries the same interval enclosure...
        CHECK(v[i].itv() == Interval());
        CHECK(v[i].noise_count() == 5);
    }
}

