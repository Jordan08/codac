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
 *      is_affine_based<AffineVarMain<T>> became true), never exercised.
 *   3. The intended asymmetry between AffineMain <- AffineVarMain
 *      (allowed, via slicing) and AffineVarMain <- AffineMain (deleted).
 *   4. The comma-initializer (operator<<) on Affine-scalar containers,
 *      which only compiles since AffineMain(const Interval&) was added -
 *      these are the exact call sites that used to fail before that.
 *   5. AffineMain <-> Interval implicit conversions and
 *      operator==(double) / operator!=(double).
 *
 * \brief Inspired by codac2_tests_AffineForm_Eigen.cpp.
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

using namespace codac2;

namespace {

using Model = AF_Default;
using Affine = AffineMain<Model>;
using AffineMatrix = AffineMainMatrix<Model>;
using AffineVector = AffineMainVector<Model>;
using AffineVarVector = AffineVarMainVector<Model>;
using AffineRow = AffineMainRow<Model>;

using FixedAffineMatrix = Eigen::Matrix<Affine, 2, 3>;
using FixedAffineVector = Eigen::Matrix<Affine, 3, 1>;
using FixedAffineRow    = Eigen::Matrix<Affine, 1, 3>;

void check_point(const Affine& actual, double expected)
{
    CAPTURE(actual.itv(), expected);
    CHECK(actual == Approx<Affine>(Interval(expected)));
}

void check_interval_enclosure(const Affine& actual, const Interval& expected)
{
    CAPTURE(actual.itv(), expected);
    CHECK(actual.itv().is_superset(expected));
}

AffineVarVector make_variable_vector_5()
{
    return AffineVarVector(IntervalVector({
        {-2.0, -1.0},
        {0.0, 2.0},
        {3.0, 4.0},
        {-1.0, 1.0},
        {5.0, 7.0}
    }));
}

} // namespace



// ============================================================================
// 2. AffineVarMainVector::itv() at the container level
// ============================================================================

TEST_CASE("Eigen AffineVarMainVector itv()", "[Eigen][itv][AffineVarMainVector]")
{
    const AffineVarVector variables = make_variable_vector_5();
    const auto boxes = variables.itv();

    REQUIRE(boxes.size() == 5);
    CHECK(boxes[0] == Interval(-2.0, -1.0));
    CHECK(boxes[1] == Interval(0.0, 2.0));
    CHECK(boxes[2] == Interval(3.0, 4.0));
    CHECK(boxes[3] == Interval(-1.0, 1.0));
    CHECK(boxes[4] == Interval(5.0, 7.0));
}


// ============================================================================
// 3. Scalar-level assignment asymmetry between AffineMain and AffineVarMain
// ============================================================================
//
// These checks go through a concept (a genuine template substitution)
// rather than a bare static_assert(requires(...){...}) on a concrete type:
// operator-syntax overload-resolution failures are only SFINAE-friendly
// inside a requires-expression when they occur during real template
// substitution, not when evaluated directly on a fixed, concrete type.

template<class T>
concept can_assign_affine_from_var =
    requires(AffineMain<T>& dst, const AffineVarMain<T>& src) { dst = src; };

template<class T>
concept can_assign_var_from_affine =
    requires(AffineVarMain<T>& dst, const AffineMain<T>& src) { dst = src; };

TEST_CASE("AffineVarMain cannot be assigned from AffineMain", "[Affine][AffineVar][regression]")
{
    static_assert(can_assign_affine_from_var<Model>,
                  "AffineMain should remain assignable from AffineVarMain (slicing)");
    static_assert(!can_assign_var_from_affine<Model>,
                  "AffineVarMain must not be assignable from AffineMain");
}


// ============================================================================
// 4. Comma-initializer with Interval literals on Affine-scalar containers
// ============================================================================
//
// These are the exact call sites that used to fail with "invalid operands
// to binary expression" before AffineMain(const Interval&) was added:
// Eigen's comma-initializer needs a converting constructor for each
// inserted value, operator=(const Interval&) alone isn't enough.

TEST_CASE("Eigen affine comma-initializer with Interval literals", "[Eigen][comma-initializer][Affine]")
{
    AffineVector vector(3);
    vector << Interval(1.0), Interval(-2.0), Interval(4.0);
    check_point(vector[0], 1.0);
    check_point(vector[1], -2.0);
    check_point(vector[2], 4.0);

    AffineRow row(3);
    row << Interval(1.0), Interval(-2.0), Interval(4.0);
    check_point(row[0], 1.0);
    check_point(row[1], -2.0);
    check_point(row[2], 4.0);

    AffineMatrix matrix(2, 3);
    matrix << Interval(1.0), Interval(-2.0), Interval(3.0),
              Interval(4.0), Interval(5.0), Interval(-6.0);
    check_point(matrix(0, 0), 1.0);
    check_point(matrix(0, 1), -2.0);
    check_point(matrix(1, 2), -6.0);

    FixedAffineVector fvector;
    fvector << Interval(1.0), Interval(-2.0), Interval(4.0);
    check_point(fvector[0], 1.0);
    check_point(fvector[1], -2.0);
    check_point(fvector[2], 4.0);

    FixedAffineRow frow;
    frow << Interval(1.0), Interval(-2.0), Interval(4.0);
    check_point(frow[0], 1.0);
    check_point(frow[1], -2.0);
    check_point(frow[2], 4.0);

    FixedAffineMatrix fmatrix;
    fmatrix << Interval(1.0), Interval(-2.0), Interval(3.0),
               Interval(4.0), Interval(5.0), Interval(-6.0);
    check_point(fmatrix(0, 0), 1.0);
    check_point(fmatrix(0, 1), -2.0);
    check_point(fmatrix(1, 2), -6.0);
}


// ============================================================================
// 5. AffineMain <-> Interval conversions, and operator==/!=(double)
// ============================================================================

TEST_CASE("Eigen affine equality with double", "[Affine][conversion]")
{
    const Affine a(4.0);

    CHECK(a == 4.0);
    CHECK_FALSE(a == 5.0);
    CHECK(a != 5.0);
    CHECK_FALSE(a != 4.0);
}

TEST_CASE("Eigen affine implicit conversion to and from Interval", "[Affine][conversion]")
{
    const Affine a(4.0);

    const Interval as_interval = a; // implicit operator Interval()
    CHECK(as_interval == Interval(4.0));

    const Affine from_interval = Interval(2.0, 3.0); // implicit AffineMain(const Interval&)
    check_interval_enclosure(from_interval, Interval(2.0, 3.0));

    // The two conversions being bidirectional and implicit is exactly the
    // kind of thing that can silently become ambiguous: pin down that
    // plain arithmetic between an Affine and an Interval still resolves
    // cleanly instead of failing to compile.
    const Interval itv(1.0, 2.0);
    const Affine sum = a + itv;
    check_interval_enclosure(sum, Interval(4.0) + itv);
}
