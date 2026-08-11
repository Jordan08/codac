/**
 * \file codac2_tests_AffineEigen.cpp
 * \brief Unit tests for Eigen operations applied to affine forms.
 */

#include <catch2/catch_test_macros.hpp>

#include <concepts>
#include <type_traits>
#include <utility>

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

using DynamicAffineMatrix =  Eigen::Matrix<Affine, Eigen::Dynamic, Eigen::Dynamic>;
using FixedAffineMatrix =    Eigen::Matrix<Affine, 2, 3>;
using FixedAffineVector =    Eigen::Matrix<Affine, 3, 1>;
using FixedAffineRow =       Eigen::Matrix<Affine, 1, 3>;

void check_point(const Affine& actual, double expected)
{
    CAPTURE(actual.itv(), expected);
    CHECK(actual == Approx<Affine>(Interval(expected)));
}

void check_interval_enclosure(
    const Affine& actual,
    const Interval& expected)
{
    CAPTURE(actual.itv(), expected);
    CHECK(actual.itv().is_superset(expected));
}

AffineMatrix make_point_matrix_3x4()
{
    AffineMatrix matrix(3, 4);
    matrix(0, 0) = Interval(1.0);  matrix(0, 1) = Interval(2.0);
    matrix(0, 2) = Interval(3.0);  matrix(0, 3) = Interval(4.0);
    matrix(1, 0) = Interval(5.0);  matrix(1, 1) = Interval(6.0);
    matrix(1, 2) = Interval(7.0);  matrix(1, 3) = Interval(8.0);
    matrix(2, 0) = Interval(9.0);  matrix(2, 1) = Interval(10.0);
    matrix(2, 2) = Interval(11.0); matrix(2, 3) = Interval(12.0);
    return matrix;
}

AffineVector make_point_vector_5()
{
    const AffineVarVector variables(
        IntervalVector({{1}, {2}, {3}, {4}, {5}}));
    return AffineVector(variables);
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

static_assert(std::same_as<typename AffineMatrix::Scalar, Affine>);
static_assert(std::same_as<typename AffineVector::Scalar, Affine>);
static_assert(std::same_as<typename AffineRow::Scalar, Affine>);
static_assert(AffineVector::ColsAtCompileTime == 1);
static_assert(AffineRow::RowsAtCompileTime == 1);

// AffineVarMain deliberately deletes scalar compound assignments.
// Keep this assertion here so that this API contract is explicit and
// cannot accidentally become enabled by the Eigen plugin.
template<class T> concept has_scalar_plus_eq  = requires(T& v, double d) { v += d; };
template<class T> concept has_scalar_minus_eq = requires(T& v, double d) { v -= d; };
template<class T> concept has_scalar_times_eq = requires(T& v, double d) { v *= d; };
template<class T> concept has_scalar_div_eq   = requires(T& v, double d) { v /= d; };

TEST_CASE("Eigen affine matrix dimensions and coefficient access")
{
    AffineMatrix matrix = make_point_matrix_3x4();

    CHECK(matrix.rows() == 3);
    CHECK(matrix.cols() == 4);
    CHECK(matrix.size() == 12);

    check_point(matrix(0, 0), 1.0);
    check_point(matrix(1, 2), 7.0);
    check_point(matrix.coeff(2, 3), 12.0);

    matrix.coeffRef(1, 1) = -6.0;
    check_point(matrix(1, 1), -6.0);
}

TEST_CASE("Eigen affine vector indexing")
{
    AffineVector vector = make_point_vector_5();

    CHECK(vector.rows() == 5);
    CHECK(vector.cols() == 1);
    CHECK(vector.size() == 5);

    check_point(vector[0], 1.0);
    check_point(vector(2), 3.0);
    check_point(vector.coeff(4), 5.0);

    vector.coeffRef(3) = -4.0;
    check_point(vector[3], -4.0);
}

TEST_CASE("Eigen affine matrix row and column views")
{
    AffineMatrix matrix = make_point_matrix_3x4();

    const AffineVector column = matrix.col(2);
    REQUIRE(column.size() == 3);
    check_point(column[0], 3.0);
    check_point(column[1], 7.0);
    check_point(column[2], 11.0);

    const AffineRow row = matrix.row(1);
    REQUIRE(row.size() == 4);
    check_point(row[0], 5.0);
    check_point(row[1], 6.0);
    check_point(row[2], 7.0);
    check_point(row[3], 8.0);
}

TEST_CASE("Eigen affine row and column assignments")
{
    AffineMatrix matrix(3, 3);
    matrix.setZero();

    const AffineVarVector column_variables(
        IntervalVector({{1}, {2}, {3}}));
    const AffineVector column(column_variables);
    matrix.col(1) = column;

    check_point(matrix(0, 1), 1.0);
    check_point(matrix(1, 1), 2.0);
    check_point(matrix(2, 1), 3.0);

    const AffineVarVector row_variables(
        IntervalVector({{4}, {5}, {6}}));
    const AffineRow row = AffineVector(row_variables).transpose();
    matrix.row(2) = row;

    check_point(matrix(2, 0), 4.0);
    check_point(matrix(2, 1), 5.0);
    check_point(matrix(2, 2), 6.0);
}

TEST_CASE("Eigen AffineVarMainVector assignment to matrix column")
{
    AffineMatrix matrix(3, 2);
    matrix.setZero();

    const AffineVarVector variables(
        IntervalVector({{-2.0, 1.0}, {3.0, 4.0}, {5.0, 7.0}}));

    matrix.col(0) = variables;

    for(Index i = 0; i < variables.size(); ++i)
        CHECK(matrix(i, 0).itv() == variables[i].itv());
}

TEST_CASE("Eigen affine block extraction")
{
    const AffineMatrix matrix = make_point_matrix_3x4();
    const AffineMatrix block = matrix.block(1, 1, 2, 3);

    REQUIRE(block.rows() == 2);
    REQUIRE(block.cols() == 3);

    check_point(block(0, 0), 6.0);
    check_point(block(0, 2), 8.0);
    check_point(block(1, 0), 10.0);
    check_point(block(1, 2), 12.0);
}

TEST_CASE("Eigen affine block assignment")
{
    AffineMatrix matrix(4, 5);
    matrix.setZero();

    AffineMatrix block(2, 3);
    block(0, 0) = Interval(1.0); block(0, 1) = Interval(2.0); block(0, 2) = Interval(3.0);
    block(1, 0) = Interval(4.0); block(1, 1) = Interval(5.0); block(1, 2) = Interval(6.0);

    matrix.block(1, 1, 2, 3) = block;

    check_point(matrix(1, 1), 1.0);
    check_point(matrix(1, 3), 3.0);
    check_point(matrix(2, 1), 4.0);
    check_point(matrix(2, 3), 6.0);
    check_point(matrix(0, 0), 0.0);
    check_point(matrix(3, 4), 0.0);
}

TEST_CASE("Eigen affine top bottom left and right corners")
{
    const AffineMatrix matrix = make_point_matrix_3x4();

    const AffineMatrix top_left = matrix.topLeftCorner(2, 2);
    const AffineMatrix bottom_right = matrix.bottomRightCorner(2, 2);

    check_point(top_left(0, 0), 1.0);
    check_point(top_left(1, 1), 6.0);
    check_point(bottom_right(0, 0), 7.0);
    check_point(bottom_right(1, 1), 12.0);
}

TEST_CASE("Eigen affine middle rows and columns")
{
    const AffineMatrix matrix = make_point_matrix_3x4();

    const AffineMatrix middle_columns = matrix.middleCols(1, 2);
    REQUIRE(middle_columns.rows() == 3);
    REQUIRE(middle_columns.cols() == 2);
    check_point(middle_columns(0, 0), 2.0);
    check_point(middle_columns(2, 1), 11.0);

    const AffineMatrix middle_rows = matrix.middleRows(1, 2);
    REQUIRE(middle_rows.rows() == 2);
    REQUIRE(middle_rows.cols() == 4);
    check_point(middle_rows(0, 0), 5.0);
    check_point(middle_rows(1, 3), 12.0);
}

TEST_CASE("Eigen affine vector segment head and tail")
{
    AffineVector vector = make_point_vector_5();

    const AffineVector segment = vector.segment(1, 3);
    const AffineVector head = vector.head(2);
    const AffineVector tail = vector.tail(2);

    REQUIRE(segment.size() == 3);
    check_point(segment[0], 2.0);
    check_point(segment[1], 3.0);
    check_point(segment[2], 4.0);

    check_point(head[0], 1.0);
    check_point(head[1], 2.0);
    check_point(tail[0], 4.0);
    check_point(tail[1], 5.0);
}

TEST_CASE("Eigen AffineVarMainVector segment conversion")
{
    const AffineVarVector variables = make_variable_vector_5();
    const AffineVector segment = variables.segment(1, 3);

    REQUIRE(segment.size() == 3);
    CHECK(segment[0].itv() == variables[1].itv());
    CHECK(segment[1].itv() == variables[2].itv());
    CHECK(segment[2].itv() == variables[3].itv());
}

TEST_CASE("Eigen affine segment assignment")
{
    AffineVector vector(6);
    vector.setZero();

    const AffineVarVector source_variables(
        IntervalVector({{2}, {4}, {6}}));
    const AffineVector source(source_variables);
    vector.segment(2, 3) = source;

    check_point(vector[0], 0.0);
    check_point(vector[1], 0.0);
    check_point(vector[2], 2.0);
    check_point(vector[3], 4.0);
    check_point(vector[4], 6.0);
    check_point(vector[5], 0.0);
}

TEST_CASE("Eigen affine transpose")
{
    const AffineMatrix matrix = make_point_matrix_3x4();
    const AffineMatrix transpose = matrix.transpose();

    REQUIRE(transpose.rows() == 4);
    REQUIRE(transpose.cols() == 3);

    check_point(transpose(0, 0), 1.0);
    check_point(transpose(2, 1), 7.0);
    check_point(transpose(3, 2), 12.0);
}

TEST_CASE("Eigen affine vector transpose to row")
{
    const AffineVector vector = make_point_vector_5();
    const AffineRow row = vector.transpose();

    REQUIRE(row.rows() == 1);
    REQUIRE(row.cols() == 5);
    check_point(row[0], 1.0);
    check_point(row[4], 5.0);
}

TEST_CASE("Eigen affine in-place transpose for square matrix")
{
    AffineMatrix matrix(3, 3);
    matrix(0, 0) = Interval(1.0); matrix(0, 1) = Interval(2.0); matrix(0, 2) = Interval(3.0);
    matrix(1, 0) = Interval(4.0); matrix(1, 1) = Interval(5.0); matrix(1, 2) = Interval(6.0);
    matrix(2, 0) = Interval(7.0); matrix(2, 1) = Interval(8.0); matrix(2, 2) = Interval(9.0);

    matrix.transposeInPlace();

    check_point(matrix(0, 1), 4.0);
    check_point(matrix(1, 0), 2.0);
    check_point(matrix(2, 1), 6.0);
}

TEST_CASE("Eigen affine diagonal extraction")
{
    AffineMatrix matrix(3, 3);
    matrix(0, 0) = Interval(1.0); matrix(0, 1) = Interval(2.0); matrix(0, 2) = Interval(3.0);
    matrix(1, 0) = Interval(4.0); matrix(1, 1) = Interval(5.0); matrix(1, 2) = Interval(6.0);
    matrix(2, 0) = Interval(7.0); matrix(2, 1) = Interval(8.0); matrix(2, 2) = Interval(9.0);

    const AffineVector diagonal = matrix.diagonal();
    REQUIRE(diagonal.size() == 3);
    check_point(diagonal[0], 1.0);
    check_point(diagonal[1], 5.0);
    check_point(diagonal[2], 9.0);
}

TEST_CASE("Eigen affine diagonal assignment")
{
    AffineMatrix matrix(3, 3);
    matrix.setZero();

    const AffineVarVector diagonal_variables(
        IntervalVector({{2}, {4}, {6}}));
    const AffineVector diagonal(diagonal_variables);
    matrix.diagonal() = diagonal;

    check_point(matrix(0, 0), 2.0);
    check_point(matrix(1, 1), 4.0);
    check_point(matrix(2, 2), 6.0);
    check_point(matrix(0, 1), 0.0);
}

TEST_CASE("Eigen affine setZero setOnes and setConstant")
{
    AffineMatrix matrix(2, 3);

    matrix.setZero();
    for(Index i = 0; i < matrix.rows(); ++i)
        for(Index j = 0; j < matrix.cols(); ++j)
            check_point(matrix(i, j), 0.0);

    matrix.setOnes();
    for(Index i = 0; i < matrix.rows(); ++i)
        for(Index j = 0; j < matrix.cols(); ++j)
            check_point(matrix(i, j), 1.0);

    matrix.setConstant(Affine(3.5));
    for(Index i = 0; i < matrix.rows(); ++i)
        for(Index j = 0; j < matrix.cols(); ++j)
            check_point(matrix(i, j), 3.5);
}

TEST_CASE("Eigen affine Zero Ones and Constant factories")
{
    const AffineMatrix zero = AffineMatrix::Zero(2, 2);
    const AffineMatrix ones = AffineMatrix::Ones(2, 2);
    const AffineMatrix constant = AffineMatrix::Constant(2, 2, Affine(-3.0));

    for(Index i = 0; i < 2; ++i)
    {
        for(Index j = 0; j < 2; ++j)
        {
            check_point(zero(i, j), 0.0);
            check_point(ones(i, j), 1.0);
            check_point(constant(i, j), -3.0);
        }
    }
}

TEST_CASE("Eigen affine Identity factory")
{
    const AffineMatrix identity = AffineMatrix::Identity(3, 3);

    for(Index i = 0; i < identity.rows(); ++i)
    {
        for(Index j = 0; j < identity.cols(); ++j)
            check_point(identity(i, j), i == j ? 1.0 : 0.0);
    }
}

TEST_CASE("Eigen affine resize and conservativeResize")
{
    AffineMatrix matrix(2, 2);
    matrix(0, 0) = Interval(1.0); matrix(0, 1) = Interval(2.0);
    matrix(1, 0) = Interval(3.0); matrix(1, 1) = Interval(4.0);

    matrix.conservativeResize(3, 3);
    REQUIRE(matrix.rows() == 3);
    REQUIRE(matrix.cols() == 3);
    check_point(matrix(0, 0), 1.0);
    check_point(matrix(0, 1), 2.0);
    check_point(matrix(1, 0), 3.0);
    check_point(matrix(1, 1), 4.0);

    matrix.resize(4, 2);
    CHECK(matrix.rows() == 4);
    CHECK(matrix.cols() == 2);
}

TEST_CASE("Eigen affine swap")
{
    const AffineVarVector lhs_variables(
        IntervalVector({{1}, {2}, {3}}));
    const AffineVarVector rhs_variables(
        IntervalVector({{4}, {5}, {6}}));
    AffineVector lhs(lhs_variables);
    AffineVector rhs(rhs_variables);

    lhs.swap(rhs);

    check_point(lhs[0], 4.0);
    check_point(lhs[2], 6.0);
    check_point(rhs[0], 1.0);
    check_point(rhs[2], 3.0);
}

TEST_CASE("Eigen affine arithmetic expressions")
{
    AffineMatrix lhs(2, 2);
    AffineMatrix rhs(2, 2);
    lhs(0, 0) = Interval(1.0); lhs(0, 1) = Interval(2.0);
    lhs(1, 0) = Interval(3.0); lhs(1, 1) = Interval(4.0);
    rhs(0, 0) = Interval(5.0); rhs(0, 1) = Interval(6.0);
    rhs(1, 0) = Interval(7.0); rhs(1, 1) = Interval(8.0);

    const AffineMatrix sum = lhs + rhs;
    const AffineMatrix difference = rhs - lhs;
    const AffineMatrix expression = 2.0 * lhs - rhs;

    check_point(sum(0, 0), 6.0);
    check_point(sum(1, 1), 12.0);
    check_point(difference(0, 1), 4.0);
    check_point(difference(1, 0), 4.0);
    check_point(expression(0, 0), -3.0);
    check_point(expression(1, 1), 0.0);
}

TEST_CASE("Eigen affine compound assignment")
{
    const AffineVarVector value_variables(
        IntervalVector({{1}, {2}, {3}}));
    const AffineVarVector other_variables(
        IntervalVector({{4}, {5}, {6}}));
    AffineVector value(value_variables);
    const AffineVector other(other_variables);

    value += other;
    check_point(value[0], 5.0);
    check_point(value[2], 9.0);

    value -= other;
    check_point(value[0], 1.0);
    check_point(value[2], 3.0);

    value += 3.0;
    check_point(value[0], 4.0);
    check_point(value[2], 6.0);

    value -= 5;
    check_point(value[0], -1.0);
    check_point(value[2], 1.0);

    value *= 2.0;
    check_point(value[0], -2.0);
    check_point(value[2], 2.0);

    value /= 2.0;
    check_point(value[0], -1.0);
    check_point(value[2], 1.0);
}


// ============================================================================
// Scalar compound assignments on all supported Eigen/CODAC containers
// ============================================================================
//
// These tests intentionally exercise +=, -=, *= and /= independently.
// The existing "Eigen affine compound assignment" test covers the basic
// AffineVector path; the tests below make sure the Eigen plugin also works
// for matrices, rows, interval containers and fixed-size Eigen types.
//

TEST_CASE("Eigen scalar compound assignments - IntervalVector",
          "[Eigen][compound-assignment][IntervalVector]")
{
    IntervalVector value({{1.0}, {-2.0}, {4.0}});

    value += 2.0;
    CHECK(value[0] == Interval(3.0));
    CHECK(value[1] == Interval(0.0));
    CHECK(value[2] == Interval(6.0));

    value -= 2.0;
    CHECK(value[0] == Interval(1.0));
    CHECK(value[1] == Interval(-2.0));
    CHECK(value[2] == Interval(4.0));

    value *= 2.0;
    CHECK(value[0] == Interval(2.0));
    CHECK(value[1] == Interval(-4.0));
    CHECK(value[2] == Interval(8.0));

    value /= 2.0;
    CHECK(value[0] == Interval(1.0));
    CHECK(value[1] == Interval(-2.0));
    CHECK(value[2] == Interval(4.0));
}

TEST_CASE("Eigen scalar compound assignments - IntervalRow",
          "[Eigen][compound-assignment][IntervalRow]")
{
    IntervalRow value(3);
    value[0] = Interval(1.0);
    value[1] = Interval(-2.0);
    value[2] = Interval(4.0);

    value += 2.0;
    CHECK(value[0] == Interval(3.0));
    CHECK(value[1] == Interval(0.0));
    CHECK(value[2] == Interval(6.0));

    value -= 2.0;
    CHECK(value[0] == Interval(1.0));
    CHECK(value[1] == Interval(-2.0));
    CHECK(value[2] == Interval(4.0));

    value *= 2.0;
    CHECK(value[0] == Interval(2.0));
    CHECK(value[1] == Interval(-4.0));
    CHECK(value[2] == Interval(8.0));

    value /= 2.0;
    CHECK(value[0] == Interval(1.0));
    CHECK(value[1] == Interval(-2.0));
    CHECK(value[2] == Interval(4.0));
}

TEST_CASE("Eigen scalar compound assignments - IntervalMatrix",
          "[Eigen][compound-assignment][IntervalMatrix]")
{
    IntervalMatrix value(2, 3);
    value << Interval(1.0), Interval(-2.0), Interval(3.0),
             Interval(4.0), Interval(5.0), Interval(-6.0);

    value += 2.0;
    CHECK(value(0, 0) == Interval(3.0));
    CHECK(value(0, 1) == Interval(0.0));
    CHECK(value(1, 2) == Interval(-4.0));

    value -= 2.0;
    CHECK(value(0, 0) == Interval(1.0));
    CHECK(value(0, 1) == Interval(-2.0));
    CHECK(value(1, 2) == Interval(-6.0));

    value *= 2.0;
    CHECK(value(0, 0) == Interval(2.0));
    CHECK(value(0, 1) == Interval(-4.0));
    CHECK(value(1, 2) == Interval(-12.0));

    value /= 2.0;
    CHECK(value(0, 0) == Interval(1.0));
    CHECK(value(0, 1) == Interval(-2.0));
    CHECK(value(1, 2) == Interval(-6.0));
}

TEST_CASE("Eigen scalar compound assignments - AffineVector",
          "[Eigen][compound-assignment][AffineVector]")
{
    AffineVector value(3);
    value[0] = Interval(1.0);
    value[1] = Interval(-2.0);
    value[2] = Interval(4.0);

    value += 2.0;
    check_point(value[0], 3.0);
    check_point(value[1], 0.0);
    check_point(value[2], 6.0);

    value -= 2.0;
    check_point(value[0], 1.0);
    check_point(value[1], -2.0);
    check_point(value[2], 4.0);

    value *= 2.0;
    check_point(value[0], 2.0);
    check_point(value[1], -4.0);
    check_point(value[2], 8.0);

    value /= 2.0;
    check_point(value[0], 1.0);
    check_point(value[1], -2.0);
    check_point(value[2], 4.0);
}

TEST_CASE("Eigen scalar compound assignments - AffineRow",
          "[Eigen][compound-assignment][AffineRow]")
{
    AffineRow value(3);
    value[0] = Interval(1.0);
    value[1] = Interval(-2.0);
    value[2] = Interval(4.0);

    value += 2.0;
    check_point(value[0], 3.0);
    check_point(value[1], 0.0);
    check_point(value[2], 6.0);

    value -= 2.0;
    check_point(value[0], 1.0);
    check_point(value[1], -2.0);
    check_point(value[2], 4.0);

    value *= 2.0;
    check_point(value[0], 2.0);
    check_point(value[1], -4.0);
    check_point(value[2], 8.0);

    value /= 2.0;
    check_point(value[0], 1.0);
    check_point(value[1], -2.0);
    check_point(value[2], 4.0);
}

TEST_CASE("Eigen scalar compound assignments - AffineMatrix",
          "[Eigen][compound-assignment][AffineMatrix]")
{
    AffineMatrix value(2, 3);
    IntervalMatrix m({{ Interval(1.0), Interval(-2.0), Interval(3.0)},{
             Interval(4.0), Interval(5.0), Interval(-6.0)}});
    value = m;

    value += 2.0;
    check_point(value(0, 0), 3.0);
    check_point(value(0, 1), 0.0);
    check_point(value(1, 2), -4.0);

    value -= 2.0;
    check_point(value(0, 0), 1.0);
    check_point(value(0, 1), -2.0);
    check_point(value(1, 2), -6.0);

    value *= 2.0;
    check_point(value(0, 0), 2.0);
    check_point(value(0, 1), -4.0);
    check_point(value(1, 2), -12.0);

    value /= 2.0;
    check_point(value(0, 0), 1.0);
    check_point(value(0, 1), -2.0);
    check_point(value(1, 2), -6.0);
}

TEST_CASE("Eigen scalar compound assignments - fixed-size Affine types",
          "[Eigen][compound-assignment][Affine][fixed-size]")
{
    FixedAffineVector vector;

    IntervalVector v({ Interval(1.0), Interval(-2.0), Interval(4.0)});
	vector = v;

    vector += 2.0;
    check_point(vector[0], 3.0);
    check_point(vector[1], 0.0);
    check_point(vector[2], 6.0);

    vector -= 2.0;
    check_point(vector[0], 1.0);
    check_point(vector[1], -2.0);
    check_point(vector[2], 4.0);

    vector *= 2.0;
    check_point(vector[0], 2.0);
    check_point(vector[1], -4.0);
    check_point(vector[2], 8.0);

    vector /= 2.0;
    check_point(vector[0], 1.0);
    check_point(vector[1], -2.0);
    check_point(vector[2], 4.0);

    FixedAffineRow row;
    IntervalVector vv({ Interval(1.0), Interval(-2.0), Interval(4.0)});
    row = vv.transpose();

    row += 2.0;
    check_point(row[0], 3.0);
    check_point(row[1], 0.0);
    check_point(row[2], 6.0);

    row -= 2.0;
    check_point(row[0], 1.0);
    check_point(row[1], -2.0);
    check_point(row[2], 4.0);

    row *= 2.0;
    check_point(row[0], 2.0);
    check_point(row[1], -4.0);
    check_point(row[2], 8.0);

    row /= 2.0;
    check_point(row[0], 1.0);
    check_point(row[1], -2.0);
    check_point(row[2], 4.0);

    FixedAffineMatrix matrix;
    IntervalMatrix mm ({{Interval(1.0), Interval(-2.0), Interval(3.0)},
              {Interval(4.0), Interval(5.0), Interval(-6.0)}});
    matrix = mm;

    matrix += 2.0;
    check_point(matrix(0, 0), 3.0);
    check_point(matrix(0, 1), 0.0);
    check_point(matrix(1, 2), -4.0);

    matrix -= 2.0;
    check_point(matrix(0, 0), 1.0);
    check_point(matrix(0, 1), -2.0);
    check_point(matrix(1, 2), -6.0);

    matrix *= 2.0;
    check_point(matrix(0, 0), 2.0);
    check_point(matrix(0, 1), -4.0);
    check_point(matrix(1, 2), -12.0);

    matrix /= 2.0;
    check_point(matrix(0, 0), 1.0);
    check_point(matrix(0, 1), -2.0);
    check_point(matrix(1, 2), -6.0);
}

TEST_CASE("Eigen scalar compound assignments - fixed-size Interval types",
          "[Eigen][compound-assignment][Interval][fixed-size]")
{
    Eigen::Matrix<Interval, 3, 1> vector;
    IntervalVector v({Interval(1.0), Interval(-2.0), Interval(4.0)});

    vector = v;

    vector += 2.0;
    CHECK(vector[0] == Interval(3.0));
    CHECK(vector[1] == Interval(0.0));
    CHECK(vector[2] == Interval(6.0));

    vector -= 2.0;
    CHECK(vector[0] == Interval(1.0));
    CHECK(vector[1] == Interval(-2.0));
    CHECK(vector[2] == Interval(4.0));

    vector *= 2.0;
    CHECK(vector[0] == Interval(2.0));
    CHECK(vector[1] == Interval(-4.0));
    CHECK(vector[2] == Interval(8.0));

    vector /= 2.0;
    CHECK(vector[0] == Interval(1.0));
    CHECK(vector[1] == Interval(-2.0));
    CHECK(vector[2] == Interval(4.0));

    Eigen::Matrix<Interval, 1, 3> row;
    IntervalVector vv({Interval(1.0), Interval(-2.0), Interval(4.0)});
    row = vv.transpose();

    row += 2.0;
    CHECK(row[0] == Interval(3.0));
    CHECK(row[1] == Interval(0.0));
    CHECK(row[2] == Interval(6.0));

    row -= 2.0;
    CHECK(row[0] == Interval(1.0));
    CHECK(row[1] == Interval(-2.0));
    CHECK(row[2] == Interval(4.0));

    row *= 2.0;
    CHECK(row[0] == Interval(2.0));
    CHECK(row[1] == Interval(-4.0));
    CHECK(row[2] == Interval(8.0));

    row /= 2.0;
    CHECK(row[0] == Interval(1.0));
    CHECK(row[1] == Interval(-2.0));
    CHECK(row[2] == Interval(4.0));

    Eigen::Matrix<Interval, 2, 3> matrix;
    IntervalMatrix mm({{ Interval(1.0), Interval(-2.0), Interval(3.0)},{
              Interval(4.0), Interval(5.0), Interval(-6.0)}});
    matrix = mm;

    matrix += 2.0;
    CHECK(matrix(0, 0) == Interval(3.0));
    CHECK(matrix(0, 1) == Interval(0.0));
    CHECK(matrix(1, 2) == Interval(-4.0));

    matrix -= 2.0;
    CHECK(matrix(0, 0) == Interval(1.0));
    CHECK(matrix(0, 1) == Interval(-2.0));
    CHECK(matrix(1, 2) == Interval(-6.0));

    matrix *= 2.0;
    CHECK(matrix(0, 0) == Interval(2.0));
    CHECK(matrix(0, 1) == Interval(-4.0));
    CHECK(matrix(1, 2) == Interval(-12.0));

    matrix /= 2.0;
    CHECK(matrix(0, 0) == Interval(1.0));
    CHECK(matrix(0, 1) == Interval(-2.0));
    CHECK(matrix(1, 2) == Interval(-6.0));
}

TEST_CASE("Eigen scalar compound assignments - AffineVarMainVector are unsupported",
          "[Eigen][compound-assignment][AffineVarMainVector]")
{


	static_assert(!has_scalar_plus_eq<AffineVarVector>);
	static_assert(!has_scalar_minus_eq<AffineVarVector>);
	static_assert(!has_scalar_times_eq<AffineVarVector>);
	static_assert(!has_scalar_div_eq<AffineVarVector>);
}

TEST_CASE("Eigen affine unary minus")
{
    const AffineVector vector = make_point_vector_5();
    const AffineVector negative = -vector;

    check_point(negative[0], -1.0);
    check_point(negative[2], -3.0);
    check_point(negative[4], -5.0);
}

TEST_CASE("Eigen affine cwise product")
{
    const AffineVarVector lhs_variables(
        IntervalVector({{1}, {-2}, {3}}));
    const AffineVarVector rhs_variables(
        IntervalVector({{4}, {5}, {-6}}));
    const AffineVector lhs(lhs_variables);
    const AffineVector rhs(rhs_variables);

    const AffineVector result = lhs.cwiseProduct(rhs);
    check_point(result[0], 4.0);
    check_point(result[1], -10.0);
    check_point(result[2], -18.0);
}

TEST_CASE("Eigen affine cwise quotient by point affine values")
{
    const AffineVarVector numerator_variables(
        IntervalVector({{4}, {-10}, {18}}));
    const AffineVarVector denominator_variables(
        IntervalVector({{2}, {5}, {-3}}));
    const AffineVector numerator(numerator_variables);
    const AffineVector denominator(denominator_variables);

    const AffineVector result = numerator.cwiseQuotient(denominator);
    check_point(result[0], 2.0);
    check_point(result[1], -2.0);
    check_point(result[2], -6.0);
}

TEST_CASE("Eigen affine matrix product")
{
    AffineMatrix lhs(2, 3);
    AffineMatrix rhs(3, 2);
    lhs(0, 0) = Interval(1.0); lhs(0, 1) = Interval(2.0); lhs(0, 2) = Interval(3.0);
    lhs(1, 0) = Interval(4.0); lhs(1, 1) = Interval(5.0); lhs(1, 2) = Interval(6.0);
    rhs(0, 0) = Interval(7.0);  rhs(0, 1) = Interval(8.0);
    rhs(1, 0) = Interval(9.0);  rhs(1, 1) = Interval(10.0);
    rhs(2, 0) = Interval(11.0); rhs(2, 1) = Interval(12.0);

    const AffineMatrix result = lhs * rhs;
    check_point(result(0, 0), 58.0);
    check_point(result(0, 1), 64.0);
    check_point(result(1, 0), 139.0);
    check_point(result(1, 1), 154.0);
}

TEST_CASE("Eigen affine matrix vector product")
{
    AffineMatrix matrix(2, 3);
    matrix(0, 0) = Interval(1.0); matrix(0, 1) = Interval(2.0); matrix(0, 2) = Interval(3.0);
    matrix(1, 0) = Interval(4.0); matrix(1, 1) = Interval(5.0); matrix(1, 2) = Interval(6.0);

    const AffineVarVector vector_variables(
        IntervalVector({{2}, {-1}, {3}}));
    const AffineVector vector(vector_variables);

    const AffineVector result = matrix * vector;
    check_point(result[0], 9.0);
    check_point(result[1], 21.0);
}

TEST_CASE("Eigen real matrix times AffineVarMainVector")
{
    IntervalMatrix matrix(2,3);
    matrix(0, 0) = Interval(1.0); matrix(0, 1) = Interval(2.0); matrix(0, 2) = Interval(3.0);
    matrix(1, 0) = Interval(4.0); matrix(1, 1) = Interval(5.0); matrix(1, 2) = Interval(6.0);

    const AffineVarVector variables(Vector({2.0, -1.0, 3.0}));
    const auto result = matrix * variables;

    static_assert(std::same_as<typename decltype(result)::Scalar, Affine>);
    check_point(result[0], 9.0);
    check_point(result[1], 21.0);
}

TEST_CASE("Eigen AffineVarMainVector outer product")
{
    const AffineVarVector variables(Vector({2.0, -1.0, 3.0}));
    IntervalRow row(2);
    row[0] = Interval(4.0);
    row[1] = Interval(-2.0);

    const auto result = variables * row;

    static_assert(std::same_as<typename decltype(result)::Scalar, Affine>);
    REQUIRE(result.rows() == 3);
    REQUIRE(result.cols() == 2);
    check_point(result(0, 0), 8.0);
    check_point(result(0, 1), -4.0);
    check_point(result(2, 0), 12.0);
    check_point(result(2, 1), -6.0);
}

TEST_CASE("Eigen affine dot product")
{
    AffineVector lhs(3);
    AffineVector rhs(3);
    lhs[0] = Interval(1.0); lhs[1] = Interval(2.0); lhs[2] = Interval(3.0);
    rhs[0] = Interval(4.0); rhs[1] = Interval(-1.0); rhs[2] = Interval(2.0);

    const Affine result = lhs.dot(rhs);
    check_point(result, 8.0);
}

TEST_CASE("Eigen affine sum")
{
    const AffineVector vector = make_point_vector_5();
    const Affine result = vector.sum();
    check_point(result, 15.0);
}

TEST_CASE("Eigen affine trace")
{
    AffineMatrix matrix(3, 3);
    matrix(0, 0) = Interval(1.0); matrix(0, 1) = Interval(2.0); matrix(0, 2) = Interval(3.0);
    matrix(1, 0) = Interval(4.0); matrix(1, 1) = Interval(5.0); matrix(1, 2) = Interval(6.0);
    matrix(2, 0) = Interval(7.0); matrix(2, 1) = Interval(8.0); matrix(2, 2) = Interval(9.0);

    const Affine result = matrix.trace();
    check_point(result, 15.0);
}

TEST_CASE("Eigen affine reverse")
{
    const AffineVector vector = make_point_vector_5();
    const AffineVector reversed = vector.reverse();

    check_point(reversed[0], 5.0);
    check_point(reversed[1], 4.0);
    check_point(reversed[4], 1.0);
}

TEST_CASE("Eigen affine replicate")
{
    const AffineVarVector row_variables(
        IntervalVector({{1}, {2}}));
    const AffineRow row = AffineVector(row_variables).transpose();

    const AffineMatrix replicated = row.replicate(3, 2);
    REQUIRE(replicated.rows() == 3);
    REQUIRE(replicated.cols() == 4);

    check_point(replicated(0, 0), 1.0);
    check_point(replicated(0, 3), 2.0);
    check_point(replicated(2, 0), 1.0);
    check_point(replicated(2, 3), 2.0);
}

TEST_CASE("Eigen affine fixed-size matrices")
{
    FixedAffineMatrix matrix;
    matrix(0, 0) = Interval(1.0); matrix(0, 1) = Interval(2.0); matrix(0, 2) = Interval(3.0);
    matrix(1, 0) = Interval(4.0); matrix(1, 1) = Interval(5.0); matrix(1, 2) = Interval(6.0);

    static_assert(FixedAffineMatrix::RowsAtCompileTime == 2);
    static_assert(FixedAffineMatrix::ColsAtCompileTime == 3);

    check_point(matrix(0, 0), 1.0);
    check_point(matrix(1, 2), 6.0);

    const FixedAffineRow row = matrix.row(0);
    check_point(row[0], 1.0);
    check_point(row[2], 3.0);

    const FixedAffineVector column = matrix.transpose().col(1);
    check_point(column[0], 4.0);
    check_point(column[2], 6.0);
}

TEST_CASE("Eigen affine interval-valued block operations")
{
    const AffineVarVector variables(
        IntervalVector({{-2.0, 1.0}, {3.0, 5.0}, {-1.0, 4.0}}));

    AffineMatrix matrix(3, 2);
    matrix.col(0) = variables;
    matrix.col(1) = -variables;

    const AffineMatrix block = matrix.block(0, 0, 2, 2);

    CHECK(block(0, 0).itv() == variables[0].itv());
    CHECK(block(1, 0).itv() == variables[1].itv());
    CHECK(block(0, 1).itv() == (-variables[0]).itv());
    CHECK(block(1, 1).itv() == (-variables[1]).itv());
}

TEST_CASE("Eigen affine interval enclosure through linear expression")
{
    const AffineVarVector variables(
        IntervalVector({{-2.0, 1.0}, {3.0, 5.0}, {-1.0, 4.0}}));
    const AffineVector affine = variables;

    const AffineVector result = 2.0 * affine - affine;

    for(Index i = 0; i < result.size(); ++i)
        check_interval_enclosure(result[i], variables[i].itv());
}

TEST_CASE("Eigen affine eval materializes expressions")
{
    const AffineVarVector lhs_variables(
        IntervalVector({{1}, {2}, {3}}));
    const AffineVarVector rhs_variables(
        IntervalVector({{4}, {5}, {6}}));
    AffineVector lhs(lhs_variables);
    AffineVector rhs(rhs_variables);

    const auto expression = lhs + rhs;
    const AffineVector evaluated = expression.eval();

    check_point(evaluated[0], 5.0);
    check_point(evaluated[1], 7.0);
    check_point(evaluated[2], 9.0);
}


// ============================================================================
// Compound assignments between AffineMain Eigen containers and all
// compatible CODAC scalar/container types.
//
// The left-hand side is always an AffineMainVector or AffineMainMatrix.
// For every compatible right-hand side, +=, -=, *= and /= are exercised
// independently.
//
// Point intervals are used deliberately so that the expected result is
// exact while still testing the Interval/Affine/AffineVar interoperability.
// ============================================================================

TEST_CASE("Eigen AffineVector compound assignments with Interval scalar",
          "[Eigen][compound-assignment][AffineVector][Interval]")
{
    AffineVector initial(3);
    initial[0] = Interval(1.0);
    initial[1] = Interval(2.0);
    initial[2] = Interval(4.0);

    SECTION("operator+=")
    {
        auto value = initial;
        value += Interval(2.0);

        check_point(value[0], 3.0);
        check_point(value[1], 4.0);
        check_point(value[2], 6.0);
    }

    SECTION("operator-=")
    {
        auto value = initial;
        value -= Interval(2.0);

        check_point(value[0], -1.0);
        check_point(value[1], 0.0);
        check_point(value[2], 2.0);
    }

    SECTION("operator*=")
    {
        auto value = initial;
        value *= Interval(2.0);

        check_point(value[0], 2.0);
        check_point(value[1], 4.0);
        check_point(value[2], 8.0);
    }

    SECTION("operator/=")
    {
        auto value = initial;
        value /= Interval(2.0);

        check_point(value[0], 0.5);
        check_point(value[1], 1.0);
        check_point(value[2], 2.0);
    }
}


TEST_CASE("Eigen AffineVector compound assignments with Affine scalar",
          "[Eigen][compound-assignment][AffineVector][Affine]")
{
    AffineVarVector initial(3);
    initial[0] = Interval(1.0);
    initial[1] = Interval(2.0);
    initial[2] = Interval(4.0);

    const Affine scalar =initial[1];

    SECTION("operator+=")
    {
        AffineVector value = initial;
        value += scalar;

        check_point(value[0], 3.0);
        check_point(value[1], 4.0);
        check_point(value[2], 6.0);
    }

    SECTION("operator-=")
    {
    	AffineVector value = initial;
        value -= scalar;

        check_point(value[0], -1.0);
        check_point(value[1], 0.0);
        check_point(value[2], 2.0);
    }

    SECTION("operator*=")
    {
    	AffineVector value = initial;
        value *= scalar;

        check_point(value[0], 2.0);
        check_point(value[1], 4.0);
        check_point(value[2], 8.0);
    }

    SECTION("operator/=")
    {
    	AffineVector value = initial;
        value /= scalar;

        check_point(value[0], 0.5);
        check_point(value[1], 1.0);
        check_point(value[2], 2.0);
    }
}


TEST_CASE("Eigen AffineVector compound assignments with IntervalVector",
          "[Eigen][compound-assignment][AffineVector][IntervalVector]")
{
    AffineVarVector initial(3);
    initial[0] = Interval(1.0);
    initial[1] = Interval(2.0);
    initial[2] = Interval(4.0);

    const IntervalVector rhs({
        {2.0},
        {3.0},
        {4.0}
    });

    SECTION("operator+=")
    {
    	AffineVector value = initial;
        value += rhs;

        check_point(value[0], 3.0);
        check_point(value[1], 5.0);
        check_point(value[2], 8.0);
    }

    SECTION("operator-=")
    {
    	AffineVector value = initial;
        value -= rhs;

        check_point(value[0], -1.0);
        check_point(value[1], -1.0);
        check_point(value[2], 0.0);
    }

//    SECTION("operator*=")
//    {
//    	AffineVector value = initial;
//        value *= rhs;
//
//        check_point(value[0], 2.0);
//        check_point(value[1], 6.0);
//        check_point(value[2], 16.0);
//    }

//    SECTION("operator/=")
//    {
//    	AffineVector value = initial;
//        value /= rhs;
//
//        check_point(value[0], 0.5);
//        check_point(value[1], 2.0 / 3.0);
//        check_point(value[2], 1.0);
//    }
}


TEST_CASE("Eigen AffineVector compound assignments with AffineVector",
          "[Eigen][compound-assignment][AffineVector][AffineVector]")
{
    AffineVarVector initial(3);
    initial[0] = Interval(1.0);
    initial[1] = Interval(2.0);
    initial[2] = Interval(4.0);

    AffineVector rhs(3);
    rhs[0] = Interval(2.0);
    rhs[1] = Interval(3.0);
    rhs[2] = Interval(4.0);

    SECTION("operator+=")
    {
    	AffineVector value = initial;
        value += rhs;

        check_point(value[0], 3.0);
        check_point(value[1], 5.0);
        check_point(value[2], 8.0);
    }

    SECTION("operator-=")
    {
    	AffineVector value = initial;
        value -= rhs;

        check_point(value[0], -1.0);
        check_point(value[1], -1.0);
        check_point(value[2], 0.0);
    }

//    SECTION("operator*=")
//    {
//    	AffineVector value = initial;
//        value *= rhs;
//
//        check_point(value[0], 2.0);
//        check_point(value[1], 6.0);
//        check_point(value[2], 16.0);
//    }

//    SECTION("operator/=")
//    {
//        auto value = initial;
//        AffineVector /= rhs;
//
//        check_point(value[0], 0.5);
//        check_point(value[1], 2.0 / 3.0);
//        check_point(value[2], 1.0);
//    }
}


TEST_CASE("Eigen AffineVector compound assignments with AffineVarVector",
          "[Eigen][compound-assignment][AffineVector][AffineVarVector]")
{
    AffineVarVector initial(3);
    initial[0] = Interval(1.0);
    initial[1] = Interval(2.0);
    initial[2] = Interval(4.0);

    const AffineVarVector rhs(
        IntervalVector({
            {2.0},
            {3.0},
            {4.0}
        }));

    SECTION("operator+=")
    {
    	AffineVector value = initial;
        value += rhs;

        check_point(value[0], 3.0);
        check_point(value[1], 5.0);
        check_point(value[2], 8.0);
    }

    SECTION("operator-=")
    {
    	AffineVector value = initial;
        value -= rhs;

        check_point(value[0], -1.0);
        check_point(value[1], -1.0);
        check_point(value[2], 0.0);
    }

//    SECTION("operator*=")
//    {
//    	AffineVector value = initial;
//        value *= rhs;
//
//        check_point(value[0], 2.0);
//        check_point(value[1], 6.0);
//        check_point(value[2], 16.0);
//    }

//    SECTION("operator/=")
//    {
//    	AffineVector value = initial;
//        value /= rhs;
//
//        check_point(value[0], 0.5);
//        check_point(value[1], 2.0 / 3.0);
//        check_point(value[2], 1.0);
//    }
}


TEST_CASE("Eigen AffineMatrix compound assignments with Interval scalar",
          "[Eigen][compound-assignment][AffineMatrix][Interval]")
{
    AffineMatrix initial(2, 2);
    initial(0, 0) = Interval(1.0);
    initial(0, 1) = Interval(2.0);
    initial(1, 0) = Interval(3.0);
    initial(1, 1) = Interval(4.0);

    SECTION("operator+=")
    {
        auto value = initial;
        value += Interval(2.0);

        check_point(value(0, 0), 3.0);
        check_point(value(0, 1), 4.0);
        check_point(value(1, 0), 5.0);
        check_point(value(1, 1), 6.0);
    }

    SECTION("operator-=")
    {
        auto value = initial;
        value -= Interval(2.0);

        check_point(value(0, 0), -1.0);
        check_point(value(0, 1), 0.0);
        check_point(value(1, 0), 1.0);
        check_point(value(1, 1), 2.0);
    }

    SECTION("operator*=")
    {
        auto value = initial;
        value *= Interval(2.0);

        check_point(value(0, 0), 2.0);
        check_point(value(0, 1), 4.0);
        check_point(value(1, 0), 6.0);
        check_point(value(1, 1), 8.0);
    }

    SECTION("operator/=")
    {
        auto value = initial;
        value /= Interval(2.0);

        check_point(value(0, 0), 0.5);
        check_point(value(0, 1), 1.0);
        check_point(value(1, 0), 1.5);
        check_point(value(1, 1), 2.0);
    }
}


TEST_CASE("Eigen AffineMatrix compound assignments with Affine scalar",
          "[Eigen][compound-assignment][AffineMatrix][Affine]")
{
    AffineMatrix initial(2, 2);
    initial(0, 0) = Interval(1.0);
    initial(0, 1) = Interval(2.0);
    initial(1, 0) = Interval(3.0);
    initial(1, 1) = Interval(4.0);

    Affine scalar;
	scalar =Interval(2.0);

    SECTION("operator+=")
    {
        auto value = initial;
        value += scalar;

        check_point(value(0, 0), 3.0);
        check_point(value(0, 1), 4.0);
        check_point(value(1, 0), 5.0);
        check_point(value(1, 1), 6.0);
    }

    SECTION("operator-=")
    {
        auto value = initial;
        value -= scalar;

        check_point(value(0, 0), -1.0);
        check_point(value(0, 1), 0.0);
        check_point(value(1, 0), 1.0);
        check_point(value(1, 1), 2.0);
    }

    SECTION("operator*=")
    {
        auto value = initial;
        value *= scalar;

        check_point(value(0, 0), 2.0);
        check_point(value(0, 1), 4.0);
        check_point(value(1, 0), 6.0);
        check_point(value(1, 1), 8.0);
    }

    SECTION("operator/=")
    {
        auto value = initial;
        value /= scalar;

        check_point(value(0, 0), 0.5);
        check_point(value(0, 1), 1.0);
        check_point(value(1, 0), 1.5);
        check_point(value(1, 1), 2.0);
    }
}

TEST_CASE("Eigen AffineMatrix compound assignments with IntervalMatrix",
          "[Eigen][compound-assignment][AffineMatrix][IntervalMatrix]")
{
    AffineMatrix initial(2, 2);
    initial(0, 0) = Interval(1.0);
    initial(0, 1) = Interval(2.0);
    initial(1, 0) = Interval(3.0);
    initial(1, 1) = Interval(4.0);

    // Use a mutable copy because IntervalMatrix's coefficient access is
    // intentionally exercised here as part of the test.
    IntervalMatrix other(2, 2);
    other(0, 0) = Interval(2.0);
    other(0, 1) = Interval(3.0);
    other(1, 0) = Interval(4.0);
    other(1, 1) = Interval(5.0);

    SECTION("operator+=")
    {
        auto value = initial;
        value += other;

        check_point(value(0, 0), 3.0);
        check_point(value(0, 1), 5.0);
        check_point(value(1, 0), 7.0);
        check_point(value(1, 1), 9.0);
    }

    SECTION("operator-=")
    {
        auto value = initial;
        value -= other;

        check_point(value(0, 0), -1.0);
        check_point(value(0, 1), -1.0);
        check_point(value(1, 0), -1.0);
        check_point(value(1, 1), -1.0);
    }

    SECTION("operator*=")
    {
        auto value = initial;
        value *= other;

        check_point(value(0, 0), 10.0);
        check_point(value(0, 1), 13.0);
        check_point(value(1, 0), 22.0);
        check_point(value(1, 1), 29.0);
    }

//    SECTION("operator/=")
//    {
//        auto value = initial;
//        value /= other;
//
//        check_point(value(0, 0), 0.5);
//        check_point(value(0, 1), 2.0 / 3.0);
//        check_point(value(1, 0), 0.75);
//        check_point(value(1, 1), 0.8);
//    }
}


TEST_CASE("Eigen IntervalMatrix compound assignments with IntervalMatrix",
          "[Eigen][compound-assignment][AffineMatrix][IntervalMatrix]")
{
    IntervalMatrix initial(2, 2);
    initial(0, 0) = Interval(1.0);
    initial(0, 1) = Interval(2.0);
    initial(1, 0) = Interval(3.0);
    initial(1, 1) = Interval(4.0);

    // Use a mutable copy because IntervalMatrix's coefficient access is
    // intentionally exercised here as part of the test.
    IntervalMatrix other(2, 2);
    other(0, 0) = Interval(2.0);
    other(0, 1) = Interval(3.0);
    other(1, 0) = Interval(4.0);
    other(1, 1) = Interval(5.0);

    SECTION("operator+=")
    {
        auto value = initial;
        value += other;

        CHECK(value(0, 0)==Interval(3.0));
        CHECK(value(0, 1)==Interval(5.0));
        CHECK(value(1, 0)==Interval(7.0));
        CHECK(value(1, 1)==Interval(9.0));
    }

    SECTION("operator-=")
    {
        auto value = initial;
        value -= other;

        CHECK(value(0, 0)==Interval(-1.0));
        CHECK(value(0, 1)==Interval(-1.0));
        CHECK(value(1, 0)==Interval(-1.0));
        CHECK(value(1, 1)==Interval(-1.0));
    }

    SECTION("operator*=")
    {
        auto value = initial;
        value *= other;

        CHECK(value(0, 0)==Interval(10.0));
        CHECK(value(0, 1)==Interval(13.0));
        CHECK(value(1, 0)==Interval(22.0));
        CHECK(value(1, 1)==Interval(29.0));
    }

//    SECTION("operator/=")
//    {
//        auto value = initial;
//        value /= other;
//
//        check_point(value(0, 0), 0.5);
//        check_point(value(0, 1), 2.0 / 3.0);
//        check_point(value(1, 0), 0.75);
//        check_point(value(1, 1), 0.8);
//    }
}


TEST_CASE("Eigen AffineMatrix compound assignments with AffineMatrix",
          "[Eigen][compound-assignment][AffineMatrix][AffineMatrix]")
{
    AffineMatrix initial(2, 2);
    initial(0, 0) = Interval(1.0);
    initial(0, 1) = Interval(2.0);
    initial(1, 0) = Interval(3.0);
    initial(1, 1) = Interval(4.0);

    AffineMatrix rhs(2, 2);
    rhs(0, 0) = Interval(2.0);
    rhs(0, 1) = Interval(3.0);
    rhs(1, 0) = Interval(4.0);
    rhs(1, 1) = Interval(5.0);

    SECTION("operator+=")
    {
        auto value = initial;
        value += rhs;

        check_point(value(0, 0), 3.0);
        check_point(value(0, 1), 5.0);
        check_point(value(1, 0), 7.0);
        check_point(value(1, 1), 9.0);
    }

    SECTION("operator-=")
    {
        auto value = initial;
        value -= rhs;

        check_point(value(0, 0), -1.0);
        check_point(value(0, 1), -1.0);
        check_point(value(1, 0), -1.0);
        check_point(value(1, 1), -1.0);
    }

    SECTION("operator*=")
    {
        auto value = initial;
        value *= rhs;

        check_point(value(0, 0), 10.0);
        check_point(value(0, 1), 13.0);
        check_point(value(1, 0), 22.0);
        check_point(value(1, 1), 29.0);
    }

//    SECTION("operator/=")
//    {
//        auto value = initial;
//        value /= rhs;
//
//        check_point(value(0, 0), 0.5);
//        check_point(value(0, 1), 2.0 / 3.0);
//        check_point(value(1, 0), 0.75);
//        check_point(value(1, 1), 0.8);
//    }
}
