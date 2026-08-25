/**
 * \file codac2_tests_AffineEigen.cpp
 * \brief Unit tests for Eigen operations applied to Affine forms.
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
using AffineT = AffineMain<Model>;
using AffineTMatrix = AffineMainMatrix<Model>;
using AffineTVector = AffineMainVector<Model>;
using AffineTVarVector = AffineVarMainVector<Model>;
using AffineTRow = AffineMainRow<Model>;

using DynamicAffineTMatrix =  Eigen::Matrix<AffineT, Eigen::Dynamic, Eigen::Dynamic>;
using FixedAffineTMatrix =    Eigen::Matrix<AffineT, 2, 3>;
using FixedAffineTVector =    Eigen::Matrix<AffineT, 3, 1>;
using FixedAffineTRow =       Eigen::Matrix<AffineT, 1, 3>;


const double ERROR = std::numeric_limits<double>::epsilon()*100; 

void check_point(const AffineT& actual, double expected)
{
    CAPTURE(actual.itv(), expected);
    CHECK(actual == Approx<AffineT>(Interval(expected),ERROR));
}

void check_interval(const AffineT& actual, const Interval expected)
{
    CAPTURE(actual.itv(), expected);
    CHECK(actual == Approx<AffineT>(expected,ERROR));
}

void check_interval_enclosure(
    const AffineT& actual,
    const Interval& expected)
{
    CAPTURE(actual.itv(), expected);
    CHECK(actual.itv().is_superset(expected));
}

AffineTMatrix make_point_matrix_3x4()
{
    AffineTMatrix matrix(3, 4);
    matrix(0, 0) = Interval(1.0);  matrix(0, 1) = Interval(2.0);
    matrix(0, 2) = Interval(3.0);  matrix(0, 3) = Interval(4.0);
    matrix(1, 0) = Interval(5.0);  matrix(1, 1) = Interval(6.0);
    matrix(1, 2) = Interval(7.0);  matrix(1, 3) = Interval(8.0);
    matrix(2, 0) = Interval(9.0);  matrix(2, 1) = Interval(10.0);
    matrix(2, 2) = Interval(11.0); matrix(2, 3) = Interval(12.0);
    return matrix;
}

AffineTVector make_point_vector_5()
{
    const AffineTVarVector variables(
        IntervalVector({{1}, {2}, {3}, {4}, {5}}));
    return AffineTVector(variables);
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

static_assert(std::same_as<typename AffineTMatrix::Scalar, AffineT>);
static_assert(std::same_as<typename AffineTVector::Scalar, AffineT>);
static_assert(std::same_as<typename AffineTRow::Scalar, AffineT>);
static_assert(AffineTVector::ColsAtCompileTime == 1);
static_assert(AffineTRow::RowsAtCompileTime == 1);

// AffineTVarMain deliberately deletes scalar compound assignments.
// Keep this assertion here so that this API contract is explicit and
// cannot accidentally become enabled by the Eigen plugin.
template<class T> concept has_scalar_plus_eq  = requires(T& v, double d) { v += d; };
template<class T> concept has_scalar_minus_eq = requires(T& v, double d) { v -= d; };
template<class T> concept has_scalar_times_eq = requires(T& v, double d) { v *= d; };
template<class T> concept has_scalar_div_eq   = requires(T& v, double d) { v /= d; };

TEST_CASE("Eigen AffineT matrix dimensions and coefficient access")
{
    AffineTMatrix matrix = make_point_matrix_3x4();

    CHECK(matrix.rows() == 3);
    CHECK(matrix.cols() == 4);
    CHECK(matrix.size() == 12);

    check_point(matrix(0, 0), 1.0);
    check_point(matrix(1, 2), 7.0);
    check_point(matrix.coeff(2, 3), 12.0);

    matrix.coeffRef(1, 1) = -6.0;
    check_point(matrix(1, 1), -6.0);
}

TEST_CASE("Eigen AffineT vector indexing")
{
    AffineTVector vector = make_point_vector_5();

    CHECK(vector.rows() == 5);
    CHECK(vector.cols() == 1);
    CHECK(vector.size() == 5);

    check_point(vector[0], 1.0);
    check_point(vector(2), 3.0);
    check_point(vector.coeff(4), 5.0);

    vector.coeffRef(3) = -4.0;
    check_point(vector[3], -4.0);
}

TEST_CASE("Eigen AffineT matrix row and column views")
{
    AffineTMatrix matrix = make_point_matrix_3x4();

    const AffineTVector column = matrix.col(2);
    REQUIRE(column.size() == 3);
    check_point(column[0], 3.0);
    check_point(column[1], 7.0);
    check_point(column[2], 11.0);

    const AffineTRow row = matrix.row(1);
    REQUIRE(row.size() == 4);
    check_point(row[0], 5.0);
    check_point(row[1], 6.0);
    check_point(row[2], 7.0);
    check_point(row[3], 8.0);
}

TEST_CASE("Eigen AffineT row and column assignments")
{
    AffineTMatrix matrix(3, 3);
    matrix.setZero();

    const AffineTVarVector column_variables(
        IntervalVector({{1}, {2}, {3}}));
    const AffineTVector column(column_variables);
    matrix.col(1) = column;

    check_point(matrix(0, 1), 1.0);
    check_point(matrix(1, 1), 2.0);
    check_point(matrix(2, 1), 3.0);

    const AffineTVarVector row_variables(
        IntervalVector({{4}, {5}, {6}}));
    const AffineTRow row = AffineTVector(row_variables).transpose();
    matrix.row(2) = row;

    check_point(matrix(2, 0), 4.0);
    check_point(matrix(2, 1), 5.0);
    check_point(matrix(2, 2), 6.0);
}

TEST_CASE("Eigen AffineTVarMainVector assignment to matrix column")
{
    AffineTMatrix matrix(3, 2);
    matrix.setZero();

    const AffineTVarVector variables(
        IntervalVector({{-2.0, 1.0}, {3.0, 4.0}, {5.0, 7.0}}));

    matrix.col(0) = variables;

    for(Index i = 0; i < variables.size(); ++i)
        CHECK(matrix(i, 0).itv() == variables[i].itv());
}

TEST_CASE("Eigen AffineT block extraction")
{
    const AffineTMatrix matrix = make_point_matrix_3x4();
    const AffineTMatrix block = matrix.block(1, 1, 2, 3);

    REQUIRE(block.rows() == 2);
    REQUIRE(block.cols() == 3);

    check_point(block(0, 0), 6.0);
    check_point(block(0, 2), 8.0);
    check_point(block(1, 0), 10.0);
    check_point(block(1, 2), 12.0);
}

TEST_CASE("Eigen AffineT block assignment")
{
    AffineTMatrix matrix(4, 5);
    matrix.setZero();

    AffineTMatrix block(2, 3);
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

TEST_CASE("Eigen AffineT top bottom left and right corners")
{
    const AffineTMatrix matrix = make_point_matrix_3x4();

    const AffineTMatrix top_left = matrix.topLeftCorner(2, 2);
    const AffineTMatrix bottom_right = matrix.bottomRightCorner(2, 2);

    check_point(top_left(0, 0), 1.0);
    check_point(top_left(1, 1), 6.0);
    check_point(bottom_right(0, 0), 7.0);
    check_point(bottom_right(1, 1), 12.0);
}

TEST_CASE("Eigen AffineT middle rows and columns")
{
    const AffineTMatrix matrix = make_point_matrix_3x4();

    const AffineTMatrix middle_columns = matrix.middleCols(1, 2);
    REQUIRE(middle_columns.rows() == 3);
    REQUIRE(middle_columns.cols() == 2);
    check_point(middle_columns(0, 0), 2.0);
    check_point(middle_columns(2, 1), 11.0);

    const AffineTMatrix middle_rows = matrix.middleRows(1, 2);
    REQUIRE(middle_rows.rows() == 2);
    REQUIRE(middle_rows.cols() == 4);
    check_point(middle_rows(0, 0), 5.0);
    check_point(middle_rows(1, 3), 12.0);
}

TEST_CASE("Eigen AffineT vector segment head and tail")
{
    AffineTVector vector = make_point_vector_5();

    const AffineTVector segment = vector.segment(1, 3);
    const AffineTVector head = vector.head(2);
    const AffineTVector tail = vector.tail(2);

    REQUIRE(segment.size() == 3);
    check_point(segment[0], 2.0);
    check_point(segment[1], 3.0);
    check_point(segment[2], 4.0);

    check_point(head[0], 1.0);
    check_point(head[1], 2.0);
    check_point(tail[0], 4.0);
    check_point(tail[1], 5.0);
}

TEST_CASE("Eigen AffineTVarMainVector segment conversion")
{
    const AffineTVarVector variables = make_variable_vector_5();
    const AffineTVector segment = variables.segment(1, 3);

    REQUIRE(segment.size() == 3);
    CHECK(segment[0].itv() == variables[1].itv());
    CHECK(segment[1].itv() == variables[2].itv());
    CHECK(segment[2].itv() == variables[3].itv());
}

TEST_CASE("Eigen AffineT segment assignment")
{
    AffineTVector vector(6);
    vector.setZero();

    const AffineTVarVector source_variables(
        IntervalVector({{2}, {4}, {6}}));
    const AffineTVector source(source_variables);
    vector.segment(2, 3) = source;

    check_point(vector[0], 0.0);
    check_point(vector[1], 0.0);
    check_point(vector[2], 2.0);
    check_point(vector[3], 4.0);
    check_point(vector[4], 6.0);
    check_point(vector[5], 0.0);
}

TEST_CASE("Eigen AffineT transpose")
{
    const AffineTMatrix matrix = make_point_matrix_3x4();
    const AffineTMatrix transpose = matrix.transpose();

    REQUIRE(transpose.rows() == 4);
    REQUIRE(transpose.cols() == 3);

    check_point(transpose(0, 0), 1.0);
    check_point(transpose(2, 1), 7.0);
    check_point(transpose(3, 2), 12.0);
}

TEST_CASE("Eigen AffineT vector transpose to row")
{
    const AffineTVector vector = make_point_vector_5();
    const AffineTRow row = vector.transpose();

    REQUIRE(row.rows() == 1);
    REQUIRE(row.cols() == 5);
    check_point(row[0], 1.0);
    check_point(row[4], 5.0);
}

TEST_CASE("Eigen AffineT in-place transpose for square matrix")
{
    AffineTMatrix matrix(3, 3);
    matrix(0, 0) = Interval(1.0); matrix(0, 1) = Interval(2.0); matrix(0, 2) = Interval(3.0);
    matrix(1, 0) = Interval(4.0); matrix(1, 1) = Interval(5.0); matrix(1, 2) = Interval(6.0);
    matrix(2, 0) = Interval(7.0); matrix(2, 1) = Interval(8.0); matrix(2, 2) = Interval(9.0);

    matrix.transposeInPlace();

    check_point(matrix(0, 1), 4.0);
    check_point(matrix(1, 0), 2.0);
    check_point(matrix(2, 1), 6.0);
}

TEST_CASE("Eigen AffineT diagonal extraction")
{
    AffineTMatrix matrix(3, 3);
    matrix(0, 0) = Interval(1.0); matrix(0, 1) = Interval(2.0); matrix(0, 2) = Interval(3.0);
    matrix(1, 0) = Interval(4.0); matrix(1, 1) = Interval(5.0); matrix(1, 2) = Interval(6.0);
    matrix(2, 0) = Interval(7.0); matrix(2, 1) = Interval(8.0); matrix(2, 2) = Interval(9.0);

    const AffineTVector diagonal = matrix.diagonal();
    REQUIRE(diagonal.size() == 3);
    check_point(diagonal[0], 1.0);
    check_point(diagonal[1], 5.0);
    check_point(diagonal[2], 9.0);
}

TEST_CASE("Eigen AffineT diagonal assignment")
{
    AffineTMatrix matrix(3, 3);
    matrix.setZero();

    const AffineTVarVector diagonal_variables(
        IntervalVector({{2}, {4}, {6}}));
    const AffineTVector diagonal(diagonal_variables);
    matrix.diagonal() = diagonal;

    check_point(matrix(0, 0), 2.0);
    check_point(matrix(1, 1), 4.0);
    check_point(matrix(2, 2), 6.0);
    check_point(matrix(0, 1), 0.0);
}

TEST_CASE("Eigen AffineT setZero setOnes and setConstant")
{
    AffineTMatrix matrix(2, 3);

    matrix.setZero();
    for(Index i = 0; i < matrix.rows(); ++i)
        for(Index j = 0; j < matrix.cols(); ++j)
            check_point(matrix(i, j), 0.0);

    matrix.setOnes();
    for(Index i = 0; i < matrix.rows(); ++i)
        for(Index j = 0; j < matrix.cols(); ++j)
            check_point(matrix(i, j), 1.0);

    matrix.setConstant(AffineT(3.5));
    for(Index i = 0; i < matrix.rows(); ++i)
        for(Index j = 0; j < matrix.cols(); ++j)
            check_point(matrix(i, j), 3.5);
}

TEST_CASE("Eigen AffineT Zero Ones and Constant factories")
{
    const AffineTMatrix zero = AffineTMatrix::Zero(2, 2);
    const AffineTMatrix ones = AffineTMatrix::Ones(2, 2);
    const AffineTMatrix constant = AffineTMatrix::Constant(2, 2, AffineT(-3.0));

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

TEST_CASE("Eigen AffineT Identity factory")
{
    const AffineTMatrix identity = AffineTMatrix::Identity(3, 3);

    for(Index i = 0; i < identity.rows(); ++i)
    {
        for(Index j = 0; j < identity.cols(); ++j)
            check_point(identity(i, j), i == j ? 1.0 : 0.0);
    }
}

TEST_CASE("Eigen AffineT resize and conservativeResize")
{
    AffineTMatrix matrix(2, 2);
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

TEST_CASE("Eigen AffineT swap")
{
    const AffineTVarVector lhs_variables(
        IntervalVector({{1}, {2}, {3}}));
    const AffineTVarVector rhs_variables(
        IntervalVector({{4}, {5}, {6}}));
    AffineTVector lhs(lhs_variables);
    AffineTVector rhs(rhs_variables);

    lhs.swap(rhs);

    check_point(lhs[0], 4.0);
    check_point(lhs[2], 6.0);
    check_point(rhs[0], 1.0);
    check_point(rhs[2], 3.0);
}

TEST_CASE("Eigen AffineT arithmetic expressions")
{
    AffineTMatrix lhs(2, 2);
    AffineTMatrix rhs(2, 2);
    lhs(0, 0) = Interval(1.0); lhs(0, 1) = Interval(2.0);
    lhs(1, 0) = Interval(3.0); lhs(1, 1) = Interval(4.0);
    rhs(0, 0) = Interval(5.0); rhs(0, 1) = Interval(6.0);
    rhs(1, 0) = Interval(7.0); rhs(1, 1) = Interval(8.0);

    const AffineTMatrix sum = lhs + rhs;
    const AffineTMatrix difference = rhs - lhs;
    const AffineTMatrix expression = 2.0 * lhs - rhs;

    check_point(sum(0, 0), 6.0);
    check_point(sum(1, 1), 12.0);
    check_point(difference(0, 1), 4.0);
    check_point(difference(1, 0), 4.0);
    check_point(expression(0, 0), -3.0);
    check_point(expression(1, 1), 0.0);
}

TEST_CASE("Eigen AffineT compound assignment")
{
    const AffineTVarVector value_variables(
        IntervalVector({{1}, {2}, {3}}));
    const AffineTVarVector other_variables(
        IntervalVector({{4}, {5}, {6}}));
    AffineTVector value(value_variables);
    const AffineTVector other(other_variables);

    value += other;
    check_point(value[0], 5.0);
    check_point(value[2], 9.0);

    value -= other;
    check_point(value[0], 1.0);
    check_point(value[2], 3.0);

    value.array() += 3.0;
    check_point(value[0], 4.0);
    check_point(value[2], 6.0);

    value.array() -= 5;
    check_point(value[0], -1.0);
    check_point(value[2], 1.0);

    value.array() *= 2.0;
    check_point(value[0], -2.0);
    check_point(value[2], 2.0);

    value.array() /= 2.0;
    check_point(value[0], -1.0);
    check_point(value[2], 1.0);
}

TEST_CASE("Eigen AffineTVar compound assignment")
{
    const AffineTVarVector value_variables(
        IntervalVector({{1,2}, {-2,5}, {3,6}}));
    AffineTVector value(value_variables);


    value -= value_variables;
    check_point(value[0], 0);
    check_point(value[1], 0);
    check_point(value[2], 0);

    value += value_variables;
    check_interval(value[0], value_variables[0].itv());
    check_interval(value[1], value_variables[1].itv());
    check_interval(value[2], value_variables[2].itv());

}

// ============================================================================
// Scalar compound assignments on all supported Eigen/CODAC containers
// ============================================================================
//
// These tests intentionally exercise +=, -=, *= and /= independently.
// The existing "Eigen AffineT compound assignment" test covers the basic
// AffineTVector path; the tests below make sure the Eigen plugin also works
// for matrices, rows, interval containers and fixed-size Eigen types.
//

TEST_CASE("Eigen scalar compound assignments - IntervalVector",
          "[Eigen][compound-assignment][IntervalVector]")
{
    IntervalVector value({{1.0}, {-2.0}, {4.0}});

    value.array() += 2.0;
    CHECK(value[0] == Interval(3.0));
    CHECK(value[1] == Interval(0.0));
    CHECK(value[2] == Interval(6.0));

    value.array() -= 2.0;
    CHECK(value[0] == Interval(1.0));
    CHECK(value[1] == Interval(-2.0));
    CHECK(value[2] == Interval(4.0));

    value.array() *= 2.0;
    CHECK(value[0] == Interval(2.0));
    CHECK(value[1] == Interval(-4.0));
    CHECK(value[2] == Interval(8.0));

    value.array() /= 2.0;
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

    value.array() += 2.0;
    CHECK(value[0] == Interval(3.0));
    CHECK(value[1] == Interval(0.0));
    CHECK(value[2] == Interval(6.0));

    value.array() -= 2.0;
    CHECK(value[0] == Interval(1.0));
    CHECK(value[1] == Interval(-2.0));
    CHECK(value[2] == Interval(4.0));

    value.array() *= 2.0;
    CHECK(value[0] == Interval(2.0));
    CHECK(value[1] == Interval(-4.0));
    CHECK(value[2] == Interval(8.0));

    value.array() /= 2.0;
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

    value.array() += 2.0;
    CHECK(value(0, 0) == Interval(3.0));
    CHECK(value(0, 1) == Interval(0.0));
    CHECK(value(1, 2) == Interval(-4.0));

    value.array() -= 2.0;
    CHECK(value(0, 0) == Interval(1.0));
    CHECK(value(0, 1) == Interval(-2.0));
    CHECK(value(1, 2) == Interval(-6.0));

    value.array() *= 2.0;
    CHECK(value(0, 0) == Interval(2.0));
    CHECK(value(0, 1) == Interval(-4.0));
    CHECK(value(1, 2) == Interval(-12.0));

    value.array() /= 2.0;
    CHECK(value(0, 0) == Interval(1.0));
    CHECK(value(0, 1) == Interval(-2.0));
    CHECK(value(1, 2) == Interval(-6.0));
}

TEST_CASE("Eigen scalar compound assignments - AffineTVector",
          "[Eigen][compound-assignment][AffineTVector]")
{
    AffineTVector value(3);
    value[0] = Interval(1.0);
    value[1] = Interval(-2.0);
    value[2] = Interval(4.0);

    value.array() += 2.0;
    check_point(value[0], 3.0);
    check_point(value[1], 0.0);
    check_point(value[2], 6.0);

    value.array() -= 2.0;
    check_point(value[0], 1.0);
    check_point(value[1], -2.0);
    check_point(value[2], 4.0);

    value.array() *= 2.0;
    check_point(value[0], 2.0);
    check_point(value[1], -4.0);
    check_point(value[2], 8.0);

    value.array() /= 2.0;
    check_point(value[0], 1.0);
    check_point(value[1], -2.0);
    check_point(value[2], 4.0);
}

TEST_CASE("Eigen scalar compound assignments - AffineTRow",
          "[Eigen][compound-assignment][AffineTRow]")
{
    AffineTRow value(3);
    value[0] = Interval(1.0);
    value[1] = Interval(-2.0);
    value[2] = Interval(4.0);

    value.array() += 2.0;
    check_point(value[0], 3.0);
    check_point(value[1], 0.0);
    check_point(value[2], 6.0);

    value.array() -= 2.0;
    check_point(value[0], 1.0);
    check_point(value[1], -2.0);
    check_point(value[2], 4.0);

    value.array() *= 2.0;
    check_point(value[0], 2.0);
    check_point(value[1], -4.0);
    check_point(value[2], 8.0);

    value.array() /= 2.0;
    check_point(value[0], 1.0);
    check_point(value[1], -2.0);
    check_point(value[2], 4.0);
}

TEST_CASE("Eigen scalar compound assignments - AffineTMatrix",
          "[Eigen][compound-assignment][AffineTMatrix]")
{
    AffineTMatrix value(2, 3);
    IntervalMatrix m({{ Interval(1.0), Interval(-2.0), Interval(3.0)},{
             Interval(4.0), Interval(5.0), Interval(-6.0)}});
    value = m;

    value.array() += 2.0;
    check_point(value(0, 0), 3.0);
    check_point(value(0, 1), 0.0);
    check_point(value(1, 2), -4.0);

    value.array() -= 2.0;
    check_point(value(0, 0), 1.0);
    check_point(value(0, 1), -2.0);
    check_point(value(1, 2), -6.0);

    value.array() *= 2.0;
    check_point(value(0, 0), 2.0);
    check_point(value(0, 1), -4.0);
    check_point(value(1, 2), -12.0);

    value.array() /= 2.0;
    check_point(value(0, 0), 1.0);
    check_point(value(0, 1), -2.0);
    check_point(value(1, 2), -6.0);
}

TEST_CASE("Eigen scalar compound assignments - fixed-size AffineT types",
          "[Eigen][compound-assignment][AffineT][fixed-size]")
{
    FixedAffineTVector vector;

    IntervalVector v({ Interval(1.0), Interval(-2.0), Interval(4.0)});
	vector = v;

    vector.array() += 2.0;
    check_point(vector[0], 3.0);
    check_point(vector[1], 0.0);
    check_point(vector[2], 6.0);

    vector.array() -= 2.0;
    check_point(vector[0], 1.0);
    check_point(vector[1], -2.0);
    check_point(vector[2], 4.0);

    vector.array() *= 2.0;
    check_point(vector[0], 2.0);
    check_point(vector[1], -4.0);
    check_point(vector[2], 8.0);

    vector.array() /= 2.0;
    check_point(vector[0], 1.0);
    check_point(vector[1], -2.0);
    check_point(vector[2], 4.0);

    FixedAffineTRow row;
    IntervalVector vv({ Interval(1.0), Interval(-2.0), Interval(4.0)});
    row = vv.transpose();

    row.array() += 2.0;
    check_point(row[0], 3.0);
    check_point(row[1], 0.0);
    check_point(row[2], 6.0);

    row.array() -= 2.0;
    check_point(row[0], 1.0);
    check_point(row[1], -2.0);
    check_point(row[2], 4.0);

    row.array() *= 2.0;
    check_point(row[0], 2.0);
    check_point(row[1], -4.0);
    check_point(row[2], 8.0);

    row.array() /= 2.0;
    check_point(row[0], 1.0);
    check_point(row[1], -2.0);
    check_point(row[2], 4.0);

    FixedAffineTMatrix matrix;
    IntervalMatrix mm ({{Interval(1.0), Interval(-2.0), Interval(3.0)},
              {Interval(4.0), Interval(5.0), Interval(-6.0)}});
    matrix = mm;

    matrix.array() += 2.0;
    check_point(matrix(0, 0), 3.0);
    check_point(matrix(0, 1), 0.0);
    check_point(matrix(1, 2), -4.0);

    matrix.array() -= 2.0;
    check_point(matrix(0, 0), 1.0);
    check_point(matrix(0, 1), -2.0);
    check_point(matrix(1, 2), -6.0);

    matrix.array() *= 2.0;
    check_point(matrix(0, 0), 2.0);
    check_point(matrix(0, 1), -4.0);
    check_point(matrix(1, 2), -12.0);

    matrix.array() /= 2.0;
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

    vector.array() += 2.0;
    CHECK(vector[0] == Interval(3.0));
    CHECK(vector[1] == Interval(0.0));
    CHECK(vector[2] == Interval(6.0));

    vector.array() -= 2.0;
    CHECK(vector[0] == Interval(1.0));
    CHECK(vector[1] == Interval(-2.0));
    CHECK(vector[2] == Interval(4.0));

    vector.array() *= 2.0;
    CHECK(vector[0] == Interval(2.0));
    CHECK(vector[1] == Interval(-4.0));
    CHECK(vector[2] == Interval(8.0));

    vector.array() /= 2.0;
    CHECK(vector[0] == Interval(1.0));
    CHECK(vector[1] == Interval(-2.0));
    CHECK(vector[2] == Interval(4.0));

    Eigen::Matrix<Interval, 1, 3> row;
    IntervalVector vv({Interval(1.0), Interval(-2.0), Interval(4.0)});
    row = vv.transpose();

    row.array() += 2.0;
    CHECK(row[0] == Interval(3.0));
    CHECK(row[1] == Interval(0.0));
    CHECK(row[2] == Interval(6.0));

    row.array() -= 2.0;
    CHECK(row[0] == Interval(1.0));
    CHECK(row[1] == Interval(-2.0));
    CHECK(row[2] == Interval(4.0));

    row.array() *= 2.0;
    CHECK(row[0] == Interval(2.0));
    CHECK(row[1] == Interval(-4.0));
    CHECK(row[2] == Interval(8.0));

    row.array() /= 2.0;
    CHECK(row[0] == Interval(1.0));
    CHECK(row[1] == Interval(-2.0));
    CHECK(row[2] == Interval(4.0));

    Eigen::Matrix<Interval, 2, 3> matrix;
    IntervalMatrix mm({{ Interval(1.0), Interval(-2.0), Interval(3.0)},{
              Interval(4.0), Interval(5.0), Interval(-6.0)}});
    matrix = mm;

    matrix.array() += 2.0;
    CHECK(matrix(0, 0) == Interval(3.0));
    CHECK(matrix(0, 1) == Interval(0.0));
    CHECK(matrix(1, 2) == Interval(-4.0));

    matrix.array() -= 2.0;
    CHECK(matrix(0, 0) == Interval(1.0));
    CHECK(matrix(0, 1) == Interval(-2.0));
    CHECK(matrix(1, 2) == Interval(-6.0));

    matrix.array() *= 2.0;
    CHECK(matrix(0, 0) == Interval(2.0));
    CHECK(matrix(0, 1) == Interval(-4.0));
    CHECK(matrix(1, 2) == Interval(-12.0));

    matrix.array() /= 2.0;
    CHECK(matrix(0, 0) == Interval(1.0));
    CHECK(matrix(0, 1) == Interval(-2.0));
    CHECK(matrix(1, 2) == Interval(-6.0));
}

TEST_CASE("Eigen scalar compound assignments - AffineTVarMainVector are unsupported",
          "[Eigen][compound-assignment][AffineTVarMainVector]")
{


	static_assert(!has_scalar_plus_eq<AffineTVarVector>);
	static_assert(!has_scalar_minus_eq<AffineTVarVector>);
	static_assert(!has_scalar_times_eq<AffineTVarVector>);
	static_assert(!has_scalar_div_eq<AffineTVarVector>);
}

TEST_CASE("Eigen AffineT unary minus")
{
    const AffineTVector vector = make_point_vector_5();
    const AffineTVector negative = -vector;

    check_point(negative[0], -1.0);
    check_point(negative[2], -3.0);
    check_point(negative[4], -5.0);
}

TEST_CASE("Eigen AffineT cwise product")
{
    const AffineTVarVector lhs_variables(
        IntervalVector({{1}, {-2}, {3}}));
    const AffineTVarVector rhs_variables(
        IntervalVector({{4}, {5}, {-6}}));
    const AffineTVector lhs(lhs_variables);
    const AffineTVector rhs(rhs_variables);

    const AffineTVector result = lhs.cwiseProduct(rhs);
    check_point(result[0], 4.0);
    check_point(result[1], -10.0);
    check_point(result[2], -18.0);
}

TEST_CASE("Eigen AffineT cwise quotient by point AffineT values")
{
    const AffineTVarVector numerator_variables(
        IntervalVector({{4}, {-10}, {18}}));
    const AffineTVarVector denominator_variables(
        IntervalVector({{2}, {5}, {-3}}));
    const AffineTVector numerator(numerator_variables);
    const AffineTVector denominator(denominator_variables);

    const AffineTVector result = numerator.cwiseQuotient(denominator);
    check_point(result[0], 2.0);
    check_point(result[1], -2.0);
    check_point(result[2], -6.0);
}

TEST_CASE("Eigen AffineT matrix product")
{
    AffineTMatrix lhs(2, 3);
    AffineTMatrix rhs(3, 2);
    lhs(0, 0) = Interval(1.0); lhs(0, 1) = Interval(2.0); lhs(0, 2) = Interval(3.0);
    lhs(1, 0) = Interval(4.0); lhs(1, 1) = Interval(5.0); lhs(1, 2) = Interval(6.0);
    rhs(0, 0) = Interval(7.0);  rhs(0, 1) = Interval(8.0);
    rhs(1, 0) = Interval(9.0);  rhs(1, 1) = Interval(10.0);
    rhs(2, 0) = Interval(11.0); rhs(2, 1) = Interval(12.0);

    const AffineTMatrix result = lhs * rhs;
    check_point(result(0, 0), 58.0);
    check_point(result(0, 1), 64.0);
    check_point(result(1, 0), 139.0);
    check_point(result(1, 1), 154.0);
}

TEST_CASE("Eigen AffineT matrix vector product")
{
    AffineTMatrix matrix(2, 3);
    matrix(0, 0) = Interval(1.0); matrix(0, 1) = Interval(2.0); matrix(0, 2) = Interval(3.0);
    matrix(1, 0) = Interval(4.0); matrix(1, 1) = Interval(5.0); matrix(1, 2) = Interval(6.0);

    const AffineTVarVector vector_variables(
        IntervalVector({{2}, {-1}, {3}}));
    const AffineTVector vector(vector_variables);

    const AffineTVector result = matrix * vector;
    check_point(result[0], 9.0);
    check_point(result[1], 21.0);
}

TEST_CASE("Eigen real matrix times AffineTVarMainVector")
{
    IntervalMatrix matrix(2,3);
    matrix(0, 0) = Interval(1.0); matrix(0, 1) = Interval(2.0); matrix(0, 2) = Interval(3.0);
    matrix(1, 0) = Interval(4.0); matrix(1, 1) = Interval(5.0); matrix(1, 2) = Interval(6.0);

    const AffineTVarVector variables(Vector({2.0, -1.0, 3.0}));
    const auto result = matrix * variables;

    static_assert(std::same_as<typename decltype(result)::Scalar, AffineT>);
    check_point(result[0], 9.0);
    check_point(result[1], 21.0);
}

TEST_CASE("Eigen AffineTVarMainVector outer product")
{
    const AffineTVarVector variables(Vector({2.0, -1.0, 3.0}));
    IntervalRow row(2);
    row[0] = Interval(4.0);
    row[1] = Interval(-2.0);

    const auto result = variables * row;

    static_assert(std::same_as<typename decltype(result)::Scalar, AffineT>);
    REQUIRE(result.rows() == 3);
    REQUIRE(result.cols() == 2);
    check_point(result(0, 0), 8.0);
    check_point(result(0, 1), -4.0);
    check_point(result(2, 0), 12.0);
    check_point(result(2, 1), -6.0);
}

TEST_CASE("Eigen AffineT dot product")
{
    AffineTVector lhs(3);
    AffineTVector rhs(3);
    lhs[0] = Interval(1.0); lhs[1] = Interval(2.0); lhs[2] = Interval(3.0);
    rhs[0] = Interval(4.0); rhs[1] = Interval(-1.0); rhs[2] = Interval(2.0);

    const AffineT result = lhs.dot(rhs);
    check_point(result, 8.0);
}

TEST_CASE("Eigen AffineT sum")
{
    const AffineTVector vector = make_point_vector_5();
    const AffineT result = vector.sum();
    check_point(result, 15.0);
}

TEST_CASE("Eigen AffineT trace")
{
    AffineTMatrix matrix(3, 3);
    matrix(0, 0) = Interval(1.0); matrix(0, 1) = Interval(2.0); matrix(0, 2) = Interval(3.0);
    matrix(1, 0) = Interval(4.0); matrix(1, 1) = Interval(5.0); matrix(1, 2) = Interval(6.0);
    matrix(2, 0) = Interval(7.0); matrix(2, 1) = Interval(8.0); matrix(2, 2) = Interval(9.0);

    const AffineT result = matrix.trace();
    check_point(result, 15.0);
}

TEST_CASE("Eigen AffineT reverse")
{
    const AffineTVector vector = make_point_vector_5();
    const AffineTVector reversed = vector.reverse();

    check_point(reversed[0], 5.0);
    check_point(reversed[1], 4.0);
    check_point(reversed[4], 1.0);
}

TEST_CASE("Eigen AffineT replicate")
{
    const AffineTVarVector row_variables(
        IntervalVector({{1}, {2}}));
    const AffineTRow row = AffineTVector(row_variables).transpose();

    const AffineTMatrix replicated = row.replicate(3, 2);
    REQUIRE(replicated.rows() == 3);
    REQUIRE(replicated.cols() == 4);

    check_point(replicated(0, 0), 1.0);
    check_point(replicated(0, 3), 2.0);
    check_point(replicated(2, 0), 1.0);
    check_point(replicated(2, 3), 2.0);
}

TEST_CASE("Eigen AffineT fixed-size matrices")
{
    FixedAffineTMatrix matrix;
    matrix(0, 0) = Interval(1.0); matrix(0, 1) = Interval(2.0); matrix(0, 2) = Interval(3.0);
    matrix(1, 0) = Interval(4.0); matrix(1, 1) = Interval(5.0); matrix(1, 2) = Interval(6.0);

    static_assert(FixedAffineTMatrix::RowsAtCompileTime == 2);
    static_assert(FixedAffineTMatrix::ColsAtCompileTime == 3);

    check_point(matrix(0, 0), 1.0);
    check_point(matrix(1, 2), 6.0);

    const FixedAffineTRow row = matrix.row(0);
    check_point(row[0], 1.0);
    check_point(row[2], 3.0);

    const FixedAffineTVector column = matrix.transpose().col(1);
    check_point(column[0], 4.0);
    check_point(column[2], 6.0);
}

TEST_CASE("Eigen AffineT interval-valued block operations")
{
    const AffineTVarVector variables(
        IntervalVector({{-2.0, 1.0}, {3.0, 5.0}, {-1.0, 4.0}}));

    AffineTMatrix matrix(3, 2);
    matrix.col(0) = variables;
    matrix.col(1) = -variables;

    const AffineTMatrix block = matrix.block(0, 0, 2, 2);

    CHECK(block(0, 0).itv() == variables[0].itv());
    CHECK(block(1, 0).itv() == variables[1].itv());
    CHECK(block(0, 1).itv() == (-variables[0]).itv());
    CHECK(block(1, 1).itv() == (-variables[1]).itv());
}

TEST_CASE("Eigen AffineT interval enclosure through linear expression")
{
    const AffineTVarVector variables(
        IntervalVector({{-2.0, 1.0}, {3.0, 5.0}, {-1.0, 4.0}}));
    const AffineTVector AffineT = variables;

    const AffineTVector result = 2.0 * AffineT - AffineT;

    for(Index i = 0; i < result.size(); ++i)
        check_interval_enclosure(result[i], variables[i].itv());
}

TEST_CASE("Eigen AffineT eval materializes expressions")
{
    const AffineTVarVector lhs_variables(
        IntervalVector({{1}, {2}, {3}}));
    const AffineTVarVector rhs_variables(
        IntervalVector({{4}, {5}, {6}}));
    AffineTVector lhs(lhs_variables);
    AffineTVector rhs(rhs_variables);

    const auto expression = lhs + rhs;
    const AffineTVector evaluated = expression.eval();

    check_point(evaluated[0], 5.0);
    check_point(evaluated[1], 7.0);
    check_point(evaluated[2], 9.0);
}


// ============================================================================
// Compound assignments between AffineTMain Eigen containers and all
// compatible CODAC scalar/container types.
//
// The left-hand side is always an AffineTMainVector or AffineTMainMatrix.
// For every compatible right-hand side, +=, -=, *= and /= are exercised
// independently.
//
// Point intervals are used deliberately so that the expected result is
// exact while still testing the Interval/AffineT/AffineTVar interoperability.
// ============================================================================

TEST_CASE("Eigen AffineTVector compound assignments with Interval scalar",
          "[Eigen][compound-assignment][AffineTVector][Interval]")
{
    AffineTVector initial(3);
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


TEST_CASE("Eigen AffineTVector compound assignments with AffineT scalar",
          "[Eigen][compound-assignment][AffineTVector][AffineT]")
{
    AffineTVarVector initial(3);
    initial[0] = Interval(1.0);
    initial[1] = Interval(2.0);
    initial[2] = Interval(4.0);

    const AffineT scalar =initial[1];

    SECTION("operator+=")
    {
        AffineTVector value = initial;
        value += scalar;

        check_point(value[0], 3.0);
        check_point(value[1], 4.0);
        check_point(value[2], 6.0);
    }

    SECTION("operator-=")
    {
    	AffineTVector value = initial;
        value -= scalar;

        check_point(value[0], -1.0);
        check_point(value[1], 0.0);
        check_point(value[2], 2.0);
    }

    SECTION("operator*=")
    {
    	AffineTVector value = initial;
        value *= scalar;

        check_point(value[0], 2.0);
        check_point(value[1], 4.0);
        check_point(value[2], 8.0);
    }

    SECTION("operator/=")
    {
    	AffineTVector value = initial;
        value /= scalar;

        check_point(value[0], 0.5);
        check_point(value[1], 1.0);
        check_point(value[2], 2.0);
    }
}


TEST_CASE("Eigen AffineTVector compound assignments with IntervalVector",
          "[Eigen][compound-assignment][AffineTVector][IntervalVector]")
{
    AffineTVarVector initial(3);
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
    	AffineTVector value = initial;
        value += rhs;

        check_point(value[0], 3.0);
        check_point(value[1], 5.0);
        check_point(value[2], 8.0);
    }

    SECTION("operator-=")
    {
    	AffineTVector value = initial;
        value -= rhs;

        check_point(value[0], -1.0);
        check_point(value[1], -1.0);
        check_point(value[2], 0.0);
    }

//    SECTION("operator*=")
//    {
//    	AffineTVector value = initial;
//        value *= rhs;
//
//        check_point(value[0], 2.0);
//        check_point(value[1], 6.0);
//        check_point(value[2], 16.0);
//    }

//    SECTION("operator/=")
//    {
//    	AffineTVector value = initial;
//        value /= rhs;
//
//        check_point(value[0], 0.5);
//        check_point(value[1], 2.0 / 3.0);
//        check_point(value[2], 1.0);
//    }
}


TEST_CASE("Eigen AffineTVector compound assignments with AffineTVector",
          "[Eigen][compound-assignment][AffineTVector][AffineTVector]")
{
    AffineTVarVector initial(3);
    initial[0] = Interval(1.0);
    initial[1] = Interval(2.0);
    initial[2] = Interval(4.0);

    AffineTVector rhs(3);
    rhs[0] = Interval(2.0);
    rhs[1] = Interval(3.0);
    rhs[2] = Interval(4.0);

    SECTION("operator+=")
    {
    	AffineTVector value = initial;
        value += rhs;

        check_point(value[0], 3.0);
        check_point(value[1], 5.0);
        check_point(value[2], 8.0);
    }

    SECTION("operator-=")
    {
    	AffineTVector value = initial;
        value -= rhs;

        check_point(value[0], -1.0);
        check_point(value[1], -1.0);
        check_point(value[2], 0.0);
    }

//    SECTION("operator*=")
//    {
//    	AffineTVector value = initial;
//        value *= rhs;
//
//        check_point(value[0], 2.0);
//        check_point(value[1], 6.0);
//        check_point(value[2], 16.0);
//    }

//    SECTION("operator/=")
//    {
//        auto value = initial;
//        AffineTVector /= rhs;
//
//        check_point(value[0], 0.5);
//        check_point(value[1], 2.0 / 3.0);
//        check_point(value[2], 1.0);
//    }
}


TEST_CASE("Eigen AffineTVector compound assignments with AffineTVarVector",
          "[Eigen][compound-assignment][AffineTVector][AffineTVarVector]")
{
    AffineTVarVector initial(3);
    initial[0] = Interval(1.0);
    initial[1] = Interval(2.0);
    initial[2] = Interval(4.0);

    const AffineTVarVector rhs(
        IntervalVector({
            {2.0},
            {3.0},
            {4.0}
        }));

    SECTION("operator+=")
    {
    	AffineTVector value = initial;
        value += rhs;

        check_point(value[0], 3.0);
        check_point(value[1], 5.0);
        check_point(value[2], 8.0);
    }

    SECTION("operator-=")
    {
    	AffineTVector value = initial;
        value -= rhs;

        check_point(value[0], -1.0);
        check_point(value[1], -1.0);
        check_point(value[2], 0.0);
    }

//    SECTION("operator*=")
//    {
//    	AffineTVector value = initial;
//        value *= rhs;
//
//        check_point(value[0], 2.0);
//        check_point(value[1], 6.0);
//        check_point(value[2], 16.0);
//    }

//    SECTION("operator/=")
//    {
//    	AffineTVector value = initial;
//        value /= rhs;
//
//        check_point(value[0], 0.5);
//        check_point(value[1], 2.0 / 3.0);
//        check_point(value[2], 1.0);
//    }
}


TEST_CASE("Eigen AffineTMatrix compound assignments with Interval scalar",
          "[Eigen][compound-assignment][AffineTMatrix][Interval]")
{
    AffineTMatrix initial(2, 2);
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


TEST_CASE("Eigen AffineTMatrix compound assignments with AffineT scalar",
          "[Eigen][compound-assignment][AffineTMatrix][AffineT]")
{
    AffineTMatrix initial(2, 2);
    initial(0, 0) = Interval(1.0);
    initial(0, 1) = Interval(2.0);
    initial(1, 0) = Interval(3.0);
    initial(1, 1) = Interval(4.0);

    AffineT scalar;
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

TEST_CASE("Eigen AffineTMatrix compound assignments with IntervalMatrix",
          "[Eigen][compound-assignment][AffineTMatrix][IntervalMatrix]")
{
    AffineTMatrix initial(2, 2);
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
          "[Eigen][compound-assignment][AffineTMatrix][IntervalMatrix]")
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


TEST_CASE("Eigen AffineTMatrix compound assignments with AffineTMatrix",
          "[Eigen][compound-assignment][AffineTMatrix][AffineTMatrix]")
{
    AffineTMatrix initial(2, 2);
    initial(0, 0) = Interval(1.0);
    initial(0, 1) = Interval(2.0);
    initial(1, 0) = Interval(3.0);
    initial(1, 1) = Interval(4.0);

    AffineTMatrix rhs(2, 2);
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
