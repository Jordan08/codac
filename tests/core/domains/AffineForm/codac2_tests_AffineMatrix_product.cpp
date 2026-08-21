/**
 * \file codac2_tests_AffineMatrix_extended.cpp
 * \brief Extensive unit tests for products involving Affine matrices,
 *        Affine vectors, Affine-variable vectors, and Affine rows.
 */

#include <catch2/catch_test_macros.hpp>

#include <concepts>
#include <sstream>
#include <type_traits>
#include <utility>

#include "codac2_Affine.h"
#include "codac2_AffineMatrix.h"
#include "codac2_Approx.h"

using namespace codac2;

namespace {

using Model = AF_Default;
using AffineT = AffineMain<Model>;
using AffineTMatrix = AffineMainMatrix<Model>;
using AffineTVector = AffineMainVector<Model>;
using AffineTVarVector = AffineVarMainVector<Model>;
using AffineTRow = AffineMainRow<Model>;
using RealMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using RealRow = Eigen::Matrix<double, 1, Eigen::Dynamic>;
using RealVector = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using IntervalMatrixEigen = Eigen::Matrix<Interval, Eigen::Dynamic, Eigen::Dynamic>;
using IntervalRowEigen = Eigen::Matrix<Interval, 1, Eigen::Dynamic>;

IntervalMatrixEigen interval_hull(const AffineTMatrix& value)
{
    IntervalMatrixEigen result(value.rows(), value.cols());
    for(Eigen::Index i = 0; i < value.rows(); ++i)
        for(Eigen::Index j = 0; j < value.cols(); ++j)
            result(i, j) = value(i, j).itv();
    return result;
}

IntervalMatrixEigen interval_hull(const AffineTRow& value)
{
    IntervalMatrixEigen result(value.rows(), value.cols());
    for(Eigen::Index j = 0; j < value.cols(); ++j)
        result(0, j) = value(0, j).itv();
    return result;
}

IntervalMatrixEigen interval_hull(const AffineTVector& value)
{
    IntervalMatrixEigen result(value.rows(), 1);
    for(Eigen::Index i = 0; i < value.rows(); ++i)
        result(i, 0) = value(i).itv();
    return result;
}

IntervalMatrixEigen interval_hull(const AffineTVarVector& value)
{
    IntervalMatrixEigen result(value.rows(), 1);
    for(Eigen::Index i = 0; i < value.rows(); ++i)
        result(i, 0) = value(i).itv();
    return result;
}

void check_encloses(
    const AffineTMatrix& actual,
    const IntervalMatrixEigen& expected)
{
    REQUIRE(actual.rows() == expected.rows());
    REQUIRE(actual.cols() == expected.cols());

    for(Eigen::Index i = 0; i < actual.rows(); ++i)
    {
        for(Eigen::Index j = 0; j < actual.cols(); ++j)
        {
            CAPTURE(i, j, actual(i, j).itv(), expected(i, j));
            CHECK(actual(i, j).itv().is_superset(expected(i, j)));
        }
    }
}

void check_encloses(
    const AffineTVector& actual,
    const IntervalMatrixEigen& expected)
{
    REQUIRE(expected.cols() == 1);
    REQUIRE(actual.rows() == expected.rows());

    for(Eigen::Index i = 0; i < actual.rows(); ++i)
    {
        CAPTURE(i, actual(i).itv(), expected(i, 0));
        CHECK(actual(i).itv().is_superset(expected(i, 0)));
    }
}

void check_encloses(
    const AffineTRow& actual,
    const IntervalMatrixEigen& expected)
{
    REQUIRE(expected.rows() == 1);
    REQUIRE(actual.cols() == expected.cols());

    for(Eigen::Index j = 0; j < actual.cols(); ++j)
    {
        CAPTURE(j, actual(0, j).itv(), expected(0, j));
        CHECK(actual(0, j).itv().is_superset(expected(0, j)));
    }
}

AffineTMatrix make_AffineT_matrix_2x3()
{
    AffineTMatrix value(2, 3);
    value(0, 0) = 1.0;
    value(0, 1) = -2.0;
    value(0, 2) = 3.0;
    value(1, 0) = 4.0;
    value(1, 1) = 0.5;
    value(1, 2) = -1.0;
    return value;
}

RealMatrix make_real_matrix_2x3()
{
    RealMatrix value(2, 3);
    value << 1.0, -2.0, 3.0,
             4.0, 0.5, -1.0;
    return value;
}

IntervalMatrixEigen make_interval_matrix_2x3()
{
    IntervalMatrixEigen value(2, 3);
    value(0, 0) = Interval(0.9, 1.1);
    value(0, 1) = Interval(-2.1, -1.9);
    value(0, 2) = Interval(2.9, 3.1);
    value(1, 0) = Interval(3.9, 4.1);
    value(1, 1) = Interval(0.4, 0.6);
    value(1, 2) = Interval(-1.1, -0.9);
    return value;
}

AffineTVarVector make_variable_vector_3()
{
    return AffineTVarVector(
        IntervalVector({{1.0, 2.0}, {-1.0, 1.0}, {3.0, 4.0}}));
}

AffineTVector make_AffineT_vector_3()
{
    const AffineTVarVector variables = make_variable_vector_3();
    return AffineTVector(variables);
}

AffineTRow make_AffineT_row_3()
{
    AffineTRow value(3);
    value(0) = Interval(1.0, 2.0);
    value(1) = Interval(-1.0, 1.0);
    value(2) = Interval(3.0, 4.0);
    return value;
}

IntervalMatrixEigen multiply_intervals(
    const IntervalMatrixEigen& lhs,
    const IntervalMatrixEigen& rhs)
{
    REQUIRE(lhs.cols() == rhs.rows());

    IntervalMatrixEigen result(lhs.rows(), rhs.cols());
    for(Eigen::Index i = 0; i < lhs.rows(); ++i)
    {
        for(Eigen::Index j = 0; j < rhs.cols(); ++j)
        {
            Interval sum(0.0);
            for(Eigen::Index k = 0; k < lhs.cols(); ++k)
                sum += lhs(i, k) * rhs(k, j);
            result(i, j) = sum;
        }
    }
    return result;
}

} // namespace

TEST_CASE("AffineTVarMainVector converts transparently to AffineTMainVector")
{
    const AffineTVarVector variables = make_variable_vector_3();

    SECTION("construction")
    {
        const AffineTVector values(variables);
        REQUIRE(values.size() == variables.size());
        for(Eigen::Index i = 0; i < values.size(); ++i)
            CHECK(values(i).itv() == variables(i).itv());
    }

    SECTION("assignment")
    {
        AffineTVector values;
        values = variables;
        REQUIRE(values.size() == variables.size());
        for(Eigen::Index i = 0; i < values.size(); ++i)
            CHECK(values(i).itv() == variables(i).itv());
    }
}

TEST_CASE("Real matrix times AffineTMainVector")
{
    const RealMatrix matrix = make_real_matrix_2x3();
    const AffineTVector vector = make_AffineT_vector_3();

    const AffineTVector result = matrix * vector;
    const IntervalMatrixEigen expected =
        multiply_intervals(matrix.cast<Interval>(), interval_hull(vector));

    check_encloses(result, expected);
}

TEST_CASE("Real matrix times AffineTVarMainVector")
{
    const RealMatrix matrix = make_real_matrix_2x3();
    const AffineTVarVector vector = make_variable_vector_3();

    const auto result = matrix * vector;
    static_assert(std::same_as<typename decltype(result)::Scalar, AffineT>);
    REQUIRE(result.rows() == 2);
    REQUIRE(result.cols() == 1);

    const IntervalMatrixEigen expected =
        multiply_intervals(matrix.cast<Interval>(), interval_hull(vector));
    check_encloses(AffineTVector(result), expected);
}

TEST_CASE("Interval matrix times AffineTMainVector")
{
    const IntervalMatrixEigen matrix = make_interval_matrix_2x3();
    const AffineTVector vector = make_AffineT_vector_3();

    const AffineTVector result = matrix * vector;
    const IntervalMatrixEigen expected =
        multiply_intervals(matrix, interval_hull(vector));

    check_encloses(result, expected);
}

TEST_CASE("Interval matrix times AffineTVarMainVector")
{
    const IntervalMatrixEigen matrix = make_interval_matrix_2x3();
    const AffineTVarVector vector = make_variable_vector_3();

    const auto result = matrix * vector;
    static_assert(std::same_as<typename decltype(result)::Scalar, AffineT>);
    REQUIRE(result.rows() == 2);
    REQUIRE(result.cols() == 1);

    const IntervalMatrixEigen expected =
        multiply_intervals(matrix, interval_hull(vector));
    check_encloses(AffineTVector(result), expected);
}

TEST_CASE("AffineT matrix times AffineTMainVector")
{
    const AffineTMatrix matrix = make_AffineT_matrix_2x3();
    const AffineTVector vector = make_AffineT_vector_3();

    const AffineTVector result = matrix * vector;
    const IntervalMatrixEigen expected =
        multiply_intervals(interval_hull(matrix), interval_hull(vector));

    check_encloses(result, expected);
}

TEST_CASE("AffineT matrix times AffineTVarMainVector")
{
    const AffineTMatrix matrix = make_AffineT_matrix_2x3();
    const AffineTVarVector vector = make_variable_vector_3();

    const auto result = matrix * vector;
    static_assert(std::same_as<typename decltype(result)::Scalar, AffineT>);
    REQUIRE(result.rows() == 2);
    REQUIRE(result.cols() == 1);

    const IntervalMatrixEigen expected =
        multiply_intervals(interval_hull(matrix), interval_hull(vector));
    check_encloses(AffineTVector(result), expected);
}

TEST_CASE("Fixed-size real matrix times AffineTVarMainVector")
{
    Eigen::Matrix<double, 2, 3> matrix;
    matrix << 1.0, 2.0, 3.0,
              -1.0, 4.0, 0.5;
    const AffineTVarVector vector = make_variable_vector_3();

    const auto result = matrix * vector;
    static_assert(std::same_as<typename decltype(result)::Scalar, AffineT>);
    static_assert(decltype(result)::RowsAtCompileTime == 2);
    static_assert(decltype(result)::ColsAtCompileTime == 1);

    const IntervalMatrixEigen expected =
        multiply_intervals(matrix.cast<Interval>(), interval_hull(vector));
    check_encloses(AffineTVector(result), expected);
}

TEST_CASE("Fixed-size interval matrix times AffineTVarMainVector")
{
    Eigen::Matrix<Interval, 2, 3> matrix;
    matrix(0, 0) = Interval(1.0);
    matrix(0, 1) = Interval(2.0);
    matrix(0, 2) = Interval(3.0);
    matrix(1, 0) = Interval(-1.0);
    matrix(1, 1) = Interval(4.0);
    matrix(1, 2) = Interval(0.5);
    const AffineTVarVector vector = make_variable_vector_3();

    const auto result = matrix * vector;
    static_assert(std::same_as<typename decltype(result)::Scalar, AffineT>);
    static_assert(decltype(result)::RowsAtCompileTime == 2);
    static_assert(decltype(result)::ColsAtCompileTime == 1);

    const IntervalMatrixEigen expected =
        multiply_intervals(IntervalMatrixEigen(matrix), interval_hull(vector));
    check_encloses(AffineTVector(result), expected);
}

TEST_CASE("Fixed-size AffineT matrix times AffineTVarMainVector")
{
    Eigen::Matrix<AffineT, 2, 3> matrix;
    matrix(0, 0) = 1.0;
    matrix(0, 1) = 2.0;
    matrix(0, 2) = 3.0;
    matrix(1, 0) = -1.0;
    matrix(1, 1) = 4.0;
    matrix(1, 2) = 0.5;
    const AffineTVarVector vector = make_variable_vector_3();

    const auto result = matrix * vector;
    static_assert(std::same_as<typename decltype(result)::Scalar, AffineT>);
    static_assert(decltype(result)::RowsAtCompileTime == 2);
    static_assert(decltype(result)::ColsAtCompileTime == 1);

    check_encloses(
        AffineTVector(result),
        multiply_intervals(interval_hull(AffineTMatrix(matrix)),
                           interval_hull(vector)));
}

TEST_CASE("Real matrix expression times AffineTVarMainVector")
{
    const RealMatrix a = make_real_matrix_2x3();
    const RealMatrix b = 2.0 * make_real_matrix_2x3();
    const AffineTVarVector vector = make_variable_vector_3();

    SECTION("sum expression")
    {
        const auto result = (a + b) * vector;
        const IntervalMatrixEigen expected = multiply_intervals(
            (a + b).eval().cast<Interval>(), interval_hull(vector));
        check_encloses(AffineTVector(result), expected);
    }

    SECTION("scaled expression")
    {
        const auto result = (3.0 * a) * vector;
        const IntervalMatrixEigen expected = multiply_intervals(
            (3.0 * a).eval().cast<Interval>(), interval_hull(vector));
        check_encloses(AffineTVector(result), expected);
    }

    SECTION("block expression")
    {
        RealMatrix large(3, 4);
        large << 1.0, -2.0, 3.0, 8.0,
                 4.0, 0.5, -1.0, 9.0,
                 7.0, 6.0, 5.0, 4.0;
        const auto block = large.block(0, 0, 2, 3);
        const auto result = block * vector;
        const IntervalMatrixEigen expected = multiply_intervals(
            block.eval().cast<Interval>(), interval_hull(vector));
        check_encloses(AffineTVector(result), expected);
    }
}

TEST_CASE("Interval matrix expression times AffineTVarMainVector")
{
    const IntervalMatrixEigen a = make_interval_matrix_2x3();
    const IntervalMatrixEigen b = make_interval_matrix_2x3();
    const AffineTVarVector vector = make_variable_vector_3();

    const auto result = (a + b) * vector;
    const IntervalMatrixEigen expected =
        multiply_intervals((a + b).eval(), interval_hull(vector));
    check_encloses(AffineTVector(result), expected);
}

TEST_CASE("AffineT matrix expression times AffineTVarMainVector")
{
    const AffineTMatrix a = make_AffineT_matrix_2x3();
    const AffineTMatrix b = make_AffineT_matrix_2x3();
    const AffineTVarVector vector = make_variable_vector_3();

    const auto result = (a + b) * vector;
    const AffineTMatrix evaluated = (a + b).eval();
    const IntervalMatrixEigen expected =
        multiply_intervals(interval_hull(evaluated), interval_hull(vector));
    check_encloses(AffineTVector(result), expected);
}

TEST_CASE("AffineTVarMainVector outer product with real row")
{
    const AffineTVarVector vector = make_variable_vector_3();
    RealRow row(4);
    row << 2.0, -1.0, 0.5, 3.0;

    const auto result = vector * row;
    static_assert(std::same_as<typename decltype(result)::Scalar, AffineT>);
    REQUIRE(result.rows() == 3);
    REQUIRE(result.cols() == 4);

    const IntervalMatrixEigen expected = multiply_intervals(
        interval_hull(vector), row.cast<Interval>());
    check_encloses(AffineTMatrix(result), expected);
}

TEST_CASE("AffineTVarMainVector outer product with interval row")
{
    const AffineTVarVector vector = make_variable_vector_3();
    IntervalRowEigen row(4);
    row(0) = Interval(1.0, 2.0);
    row(1) = Interval(-1.0, 1.0);
    row(2) = Interval(0.5);
    row(3) = Interval(2.0, 3.0);

    const auto result = vector * row;
    static_assert(std::same_as<typename decltype(result)::Scalar, AffineT>);
    REQUIRE(result.rows() == 3);
    REQUIRE(result.cols() == 4);

    const IntervalMatrixEigen expected =
        multiply_intervals(interval_hull(vector), IntervalMatrixEigen(row));
    check_encloses(AffineTMatrix(result), expected);
}

TEST_CASE("AffineTVarMainVector outer product with AffineTRow")
{
    const AffineTVarVector vector = make_variable_vector_3();
    AffineTRow row(4);
    row(0) = Interval(1.0, 2.0);
    row(1) = Interval(-1.0, 1.0);
    row(2) = 0.5;
    row(3) = Interval(2.0, 3.0);

    const auto result = vector * row;
    static_assert(std::same_as<typename decltype(result)::Scalar, AffineT>);
    REQUIRE(result.rows() == 3);
    REQUIRE(result.cols() == 4);

    const IntervalMatrixEigen expected =
        multiply_intervals(interval_hull(vector), interval_hull(row));
    check_encloses(AffineTMatrix(result), expected);
}

TEST_CASE("AffineTMainVector outer product with real row")
{
    const AffineTVector vector = make_AffineT_vector_3();
    RealRow row(2);
    row << -2.0, 4.0;

    const AffineTMatrix result = vector * row;
    const IntervalMatrixEigen expected = multiply_intervals(
        interval_hull(vector), row.cast<Interval>());
    check_encloses(result, expected);
}

TEST_CASE("AffineTMainVector outer product with interval row")
{
    const AffineTVector vector = make_AffineT_vector_3();
    IntervalRowEigen row(2);
    row(0) = Interval(-2.0, -1.0);
    row(1) = Interval(3.0, 4.0);

    const AffineTMatrix result = vector * row;
    const IntervalMatrixEigen expected = multiply_intervals(
        interval_hull(vector), IntervalMatrixEigen(row));
    check_encloses(result, expected);
}

TEST_CASE("AffineTMainVector outer product with AffineTRow")
{
    const AffineTVector vector = make_AffineT_vector_3();
    AffineTRow row(2);
    row(0) = Interval(-2.0, -1.0);
    row(1) = Interval(3.0, 4.0);

    const AffineTMatrix result = vector * row;
    const IntervalMatrixEigen expected = multiply_intervals(
        interval_hull(vector), interval_hull(row));
    check_encloses(result, expected);
}

TEST_CASE("AffineTRow times real matrix")
{
    const AffineTRow row = make_AffineT_row_3();
    RealMatrix matrix(3, 2);
    matrix << 1.0, -1.0,
              2.0, 3.0,
              -2.0, 4.0;

    const AffineTRow result = row * matrix;
    const IntervalMatrixEigen expected = multiply_intervals(
        interval_hull(row), matrix.cast<Interval>());
    check_encloses(result, expected);
}

TEST_CASE("AffineTRow times interval matrix")
{
    const AffineTRow row = make_AffineT_row_3();
    IntervalMatrixEigen matrix(3, 2);
    matrix(0, 0) = Interval(1.0, 1.1);
    matrix(0, 1) = Interval(-1.0, -0.9);
    matrix(1, 0) = Interval(2.0, 2.1);
    matrix(1, 1) = Interval(3.0, 3.1);
    matrix(2, 0) = Interval(-2.0, -1.9);
    matrix(2, 1) = Interval(4.0, 4.1);

    const AffineTRow result = row * matrix;
    const IntervalMatrixEigen expected =
        multiply_intervals(interval_hull(row), matrix);
    check_encloses(result, expected);
}

TEST_CASE("AffineTRow times AffineT matrix")
{
    const AffineTRow row = make_AffineT_row_3();
    AffineTMatrix matrix(3, 2);
    matrix(0, 0) = 1.0;
    matrix(0, 1) = -1.0;
    matrix(1, 0) = 2.0;
    matrix(1, 1) = 3.0;
    matrix(2, 0) = -2.0;
    matrix(2, 1) = 4.0;

    const AffineTRow result = row * matrix;
    const IntervalMatrixEigen expected = multiply_intervals(
        interval_hull(row), interval_hull(matrix));
    check_encloses(result, expected);
}

TEST_CASE("Real column vector times AffineTRow")
{
    RealVector vector(3);
    vector << 1.0, -2.0, 4.0;
    const AffineTRow row = make_AffineT_row_3();

    const AffineTMatrix result = vector * row;
    const IntervalMatrixEigen expected = multiply_intervals(
        vector.cast<Interval>(), interval_hull(row));
    check_encloses(result, expected);
}

TEST_CASE("AffineT matrix products preserve dynamic dimensions")
{
    AffineTMatrix lhs(4, 3);
    AffineTMatrix rhs(3, 5);

    for(Eigen::Index i = 0; i < lhs.rows(); ++i)
        for(Eigen::Index j = 0; j < lhs.cols(); ++j)
            lhs(i, j) = static_cast<double>(1 + i + 2 * j);

    for(Eigen::Index i = 0; i < rhs.rows(); ++i)
        for(Eigen::Index j = 0; j < rhs.cols(); ++j)
            rhs(i, j) = static_cast<double>(2 - i + j);

    const AffineTMatrix result = lhs * rhs;
    REQUIRE(result.rows() == 4);
    REQUIRE(result.cols() == 5);

    const IntervalMatrixEigen expected = multiply_intervals(
        interval_hull(lhs), interval_hull(rhs));
    check_encloses(result, expected);
}

TEST_CASE("Degenerate products produce degenerate AffineT results")
{
    RealMatrix matrix(2, 3);
    matrix << 1.0, 2.0, 3.0,
              -1.0, 4.0, 0.5;

    AffineTVarVector variables(Vector({2.0, -1.0, 3.0}));
    const auto result = matrix * variables;

    REQUIRE(result.rows() == 2);
    CHECK(result(0).itv() == Interval(9.0));
    CHECK(result(1).itv() == Interval(-4.5));
}

TEST_CASE("Identity matrix preserves AffineTVarMainVector hulls")
{
    const AffineTVarVector variables = make_variable_vector_3();
    const RealMatrix identity = RealMatrix::Identity(3, 3);

    const auto result = identity * variables;
    REQUIRE(result.size() == variables.size());

    for(Eigen::Index i = 0; i < variables.size(); ++i)
        CHECK(result(i).itv().is_superset(variables(i).itv()));
}

TEST_CASE("Zero real matrix annihilates AffineTVarMainVector")
{
    const AffineTVarVector variables = make_variable_vector_3();
    const RealMatrix zero = RealMatrix::Zero(2, 3);

    const auto result = zero * variables;
    REQUIRE(result.rows() == 2);
    REQUIRE(result.cols() == 1);
    CHECK(result(0).itv() == Interval(0.0));
    CHECK(result(1).itv() == Interval(0.0));
}

TEST_CASE("One-by-one products retain AffineTMain scalar")
{
    const AffineTVarVector variables(IntervalVector({{-2.0, 3.0}}));

    SECTION("real matrix")
    {
        RealMatrix matrix(1, 1);
        matrix(0, 0) = -2.0;
        const auto result = matrix * variables;
        static_assert(std::same_as<typename decltype(result)::Scalar, AffineT>);
        CHECK(result(0).itv().is_superset(Interval(-6.0, 4.0)));
    }

    SECTION("interval matrix")
    {
        IntervalMatrixEigen matrix(1, 1);
        matrix(0, 0) = Interval(-2.0, -1.0);
        const auto result = matrix * variables;
        static_assert(std::same_as<typename decltype(result)::Scalar, AffineT>);
        CHECK(result(0).itv().is_superset(
            matrix(0, 0) * variables(0).itv()));
    }

    SECTION("AffineT matrix")
    {
        AffineTMatrix matrix(1, 1);
        matrix(0, 0) = Interval(-2.0, -1.0);
        const auto result = matrix * variables;
        static_assert(std::same_as<typename decltype(result)::Scalar, AffineT>);
        CHECK(result(0).itv().is_superset(
            matrix(0, 0).itv() * variables(0).itv()));
    }
}

TEST_CASE("Empty coefficient propagation through products")
{
    AffineTMatrix matrix(1, 2);
    matrix(0, 0) = 1.0;
    matrix(0, 1) = 2.0;
    matrix(0, 1).set_empty();

    const AffineTVarVector variables(
        IntervalVector({{1.0, 2.0}, {3.0, 4.0}}));
    const auto result = matrix * variables;

    REQUIRE(result.size() == 1);
    CHECK(result(0).is_empty());
}


TEST_CASE("Matrix product preserves dependency: A*x - A*x collapses to zero")
{
    const RealMatrix matrix = make_real_matrix_2x3();
    const AffineTVarVector x = make_variable_vector_3();

    const AffineTVector result = (matrix * x) - (matrix * x);

    for(Eigen::Index i = 0; i < result.size(); ++i)
        CHECK(result(i).itv() == Interval(0.0));
}


