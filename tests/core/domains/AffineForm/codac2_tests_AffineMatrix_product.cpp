/**
 * \file codac2_tests_AffineMatrix_extended.cpp
 * \brief Extensive unit tests for products involving affine matrices,
 *        affine vectors, affine-variable vectors, and affine rows.
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
using Affine = AffineMain<Model>;
using AffineMatrix = AffineMainMatrix<Model>;
using AffineVector = AffineMainVector<Model>;
using AffineVarVector = AffineVarMainVector<Model>;
using AffineRow = AffineMainRow<Model>;
using RealMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using RealRow = Eigen::Matrix<double, 1, Eigen::Dynamic>;
using RealVector = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using IntervalMatrixEigen = Eigen::Matrix<Interval, Eigen::Dynamic, Eigen::Dynamic>;
using IntervalRowEigen = Eigen::Matrix<Interval, 1, Eigen::Dynamic>;

IntervalMatrixEigen interval_hull(const AffineMatrix& value)
{
    IntervalMatrixEigen result(value.rows(), value.cols());
    for(Eigen::Index i = 0; i < value.rows(); ++i)
        for(Eigen::Index j = 0; j < value.cols(); ++j)
            result(i, j) = value(i, j).itv();
    return result;
}

IntervalMatrixEigen interval_hull(const AffineRow& value)
{
    IntervalMatrixEigen result(value.rows(), value.cols());
    for(Eigen::Index j = 0; j < value.cols(); ++j)
        result(0, j) = value(0, j).itv();
    return result;
}

IntervalMatrixEigen interval_hull(const AffineVector& value)
{
    IntervalMatrixEigen result(value.rows(), 1);
    for(Eigen::Index i = 0; i < value.rows(); ++i)
        result(i, 0) = value(i).itv();
    return result;
}

IntervalMatrixEigen interval_hull(const AffineVarVector& value)
{
    IntervalMatrixEigen result(value.rows(), 1);
    for(Eigen::Index i = 0; i < value.rows(); ++i)
        result(i, 0) = value(i).itv();
    return result;
}

void check_encloses(
    const AffineMatrix& actual,
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
    const AffineVector& actual,
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
    const AffineRow& actual,
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

AffineMatrix make_affine_matrix_2x3()
{
    AffineMatrix value(2, 3);
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

AffineVarVector make_variable_vector_3()
{
    return AffineVarVector(
        IntervalVector({{1.0, 2.0}, {-1.0, 1.0}, {3.0, 4.0}}));
}

AffineVector make_affine_vector_3()
{
    const AffineVarVector variables = make_variable_vector_3();
    return AffineVector(variables);
}

AffineRow make_affine_row_3()
{
    AffineRow value(3);
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

TEST_CASE("AffineVarMainVector converts transparently to AffineMainVector")
{
    const AffineVarVector variables = make_variable_vector_3();

    SECTION("construction")
    {
        const AffineVector values(variables);
        REQUIRE(values.size() == variables.size());
        for(Eigen::Index i = 0; i < values.size(); ++i)
            CHECK(values(i).itv() == variables(i).itv());
    }

    SECTION("assignment")
    {
        AffineVector values;
        values = variables;
        REQUIRE(values.size() == variables.size());
        for(Eigen::Index i = 0; i < values.size(); ++i)
            CHECK(values(i).itv() == variables(i).itv());
    }
}

TEST_CASE("Real matrix times AffineMainVector")
{
    const RealMatrix matrix = make_real_matrix_2x3();
    const AffineVector vector = make_affine_vector_3();

    const AffineVector result = matrix * vector;
    const IntervalMatrixEigen expected =
        multiply_intervals(matrix.cast<Interval>(), interval_hull(vector));

    check_encloses(result, expected);
}

TEST_CASE("Real matrix times AffineVarMainVector")
{
    const RealMatrix matrix = make_real_matrix_2x3();
    const AffineVarVector vector = make_variable_vector_3();

    const auto result = matrix * vector;
    static_assert(std::same_as<typename decltype(result)::Scalar, Affine>);
    REQUIRE(result.rows() == 2);
    REQUIRE(result.cols() == 1);

    const IntervalMatrixEigen expected =
        multiply_intervals(matrix.cast<Interval>(), interval_hull(vector));
    check_encloses(AffineVector(result), expected);
}

TEST_CASE("Interval matrix times AffineMainVector")
{
    const IntervalMatrixEigen matrix = make_interval_matrix_2x3();
    const AffineVector vector = make_affine_vector_3();

    const AffineVector result = matrix * vector;
    const IntervalMatrixEigen expected =
        multiply_intervals(matrix, interval_hull(vector));

    check_encloses(result, expected);
}

TEST_CASE("Interval matrix times AffineVarMainVector")
{
    const IntervalMatrixEigen matrix = make_interval_matrix_2x3();
    const AffineVarVector vector = make_variable_vector_3();

    const auto result = matrix * vector;
    static_assert(std::same_as<typename decltype(result)::Scalar, Affine>);
    REQUIRE(result.rows() == 2);
    REQUIRE(result.cols() == 1);

    const IntervalMatrixEigen expected =
        multiply_intervals(matrix, interval_hull(vector));
    check_encloses(AffineVector(result), expected);
}

TEST_CASE("Affine matrix times AffineMainVector")
{
    const AffineMatrix matrix = make_affine_matrix_2x3();
    const AffineVector vector = make_affine_vector_3();

    const AffineVector result = matrix * vector;
    const IntervalMatrixEigen expected =
        multiply_intervals(interval_hull(matrix), interval_hull(vector));

    check_encloses(result, expected);
}

TEST_CASE("Affine matrix times AffineVarMainVector")
{
    const AffineMatrix matrix = make_affine_matrix_2x3();
    const AffineVarVector vector = make_variable_vector_3();

    const auto result = matrix * vector;
    static_assert(std::same_as<typename decltype(result)::Scalar, Affine>);
    REQUIRE(result.rows() == 2);
    REQUIRE(result.cols() == 1);

    const IntervalMatrixEigen expected =
        multiply_intervals(interval_hull(matrix), interval_hull(vector));
    check_encloses(AffineVector(result), expected);
}

TEST_CASE("Fixed-size real matrix times AffineVarMainVector")
{
    Eigen::Matrix<double, 2, 3> matrix;
    matrix << 1.0, 2.0, 3.0,
              -1.0, 4.0, 0.5;
    const AffineVarVector vector = make_variable_vector_3();

    const auto result = matrix * vector;
    static_assert(std::same_as<typename decltype(result)::Scalar, Affine>);
    static_assert(decltype(result)::RowsAtCompileTime == 2);
    static_assert(decltype(result)::ColsAtCompileTime == 1);

    const IntervalMatrixEigen expected =
        multiply_intervals(matrix.cast<Interval>(), interval_hull(vector));
    check_encloses(AffineVector(result), expected);
}

TEST_CASE("Fixed-size interval matrix times AffineVarMainVector")
{
    Eigen::Matrix<Interval, 2, 3> matrix;
    matrix(0, 0) = Interval(1.0);
    matrix(0, 1) = Interval(2.0);
    matrix(0, 2) = Interval(3.0);
    matrix(1, 0) = Interval(-1.0);
    matrix(1, 1) = Interval(4.0);
    matrix(1, 2) = Interval(0.5);
    const AffineVarVector vector = make_variable_vector_3();

    const auto result = matrix * vector;
    static_assert(std::same_as<typename decltype(result)::Scalar, Affine>);
    static_assert(decltype(result)::RowsAtCompileTime == 2);
    static_assert(decltype(result)::ColsAtCompileTime == 1);

    const IntervalMatrixEigen expected =
        multiply_intervals(IntervalMatrixEigen(matrix), interval_hull(vector));
    check_encloses(AffineVector(result), expected);
}

TEST_CASE("Fixed-size affine matrix times AffineVarMainVector")
{
    Eigen::Matrix<Affine, 2, 3> matrix;
    matrix(0, 0) = 1.0;
    matrix(0, 1) = 2.0;
    matrix(0, 2) = 3.0;
    matrix(1, 0) = -1.0;
    matrix(1, 1) = 4.0;
    matrix(1, 2) = 0.5;
    const AffineVarVector vector = make_variable_vector_3();

    const auto result = matrix * vector;
    static_assert(std::same_as<typename decltype(result)::Scalar, Affine>);
    static_assert(decltype(result)::RowsAtCompileTime == 2);
    static_assert(decltype(result)::ColsAtCompileTime == 1);

    check_encloses(
        AffineVector(result),
        multiply_intervals(interval_hull(AffineMatrix(matrix)),
                           interval_hull(vector)));
}

TEST_CASE("Real matrix expression times AffineVarMainVector")
{
    const RealMatrix a = make_real_matrix_2x3();
    const RealMatrix b = 2.0 * make_real_matrix_2x3();
    const AffineVarVector vector = make_variable_vector_3();

    SECTION("sum expression")
    {
        const auto result = (a + b) * vector;
        const IntervalMatrixEigen expected = multiply_intervals(
            (a + b).eval().cast<Interval>(), interval_hull(vector));
        check_encloses(AffineVector(result), expected);
    }

    SECTION("scaled expression")
    {
        const auto result = (3.0 * a) * vector;
        const IntervalMatrixEigen expected = multiply_intervals(
            (3.0 * a).eval().cast<Interval>(), interval_hull(vector));
        check_encloses(AffineVector(result), expected);
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
        check_encloses(AffineVector(result), expected);
    }
}

TEST_CASE("Interval matrix expression times AffineVarMainVector")
{
    const IntervalMatrixEigen a = make_interval_matrix_2x3();
    const IntervalMatrixEigen b = make_interval_matrix_2x3();
    const AffineVarVector vector = make_variable_vector_3();

    const auto result = (a + b) * vector;
    const IntervalMatrixEigen expected =
        multiply_intervals((a + b).eval(), interval_hull(vector));
    check_encloses(AffineVector(result), expected);
}

TEST_CASE("Affine matrix expression times AffineVarMainVector")
{
    const AffineMatrix a = make_affine_matrix_2x3();
    const AffineMatrix b = make_affine_matrix_2x3();
    const AffineVarVector vector = make_variable_vector_3();

    const auto result = (a + b) * vector;
    const AffineMatrix evaluated = (a + b).eval();
    const IntervalMatrixEigen expected =
        multiply_intervals(interval_hull(evaluated), interval_hull(vector));
    check_encloses(AffineVector(result), expected);
}

TEST_CASE("AffineVarMainVector outer product with real row")
{
    const AffineVarVector vector = make_variable_vector_3();
    RealRow row(4);
    row << 2.0, -1.0, 0.5, 3.0;

    const auto result = vector * row;
    static_assert(std::same_as<typename decltype(result)::Scalar, Affine>);
    REQUIRE(result.rows() == 3);
    REQUIRE(result.cols() == 4);

    const IntervalMatrixEigen expected = multiply_intervals(
        interval_hull(vector), row.cast<Interval>());
    check_encloses(AffineMatrix(result), expected);
}

TEST_CASE("AffineVarMainVector outer product with interval row")
{
    const AffineVarVector vector = make_variable_vector_3();
    IntervalRowEigen row(4);
    row(0) = Interval(1.0, 2.0);
    row(1) = Interval(-1.0, 1.0);
    row(2) = Interval(0.5);
    row(3) = Interval(2.0, 3.0);

    const auto result = vector * row;
    static_assert(std::same_as<typename decltype(result)::Scalar, Affine>);
    REQUIRE(result.rows() == 3);
    REQUIRE(result.cols() == 4);

    const IntervalMatrixEigen expected =
        multiply_intervals(interval_hull(vector), IntervalMatrixEigen(row));
    check_encloses(AffineMatrix(result), expected);
}

TEST_CASE("AffineVarMainVector outer product with AffineRow")
{
    const AffineVarVector vector = make_variable_vector_3();
    AffineRow row(4);
    row(0) = Interval(1.0, 2.0);
    row(1) = Interval(-1.0, 1.0);
    row(2) = 0.5;
    row(3) = Interval(2.0, 3.0);

    const auto result = vector * row;
    static_assert(std::same_as<typename decltype(result)::Scalar, Affine>);
    REQUIRE(result.rows() == 3);
    REQUIRE(result.cols() == 4);

    const IntervalMatrixEigen expected =
        multiply_intervals(interval_hull(vector), interval_hull(row));
    check_encloses(AffineMatrix(result), expected);
}

TEST_CASE("AffineMainVector outer product with real row")
{
    const AffineVector vector = make_affine_vector_3();
    RealRow row(2);
    row << -2.0, 4.0;

    const AffineMatrix result = vector * row;
    const IntervalMatrixEigen expected = multiply_intervals(
        interval_hull(vector), row.cast<Interval>());
    check_encloses(result, expected);
}

TEST_CASE("AffineMainVector outer product with interval row")
{
    const AffineVector vector = make_affine_vector_3();
    IntervalRowEigen row(2);
    row(0) = Interval(-2.0, -1.0);
    row(1) = Interval(3.0, 4.0);

    const AffineMatrix result = vector * row;
    const IntervalMatrixEigen expected = multiply_intervals(
        interval_hull(vector), IntervalMatrixEigen(row));
    check_encloses(result, expected);
}

TEST_CASE("AffineMainVector outer product with AffineRow")
{
    const AffineVector vector = make_affine_vector_3();
    AffineRow row(2);
    row(0) = Interval(-2.0, -1.0);
    row(1) = Interval(3.0, 4.0);

    const AffineMatrix result = vector * row;
    const IntervalMatrixEigen expected = multiply_intervals(
        interval_hull(vector), interval_hull(row));
    check_encloses(result, expected);
}

TEST_CASE("AffineRow times real matrix")
{
    const AffineRow row = make_affine_row_3();
    RealMatrix matrix(3, 2);
    matrix << 1.0, -1.0,
              2.0, 3.0,
              -2.0, 4.0;

    const AffineRow result = row * matrix;
    const IntervalMatrixEigen expected = multiply_intervals(
        interval_hull(row), matrix.cast<Interval>());
    check_encloses(result, expected);
}

TEST_CASE("AffineRow times interval matrix")
{
    const AffineRow row = make_affine_row_3();
    IntervalMatrixEigen matrix(3, 2);
    matrix(0, 0) = Interval(1.0, 1.1);
    matrix(0, 1) = Interval(-1.0, -0.9);
    matrix(1, 0) = Interval(2.0, 2.1);
    matrix(1, 1) = Interval(3.0, 3.1);
    matrix(2, 0) = Interval(-2.0, -1.9);
    matrix(2, 1) = Interval(4.0, 4.1);

    const AffineRow result = row * matrix;
    const IntervalMatrixEigen expected =
        multiply_intervals(interval_hull(row), matrix);
    check_encloses(result, expected);
}

TEST_CASE("AffineRow times affine matrix")
{
    const AffineRow row = make_affine_row_3();
    AffineMatrix matrix(3, 2);
    matrix(0, 0) = 1.0;
    matrix(0, 1) = -1.0;
    matrix(1, 0) = 2.0;
    matrix(1, 1) = 3.0;
    matrix(2, 0) = -2.0;
    matrix(2, 1) = 4.0;

    const AffineRow result = row * matrix;
    const IntervalMatrixEigen expected = multiply_intervals(
        interval_hull(row), interval_hull(matrix));
    check_encloses(result, expected);
}

TEST_CASE("Real column vector times AffineRow")
{
    RealVector vector(3);
    vector << 1.0, -2.0, 4.0;
    const AffineRow row = make_affine_row_3();

    const AffineMatrix result = vector * row;
    const IntervalMatrixEigen expected = multiply_intervals(
        vector.cast<Interval>(), interval_hull(row));
    check_encloses(result, expected);
}

TEST_CASE("Affine matrix products preserve dynamic dimensions")
{
    AffineMatrix lhs(4, 3);
    AffineMatrix rhs(3, 5);

    for(Eigen::Index i = 0; i < lhs.rows(); ++i)
        for(Eigen::Index j = 0; j < lhs.cols(); ++j)
            lhs(i, j) = static_cast<double>(1 + i + 2 * j);

    for(Eigen::Index i = 0; i < rhs.rows(); ++i)
        for(Eigen::Index j = 0; j < rhs.cols(); ++j)
            rhs(i, j) = static_cast<double>(2 - i + j);

    const AffineMatrix result = lhs * rhs;
    REQUIRE(result.rows() == 4);
    REQUIRE(result.cols() == 5);

    const IntervalMatrixEigen expected = multiply_intervals(
        interval_hull(lhs), interval_hull(rhs));
    check_encloses(result, expected);
}

TEST_CASE("Degenerate products produce degenerate affine results")
{
    RealMatrix matrix(2, 3);
    matrix << 1.0, 2.0, 3.0,
              -1.0, 4.0, 0.5;

    AffineVarVector variables(Vector({2.0, -1.0, 3.0}));
    const auto result = matrix * variables;

    REQUIRE(result.rows() == 2);
    CHECK(result(0).itv() == Interval(9.0));
    CHECK(result(1).itv() == Interval(-4.5));
}

TEST_CASE("Identity matrix preserves AffineVarMainVector hulls")
{
    const AffineVarVector variables = make_variable_vector_3();
    const RealMatrix identity = RealMatrix::Identity(3, 3);

    const auto result = identity * variables;
    REQUIRE(result.size() == variables.size());

    for(Eigen::Index i = 0; i < variables.size(); ++i)
        CHECK(result(i).itv().is_superset(variables(i).itv()));
}

TEST_CASE("Zero real matrix annihilates AffineVarMainVector")
{
    const AffineVarVector variables = make_variable_vector_3();
    const RealMatrix zero = RealMatrix::Zero(2, 3);

    const auto result = zero * variables;
    REQUIRE(result.rows() == 2);
    REQUIRE(result.cols() == 1);
    CHECK(result(0).itv() == Interval(0.0));
    CHECK(result(1).itv() == Interval(0.0));
}

TEST_CASE("One-by-one products retain AffineMain scalar")
{
    const AffineVarVector variables(IntervalVector({{-2.0, 3.0}}));

    SECTION("real matrix")
    {
        RealMatrix matrix(1, 1);
        matrix(0, 0) = -2.0;
        const auto result = matrix * variables;
        static_assert(std::same_as<typename decltype(result)::Scalar, Affine>);
        CHECK(result(0).itv().is_superset(Interval(-6.0, 4.0)));
    }

    SECTION("interval matrix")
    {
        IntervalMatrixEigen matrix(1, 1);
        matrix(0, 0) = Interval(-2.0, -1.0);
        const auto result = matrix * variables;
        static_assert(std::same_as<typename decltype(result)::Scalar, Affine>);
        CHECK(result(0).itv().is_superset(
            matrix(0, 0) * variables(0).itv()));
    }

    SECTION("affine matrix")
    {
        AffineMatrix matrix(1, 1);
        matrix(0, 0) = Interval(-2.0, -1.0);
        const auto result = matrix * variables;
        static_assert(std::same_as<typename decltype(result)::Scalar, Affine>);
        CHECK(result(0).itv().is_superset(
            matrix(0, 0).itv() * variables(0).itv()));
    }
}

TEST_CASE("Empty coefficient propagation through products")
{
    AffineMatrix matrix(1, 2);
    matrix(0, 0) = 1.0;
    matrix(0, 1) = 2.0;
    matrix(0, 1).set_empty();

    const AffineVarVector variables(
        IntervalVector({{1.0, 2.0}, {3.0, 4.0}}));
    const auto result = matrix * variables;

    REQUIRE(result.size() == 1);
    CHECK(result(0).is_empty());
}
