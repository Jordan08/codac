/**
 * \file codac2_tests_AffineMatrix.cpp
 * \brief Unit tests for AffineMainMatrix: construction/access/resize,
 *        block/row/column/diagonal views and assignment, transpose,
 *        factories (Zero/Ones/Constant/Identity), compound assignments,
 *        and every product involving a matrix or a vector-by-row outer
 *        product (real, Interval and AffineT operands, dynamic and
 *        fixed-size).
 */

#include <catch2/catch_test_macros.hpp>

#include <concepts>
#include <sstream>
#include <iostream>
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

using DynamicAffineTMatrix = Eigen::Matrix<AffineT, Eigen::Dynamic, Eigen::Dynamic>;
using FixedAffineTMatrix =   Eigen::Matrix<AffineT, 2, 3>;
using FixedAffineTVector =   Eigen::Matrix<AffineT, 3, 1>;
using FixedAffineTRow =      Eigen::Matrix<AffineT, 1, 3>;

using RealMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using RealRow = Eigen::Matrix<double, 1, Eigen::Dynamic>;
using RealVector = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using IntervalMatrixEigen = Eigen::Matrix<Interval, Eigen::Dynamic, Eigen::Dynamic>;
using IntervalRowEigen = Eigen::Matrix<Interval, 1, Eigen::Dynamic>;

const double ERROR = std::numeric_limits<double>::epsilon()*100;

void check_point(const AffineT& actual, double expected)
{
    CAPTURE(actual.itv(), expected);
    CHECK(actual == Approx<AffineT>(Interval(expected),ERROR));
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

static_assert(std::same_as<typename AffineTMatrix::Scalar, AffineT>);

TEST_CASE("AffineMainMatrix construction, access and assignment")
{
  AffineTMatrix a(2, 3);
  REQUIRE(a.rows() == 2);
  REQUIRE(a.cols() == 3);

  a(0,0) = Interval(1.,2.); a(0,1) = 3.; a(0,2) = Interval(-1.,1.);
  a(1,0) = -2.;            a(1,1) = Interval(4.,5.); a(1,2) = 0.;

  CHECK(a(0,0).itv() == Interval(1.,2.));
  CHECK(a(0,1).itv() == Interval(3.));
  CHECK(a(1,1).itv() == Interval(4.,5.));

  AffineTMatrix b = a;
  CHECK(b.rows() == 2);
  CHECK(b.cols() == 3);
  CHECK(b.itv() == a.itv());

  b(0,0) = Interval(10.,11.);
  CHECK(a(0,0).itv() == Interval(1.,2.));
  CHECK(b(0,0).itv() == Interval(10.,11.));

  AffineTMatrix c(1,1);
  c = a;
  CHECK(c.itv() == a.itv());
}

TEST_CASE("AffineMainMatrix resize preserves Eigen coefficients")
{
  AffineTMatrix a(2,2);
  a(0,0)=1.; a(0,1)=2.; a(1,0)=3.; a(1,1)=4.;

  a.conservativeResize(3,3);
  CHECK(a.rows() == 3);
  CHECK(a.cols() == 3);
  CHECK(a(0,0).itv() == Interval(1.));
  CHECK(a(0,1).itv() == Interval(2.));
  CHECK(a(1,0).itv() == Interval(3.));
  CHECK(a(1,1).itv() == Interval(4.));

  a(2,0)=5.; a(2,1)=6.; a(2,2)=7.; a(0,2)=8.; a(1,2)=9.;
  a.conservativeResize(2,1);
  CHECK(a.rows() == 2);
  CHECK(a.cols() == 1);
  CHECK(a(0,0).itv() == Interval(1.));
  CHECK(a(1,0).itv() == Interval(3.));

  a.resize(4,2);
  CHECK(a.rows() == 4);
  CHECK(a.cols() == 2);
}

TEST_CASE("AffineMainMatrix empty coefficients and stream output")
{
  AffineTMatrix a(2,2);
  a(0,0)=1.; a(0,1)=2.; a(1,0)=3.; a(1,1)=4.;
  a(0,1).set_empty();

  CHECK(a.is_empty());
  CHECK(a.itv()(0,0) == Interval(1.));
  CHECK(a.itv()(0,1).is_empty());

  std::ostringstream stream;
  stream << a;
  CHECK(stream.str() == "[ empty 2x2 mat ]");

  a.set_empty();
  for (Eigen::Index i=0; i<a.rows(); ++i)
    for (Eigen::Index j=0; j<a.cols(); ++j)
      CHECK(a(i,j).is_empty());
}

TEST_CASE("AffineMainMatrix linear addition and subtraction")
{
  AffineTMatrix a(2,2), b(2,2);
  a(0,0)=Interval(1.,2.); a(0,1)=Interval(-1.,1.);
  a(1,0)=3.;              a(1,1)=Interval(4.,5.);
  b(0,0)=Interval(2.,3.); b(0,1)=2.;
  b(1,0)=Interval(-2.,0.); b(1,1)=Interval(1.,2.);

  AffineTMatrix sum = a + b;
  AffineTMatrix difference = a - b;
  for (Eigen::Index i=0; i<2; ++i) {
    for (Eigen::Index j=0; j<2; ++j) {
      CHECK(sum(i,j).itv().is_superset(a(i,j).itv()+b(i,j).itv()));
      CHECK(difference(i,j).itv().is_superset(a(i,j).itv()-b(i,j).itv()));
    }
  }
  AffineTVarVector variables(IntervalVector({{1.,2.},{-1.,1.}}));
  a.col(0) = variables;
  a.col(1) = variables;
  AffineTMatrix zero = a - a;

  for (Eigen::Index i=0; i<2; ++i)
    for (Eigen::Index j=0; j<2; ++j)
      CHECK(zero(i,j).itv() == Interval(0.));


}

TEST_CASE("AffineMainMatrix linear addition with row and segment assignment")
{
  AffineTMatrix a(2,2), b(2,2);
  AffineTVarVector v(IntervalVector({{1.,2.},{-1.,1.},{3.,4.},{5.,6.}}));
  a.row(0) = v.segment(0,2);
  a.row(1) = v.segment(2,2);
  b(0,0) = v[0];
  b(0,1) = v[1];
  b(1,0) = v[2];
  b(1,1) = v[3];

  AffineTMatrix zero = a - b;

  for (Eigen::Index i=0; i<2; ++i)
    for (Eigen::Index j=0; j<2; ++j)
      CHECK(zero(i,j).itv() == Interval(0.));
}

TEST_CASE("AffineMainMatrix scalar operations")
{
  AffineTMatrix a(2,2);
  a(0,0)=Interval(1.,2.); a(0,1)=Interval(-3.,-1.);
  a(1,0)=Interval(0.,4.); a(1,1)=5.;

  AffineTMatrix twice = 2. * a;
  AffineTMatrix negative = a * -3.;
  AffineTMatrix half = a / 2.;
  for (Eigen::Index i=0; i<2; ++i) {
    for (Eigen::Index j=0; j<2; ++j) {
      CHECK(twice(i,j).itv().is_superset(2.*a(i,j).itv()));
      CHECK(negative(i,j).itv().is_superset(-3.*a(i,j).itv()));
      CHECK(half(i,j).itv().is_superset(a(i,j).itv()/2.));
    }
  }
}

TEST_CASE("AffineMainMatrix matrix product")
{
  AffineTMatrix a(2,3), b(3,2);
  const double av[2][3]={{1.,2.,3.},{-1.,4.,2.}};
  const double bv[3][2]={{2.,1.},{0.,-1.},{3.,2.}};
  for (Eigen::Index i=0;i<2;++i) for(Eigen::Index j=0;j<3;++j) a(i,j)=av[i][j];
  for (Eigen::Index i=0;i<3;++i) for(Eigen::Index j=0;j<2;++j) b(i,j)=bv[i][j];

  AffineTMatrix result = a * b;
  REQUIRE(result.rows()==2);
  REQUIRE(result.cols()==2);
  CHECK(result(0,0).itv() == Interval(11.));
  CHECK(result(0,1).itv() == Interval(5.));
  CHECK(result(1,0).itv() == Interval(4.));
  CHECK(result(1,1).itv() == Interval(-1.));
}

TEST_CASE("AffineMainMatrix times AffineMainVector")
{
  AffineTMatrix a(2,3);
  a(0,0)=1.; a(0,1)=2.;  a(0,2)=-1.;
  a(1,0)=3.; a(1,1)=-2.; a(1,2)=4.;

  AffineTVarVector variables(IntervalVector({{1.,2.},{-1.,1.},{3.,4.}}));
  AffineTVector x = variables;
  AffineTVector y = a * x;

  REQUIRE(y.size()==2);
  Interval expected0 = x[0].itv()+2.*x[1].itv()-x[2].itv();
  Interval expected1 = 3.*x[0].itv()-2.*x[1].itv()+4.*x[2].itv();
  CHECK(y[0].itv().is_superset(expected0));
  CHECK(y[1].itv().is_superset(expected1));
}

TEST_CASE("AffineMainMatrix times AffineVarMainVector")
{
  AffineTMatrix a(2,2);
  a(0,0)=1.; a(0,1)=-1.;
  a(1,0)=2.; a(1,1)=3.;

  AffineTVarVector x(IntervalVector({{-2.,5.},{1.,4.}}));
  AffineTVector y = a * x;

  REQUIRE(y.size()==2);
  CHECK(y[0].itv().is_superset(x[0].itv()-x[1].itv()));
  CHECK(y[1].itv().is_superset(2.*x[0].itv()+3.*x[1].itv()));
}

TEST_CASE("AffineMainVector row times AffineMainMatrix")
{
  AffineTMatrix a(2,3);
  a(0,0)=1.; a(0,1)=2.; a(0,2)=3.;
  a(1,0)=4.; a(1,1)=5.; a(1,2)=6.;

  AffineTVarVector variables(Vector({2.,-1.}));
  AffineTVector x = variables;
  Eigen::Matrix<AffineMain<AF_Default>,1,-1> y = x.transpose() * a;

  CHECK(y(0,0).itv() == Interval(-2.));
  CHECK(y(0,1).itv() == Interval(-1.));
  CHECK(y(0,2).itv() == Interval(0.));
}

TEST_CASE("AffineMainMatrix chained linear expression")
{
  AffineTMatrix a(2,2), b(2,2), c(2,2);
  a(0,0)=1.; a(0,1)=2.; a(1,0)=3.; a(1,1)=4.;
  b(0,0)=2.; b(0,1)=-1.; b(1,0)=0.; b(1,1)=5.;
  c(0,0)=-1.; c(0,1)=3.; c(1,0)=2.; c(1,1)=1.;

  AffineTMatrix result = 2.*a - 3.*b + c;
  CHECK(result(0,0).itv() == Interval(-5.));
  CHECK(result(0,1).itv() == Interval(10.));
  CHECK(result(1,0).itv() == Interval(8.));
  CHECK(result(1,1).itv() == Interval(-6.));
}


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



TEST_CASE("Eigen scalar compound assignments - fixed-size AffineTMatrix",
          "[Eigen][compound-assignment][AffineT][fixed-size]")
{
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


TEST_CASE("Eigen scalar compound assignments - fixed-size IntervalMatrix",
          "[Eigen][compound-assignment][Interval][fixed-size]")
{
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



TEST_CASE("Eigen AffineT trace")
{
    AffineTMatrix matrix(3, 3);
    matrix(0, 0) = Interval(1.0); matrix(0, 1) = Interval(2.0); matrix(0, 2) = Interval(3.0);
    matrix(1, 0) = Interval(4.0); matrix(1, 1) = Interval(5.0); matrix(1, 2) = Interval(6.0);
    matrix(2, 0) = Interval(7.0); matrix(2, 1) = Interval(8.0); matrix(2, 2) = Interval(9.0);

    const AffineT result = matrix.trace();
    check_point(result, 15.0);
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


TEST_CASE("AffineMainMatrix stream output for a non-empty matrix")
{
    const AffineTMatrix matrix = make_point_matrix_3x4();
    std::ostringstream stream;
    stream << matrix;

    CAPTURE(stream.str());
    CHECK_FALSE(stream.str().empty());
    CHECK(stream.str().find("empty") == std::string::npos);
}


