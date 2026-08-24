/**
 * \file codac2_tests_AffineMatrix.h
 * Unit tests for AffineMainMatrix.
 */

#include <catch2/catch_test_macros.hpp>
#include <sstream>
#include <iostream>
#include "codac2_Affine.h"
#include "codac2_AffineMatrix.h"
#include "codac2_Approx.h"

using namespace codac2;

namespace {
using MatrixAA = AffineMainMatrix<AF_Default>;
using VectorAA = AffineMainVector<AF_Default>;
using VarVectorAA = AffineVarMainVector<AF_Default>;


} // namespace

TEST_CASE("AffineMainMatrix construction, access and assignment")
{
  MatrixAA a(2, 3);
  REQUIRE(a.rows() == 2);
  REQUIRE(a.cols() == 3);

  a(0,0) = Interval(1.,2.); a(0,1) = 3.; a(0,2) = Interval(-1.,1.);
  a(1,0) = -2.;            a(1,1) = Interval(4.,5.); a(1,2) = 0.;

  CHECK(a(0,0).itv() == Interval(1.,2.));
  CHECK(a(0,1).itv() == Interval(3.));
  CHECK(a(1,1).itv() == Interval(4.,5.));

  MatrixAA b = a;
  CHECK(b.rows() == 2);
  CHECK(b.cols() == 3);
  CHECK(b.itv() == a.itv());

  b(0,0) = Interval(10.,11.);
  CHECK(a(0,0).itv() == Interval(1.,2.));
  CHECK(b(0,0).itv() == Interval(10.,11.));

  MatrixAA c(1,1);
  c = a;
  CHECK(c.itv() == a.itv());
}

TEST_CASE("AffineMainMatrix resize preserves Eigen coefficients")
{
  MatrixAA a(2,2);
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
  MatrixAA a(2,2);
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
  MatrixAA a(2,2), b(2,2);
  a(0,0)=Interval(1.,2.); a(0,1)=Interval(-1.,1.);
  a(1,0)=3.;              a(1,1)=Interval(4.,5.);
  b(0,0)=Interval(2.,3.); b(0,1)=2.;
  b(1,0)=Interval(-2.,0.); b(1,1)=Interval(1.,2.);

  MatrixAA sum = a + b;
  MatrixAA difference = a - b;
  for (Eigen::Index i=0; i<2; ++i) {
    for (Eigen::Index j=0; j<2; ++j) {
      CHECK(sum(i,j).itv().is_superset(a(i,j).itv()+b(i,j).itv()));
      CHECK(difference(i,j).itv().is_superset(a(i,j).itv()-b(i,j).itv()));
    }
  }
  VarVectorAA variables(IntervalVector({{1.,2.},{-1.,1.}}));
  a.col(0) = variables;
  a.col(1) = variables;
  MatrixAA zero = a - a;

  for (Eigen::Index i=0; i<2; ++i)
    for (Eigen::Index j=0; j<2; ++j)
      CHECK(zero(i,j).itv() == Interval(0.));


}

TEST_CASE("AffineMainMatrix linear addition with row and segment assignment")
{
  MatrixAA a(2,2), b(2,2);
  VarVectorAA v(IntervalVector({{1.,2.},{-1.,1.},{3.,4.},{5.,6.}}));
  a.row(0) = v.segment(0,2);
  a.row(1) = v.segment(2,2);
  b(0,0) = v[0];
  b(0,1) = v[1];
  b(1,0) = v[2];
  b(1,1) = v[3];

  MatrixAA zero = a - b;

  for (Eigen::Index i=0; i<2; ++i)
    for (Eigen::Index j=0; j<2; ++j)
      CHECK(zero(i,j).itv() == Interval(0.));
}

TEST_CASE("AffineMainMatrix scalar operations")
{
  MatrixAA a(2,2);
  a(0,0)=Interval(1.,2.); a(0,1)=Interval(-3.,-1.);
  a(1,0)=Interval(0.,4.); a(1,1)=5.;

  MatrixAA twice = 2. * a;
  MatrixAA negative = a * -3.;
  MatrixAA half = a / 2.;
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
  MatrixAA a(2,3), b(3,2);
  const double av[2][3]={{1.,2.,3.},{-1.,4.,2.}};
  const double bv[3][2]={{2.,1.},{0.,-1.},{3.,2.}};
  for (Eigen::Index i=0;i<2;++i) for(Eigen::Index j=0;j<3;++j) a(i,j)=av[i][j];
  for (Eigen::Index i=0;i<3;++i) for(Eigen::Index j=0;j<2;++j) b(i,j)=bv[i][j];

  MatrixAA result = a * b;
  REQUIRE(result.rows()==2);
  REQUIRE(result.cols()==2);
  CHECK(result(0,0).itv() == Interval(11.));
  CHECK(result(0,1).itv() == Interval(5.));
  CHECK(result(1,0).itv() == Interval(4.));
  CHECK(result(1,1).itv() == Interval(-1.));
}

TEST_CASE("AffineMainMatrix times AffineMainVector")
{
  MatrixAA a(2,3);
  a(0,0)=1.; a(0,1)=2.;  a(0,2)=-1.;
  a(1,0)=3.; a(1,1)=-2.; a(1,2)=4.;

  VarVectorAA variables(IntervalVector({{1.,2.},{-1.,1.},{3.,4.}}));
  VectorAA x = variables;
  VectorAA y = a * x;

  REQUIRE(y.size()==2);
  Interval expected0 = x[0].itv()+2.*x[1].itv()-x[2].itv();
  Interval expected1 = 3.*x[0].itv()-2.*x[1].itv()+4.*x[2].itv();
  CHECK(y[0].itv().is_superset(expected0));
  CHECK(y[1].itv().is_superset(expected1));
}

TEST_CASE("AffineMainMatrix times AffineVarMainVector")
{
  MatrixAA a(2,2);
  a(0,0)=1.; a(0,1)=-1.;
  a(1,0)=2.; a(1,1)=3.;

  VarVectorAA x(IntervalVector({{-2.,5.},{1.,4.}}));
  VectorAA y = a * x;

  REQUIRE(y.size()==2);
  CHECK(y[0].itv().is_superset(x[0].itv()-x[1].itv()));
  CHECK(y[1].itv().is_superset(2.*x[0].itv()+3.*x[1].itv()));
}

TEST_CASE("AffineMainVector row times AffineMainMatrix")
{
  MatrixAA a(2,3);
  a(0,0)=1.; a(0,1)=2.; a(0,2)=3.;
  a(1,0)=4.; a(1,1)=5.; a(1,2)=6.;

  VarVectorAA variables(Vector({2.,-1.}));
  VectorAA x = variables;
  Eigen::Matrix<AffineMain<AF_Default>,1,-1> y = x.transpose() * a;

  CHECK(y(0,0).itv() == Interval(-2.));
  CHECK(y(0,1).itv() == Interval(-1.));
  CHECK(y(0,2).itv() == Interval(0.));
}

TEST_CASE("AffineMainMatrix chained linear expression")
{
  MatrixAA a(2,2), b(2,2), c(2,2);
  a(0,0)=1.; a(0,1)=2.; a(1,0)=3.; a(1,1)=4.;
  b(0,0)=2.; b(0,1)=-1.; b(1,0)=0.; b(1,1)=5.;
  c(0,0)=-1.; c(0,1)=3.; c(1,0)=2.; c(1,1)=1.;

  MatrixAA result = 2.*a - 3.*b + c;
  CHECK(result(0,0).itv() == Interval(-5.));
  CHECK(result(0,1).itv() == Interval(10.));
  CHECK(result(1,0).itv() == Interval(8.));
  CHECK(result(1,1).itv() == Interval(-6.));
}
