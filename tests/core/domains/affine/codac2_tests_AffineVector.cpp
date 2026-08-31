/**
 * \file codac2_tests_AffineVector.cpp
 * \brief Unit tests for the "value" containers AffineMainVector and
 *        AffineMainRow: construction, indexing, segment/head/tail/block
 *        views, transpose, dot product, sum/reverse, cwiseProduct/
 *        cwiseQuotient, compound assignments (scalar, IntervalVector,
 *        AffineTVector, AffineTVarVector), and the arithmetic of
 *        AffineMainVector (linear combinations, empty propagation,
 *        component-wise nonlinear expressions).
 * \date       2026
 * \author     Jordan Ninin
 * \copyright  Copyright 2026 Codac Team
 * \license    GNU Lesser General Public License (LGPL)
 */

#include <catch2/catch_test_macros.hpp>

#include <concepts>
#include <sstream>
#include <type_traits>
#include <utility>

#include "codac2_Affine.h"
#include "codac2_AffineMatrix.h"
#include "codac2_AffineRow.h"
#include "codac2_AffineVector.h"
#include "codac2_IntervalMatrix.h"
#include "codac2_IntervalRow.h"
#include "codac2_IntervalVector.h"
#include "codac2_Approx.h"

using namespace std;
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

AffineTVector make_point_vector_5()
{
    const AffineTVarVector variables(
        IntervalVector({{1}, {2}, {3}, {4}, {5}}));
    return AffineTVector(variables);
}

} // namespace

typedef AF_Default AA;

static_assert(std::same_as<typename AffineTVector::Scalar, AffineT>);
static_assert(std::same_as<typename AffineTRow::Scalar, AffineT>);
static_assert(AffineTVector::ColsAtCompileTime == 1);
static_assert(AffineTRow::RowsAtCompileTime == 1);

TEST_CASE("AffineMainVector construction, copy and access")
{
SECTION("Construction from an AffineVarMainVector and copy")
  {
    AffineVarMainVector<AA> x(3);
    x = IntervalVector({{1.,2.},{3.,4.},{5.,6.}});

    AffineMainVector<AA> y = x;

    CHECK(y.size() == 3);
    CHECK(y.itv() == IntervalVector({{1.,2.},{3.,4.},{5.,6.}}));
    CHECK(y[0].itv() == Interval(1.,2.));
    CHECK(y[1].itv() == Interval(3.,4.));
    CHECK(y[2].itv() == Interval(5.,6.));
  }

SECTION("Copy of a copy")
  {
    AffineVarMainVector<AA> x(2);
    x = IntervalVector({{-2.,1.},{4.,7.}});

    AffineMainVector<AA> y = x;
    AffineMainVector<AA> z = y;

    CHECK(z.size() == 2);
    CHECK(z.itv() == y.itv());
    CHECK(z.itv() == IntervalVector({{-2.,1.},{4.,7.}}));
  }

SECTION("Component assignment on the destination leaves the source vector unchanged")
  {
    AffineVarMainVector<AA> x(2);
    x = IntervalVector({{1.,2.},{3.,4.}});

    AffineMainVector<AA> y = x;
    y[0] = Interval(-1.,1.);

    CHECK(y.itv() == IntervalVector({{-1.,1.},{3.,4.}}));
    CHECK(x.itv() == IntervalVector({{1.,2.},{3.,4.}}));
  }

SECTION("Direct component assignment via operator[]")
  {
	AffineVarMainVector<AA> x(3);
	x[2] = Interval(-1., 1.);

	CHECK(x[2].itv() == Interval(-1., 1.));
  }

SECTION("conservativeResize grows and then shrinks the vector")
  {
	  AffineVarMainVector<AA> x( IntervalVector({ {-1.,2.}, {3.,4.} }));

	  x.conservativeResize(3);
	  CHECK(x.size() == 3);
	  CHECK(x[0].itv() == Interval(-1.,2.));
	  CHECK(x[1].itv() == Interval(3.,4.));
	  CHECK(x[2].itv() == Interval());
	  x.conservativeResize(1);

	  CHECK(x.size() == 1);
	  CHECK(x[0].itv() == Interval(-1.,2.));
  }

SECTION("conservativeResize to the same size is a no-op")
  {
    AffineVarMainVector<AA> x(IntervalVector({{1.,2.},{3.,4.}}));
    x.conservativeResize(2); // même taille
    CHECK(x[0].itv() == Interval(1.,2.));
    CHECK(x[1].itv() == Interval(3.,4.));
  }

SECTION("Construction through AffineVarMainVector from an IntervalVector.")
  {
	IntervalVector box({{-5.,-2.},{0.},{7.,9.},{-oo,oo}});
    AffineVarMainVector<AA> x(box);
    AffineMainVector<AA> y = x;

    CHECK(y.size() == 4);
    CHECK(y.rows() == 4);
    CHECK(y.cols() == 1);
    CHECK(y.itv() == box);
  }

SECTION("Construction through AffineVarMainVector from a real Vector.")
  {
    Vector values({-2.,0.,3.5,8.});
    AffineVarMainVector<AA> x(values);
    AffineMainVector<AA> y = x;

    CHECK(y.size() == 4);
    CHECK(y.itv() == IntervalVector({{-2.},{0.},{3.5},{8.}}));
  }

SECTION("The default affine variables are unbounded.")
  {
    AffineVarMainVector<AA> x(3);
    AffineMainVector<AA> y = x;

    CHECK(y.size() == 3);
    CHECK(y.itv() == IntervalVector(3));
    CHECK(!y.is_empty());
  }

SECTION("Parentheses, brackets and Eigen's coeff() address the same elements.")
  {
    AffineVarMainVector<AA> x(3);
    x = IntervalVector({{-1.,1.},{2.,4.},{6.}});
    AffineMainVector<AA> y = x;

    CHECK(y[0].itv() == y(0).itv());
    CHECK(y[1].itv() == y.coeff(1).itv());
    CHECK(y(2).itv() == Interval(6.));
  }

SECTION("Component assignment does not change the vector dimension.")
  {
    AffineVarMainVector<AA> x(4);
    x = IntervalVector({{1.,2.},{3.,4.},{5.,6.},{7.,8.}});
    AffineMainVector<AA> y = x;

    y[1] = Interval(-10.,-8.);
    y[3] = Interval(0.);

    CHECK(y.size() == 4);
    CHECK(y.itv() == IntervalVector({{1.,2.},{-10.,-8.},{5.,6.},{0.}}));
  }

SECTION("Unbounded components are preserved by itv().")
  {
    AffineVarMainVector<AA> x(3);
    x = IntervalVector({{-oo,0.},{1.,oo},{-oo,oo}});
    AffineMainVector<AA> y = x;

    CHECK(y.itv() == IntervalVector({{-oo,0.},{1.,oo},{-oo,oo}}));
    CHECK(y.itv().is_unbounded());
  }

SECTION("AffineVarMainVector must be built directly from the IntervalVector.; Assigning an IntervalVector afterwards instantiates resize(), which is; currently incompatible with the protected copy operations of AffineVarMain.")
  {
    AffineVarMainVector<AA> variables(
      IntervalVector({{1.,2.},{3.,4.},{5.,6.}})
    );

    CHECK(variables.size() == 3);
    CHECK(variables == Approx<AffineMainVector<AA>>(IntervalVector({{1.,2.},{3.,4.},{5.,6.}}),ERROR));
  }

SECTION("Explicit conversion from affine variables to ordinary affine forms.; Eigen cannot assign AffineVarMain directly to AffineMain without a cast.")
  {
    AffineVarMainVector<AA> variables(
      IntervalVector({{-2.,1.},{4.,7.}})
    );

    AffineMainVector<AA> x =
      variables.template cast<AffineMain<AA>>();

    CHECK(x.size() == 2);
    CHECK(x.itv() == IntervalVector({{-2.,1.},{4.,7.}}));
  }

SECTION("Copy construction of an AffineMainVector.")
  {
    AffineVarMainVector<AA> variables(
      IntervalVector({{1.,2.},{3.,4.},{5.,6.}})
    );
    AffineMainVector<AA> x = variables;
    AffineMainVector<AA> y = x;

    CHECK(y.size() == 3);
    CHECK(y.itv() == x.itv());
  }

SECTION("Access through [], () and coeff().")
  {
    AffineVarMainVector<AA> variables(
      IntervalVector({{-1.,1.},{2.,4.},{6.}})
    );
    AffineMainVector<AA> x =
      variables.template cast<AffineMain<AA>>();

    CHECK(x[0].itv() == Interval(-1.,1.));
    CHECK(x(1).itv() == Interval(2.,4.));
    CHECK(x.coeff(2).itv() == Interval(6.));
  }

SECTION("Modification of components and independence after copy.")
  {
    AffineVarMainVector<AA> variables(
      IntervalVector({{1.,2.},{3.,4.},{5.,6.}})
    );
    AffineMainVector<AA> x =
      variables.template cast<AffineMain<AA>>();
    AffineMainVector<AA> y = x;

    y[0] = Interval(-1.,1.);
    y[2] = Interval(0.);

    CHECK(x.itv() == IntervalVector({{1.,2.},{3.,4.},{5.,6.}}));
    CHECK(y.itv() == IntervalVector({{-1.,1.},{3.,4.},{0.}}));
  }

}


TEST_CASE("AffineMainVector empty propagation")
{
SECTION("set_empty() on a copy leaves the source unaffected")
  {
    AffineVarMainVector<AA> x(3);
    x = IntervalVector({{1.,2.},{3.,4.},{5.,6.}});

    AffineMainVector<AA> y = x;

    CHECK(!y.is_empty());
    y.set_empty();
    CHECK(y.is_empty());
    CHECK(y.size() == 3);
  }

SECTION("A single empty component makes the whole vector empty")
  {
    AffineVarMainVector<AA> x(3);
    x = IntervalVector({{1.,2.},Interval::empty(),{5.,6.}});

    AffineMainVector<AA> y = x;

    CHECK(y.is_empty());
    CHECK(y.size() == 3);
    CHECK(y.itv().is_empty());
  }

SECTION("An empty component survives being copied twice")
  {
	  AffineVarMainVector<AA> x(
		IntervalVector({
		  {1.,2.},
		  Interval::empty(),
		  {3.,4.}
		})
	  );

	  AffineMainVector<AA> y = x;
	  AffineMainVector<AA> z = y;

	  CHECK(y.size() == 3);
	  CHECK(y.is_empty());
	  CHECK(z.size() == 3);
	  CHECK(z.is_empty());
  }

SECTION("Setting one component to empty makes the whole vector empty.")
  {
    AffineVarMainVector<AA> x(3);
    x = IntervalVector({{1.,2.},{3.,4.},{5.,6.}});
    AffineMainVector<AA> y = x;

    y[1].set_empty();

    CHECK(y.is_empty());
    CHECK(y.itv().is_empty());
    CHECK(y.size() == 3);
  }

SECTION("Empty vectors keep their dimension through a copy.")
  {
    AffineVarMainVector<AA> x(2);
    x = IntervalVector({{1.,2.},{3.,4.}});
    AffineMainVector<AA> y = x;
    y.set_empty();
    AffineMainVector<AA> z = y;

    CHECK(z.size() == 2);
    CHECK(z.is_empty());
    CHECK(z.itv().is_empty());
  }

SECTION("Empty values propagate through vector addition and subtraction.")
  {
    AffineVarMainVector<AA> x(2);
    AffineVarMainVector<AA> e(2);
    x = IntervalVector({{1.,2.},{3.,4.}});
    e = IntervalVector({Interval::empty(),Interval::empty()});

    AffineMainVector<AA> sum = x + e;
    AffineMainVector<AA> difference = e - x;

    CHECK(sum.is_empty());
    CHECK(difference.is_empty());
    CHECK(sum.size() == 2);
    CHECK(difference.size() == 2);
  }

SECTION("Empty vector propagation through unary minus and scalar operations.")
  {
    AffineVarMainVector<AA> e(3);
    e = IntervalVector({Interval::empty(),Interval::empty(),Interval::empty()});

    AffineMainVector<AA> opposite = -e;
    AffineMainVector<AA> product_left = 2. * e;
    AffineMainVector<AA> product_zero = 0. * e;
    AffineMainVector<AA> quotient = e / 2.;

    CHECK(opposite.is_empty());
    CHECK(product_left.is_empty());
    CHECK(product_zero.is_empty());
    CHECK(quotient.is_empty());
    CHECK(opposite.size() == 3);
    CHECK(product_left.size() == 3);
    CHECK(product_zero.size() == 3);
    CHECK(quotient.size() == 3);
  }

SECTION("Empty vector propagation for both operand orders.")
  {
    AffineVarMainVector<AA> x(3);
    AffineVarMainVector<AA> e(3);
    x = IntervalVector({{1.,2.},{-3.,4.},{5.,6.}});
    e = IntervalVector({Interval::empty(),Interval::empty(),Interval::empty()});

    AffineMainVector<AA> sum_left = e + x;
    AffineMainVector<AA> sum_right = x + e;
    AffineMainVector<AA> difference_left = e - x;
    AffineMainVector<AA> difference_right = x - e;

    CHECK(sum_left.is_empty());
    CHECK(sum_right.is_empty());
    CHECK(difference_left.is_empty());
    CHECK(difference_right.is_empty());
    CHECK(sum_left.size() == 3);
    CHECK(sum_right.size() == 3);
    CHECK(difference_left.size() == 3);
    CHECK(difference_right.size() == 3);
  }

SECTION("Empty propagation through a nested expression.")
  {
    AffineVarMainVector<AA> x(2);
    AffineVarMainVector<AA> y(2);
    AffineVarMainVector<AA> e(2);
    x = IntervalVector({{1.,2.},{3.,4.}});
    y = IntervalVector({{-2.,1.},{5.,7.}});
    e = IntervalVector({Interval::empty(),Interval::empty()});

    AffineMainVector<AA> result = 2. * (x + y) - (3. * e - x) / 2.;

    CHECK(result.is_empty());
    CHECK(result.size() == 2);
    CHECK(result.itv().is_empty());
  }

SECTION("Empty propagation through intermediate AffineMainVector results.")
  {
    AffineVarMainVector<AA> x(2);
    AffineVarMainVector<AA> e(2);
    x = IntervalVector({{1.,3.},{-2.,4.}});
    e = IntervalVector({Interval::empty(),Interval::empty()});

    AffineMainVector<AA> first = x + e;
    AffineMainVector<AA> second = 3. * first - x;
    AffineMainVector<AA> third = (second + x) / 4.;

    CHECK(first.is_empty());
    CHECK(second.is_empty());
    CHECK(third.is_empty());
    CHECK(first.size() == 2);
    CHECK(second.size() == 2);
    CHECK(third.size() == 2);
  }

SECTION("A single empty component propagates through every resulting vector.")
  {
    AffineVarMainVector<AA> x(3);
    AffineVarMainVector<AA> y(3);
    x = IntervalVector({{1.,2.},Interval::empty(),{5.,6.}});
    y = IntervalVector({{-1.,1.},{3.,4.},{2.,3.}});

    AffineMainVector<AA> sum = x + y;
    AffineMainVector<AA> difference = y - x;
    AffineMainVector<AA> expression = 2. * x - 3. * y;

    CHECK(sum.is_empty());
    CHECK(difference.is_empty());
    CHECK(expression.is_empty());
    CHECK(sum.size() == 3);
    CHECK(difference.size() == 3);
    CHECK(expression.size() == 3);
    CHECK(sum.itv().is_empty());
    CHECK(difference.itv().is_empty());
    CHECK(expression.itv().is_empty());
  }

SECTION("Assigning an empty expression to a non-empty vector propagates emptiness.")
  {
    AffineVarMainVector<AA> x(2);
    AffineVarMainVector<AA> e(2);
    x = IntervalVector({{1.,2.},{3.,4.}});
    e = IntervalVector({Interval::empty(),Interval::empty()});

    AffineMainVector<AA> result = x;
    CHECK(!result.is_empty());

    result = x + e;

    CHECK(result.is_empty());
    CHECK(result.size() == 2);
    CHECK(result.itv().is_empty());
  }

SECTION("Empty vector and dimension preservation.")
  {
    AffineVarMainVector<AA> variables(
      IntervalVector({{1.,2.},{3.,4.},{5.,6.}})
    );
    AffineMainVector<AA> x =
      variables.template cast<AffineMain<AA>>();

    CHECK(!x.is_empty());
    x.set_empty();
    CHECK(x.is_empty());
    CHECK(x.size() == 3);
    CHECK(x.itv().is_empty());
  }

SECTION("Propagation of one empty component during conversion.")
  {
    AffineVarMainVector<AA> variables(
      IntervalVector({{1.,2.},Interval::empty(),{5.,6.}})
    );
    AffineMainVector<AA> x =
      variables.template cast<AffineMain<AA>>();

    CHECK(x.is_empty());
    CHECK(x.size() == 3);
    CHECK(x.itv().is_empty());
  }

SECTION("Propagation of emptiness through addition, subtraction and unary minus.")
  {
    AffineVarMainVector<AA> vx(
      IntervalVector({{1.,2.},{3.,4.}})
    );
    AffineVarMainVector<AA> ve(
      IntervalVector({Interval::empty(),Interval::empty()})
    );
    AffineMainVector<AA> x = vx.template cast<AffineMain<AA>>();
    AffineMainVector<AA> empty = ve.template cast<AffineMain<AA>>();

    AffineMainVector<AA> sum_left = empty + x;
    AffineMainVector<AA> sum_right = x + empty;
    AffineMainVector<AA> difference_left = empty - x;
    AffineMainVector<AA> difference_right = x - empty;
    AffineMainVector<AA> opposite = -empty;

    CHECK(sum_left.is_empty());
    CHECK(sum_right.is_empty());
    CHECK(difference_left.is_empty());
    CHECK(difference_right.is_empty());
    CHECK(opposite.is_empty());
    CHECK(sum_left.size() == 2);
    CHECK(difference_right.size() == 2);
  }

SECTION("Propagation of emptiness through scalar and nested operations.")
  {
    AffineVarMainVector<AA> vx(
      IntervalVector({{1.,2.},{3.,4.}})
    );
    AffineVarMainVector<AA> ve(
      IntervalVector({Interval::empty(),Interval::empty()})
    );
    AffineMainVector<AA> x = vx.template cast<AffineMain<AA>>();
    AffineMainVector<AA> empty = ve.template cast<AffineMain<AA>>();

    AffineMainVector<AA> product = 2. * empty;
    AffineMainVector<AA> zero_product = 0. * empty;
    AffineMainVector<AA> nested = 2. * (x + empty) - 0.5 * x;

    CHECK(product.is_empty());
    CHECK(zero_product.is_empty());
    CHECK(nested.is_empty());
    CHECK(product.size() == 2);
    CHECK(nested.itv().is_empty());
  }

}


TEST_CASE("AffineMainVector linear arithmetic")
{
SECTION("Addition of two affine-variable vectors")
  {
    AffineVarMainVector<AA> x(2);
    AffineVarMainVector<AA> y(2);
    x = IntervalVector({{1.,2.},{-3.,-1.}});
    y = IntervalVector({{4.,5.},{2.,6.}});

    AffineMainVector<AA> z = x + y;

    CHECK(z.size() == 2);
    CHECK(z.itv() == IntervalVector({{5.,7.},{-1.,5.}}));
  }

SECTION("Subtraction of an affine-variable vector from an equal copy is exactly zero")
  {
    AffineVarMainVector<AA> x(2);
    AffineVarMainVector<AA> y(2);
    x = IntervalVector({{1.,3.},{-3.,1.}});
    y = x;

    AffineMainVector<AA> z = x - y;

    CHECK(z.size() == 2);
    CHECK(z.itv() == IntervalVector({{0.,0.},{0.,0.}}));
  }

SECTION("Subtraction of two independent affine-variable vectors")
  {
    AffineVarMainVector<AA> x(2);
    AffineVarMainVector<AA> y(2);
    x = IntervalVector({{1.,3.},{-3.,-1.}});
    y = IntervalVector({{4.,6.},{2.,6.}});

    AffineMainVector<AA> z = x - y;

    CHECK(z.size() == 2);
    CHECK(z.itv() == IntervalVector({{-3.,-3.},{-7.,-5.}}));
  }

SECTION("Unary minus")
  {
    AffineVarMainVector<AA> x(3);
    x = IntervalVector({{1.,2.},{-3.,-1.},{0.,4.}});

    AffineMainVector<AA> y = -x;

    CHECK(y.size() == 3);
    CHECK(y.itv() == IntervalVector({{-2.,-1.},{1.,3.},{-4.,0.}}));
  }

SECTION("Scalar multiplication by a positive constant")
  {
    AffineVarMainVector<AA> x(3);
    x = IntervalVector({{1.,2.},{-3.,-1.},{0.,4.}});

    AffineMainVector<AA> y = 2. * x;

    CHECK(y.size() == 3);
    CHECK(y.itv() == IntervalVector({{2.,4.},{-6.,-2.},{0.,8.}}));
  }

SECTION("Scalar division by a positive constant")
  {
    AffineVarMainVector<AA> x(3);
    x = IntervalVector({{1.,2.},{-3.,-1.},{0.,4.}});

    AffineMainVector<AA> y = x / 2.;

    CHECK(y.size() == 3);
    CHECK(y.itv() == IntervalVector({{0.5,1.},{-1.5,-0.5},{0.,2.}}));
  }

SECTION("Addition with a zero vector.")
  {
    AffineVarMainVector<AA> x(3);
    AffineVarMainVector<AA> zero(3);
    x = IntervalVector({{-2.,3.},{4.,5.},{-7.,-1.}});
    zero = IntervalVector({{0.},{0.},{0.}});

    AffineMainVector<AA> y = x + zero;
    AffineMainVector<AA> z = zero + x;

    CHECK(y.itv() == x.itv());
    CHECK(z.itv() == x.itv());
  }

SECTION("Subtraction from a zero vector.")
  {
    AffineVarMainVector<AA> x(3);
    AffineVarMainVector<AA> zero(3);
    x = IntervalVector({{-2.,3.},{4.,5.},{-7.,-1.}});
    zero = IntervalVector({{0.},{0.},{0.}});

    AffineMainVector<AA> y = zero - x;

    CHECK(y.itv() == IntervalVector({{-3.,2.},{-5.,-4.},{1.,7.}}));
  }

SECTION("Dependency is preserved: an affine variable minus itself is exactly zero.")
  {
    AffineVarMainVector<AA> x(3);
    x = IntervalVector({{-2.,3.},{4.,5.},{-7.,-1.}});

    AffineMainVector<AA> y = x - x;

    CHECK(y.itv() == IntervalVector({{0.},{0.},{0.}}));
  }

SECTION("Dependency is preserved in a linear expression.")
  {
    AffineVarMainVector<AA> x(2);
    x = IntervalVector({{1.,2.},{-3.,4.}});

    AffineMainVector<AA> y = x + 2. * x;

    CHECK(y.itv() == IntervalVector({{3.,6.},{-9.,12.}}));
  }

SECTION("Difference of two independent affine-variable vectors.")
  {
    AffineVarMainVector<AA> x(2);
    AffineVarMainVector<AA> y(2);
    x = IntervalVector({{0.,2.},{4.,6.}});
    y = IntervalVector({{-1.,1.},{2.,3.}});

    AffineMainVector<AA> z1 = x - y;
    AffineMainVector<AA> z2 = y - x;

    CHECK(z1.itv() == IntervalVector({{1.,1.},{2.,3.}}));
    CHECK(z2.itv() == -z1.itv());
  }

SECTION("Multiplication by zero.")
  {
    AffineVarMainVector<AA> x(3);
    x = IntervalVector({{-8.,10.},{-3.,4.},{5.}});

    AffineMainVector<AA> y = 0. * x;

    CHECK(y.itv() == IntervalVector({{0.},{0.},{0.}}));
  }

SECTION("Multiplication by a negative scalar reverses interval bounds.")
  {
    AffineVarMainVector<AA> x(3);
    x = IntervalVector({{1.,2.},{-3.,-1.},{0.,4.}});

    AffineMainVector<AA> y = -2. * x;

    CHECK(y.itv() == IntervalVector({{-4.,-2.},{2.,6.},{-8.,0.}}));
  }

SECTION("Division by a negative scalar.")
  {
    AffineVarMainVector<AA> x(3);
    x = IntervalVector({{2.,4.},{-6.,-2.},{0.,8.}});

    AffineMainVector<AA> y = x / -2.;

    CHECK(y.itv() == IntervalVector({{-2.,-1.},{1.,3.},{-4.,0.}}));
  }

SECTION("Chained linear expression inspired by ex_affineform.cpp.")
  {
    AffineVarMainVector<AA> x(3);
    x = IntervalVector({{1.,2.},{1.,3.},{1.,2.}});

    AffineMainVector<AA> y = x + 2. * x;
    AffineMainVector<AA> z = y - x;

    CHECK(y.itv() == IntervalVector({{3.,6.},{3.,9.},{3.,6.}}));
    CHECK(z.itv() == IntervalVector({{2.,4.},{2.,6.},{2.,4.}}));
  }

SECTION("Assignment of an Eigen expression to an existing AffineMainVector.")
  {
    AffineVarMainVector<AA> x(2);
    AffineVarMainVector<AA> y(2);
    x = IntervalVector({{1.,2.},{3.,4.}});
    y = IntervalVector({{-1.,1.},{5.,6.}});

    AffineMainVector<AA> z = x;
    z = x + y;

    CHECK(z.size() == 2);
    CHECK(z.itv() == IntervalVector({{0.,3.},{8.,10.}}));
  }

SECTION("Aliasing during assignment must not corrupt the expression.")
  {
    AffineVarMainVector<AA> x(2);
    x = IntervalVector({{1.,2.},{-3.,4.}});
    AffineMainVector<AA> y = x;

    y = y + x;

    CHECK(y.itv() == IntervalVector({{2.,4.},{-6.,8.}}));
  }

SECTION("A vector of degenerate intervals behaves as a real vector.")
  {
    AffineVarMainVector<AA> x(Vector({1.,-2.,3.}));
    AffineVarMainVector<AA> y(Vector({4.,5.,-6.}));

    AffineMainVector<AA> sum = x + y;
    AffineMainVector<AA> difference = x - y;

    CHECK(sum.itv() == IntervalVector({{5.},{3.},{-3.}}));
    CHECK(difference.itv() == IntervalVector({{-3.},{-7.},{9.}}));
  }

SECTION("Complex linear expression with three independent affine vectors.")
  {
    AffineVarMainVector<AA> x(3);
    AffineVarMainVector<AA> y(3);
    AffineVarMainVector<AA> z(3);
    x = IntervalVector({{1.,2.},{-2.,1.},{4.,5.}});
    y = IntervalVector({{-1.,1.},{3.,4.},{-2.,0.}});
    z = IntervalVector({{2.,6.},{-4.,2.},{1.,3.}});

    AffineMainVector<AA> result = 2. * x - 3. * y + z / 2.;

    CHECK(result.size() == 3);
    CHECK(result.itv() == IntervalVector({{4.,6.},{-15.,-9.},{11.5,14.5}}));
  }

SECTION("Nested expression combining sums, differences and scalar operations.")
  {
    AffineVarMainVector<AA> x(2);
    AffineVarMainVector<AA> y(2);
    x = IntervalVector({{-2.,4.},{1.,3.}});
    y = IntervalVector({{5.,7.},{-4.,2.}});

    AffineMainVector<AA> result = (2. * (x + y) - (x - y)) / 3.;
    //std::cout.precision(10);

    // Algebraically, the expression is (x + 3*y) / 3.
    CHECK(result == Approx<AffineMainVector<AA>>(IntervalVector({{13./3.,25./3.},{-11./3.,3.}}),ERROR));
  }

SECTION("A complex expression must preserve affine dependencies.")
  {
    AffineVarMainVector<AA> x(3);
    x = IntervalVector({{-5.,2.},{1.,7.},{-3.,4.}});

    AffineMainVector<AA> result = 4. * x - 2. * (x + x) + (x - x);

    CHECK(result.itv() == IntervalVector({{0.},{0.},{0.}}));
  }

SECTION("Several intermediate AffineMainVector expressions can be chained.")
  {
    AffineVarMainVector<AA> x(3);
    AffineVarMainVector<AA> y(3);
    x = IntervalVector({{1.,2.},{-1.,3.},{4.,6.}});
    y = IntervalVector({{-2.,0.},{2.,5.},{-3.,1.}});

    AffineMainVector<AA> u = 2. * x - y;
    AffineMainVector<AA> v = x + 3. * y;
    AffineMainVector<AA> result = (u + v) / 3.;

    // (2*x-y+x+3*y)/3 = x + 2*y/3.
    CHECK(result == Approx<AffineMainVector<AA>>(IntervalVector({{-1./3.,2.},{1./3.,19./3.},{2.,20./3.}}),ERROR));
  }

SECTION("Addition between vectors having the same Eigen scalar type.")
  {
    AffineVarMainVector<AA> vx(
      IntervalVector({{1.,2.},{-3.,-1.}})
    );
    AffineVarMainVector<AA> vy(
      IntervalVector({{4.,5.},{2.,6.}})
    );
    AffineMainVector<AA> x = vx;
    AffineMainVector<AA> y = vy;

    AffineMainVector<AA> result = x + y;

    CHECK(result.itv() == IntervalVector({{5.,7.},{-1.,5.}}));
  }

SECTION("Subtraction and unary minus.")
  {
    AffineVarMainVector<AA> vx(
      IntervalVector({{1.,2.},{-3.,-1.}})
    );
    AffineVarMainVector<AA> vy(
      IntervalVector({{4.,5.},{2.,6.}})
    );
    AffineMainVector<AA> x = vx;
    AffineMainVector<AA> y = vy;

    AffineMainVector<AA> difference = x - y;
    AffineMainVector<AA> opposite = -x;

    CHECK(difference.itv() == IntervalVector({{-3.,-3.},{-7.,-5.}}));
    CHECK(opposite.itv() == IntervalVector({{-2.,-1.},{1.,3.}}));
  }

SECTION("Scalar multiplication. Multiplication by the inverse is used instead of; vector-expression division, which is not defined by the current API.")
  {
    AffineVarMainVector<AA> variables(
    		IntervalVector({{1.,2.},{-3.,-1.},{0.,4.}})
    );
    AffineMainVector<AA> x = variables;

    AffineMainVector<AA> twice = 2. * x;
    AffineMainVector<AA> half = 0.5 * x;
    AffineMainVector<AA> negative = -2. * x;

    CHECK(twice.itv() == IntervalVector({{2.,4.},{-6.,-2.},{0.,8.}}));
    CHECK(half.itv() == IntervalVector({{0.5,1.},{-1.5,-0.5},{0.,2.}}));
    CHECK(negative.itv() == IntervalVector({{-4.,-2.},{2.,6.},{-8.,0.}}));
  }

SECTION("Affine dependency: x-x is exactly zero.")
  {
    AffineVarMainVector<AA> variables(
      IntervalVector({{-2.,3.},{4.,5.},{-7.,-1.}})
    );
    AffineMainVector<AA> x =
      variables.template cast<AffineMain<AA>>();

    AffineMainVector<AA> result = x - x;

    CHECK(result.itv() == IntervalVector({{0.},{0.},{0.}}));
  }

SECTION("Complex linear expression.")
  {
    AffineVarMainVector<AA> vx(
      IntervalVector({{1.,2.},{-2.,1.},{4.,5.}})
    );
    AffineVarMainVector<AA> vy(
      IntervalVector({{-1.,1.},{3.,4.},{-2.,0.}})
    );
    AffineVarMainVector<AA> vz(
      IntervalVector({{2.,6.},{-4.,2.},{1.,3.}})
    );
    AffineMainVector<AA> x = vx.template cast<AffineMain<AA>>();
    AffineMainVector<AA> y = vy.template cast<AffineMain<AA>>();
    AffineMainVector<AA> z = vz.template cast<AffineMain<AA>>();

    AffineMainVector<AA> result = 2. * x - 3. * y + 0.5 * z;

    CHECK(result.size() == 3);
    CHECK(result.itv() == IntervalVector({{4.,6.},{-15.,-9.},{11.5,14.5}}));
  }

SECTION("Complex expression preserving dependencies.")
  {
    AffineVarMainVector<AA> variables(
      IntervalVector({{-5.,2.},{1.,7.},{-3.,4.}})
    );
    AffineMainVector<AA> x =
      variables.template cast<AffineMain<AA>>();

    AffineMainVector<AA> result =
      4. * x - 2. * (x + x) + (x - x);

    CHECK(result.itv() == IntervalVector({{0.},{0.},{0.}}));
  }

}


TEST_CASE("AffineMainVector product by an Interval scalar")
{
SECTION("Product of an affine vector by a degenerate interval.")
  {
    AffineVarMainVector<AA> vx(
      IntervalVector({{1.,2.},{-3.,-1.},{0.,4.}})
    );
    AffineMainVector<AA> x = vx;
    const Interval a(2.);

    AffineMainVector<AA> result = x * a;

    CHECK(result.size() == 3);
    CHECK(result.itv() == IntervalVector({{2.,4.},{-6.,-2.},{0.,8.}}));
  }

SECTION("Product with a degenerate interval in the opposite operand order.")
  {
    AffineVarMainVector<AA> vx(
      IntervalVector({{1.,2.},{-3.,-1.},{0.,4.}})
    );
    AffineMainVector<AA> x = vx;
    const Interval a(-2.);

    AffineMainVector<AA> result = a * x;

    CHECK(result.size() == 3);
    CHECK(result.itv() == IntervalVector({{-4.,-2.},{2.,6.},{-8.,0.}}));
  }

SECTION("Product by a non-degenerate positive interval encloses interval arithmetic.")
  {
    AffineVarMainVector<AA> vx(
      IntervalVector({{1.,2.},{-3.,-1.},{0.,4.}})
    );
    AffineMainVector<AA> x = vx;
    const Interval a(2.,3.);

    AffineMainVector<AA> result = x * a;
    const IntervalVector expected({{2.,6.},{-9.,-2.},{0.,12.}});

    CHECK(result.size() == 3);
    for(Index i = 0; i < result.size(); ++i)
      CHECK(result[i].itv().is_superset(expected[i]));
  }

SECTION("Product by an interval containing zero.")
  {
    AffineVarMainVector<AA> vx(
      IntervalVector({{1.,2.},{-3.,-1.},{0.,4.}})
    );
    AffineMainVector<AA> x = vx;
    const Interval a(-1.,2.);

    AffineMainVector<AA> result = a * x;
    const IntervalVector expected({{-2.,4.},{-6.,3.},{-4.,8.}});

    for(Index i = 0; i < result.size(); ++i)
      CHECK(result[i].itv().is_superset(expected[i]));
  }

SECTION("Empty interval propagation through a vector product.")
  {
    AffineVarMainVector<AA> vx(
      IntervalVector({{1.,2.},{-3.,-1.},{0.,4.}})
    );
    AffineMainVector<AA> x = vx;

    AffineMainVector<AA> result = x * Interval::empty();

    CHECK(result.size() == 3);
    CHECK(result.is_empty());
    CHECK(result.itv().is_empty());
  }

}


TEST_CASE("AffineMainVector dot product")
{
SECTION("Dot product with a point IntervalVector, via an AffineMainVector")
  {
    AffineVarMainVector<AA> vx(
      IntervalVector({{1.,2.},{-3.,-1.},{4.,5.}})
    );
    AffineMainVector<AA> x (3);
    x = vx;
    const IntervalVector weights({{2.},{-1.},{3.}});

    const AffineMain<AA> result = x.dot(weights);
    CHECK(result.itv() == Interval(15.,22.));
  }

SECTION("Dot product with a point IntervalVector, called directly on the AffineVarMainVector")
  {
    AffineVarMainVector<AA> vx(
      IntervalVector({{1.,2.},{-3.,-1.},{4.,5.}})
    );
    const IntervalVector weights({{2.},{-1.},{3.}});

    const AffineMain<AA> result = vx.dot(weights);
    CHECK(result.itv() == Interval(15.,22.));
  }

SECTION("Dot product with the IntervalVector on the left.")
  {
    AffineVarMainVector<AA> vx(Vector({2.,-3.,4.}));
    AffineMainVector<AA> x = vx;
    const IntervalVector weights({{5.},{2.},{-1.}});

    const AffineMain<AA> result = weights.dot(x);

    CHECK(result.itv() == Interval(0.));
  }

SECTION("Dot product with non-degenerate interval coefficients.")
  {
    AffineVarMainVector<AA> vx(
      IntervalVector({{1.,2.},{-3.,-1.}})
    );
    AffineMainVector<AA> x = vx;
    const IntervalVector weights({{2.,3.},{-1.,2.}});

    const AffineMain<AA> result = x.dot(weights);
    const Interval expected(-4.,9.);

    CHECK(result.itv().is_superset(expected));
  }

SECTION("Dependency in a dot product: [1,-1] dot [x,x] is exactly zero.")
  {
    AffineVarMainVector<AA> variable(1);
    variable = IntervalVector({{-2.,5.}});
    AffineMainVector<AA> repeated(2);
    repeated[0] = variable[0];
    repeated[1] = variable[0];
    const IntervalVector weights({{1.},{-1.}});

    const AffineMain<AA> result = repeated.dot(weights);

    CHECK(result.itv() == Interval(0.));
  }

SECTION("Empty coefficient propagation through a dot product.")
  {
    AffineVarMainVector<AA> vx(Vector({1.,2.,3.}));
    AffineMainVector<AA> x = vx;
    const IntervalVector weights(
      {Interval(1.),Interval::empty(),Interval(2.)}
    );

    const AffineMain<AA> result = x.dot(weights);

    CHECK(result.is_empty());
    CHECK(result.itv().is_empty());
  }

SECTION("Dot product between two affine vectors.")
  {
    AffineVarMainVector<AA> vx(Vector({2.,-3.,4.}));
    AffineVarMainVector<AA> vy(Vector({5.,6.,-2.}));
    AffineMainVector<AA> x = vx;
    AffineMainVector<AA> y = vy;

    const AffineMain<AA> result = x.dot(y);

    // 2*5 + (-3)*6 + 4*(-2) = -16.
    CHECK(result.itv() == Interval(-16.));
  }

SECTION("Dot product in reverse order between two affine vectors.")
  {
    AffineVarMainVector<AA> vx(Vector({2.,-3.,4.}));
    AffineVarMainVector<AA> vy(Vector({5.,6.,-2.}));
    AffineMainVector<AA> x = vx;
    AffineMainVector<AA> y = vy;

    const AffineMain<AA> direct = x.dot(y);
    const AffineMain<AA> reverse = y.dot(x);

    CHECK(direct.itv() == Interval(-16.));
    CHECK(reverse.itv() == Interval(-16.));
    CHECK(reverse.itv() == direct.itv());
  }

SECTION("Dot product between affine-variable vectors with disjoint per-component ranges")
  {
    AffineVarMainVector<AA> x(
      IntervalVector({{1.,2.},{-3.,-1.},{4.,5.}})
    );
    AffineVarMainVector<AA> y(
      IntervalVector({{2.,3.},{-2.,1.},{-1.,2.}})
    );

    const AffineMain<AA> direct = x.dot(y);
    const AffineMain<AA> reverse = y.dot(x);
    const Interval expected(-5.5,22.);
    CHECK(direct.itv().is_superset(expected));
    CHECK(reverse.itv().is_superset(expected));
    CHECK(direct == reverse);
  }

SECTION("Dot product between affine-variable vectors with overlapping per-component ranges")
  {
    AffineVarMainVector<AA> x(
      IntervalVector({{1.,2.},{3.,10.},{4.,5.}})
    );
    AffineVarMainVector<AA> y(
      IntervalVector({{2.,3.},{2.,10.},{1.,2.}})
    );

    const AffineMain<AA> direct = x.dot(y);
    const AffineMain<AA> reverse = y.dot(x);
    const Interval expected(12,116.);
    CHECK(direct.is_superset(expected));
    CHECK(reverse.is_superset(expected));
    CHECK(reverse == direct);
  }

SECTION("Dependency is preserved in both orders for IntervalVector initialization.")
  {
    AffineVarMainVector<AA> x(
      IntervalVector({{-2.,5.},{-2.,5.}})
    );
    AffineMainVector<AA> repeated(2);
    repeated[0] = x[0];
    repeated[1] = x[0];

    AffineVarMainVector<AA> weights(
      IntervalVector({{1.},{-1.}})
    );
    AffineMainVector<AA> affine_weights = weights;

    const AffineMain<AA> direct = repeated.dot(affine_weights);
    const AffineMain<AA> reverse = affine_weights.dot(repeated);

    CHECK(direct.itv() == Interval(0.));
    CHECK(reverse.itv() == Interval(0.));
  }

}


TEST_CASE("AffineMainVector component-wise nonlinear expressions")
{
SECTION("Component-wise nonlinear expressions on point affine variables.")
  {
    AffineVarMainVector<AA> x(Vector({2.,3.,4.}));
    AffineVarMainVector<AA> y(Vector({5.,6.,8.}));
    AffineMainVector<AA> result = x;

    result[0] = (x[0] + y[0]) * (y[0] - x[0]);
    result[1] = x[1] * y[1] + x[1] / y[1];
    result[2] = (x[2] * x[2] + y[2]) / x[2];

    CHECK(result == Approx<AffineMainVector<AA>>( IntervalVector({{21.},{18.5},{6.}}),ERROR));
  }

SECTION("Nonlinear expression 1: component-wise quadratic polynomial (implicit conversion from AffineVarMainVector)")
  {
    AffineVarMainVector<AA> x(Vector({-2.,3.,4.}));
    AffineMainVector<AA> result = x;

    for(Index i = 0 ; i < x.size() ; i++)
      result[i] = x[i] * x[i] + 2. * x[i] + 1.;

    CHECK(result.size() == 3);
    CHECK(result.itv() == IntervalVector({{1.},{16.},{25.}}));
  }

SECTION("Nonlinear expression 2: products between vector components (implicit conversion from AffineVarMainVector)")
  {
    AffineVarMainVector<AA> x(Vector({2.,-3.,4.}));
    AffineVarMainVector<AA> y(Vector({5.,6.,-2.}));
    AffineMainVector<AA> result = x;

    result[0] = x[0] * y[1] + x[1] * y[0];
    result[1] = (x[0] + y[0]) * (x[1] - y[1]);
    result[2] = x[2] * y[2] - x[0] * x[1];

    CHECK(result.size() == 3);
    CHECK(result.itv() == IntervalVector({{-3.},{-63.},{-2.}}));
  }

SECTION("Nonlinear expression 3: rational expressions with point denominators (implicit conversion from AffineVarMainVector)")
  {
    AffineVarMainVector<AA> x(Vector({2.,3.,4.}));
    AffineVarMainVector<AA> y(Vector({5.,6.,8.}));
    AffineMainVector<AA> result = x;

    result[0] = (x[0] * x[0] + y[0]) / x[0];
    result[1] = (x[1] + y[1]) / (y[1] - x[1]);
    result[2] = (x[2] * y[2]) / (x[2] + y[2]);

    CHECK(result.size() == 3);
    CHECK(result ==  Approx<AffineMainVector<AA>>(IntervalVector({{4.5},{3.},{8./3.}}),ERROR));
  }

SECTION("Nonlinear expression 4: cubic and mixed polynomial terms (implicit conversion from AffineVarMainVector)")
  {
    AffineVarMainVector<AA> x(Vector({-2.,3.,4.}));
    AffineVarMainVector<AA> y(Vector({1.,-1.,2.}));
    AffineMainVector<AA> result = x;

    result[0] = x[0] * x[0] * x[0] - y[0];
    result[1] = x[1] * x[1] * y[1] + 2. * x[1];
    result[2] = (x[2] + y[2]) * (x[2] - y[2]) + y[2] * y[2];

    CHECK(result.size() == 3);
    CHECK(result == Approx<AffineMainVector<AA>>(IntervalVector({{-9.},{-3.},{16.}}),ERROR));
  }

SECTION("Nonlinear expression 1: quadratic polynomial on point variables (explicit cast from AffineVarMainVector)")
  {
    AffineVarMainVector<AA> variables(Vector({-2.,3.,4.}));
    AffineMainVector<AA> x =
      variables.template cast<AffineMain<AA>>();
    AffineMainVector<AA> result = x;

    for(Index i = 0 ; i < x.size() ; i++)
      result[i] = x[i] * x[i] + 2. * x[i] + 1.;

    CHECK(result.itv() == IntervalVector({{1.},{16.},{25.}}));
  }

SECTION("Nonlinear expression 2: crossed products on point variables (explicit cast from AffineVarMainVector)")
  {
    AffineVarMainVector<AA> vx(Vector({2.,-3.,4.}));
    AffineVarMainVector<AA> vy(Vector({5.,6.,-2.}));
    AffineMainVector<AA> x = vx.template cast<AffineMain<AA>>();
    AffineMainVector<AA> y = vy.template cast<AffineMain<AA>>();
    AffineMainVector<AA> result = x;

    result[0] = x[0] * y[1] + x[1] * y[0];
    result[1] = (x[0] + y[0]) * (x[1] - y[1]);
    result[2] = x[2] * y[2] - x[0] * x[1];

    CHECK(result.itv() == IntervalVector({{-3.},{-63.},{-2.}}));
  }

SECTION("Nonlinear expression 3: rational component expressions (explicit cast from AffineVarMainVector)")
  {
    AffineVarMainVector<AA> vx(Vector({2.,3.,4.}));
    AffineVarMainVector<AA> vy(Vector({5.,6.,8.}));
    AffineMainVector<AA> x = vx.template cast<AffineMain<AA>>();
    AffineMainVector<AA> y = vy.template cast<AffineMain<AA>>();
    AffineMainVector<AA> result = x;

    result[0] = (x[0] * x[0] + y[0]) / x[0];
    result[1] = (x[1] + y[1]) / (y[1] - x[1]);
    result[2] = (x[2] * y[2]) / (x[2] + y[2]);

    CHECK(result == Approx<AffineMainVector<AA>>(IntervalVector({{4.5},{3.},{8./3.}}),ERROR));
  }

SECTION("Nonlinear expression 4: cubic and mixed polynomial terms (explicit cast from AffineVarMainVector)")
  {
    AffineVarMainVector<AA> vx(Vector({-2.,3.,4.}));
    AffineVarMainVector<AA> vy(Vector({1.,-1.,2.}));
    AffineMainVector<AA> x = vx.template cast<AffineMain<AA>>();
    AffineMainVector<AA> y = vy.template cast<AffineMain<AA>>();
    AffineMainVector<AA> result = x;

    result[0] = x[0] * x[0] * x[0] - y[0];
    result[1] = x[1] * x[1] * y[1] + 2. * x[1];
    result[2] = (x[2] + y[2]) * (x[2] - y[2]) + y[2] * y[2];

    CHECK(result.itv() == IntervalVector({{-9.},{-3.},{16.}}));
  }

}


TEST_CASE("AffineMainVector component-wise elementary functions")
{
SECTION("Component-wise sqr() and abs() on point (degenerate) affine variables.")
  {
    AffineVarMainVector<AA> vx(Vector({2.,-3.,-0.5}));
    AffineMainVector<AA> x = vx;
    AffineMainVector<AA> result = x;

    result[0] = sqr(x[0]);
    result[1] = sqr(x[1]);
    result[2] = abs(x[2]);

    CHECK(result == Approx<AffineMainVector<AA>>(IntervalVector({{4.},{9.},{0.5}}),ERROR));
  }

SECTION("Component-wise sqrt() encloses the exact interval-arithmetic result.")
  {
    AffineVarMainVector<AA> vx(
      IntervalVector({{4.,9.},{1.,4.}})
    );
    AffineMainVector<AA> x = vx;
    AffineMainVector<AA> result = x;

    result[0] = sqrt(x[0]);
    result[1] = sqrt(x[1]);

    CHECK(result[0].itv().is_superset(Interval(2.,3.)));
    CHECK(result[1].itv().is_superset(Interval(1.,2.)));
  }

SECTION("Component-wise exp() and log() enclose the exact interval-arithmetic result.")
  {
    AffineVarMainVector<AA> vx(
      IntervalVector({{0.,1.},{1.,2.}})
    );
    AffineMainVector<AA> x = vx;
    AffineMainVector<AA> result = x;

    result[0] = exp(x[0]);
    result[1] = log(x[1]);

    CHECK(result[0].itv().is_superset(Interval(std::exp(0.),std::exp(1.))));
    CHECK(result[1].itv().is_superset(Interval(std::log(1.),std::log(2.))));
  }

SECTION("Component-wise sin() and cos() enclose the exact interval-arithmetic; result, on a range where both functions are monotonic.")
  {
    AffineVarMainVector<AA> vx(
      IntervalVector({{0.,0.5},{0.,0.5}})
    );
    AffineMainVector<AA> x = vx;
    AffineMainVector<AA> result = x;

    result[0] = sin(x[0]);
    result[1] = cos(x[1]);

    CHECK(result[0].itv().is_superset(Interval(std::sin(0.),std::sin(0.5))));
    CHECK(result[1].itv().is_superset(Interval(std::cos(0.5),std::cos(0.))));
  }

}


TEST_CASE("AffineMainVector component-wise min, max and set operations")
{
SECTION("Component-wise min(), max(), intersection (&) and union (|) between two; point affine-variable vectors.")
  {
    AffineVarMainVector<AA> va(Vector({1.,5.,3.}));
    AffineVarMainVector<AA> vb(Vector({4.,2.,3.}));
    AffineMainVector<AA> a = va;
    AffineMainVector<AA> b = vb;

    IntervalVector result_min(3), result_max(3), result_union(3);
    for(Index i = 0 ; i < a.size() ; i++)
    {
      result_min[i] = min(a[i], b[i]);
      result_max[i] = max(a[i], b[i]);
      result_union[i] = a[i] | b[i];
    }
    const Interval inter_0 = a[0] & b[0]; // disjoint point values
    const Interval inter_2 = a[2] & b[2]; // equal point values

    CHECK(result_min == IntervalVector({{1.},{2.},{3.}}));
    CHECK(result_max == IntervalVector({{4.},{5.},{3.}}));
    CHECK(result_union == IntervalVector({{1.,4.},{2.,5.},{3.}}));
    CHECK(inter_0.is_empty());
    CHECK(inter_2 == Interval(3.));
  }

}


TEST_CASE("AffineMainVector: test des addon Eigen (lb, ub, mid)")
{
	AffineVarMainVector<AA> x(2);
	const IntervalVector values({
	    {-1.0, 1.0},
	    {2.0, 3.0},
	    {4.0, 5.0}
	});

	x = values;
	CHECK(x.lb() == Vector({-1.0, 2.0, 4.0}));
	CHECK(x.ub() == Vector({1.0, 3.0, 5.0}));
	CHECK(x.size() == 3);
	CHECK(x.mid() == Vector({0.0, 2.5, 4.5}));

	AffineMainVector<AA> y = x;
	CHECK(y.lb() == Vector({-1.0, 2.0, 4.0}));
	CHECK(y.ub() == Vector({1.0, 3.0, 5.0}));
	CHECK(y.size() == 3);
	CHECK(y.mid() == Vector({0.0, 2.5, 4.5}));

	y.set_empty();
	CHECK(y.is_empty());

	x.set_empty();
	CHECK(x.is_empty());

	x= values;
	CHECK(x.itv() == values);

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



TEST_CASE("Eigen AffineT vector transpose to row")
{
    const AffineTVector vector = make_point_vector_5();
    const AffineTRow row = vector.transpose();

    REQUIRE(row.rows() == 1);
    REQUIRE(row.cols() == 5);
    check_point(row[0], 1.0);
    check_point(row[4], 5.0);
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



TEST_CASE("Eigen scalar compound assignments - fixed-size AffineTVector and AffineTRow",
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
}


TEST_CASE("Eigen scalar compound assignments - fixed-size IntervalVector and IntervalRow",
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



TEST_CASE("Eigen AffineT reverse")
{
    const AffineTVector vector = make_point_vector_5();
    const AffineTVector reversed = vector.reverse();

    check_point(reversed[0], 5.0);
    check_point(reversed[1], 4.0);
    check_point(reversed[4], 1.0);
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


TEST_CASE("AffineMainVector and AffineMainRow stream output")
{
    SECTION("AffineMainVector, non-empty")
    {
        AffineTVector v = make_point_vector_5();
        std::ostringstream stream;
        stream << v;
        CAPTURE(stream.str());
        CHECK_FALSE(stream.str().empty());
        CHECK(stream.str().find("empty") == std::string::npos);
    }

    SECTION("AffineMainVector, empty")
    {
        AffineTVector v(2);
        v[0] = Interval::empty();
        std::ostringstream stream;
        stream << v;
        CHECK(stream.str() == "[ empty 2d box ]");
    }

    SECTION("AffineMainRow, non-empty")
    {
        AffineTRow row(3);
        row[0] = Interval(1.0);
        row[1] = Interval(2.0);
        row[2] = Interval(3.0);
        std::ostringstream stream;
        stream << row;
        CAPTURE(stream.str());
        CHECK_FALSE(stream.str().empty());
        CHECK(stream.str().find("empty") == std::string::npos);
    }

    SECTION("AffineMainRow, empty")
    {
        AffineTRow row(2);
        row[0] = Interval::empty();
        std::ostringstream stream;
        stream << row;
        CHECK(stream.str() == "[ empty row ]");
    }
}


