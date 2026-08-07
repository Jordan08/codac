/**
 * \file codac2_tests_AffineVector.cpp
 * Unit tests for AffineMainVector.
 * \date       2026
 * \author     Jordan Ninin
 * \copyright  Copyright 2026 Codac Team
 * \license    GNU Lesser General Public License (LGPL)
 */

#include <catch2/catch_test_macros.hpp>

#include "codac2_Affine.h"
#include "codac2_IntervalVector.h"
#include "codac2_Approx.h"

using namespace std;
using namespace codac2;

const double ERROR = std::numeric_limits<double>::epsilon()*10; //__builtin_powi(2.0, -50);

TEST_CASE("AffineMainVector")
{
	  typedef AF_Default AA;
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

  {
    AffineVarMainVector<AA> x(2);
    x = IntervalVector({{-2.,1.},{4.,7.}});

    AffineMainVector<AA> y = x;
    AffineMainVector<AA> z = y;

    CHECK(z.size() == 2);
    CHECK(z.itv() == y.itv());
    CHECK(z.itv() == IntervalVector({{-2.,1.},{4.,7.}}));
  }

  {
    AffineVarMainVector<AA> x(2);
    AffineVarMainVector<AA> y(2);
    x = IntervalVector({{1.,2.},{-3.,-1.}});
    y = IntervalVector({{4.,5.},{2.,6.}});

    AffineMainVector<AA> z = x + y;

    CHECK(z.size() == 2);
    CHECK(z.itv() == IntervalVector({{5.,7.},{-1.,5.}}));
  }

  {
    AffineVarMainVector<AA> x(2);
    AffineVarMainVector<AA> y(2);
    x = IntervalVector({{1.,3.},{-3.,1.}});
    y = x;

    AffineMainVector<AA> z = x - y;

    CHECK(z.size() == 2);
    CHECK(z.itv() == IntervalVector({{0.,0.},{0.,0.}}));
  }

  {
    AffineVarMainVector<AA> x(2);
    AffineVarMainVector<AA> y(2);
    x = IntervalVector({{1.,3.},{-3.,-1.}});
    y = IntervalVector({{4.,6.},{2.,6.}});

    AffineMainVector<AA> z = x - y;

    CHECK(z.size() == 2);
    CHECK(z.itv() == IntervalVector({{-3.,-3.},{-7.,-5.}}));
  }
  {
    AffineVarMainVector<AA> x(3);
    x = IntervalVector({{1.,2.},{-3.,-1.},{0.,4.}});

    AffineMainVector<AA> y = -x;

    CHECK(y.size() == 3);
    CHECK(y.itv() == IntervalVector({{-2.,-1.},{1.,3.},{-4.,0.}}));
  }

  {
    AffineVarMainVector<AA> x(3);
    x = IntervalVector({{1.,2.},{-3.,-1.},{0.,4.}});

    AffineMainVector<AA> y = 2. * x;

    CHECK(y.size() == 3);
    CHECK(y.itv() == IntervalVector({{2.,4.},{-6.,-2.},{0.,8.}}));
  }

  {
    AffineVarMainVector<AA> x(3);
    x = IntervalVector({{1.,2.},{-3.,-1.},{0.,4.}});

    AffineMainVector<AA> y = x / 2.;

    CHECK(y.size() == 3);
    CHECK(y.itv() == IntervalVector({{0.5,1.},{-1.5,-0.5},{0.,2.}}));
  }

  {
    AffineVarMainVector<AA> x(2);
    x = IntervalVector({{1.,2.},{3.,4.}});

    AffineMainVector<AA> y = x;
    y[0] = Interval(-1.,1.);

    CHECK(y.itv() == IntervalVector({{-1.,1.},{3.,4.}}));
    CHECK(x.itv() == IntervalVector({{1.,2.},{3.,4.}}));
  }

  {
    AffineVarMainVector<AA> x(3);
    x = IntervalVector({{1.,2.},{3.,4.},{5.,6.}});

    AffineMainVector<AA> y = x;

    CHECK(!y.is_empty());
    y.set_empty();
    CHECK(y.is_empty());
    CHECK(y.size() == 3);
  }

  {
    AffineVarMainVector<AA> x(3);
    x = IntervalVector({{1.,2.},Interval::empty(),{5.,6.}});

    AffineMainVector<AA> y = x;

    CHECK(y.is_empty());
    CHECK(y.size() == 3);
    CHECK(y.itv().is_empty());
  }

  // Construction through AffineVarMainVector from an IntervalVector.
  {
    IntervalVector box({{-5.,-2.},{0.},{7.,9.},{-oo,oo}});
    AffineVarMainVector<AA> x(box);
    AffineMainVector<AA> y = x;

    CHECK(y.size() == 4);
    CHECK(y.rows() == 4);
    CHECK(y.cols() == 1);
    CHECK(y.itv() == box);
  }

  // Construction through AffineVarMainVector from a real Vector.
  {
    Vector values({-2.,0.,3.5,8.});
    AffineVarMainVector<AA> x(values);
    AffineMainVector<AA> y = x;

    CHECK(y.size() == 4);
    CHECK(y.itv() == IntervalVector({{-2.},{0.},{3.5},{8.}}));
  }

  // The default affine variables are unbounded.
  {
    AffineVarMainVector<AA> x(3);
    AffineMainVector<AA> y = x;

    CHECK(y.size() == 3);
    CHECK(y.itv() == IntervalVector(3));
    CHECK(!y.is_empty());
  }

  // Parentheses, brackets and Eigen's coeff() address the same elements.
  {
    AffineVarMainVector<AA> x(3);
    x = IntervalVector({{-1.,1.},{2.,4.},{6.}});
    AffineMainVector<AA> y = x;

    CHECK(y[0].itv() == y(0).itv());
    CHECK(y[1].itv() == y.coeff(1).itv());
    CHECK(y(2).itv() == Interval(6.));
  }

  // Addition with a zero vector.
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

  // Subtraction from a zero vector.
  {
    AffineVarMainVector<AA> x(3);
    AffineVarMainVector<AA> zero(3);
    x = IntervalVector({{-2.,3.},{4.,5.},{-7.,-1.}});
    zero = IntervalVector({{0.},{0.},{0.}});

    AffineMainVector<AA> y = zero - x;

    CHECK(y.itv() == IntervalVector({{-3.,2.},{-5.,-4.},{1.,7.}}));
  }

  // Dependency is preserved: an affine variable minus itself is exactly zero.
  {
    AffineVarMainVector<AA> x(3);
    x = IntervalVector({{-2.,3.},{4.,5.},{-7.,-1.}});

    AffineMainVector<AA> y = x - x;

    CHECK(y.itv() == IntervalVector({{0.},{0.},{0.}}));
  }

  // Dependency is preserved in a linear expression.
  {
    AffineVarMainVector<AA> x(2);
    x = IntervalVector({{1.,2.},{-3.,4.}});

    AffineMainVector<AA> y = x + 2. * x;

    CHECK(y.itv() == IntervalVector({{3.,6.},{-9.,12.}}));
  }

  // Difference of two independent affine-variable vectors.
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

  // Multiplication by zero.
  {
    AffineVarMainVector<AA> x(3);
    x = IntervalVector({{-8.,10.},{-3.,4.},{5.}});

    AffineMainVector<AA> y = 0. * x;

    CHECK(y.itv() == IntervalVector({{0.},{0.},{0.}}));
  }

  // Multiplication by a negative scalar reverses interval bounds.
  {
    AffineVarMainVector<AA> x(3);
    x = IntervalVector({{1.,2.},{-3.,-1.},{0.,4.}});

    AffineMainVector<AA> y = -2. * x;

    CHECK(y.itv() == IntervalVector({{-4.,-2.},{2.,6.},{-8.,0.}}));
  }

  // Division by a negative scalar.
  {
    AffineVarMainVector<AA> x(3);
    x = IntervalVector({{2.,4.},{-6.,-2.},{0.,8.}});

    AffineMainVector<AA> y = x / -2.;

    CHECK(y.itv() == IntervalVector({{-2.,-1.},{1.,3.},{-4.,0.}}));
  }

  // Chained linear expression inspired by ex_affineform.cpp.
  {
    AffineVarMainVector<AA> x(3);
    x = IntervalVector({{1.,2.},{1.,3.},{1.,2.}});

    AffineMainVector<AA> y = x + 2. * x;
    AffineMainVector<AA> z = y - x;

    CHECK(y.itv() == IntervalVector({{3.,6.},{3.,9.},{3.,6.}}));
    CHECK(z.itv() == IntervalVector({{2.,4.},{2.,6.},{2.,4.}}));
  }

  // Assignment of an Eigen expression to an existing AffineMainVector.
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

  // Aliasing during assignment must not corrupt the expression.
  {
    AffineVarMainVector<AA> x(2);
    x = IntervalVector({{1.,2.},{-3.,4.}});
    AffineMainVector<AA> y = x;

    y = y + x;

    CHECK(y.itv() == IntervalVector({{2.,4.},{-6.,8.}}));
  }

  // Component assignment does not change the vector dimension.
  {
    AffineVarMainVector<AA> x(4);
    x = IntervalVector({{1.,2.},{3.,4.},{5.,6.},{7.,8.}});
    AffineMainVector<AA> y = x;

    y[1] = Interval(-10.,-8.);
    y[3] = Interval(0.);

    CHECK(y.size() == 4);
    CHECK(y.itv() == IntervalVector({{1.,2.},{-10.,-8.},{5.,6.},{0.}}));
  }

  // Setting one component to empty makes the whole vector empty.
  {
    AffineVarMainVector<AA> x(3);
    x = IntervalVector({{1.,2.},{3.,4.},{5.,6.}});
    AffineMainVector<AA> y = x;

    y[1].set_empty();

    CHECK(y.is_empty());
    CHECK(y.itv().is_empty());
    CHECK(y.size() == 3);
  }

  // Empty vectors keep their dimension through a copy.
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

  // Empty values propagate through vector addition and subtraction.
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

  // Unbounded components are preserved by itv().
  {
    AffineVarMainVector<AA> x(3);
    x = IntervalVector({{-oo,0.},{1.,oo},{-oo,oo}});
    AffineMainVector<AA> y = x;

    CHECK(y.itv() == IntervalVector({{-oo,0.},{1.,oo},{-oo,oo}}));
    CHECK(y.itv().is_unbounded());
  }

  // A vector of degenerate intervals behaves as a real vector.
  {
    AffineVarMainVector<AA> x(Vector({1.,-2.,3.}));
    AffineVarMainVector<AA> y(Vector({4.,5.,-6.}));

    AffineMainVector<AA> sum = x + y;
    AffineMainVector<AA> difference = x - y;

    CHECK(sum.itv() == IntervalVector({{5.},{3.},{-3.}}));
    CHECK(difference.itv() == IntervalVector({{-3.},{-7.},{9.}}));
  }


  // Complex linear expression with three independent affine vectors.
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

  // Nested expression combining sums, differences and scalar operations.
  {
    AffineVarMainVector<AA> x(2);
    AffineVarMainVector<AA> y(2);
    x = IntervalVector({{-2.,4.},{1.,3.}});
    y = IntervalVector({{5.,7.},{-4.,2.}});

    AffineMainVector<AA> result = (2. * (x + y) - (x - y)) / 3.;
    std::cout.precision(10);

    // Algebraically, the expression is (x + 3*y) / 3.
    CHECK(result == Approx<AffineMainVector<AA>>(IntervalVector({{13./3.,25./3.},{-11./3.,3.}}),ERROR));
  }

  // A complex expression must preserve affine dependencies.
  {
    AffineVarMainVector<AA> x(3);
    x = IntervalVector({{-5.,2.},{1.,7.},{-3.,4.}});

    AffineMainVector<AA> result = 4. * x - 2. * (x + x) + (x - x);

    CHECK(result.itv() == IntervalVector({{0.},{0.},{0.}}));
  }

  // Several intermediate AffineMainVector expressions can be chained.
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


  // Component-wise nonlinear expressions on point affine variables.
  {
    AffineVarMainVector<AA> x(Vector({2.,3.,4.}));
    AffineVarMainVector<AA> y(Vector({5.,6.,8.}));
    AffineMainVector<AA> result = x;

    result[0] = (x[0] + y[0]) * (y[0] - x[0]);
    result[1] = x[1] * y[1] + x[1] / y[1];
    result[2] = (x[2] * x[2] + y[2]) / x[2];

    CHECK(result == Approx<AffineMainVector<AA>>( IntervalVector({{21.},{18.5},{6.}}),ERROR));
  }

  // Empty vector propagation through unary minus and scalar operations.
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

  // Empty vector propagation for both operand orders.
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

  // Empty propagation through a nested expression.
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

  // Empty propagation through intermediate AffineMainVector results.
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

  // A single empty component propagates through every resulting vector.
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

  // Assigning an empty expression to a non-empty vector propagates emptiness.
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


  // Nonlinear expression 1: component-wise quadratic polynomial.
  {
    AffineVarMainVector<AA> x(Vector({-2.,3.,4.}));
    AffineMainVector<AA> result = x;

    for(Index i = 0 ; i < x.size() ; i++)
      result[i] = x[i] * x[i] + 2. * x[i] + 1.;

    CHECK(result.size() == 3);
    CHECK(result.itv() == IntervalVector({{1.},{16.},{25.}}));
  }

  // Nonlinear expression 2: products between vector components.
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

  // Nonlinear expression 3: rational expressions with point denominators.
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

  // Nonlinear expression 4: cubic and mixed polynomial terms.
  {
    AffineVarMainVector<AA> x(Vector({-2.,3.,4.}));
    AffineVarMainVector<AA> y(Vector({1.,-1.,2.}));
    AffineMainVector<AA> result = x;

    result[0] = x[0] * x[0] * x[0] - y[0];
    result[1] = x[1] * x[1] * y[1] + 2. * x[1];
    result[2] = (x[2] + y[2]) * (x[2] - y[2]) + y[2] * y[2];

    CHECK(result.size() == 3);
    CHECK(result.itv() == IntervalVector({{-9.},{-3.},{16.}}));
  }

}
