/**
 * \file codac2_tests_AffineForm_arithmetic.cpp
 * \brief Tests for AffineMain<T> arithmetic operators (+,-,*,/) and
 *        elementary/transcendental functions (sqrt, log, exp, pow, root,
 *        trigonometric and hyperbolic functions, floor/ceil/sign,
 *        min/max/intersection/union), including a sampling-based
 *        regression suite checking that every nonlinear function returns
 *        a sound affine enclosure of its reference interval image.
 */

#include <catch2/catch_test_macros.hpp>
#include <algorithm>
#include <cmath>
#include "codac2_Affine.h"
#include "codac2_Approx.h"
#include <limits>
#include <vector>
#include "codac2_Matrix.h"
#include <iostream>
using namespace codac2;


using AA = AF_Default;
using AffineT = AffineMain<AA>;
using AffineTMatrix = AffineMainMatrix<AA>;
using AffineTVector = AffineMainVector<AA>;
using AffineTVarVector = AffineVarMainVector<AA>;

const double MAX_DOUBLE = std::numeric_limits<double>::max();
const int SAMPLE_SIZE = 300;
const double ERROR = std::numeric_limits<double>::epsilon()*10; //__builtin_powi(2.0, -50);

#define piL Interval::pi().lb()
#define piU Interval::pi().ub()

template<class T>
void CHECK_affine_inclu(const AffineMain<T>& y_actual, const Interval& y_expected)
{
  if (y_expected.is_empty()) { CHECK(y_actual.is_empty()); return; }

  CHECK_FALSE(y_actual.is_empty());
  CHECK(y_expected.lb()>=y_actual.itv().lb());
  CHECK(y_expected.ub()<=y_actual.itv().ub());
}

template<class T>
void  CHECK_add(const Interval& xi, double z, const Interval& y_expected) {
	AffineVarMainVector<T> x(1);
	x[0] = xi;
	AffineMain<T>  y_actual=x[0]+z;
    CHECK(y_actual == Approx<AffineMain<T>>(y_expected, ERROR));

	// test the symmetrical case
	y_actual = z+x[0];
    CHECK(y_actual == Approx<AffineMain<T>>(y_expected, ERROR));

	// test the +=operator
	y_actual = x[0];
	y_actual += z;
    CHECK(y_actual == Approx<AffineMain<T>>(y_expected, ERROR));

	// test the +=operator in the other direction
	y_actual = z;
	y_actual += x[0];
    CHECK(y_actual == Approx<AffineMain<T>>(y_expected, ERROR));

	// test subtraction
	y_actual=(-x[0])-z;
    CHECK(y_actual == Approx<AffineMain<T>>(-y_expected, ERROR));

	// test the symmetrical case
	y_actual = (-z)-x[0];
    CHECK(y_actual == Approx<AffineMain<T>>(-y_expected, ERROR));

	// test the -=operator
	y_actual = -x[0];
	y_actual -= z;
    CHECK(y_actual == Approx<AffineMain<T>>(-y_expected, ERROR));

	// test the -=operator in the other direction
	y_actual = -z;
	y_actual -= x[0];
    CHECK(y_actual == Approx<AffineMain<T>>(-y_expected, ERROR));

	// test the linear simplification
	if (x[0].is_active() && std::fabs(z)<1.e100) {
		y_actual = z ;
		y_actual += x[0] ;
		y_actual -= x[0] ;
		y_actual += x[0] ;
	    CHECK(y_actual == Approx<AffineMain<T>>(y_expected, ERROR));

		y_actual = x[0] ;
		y_actual += z ;
		y_actual -= z ;
		y_actual += z ;
	    CHECK(y_actual == Approx<AffineMain<T>>(y_expected, ERROR));
	}
}

template<class T>
void  CHECK_add(const AffineMain<T>& x, const AffineMain<T>& z, const Interval& y_expected) {
	AffineMain<T>  y_actual=x+z;
    CHECK(y_actual == Approx<AffineMain<T>>(y_expected, ERROR));

	// test the symmetrical case
	y_actual = z+x;
    CHECK(y_actual == Approx<AffineMain<T>>(y_expected, ERROR));

	// test the +=operator
	y_actual = x;
	y_actual += z;
    CHECK(y_actual == Approx<AffineMain<T>>(y_expected, ERROR));

	// test the +=operator in the other direction
	y_actual = z;
	y_actual += x;
    CHECK(y_actual == Approx<AffineMain<T>>(y_expected, ERROR));

	// test subtraction
	y_actual=(-x)-z;
    CHECK(y_actual == Approx<AffineMain<T>>(-y_expected, ERROR));

	// test the symmetrical case
	y_actual = (-z)-x;
    CHECK(y_actual == Approx<AffineMain<T>>(-y_expected, ERROR));

	// test the -=operator
	y_actual = -x;
	y_actual -= z;
    CHECK(y_actual == Approx<AffineMain<T>>(-y_expected, ERROR));

	// test the -=operator in the other direction
	y_actual = -z;
	y_actual -= x;
    CHECK(y_actual == Approx<AffineMain<T>>(-y_expected, ERROR));

	// test the linear simplification
	if (x.is_active() && z.is_active()) {
		y_actual = z ;
		y_actual += x ;
		y_actual -= x ;
		y_actual += x ;
	    CHECK(y_actual == Approx<AffineMain<T>>(y_expected, ERROR));

		y_actual = x ;
		y_actual += z ;
		y_actual -= z ;
		y_actual += z ;
	    CHECK(y_actual == Approx<AffineMain<T>>(y_expected, ERROR));
	}
}

template<class T>
void  CHECK_add(const AffineMain<T>& x, const Interval& z, const Interval& y_expected) {
	AffineMain<T>  y_actual=x+z;
    CHECK(y_actual == Approx<AffineMain<T>>(y_expected, ERROR));

	// test the symmetrical case
	y_actual = z+x;
    CHECK(y_actual == Approx<AffineMain<T>>(y_expected, ERROR));

	// test the +=operator
	y_actual = x;
	y_actual += z;
    CHECK(y_actual == Approx<AffineMain<T>>(y_expected, ERROR));

	// test the +=operator in the other direction
	y_actual = z;
	y_actual += x;
    CHECK(y_actual == Approx<AffineMain<T>>(y_expected, ERROR));

	// test subtraction
	y_actual=(-x)-z;
    CHECK(y_actual == Approx<AffineMain<T>>(-y_expected, ERROR));

	// test the symmetrical case
	y_actual = (-z)-x;
    CHECK(y_actual == Approx<AffineMain<T>>(-y_expected, ERROR));

	// test the -=operator
	y_actual = -x;
	y_actual -= z;
    CHECK(y_actual == Approx<AffineMain<T>>(-y_expected, ERROR));

	// test the -=operator in the other direction
	y_actual = -z;
	y_actual -= x;
    CHECK(y_actual == Approx<AffineMain<T>>(-y_expected, ERROR));

	// test the linear simplification
	if (x.is_active() && !z.is_unbounded() && !z.is_empty()) {
		y_actual = z ;
		y_actual += x ;
		y_actual -= x ;
		y_actual += x ;
	    CHECK(y_actual == Approx<AffineMain<T>>(-y_expected, ERROR));

		y_actual = x ;
		y_actual += z ;
		y_actual -= z ;
		y_actual += z ;
	    CHECK(y_actual == Approx<AffineMain<T>>(-y_expected, ERROR));
	}
}


template<class T>
void  CHECK_add(const Interval& x, const Interval& z, const Interval& y_expected) {
	AffineVarMainVector<T> ax(2);
	ax[0]=x;
	ax[1]=z;
	CHECK_add<T>(ax[0],ax[1],y_expected);
}

template<class T>
void  CHECK_add_scalar(const Interval& x, double z, const Interval& y_expected) {
	CHECK_add<T>(x,z, y_expected);
	CHECK_add<T>((-x),-z, -y_expected);
}

template<class T>
void  CHECK_add_scalar(const Interval& x, const Interval& z, const Interval& y_expected) {
	AffineVarMainVector<T> xa(1,x);
	CHECK_add<T>(xa[0],z,y_expected);
	CHECK_add<T>((-xa[0]),-z, -y_expected);
}



template<class T>
void CHECK_mul(const AffineMain<T>& x, const AffineMain<T>& z, const Interval& y_expected)
{
  AffineMain<T> y_actual=x*z;
  CHECK_affine_inclu<T>(y_actual, y_expected);

  y_actual = z*x; // symmetrical case
  CHECK_affine_inclu<T>(y_actual, y_expected);

  y_actual = x;
  y_actual *= z; // *=operator
  CHECK_affine_inclu<T>(y_actual, y_expected);
}

template<class T>
void CHECK_mul(const Interval& x, const Interval& z, const Interval& y_expected)
{
  AffineVarMainVector<T> ax(2);
  ax[0]=x;
  ax[1]=z;
  CHECK_mul<T>(ax[0],ax[1],y_expected);
}

template<class T>
void CHECK_mul_scal(const Interval& x, double z, const Interval& y_expected)
{
  AffineVarMainVector<T> xa(1);
  xa[0] = x;
  CHECK(xa[0]*z == Approx<AffineMain<T>>(y_expected, ERROR));
  CHECK(z*xa[0] == Approx<AffineMain<T>>(y_expected, ERROR));
}

template<class T>
void CHECK_div(const AffineMain<T>& x, const AffineMain<T>& z, const Interval& y_expected)
{
  AffineMain<T> y_actual=x/z;
  CHECK_affine_inclu<T>(y_actual, y_expected);

  y_actual = x;
  y_actual /= z; // /=operator
  CHECK_affine_inclu<T>(y_actual, y_expected);
}

template<class T>
void CHECK_div(const Interval& x, const Interval& z, const Interval& y_expected)
{
  AffineVarMainVector<T> ax(2);
  ax[0]=x;
  ax[1]=z;
  CHECK_div<T>(ax[0],ax[1],y_expected);
}

template<class T>
void CHECK_div_scal(const Interval& x, double z, const Interval& y_expected)
{
  AffineVarMainVector<T> xa(1);
  xa[0] = x;
  AffineMain<T> y_actual=xa[0]/z;
  CHECK_affine_inclu<T>(y_actual, y_expected);

  y_actual = xa[0];
  y_actual /= z; // /=operator
  CHECK_affine_inclu<T>(y_actual, y_expected);
}

template<class T>
void CHECK_sqrt(const Interval& x, const Interval& y_expected)
{
  AffineVarMainVector<T> ax(1);
  ax[0] = x;
  CHECK_affine_inclu<T>(sqrt(ax[0]), y_expected);
}

template<class T>
void CHECK_log(const Interval& x, const Interval& y_expected)
{
  AffineVarMainVector<T> ax(1);
  ax[0] = x;
  CHECK_affine_inclu<T>(log(ax[0]), y_expected);
}

template<class T>
void CHECK_exp(const Interval& x, const Interval& y_expected)
{
  AffineVarMainVector<T> ax(1);
  ax[0] = x;
  CHECK_affine_inclu<T>(exp(ax[0]), y_expected);
}

template<class T>
void CHECK_trigo(const Interval& x, const Interval& expected)
{
  AffineVarMainVector<T> ax(1);
  ax[0] = x;
  CHECK_affine_inclu<T>(sin(ax[0]), expected);
  CHECK_affine_inclu<T>(sin(Interval::pi()-ax[0]), expected);
  CHECK_affine_inclu<T>(sin(ax[0]+Interval::two_pi()), expected);
  CHECK_affine_inclu<T>(sin(-ax[0]), -expected);
  CHECK_affine_inclu<T>(cos(ax[0]-Interval::half_pi()), expected);
  CHECK_affine_inclu<T>(cos(Interval::half_pi()-ax[0]), expected);
  CHECK_affine_inclu<T>(cos(ax[0]+Interval::half_pi()), -expected);
  CHECK_affine_inclu<T>(cos(ax[0]+Interval::two_pi()-Interval::half_pi()), expected);
}

template<class T>
void CHECK_pow(const Interval& x, int p, const Interval& y_expected)
{
  AffineVarMainVector<T> ax(1);
  ax[0] = x;
  CHECK_affine_inclu<T>(pow(ax[0],p),y_expected);
  CHECK_affine_inclu<T>(pow(-ax[0],p),(p%2==0)? y_expected : -y_expected);
}

template<class T>
void CHECK_root(const Interval& x, int p, const Interval& y_expected)
{
  AffineVarMainVector<T> ax(1);
  ax[0] = x;
  CHECK_affine_inclu<T>(root(ax[0],p), y_expected);
}

template<class T>
void CHECK_sqr(const Interval& x, const Interval& y_expected)
{
  AffineVarMainVector<T> ax(1);
  ax[0] = x;
  CHECK_affine_inclu<T>(sqr(ax[0]), y_expected);


}

template<class T>
void CHECK_asin(const Interval& x, const Interval& y_expected)
{
	AffineVarMainVector<T> ax(1);
	ax[0] = x;
	CHECK_affine_inclu<T>(asin(ax[0]), y_expected);

	if (!x.is_empty() && !x.is_unbounded() && !y_expected.is_empty() && !y_expected.is_unbounded()) {
		AffineVarMainVector<T> variables(IntervalVector({x}));
		const AffineMain<T> y = asin(variables[0]);
		for (int k = 0; k < 10000; ++k) {
			const double d = x.lb() + (x.ub() - x.lb()) * k / 10000.0;
			double dd = std::asin(d);
			if (!std::isnan(dd)) {
				CAPTURE(x, d,dd, y.itv());
				CHECK(y.contains(dd));
			}
		}
	}
}

template<class T>
void CHECK_acos(const Interval& x, const Interval& y_expected)
{
	AffineVarMainVector<T> ax(1);
	ax[0] = x;
	CHECK_affine_inclu<T>(acos(ax[0]), y_expected);

	if (!x.is_empty() && !x.is_unbounded() && !y_expected.is_empty() && !y_expected.is_unbounded()) {
		AffineVarMainVector<T> variables(IntervalVector({x}));
		const AffineMain<T> y = acos(variables[0]);
		for (int k = 0; k < 10000; ++k) {
			const double d = x.lb() + (x.ub() - x.lb()) * k / 10000.0;
			double dd = std::acos(d);
			if (!std::isnan(dd)) {
				CAPTURE(x, d,dd, y.itv());
				CHECK(y.contains(dd));
			}
		}
	}
}

template<class T>
void CHECK_atan(const Interval& x, const Interval& y_expected)
{
	AffineVarMainVector<T> ax(1);
	ax[0] = x;
	CAPTURE(x, atan(ax[0]).itv(),y_expected);
	CHECK_affine_inclu<T>(atan(ax[0]), y_expected);

	if (!x.is_empty() && !x.is_unbounded() && !y_expected.is_empty() && !y_expected.is_unbounded()) {
		AffineVarMainVector<T> variables(IntervalVector({x}));
		const AffineMain<T> y = atan(variables[0]);
		for (int k = 0; k < 10000; ++k) {
			const double d = x.lb() + (x.ub() - x.lb()) * k / 10000.0;
			double dd = std::atan(d);
			if (!std::isnan(dd)) {
				CAPTURE(x, d,dd, y.itv());
				CHECK(y.contains(dd));
			}
		}
	}
}

template<class T>
void CHECK_tanh(const Interval& x)
{
	const Interval y_expected = tanh(x);
	AffineVarMainVector<T> ax(1);
	ax[0] = x;
	CAPTURE(x, tanh(ax[0]).itv(),y_expected);
	CHECK_affine_inclu<T>(tanh(ax[0]), y_expected);

	if (!x.is_empty() && !x.is_unbounded() && !y_expected.is_empty() && !y_expected.is_unbounded()) {
		AffineVarMainVector<T> variables(IntervalVector({x}));
		const AffineMain<T> y = tanh(variables[0]);
		for (int k = 0; k < 10000; ++k) {
			const double d = x.lb() + (x.ub() - x.lb()) * k / 10000.0;
			double dd = std::tanh(d);
			if (!std::isnan(dd)) {
				CAPTURE(x, d,dd, y.itv());
				CHECK(y.contains(dd));
			}
		}
	}
}

template<class T>
void CHECK_sinh(const Interval& x)
{
	const Interval y_expected = sinh(x);
	AffineVarMainVector<T> ax(1);
	ax[0] = x;
	CAPTURE(x, sinh(ax[0]).itv(),y_expected);
	CHECK_affine_inclu<T>(sinh(ax[0]), y_expected);


	if (!x.is_empty() && !x.is_unbounded() && !y_expected.is_empty() && !y_expected.is_unbounded()) {
		AffineVarMainVector<T> variables(IntervalVector({x}));
		const AffineMain<T> y = sinh(variables[0]);
		for (int k = 0; k < 10000; ++k) {
			const double d = x.lb() + (x.ub() - x.lb()) * k / 10000.0;
			double dd = std::sinh(d);
			if (!std::isnan(dd)) {
				CAPTURE(x, d,dd, y.itv());
				CHECK(y.contains(dd));
			}
		}
	}
}

template<class T>
void CHECK_cosh(const Interval& x)
{
	const Interval y_expected = cosh(x);
	AffineVarMainVector<T> ax(1);
	ax[0] = x;
	CAPTURE(x, cosh(ax[0]).itv(),y_expected);
	CHECK_affine_inclu<T>(cosh(ax[0]), y_expected);


	if (!x.is_empty() && !x.is_unbounded() && !y_expected.is_empty() && !y_expected.is_unbounded()) {
		AffineVarMainVector<T> variables(IntervalVector({x}));
		const AffineMain<T> y = cosh(variables[0]);
		for (int k = 0; k < 10000; ++k) {
			const double d = x.lb() + (x.ub() - x.lb()) * k / 10000.0;
			double dd = std::cosh(d);
			if (!std::isnan(dd)) {
				CAPTURE(x, d,dd, y.itv());
				CHECK(y.contains(dd));
			}
		}
	}
}


TEST_CASE("AffineForm arithmetic operators and elementary functions")
{
  /* test:
   * =======
   *   operator+=(double d)
   *   operator+=(const Interval& x)
   *   operator+(const Interval& x, double d)
   *   operator+(double d,const Interval& x)
   *   operator+(const Interval& x1, const Interval& x2)
   *   operator-=(double d)
   *   operator-=(const Interval& x)
   *   operator-(const Interval& x, double d)
   *   operator-(double d,const Interval& x)
   *   operator-(const Interval& x1, const Interval& x2)
   */

  CHECK_add<AA>(Interval::empty(), Interval(0,1), Interval::empty());
  CHECK_add<AA>(Interval(0,1), Interval::empty(), Interval::empty());
  CHECK_add<AA>(Interval(-oo,1), Interval(0,1), Interval(-oo, 2));
  CHECK_add<AA>(Interval(1,oo), Interval(0,1), Interval(1,oo));
  CHECK_add<AA>(Interval(-oo,oo), Interval(0,1), Interval(-oo,oo));
  CHECK_add<AA>(Interval(MAX_DOUBLE,oo), 1, Interval(MAX_DOUBLE,oo));
  CHECK_add<AA>(Interval(MAX_DOUBLE,oo), -1, Interval(codac2::prev_float(MAX_DOUBLE),oo));
  CHECK_add<AA>(Interval(MAX_DOUBLE,oo), Interval(MAX_DOUBLE,oo), Interval(MAX_DOUBLE,oo));
  CHECK_add<AA>(Interval(MAX_DOUBLE,oo), -oo, Interval::empty());
  CHECK_add<AA>(Interval(MAX_DOUBLE,oo), oo, Interval::empty());
  CHECK_add<AA>(Interval(-oo,-MAX_DOUBLE), 1, Interval(-oo,codac2::next_float(-MAX_DOUBLE)));
  CHECK_add<AA>(Interval(-oo,-MAX_DOUBLE), -1, Interval(-oo,-MAX_DOUBLE));
  CHECK_add<AA>(Interval(-oo,-MAX_DOUBLE), Interval(-oo,-MAX_DOUBLE), Interval(-oo,-MAX_DOUBLE));

  CHECK_add_scalar<AA>(Interval::empty(), oo, Interval::empty());
  CHECK_add_scalar<AA>(Interval::empty(), 0, Interval::empty());
  CHECK_add_scalar<AA>(Interval(0,1), 1, Interval(1,2));
  CHECK_add_scalar<AA>(Interval(0,1), -oo, Interval::empty());
  CHECK_add_scalar<AA>(Interval(0,1), oo, Interval::empty());
  CHECK_add_scalar<AA>(Interval(-oo,1), 1, Interval(-oo,2));

  /* test:
   * =======
   *   operator*=(double d)
   *   operator*=(const Interval& x)
   *   operator*(const Interval& x, double d)
   *   operator*(double d,const Interval& x)
   *   operator*(const Interval& x1, const Interval& x2)
   */
  CHECK_mul<AA>(Interval::empty(), Interval(0,1), Interval::empty());
  CHECK_mul<AA>(Interval::zero(), Interval(), Interval::zero());
  CHECK_mul<AA>(Interval(-1,1), Interval(-oo,0), Interval());
  CHECK_mul<AA>(Interval(-oo,-1), Interval(-1,0), Interval(0,oo));
  CHECK_mul<AA>(Interval(-oo, 1), Interval(-1,0), Interval(-1,oo));
  CHECK_mul<AA>(Interval(0, 1), Interval(1,oo), Interval(0,oo));
  CHECK_mul<AA>(Interval(0, 1), Interval(-1,oo), Interval(-1,oo));
  CHECK_mul<AA>(Interval(-oo,-1), Interval(0,1), Interval(-oo,0));
  CHECK_mul<AA>(Interval(-oo, 1), Interval(0,1), Interval(-oo,1));
  CHECK_mul<AA>(Interval(0, 1), Interval(-oo,-1), Interval(-oo,0));
  CHECK_mul<AA>(Interval(0, 1), Interval(-oo,1), Interval(-oo,1));
  CHECK_mul<AA>(Interval(1,oo), Interval(0,1), Interval(0,oo));
  CHECK_mul<AA>(Interval(-1,oo), Interval(0,1), Interval(-1,oo));
  CHECK_mul<AA>(Interval(1,2), Interval(1,2), Interval(1,4));
  CHECK_mul<AA>(Interval(1,2), Interval(-2,3), Interval(-4,6));
  CHECK_mul<AA>(Interval(-1,1), Interval(-2,3), Interval(-3,3));
  CHECK_mul_scal<AA>(Interval(1,2), -oo, Interval::empty());
  CHECK_mul_scal<AA>(Interval(1,2), oo, Interval::empty());
  CHECK_mul_scal<AA>(Interval(1,2), -1, Interval(-2,-1));

  {
    AffineTVarVector a(2);   // contexte à 2 variables
    a[0] = Interval(1.,2.);
    AffineTVarVector b(4);   // contexte à 4 variables
    b[3] = Interval(3.,4.);

    AffineT p = a[0] * b[3];
    CHECK(p.itv().is_superset(Interval(1.,2.)*Interval(3.,4.)));
  }

  /* test:
   * =======
   *   operator/=(double d)
   *   operator/=(const Interval& x)
   *   operator/(const Interval& x, double d)
   *   operator/(double d,const Interval& x)
   *   operator/(const Interval& x1, const Interval& x2)
   */
  CHECK_div<AA>(Interval::empty(), Interval(0,1), Interval::empty());
  CHECK_div<AA>(Interval::zero(), Interval::zero(), Interval::empty());
  CHECK_div<AA>(Interval(1,2), Interval::zero(), Interval::empty());
  CHECK_div<AA>(Interval(), Interval::zero(), Interval::empty());
  CHECK_div<AA>(Interval::zero(), Interval(0,1), Interval::zero());
  CHECK_div<AA>(Interval::zero(), Interval(), Interval::zero());
  CHECK_div<AA>(Interval(6,12), Interval(2,3), Interval(2,6));
  CHECK_div<AA>(Interval(6,12), Interval(-3,-2), Interval(-6,-2));
  CHECK_div<AA>(Interval(6,12), Interval(-2,3), Interval());
  CHECK_div<AA>(Interval(-oo,-1), Interval(-1,0), Interval(1,oo));
  CHECK_div<AA>(Interval(-oo,-1), Interval(0,1), Interval(-oo,-1));
  CHECK_div<AA>(Interval(1,oo), Interval(-1,0), Interval(-oo,-1));
  CHECK_div<AA>(Interval(1,oo), Interval(0,1), Interval(1,oo));
  CHECK_div<AA>(Interval(-1,1), Interval(-1,1), Interval());
  CHECK_div_scal<AA>(Interval(1,2), -oo, Interval::empty());
  CHECK_div_scal<AA>(Interval(1,2), oo, Interval::empty());
  CHECK_div_scal<AA>(Interval(1,2), -1, Interval(-2,-1));


  CHECK_log<AA>(Interval::empty(), Interval::empty());
  CHECK_log<AA>(Interval(), Interval());
  CHECK_log<AA>(Interval(0,oo), Interval());
  CHECK_log<AA>(Interval(-oo,0), Interval::empty());
  CHECK_log<AA>(Interval(1,2), Interval(0,std::log(2)));
  CHECK_log<AA>(Interval(-1,1), Interval(-oo,0));
  {
    AffineTVarVector ax(1);
    ax[0] = Interval(0,codac2::next_float(0));
    CHECK(log(ax[0]).itv().ub() > -744.5);
  }
  CHECK_log<AA>(Interval(0,1), Interval(-oo,0));
  CHECK_log<AA>(Interval(1,oo), Interval(0,oo));
  CHECK_log<AA>(Interval(0), Interval::empty());
  CHECK_log<AA>(Interval(-2,-1), Interval::empty());

  CHECK_exp<AA>(Interval::empty(), Interval::empty());
  CHECK_exp<AA>(Interval(), Interval(0,oo));
  CHECK_exp<AA>(Interval(0,oo), Interval(1,oo));
  CHECK_exp<AA>(Interval(-oo,0), Interval(0,1));
  CHECK_exp<AA>(Interval(0,2), Interval(1,std::exp(2)));
  CHECK_exp<AA>(Interval(-1,1), Interval(std::exp(-1),std::exp(1)));
  CHECK_exp<AA>(Interval(1.e100,1.e111), Interval(MAX_DOUBLE,oo));
  CHECK_exp<AA>(Interval(MAX_DOUBLE,oo), Interval(MAX_DOUBLE,oo));
  CHECK_exp<AA>(Interval(0, MAX_DOUBLE), Interval(1,oo));

  {
    AffineTVarVector av(IntervalVector({{-8.}}) );

    CHECK(root(av[0], 3).itv() == Approx(Interval(-2.),ERROR));
    CHECK(root(av[0], 5).itv() == Approx(Interval(-std::pow(8., 0.2)),ERROR));

    av[0] = Interval(1., 5.);
    CHECK(root(av[0], 0).itv() ==	Interval::empty());
    CHECK(root(av[0].itv(), 0) ==	Interval::empty());
  }
  {
    AffineTVarVector ax(1);
    ax[0] = Interval(-8., 8.);
    AffineT y = root(ax[0], 3);
    CHECK(y.itv().is_superset(Interval(-2., 2.))); // root(-8,3)=-2, root(8,3)=2
    // vérifie que la forme reste affine (pas juste un intervalle plat) :
    CHECK_FALSE(y.is_degenerated());
  }

  /* test:
   * sin
   * cos
   */
  CHECK_trigo<AA>(Interval(), Interval(-1,1));
  CHECK_trigo<AA>(Interval::empty(), Interval::empty());
  CHECK_trigo<AA>(Interval(0,piU/2.0), Interval(0,1));
  CHECK_trigo<AA>(Interval(0,piU), Interval(0,1));
  CHECK_trigo<AA>(Interval(0,3*piU/2.0), Interval(-1,1));
  CHECK_trigo<AA>(Interval(piL,3*piU/2.0), Interval(-1,0));
  CHECK_trigo<AA>(Interval(0.5,1.5), Interval(std::sin(0.5),std::sin(1.5)));
  CHECK_trigo<AA>(Interval(1.5,3), Interval(std::sin(3.0),1));
  CHECK_trigo<AA>(Interval(3,4), Interval(std::sin(4.0),std::sin(3.0)));
  CHECK_trigo<AA>(Interval(3,5), Interval(-1,std::sin(3.0)));
  CHECK_trigo<AA>(Interval(3,2*piU+1.5), Interval(-1,std::sin(1.5)));
  CHECK_trigo<AA>(Interval(5,2*piU+1.5), Interval(std::sin(5.0),std::sin(1.5)));
  CHECK_trigo<AA>(Interval(5,2*piU+3), Interval(std::sin(5.0),1));


  {
	  AffineTVarVector ax(1);
	  CHECK_affine_inclu<AA>(tan(ax[0]), Interval());
  }
  {
	  AffineTVarVector ax(1);
	  ax[0] = (-Interval::pi()/4.0 | Interval::pi()/4.0);
	  CHECK_affine_inclu<AA>(tan(ax[0]), Interval(-1,1));
  }
  { // tan(pi/4,pi/2)=[1,+oo)
	  AffineTVarVector x(1); // upper bound of x is close to pi/2
	  x[0] = Interval(piL/4.0,(1-1e-10)*piL/2.0);
	  AffineT y=tan(x[0]);
	  CHECK(y.itv().lb()<=1.0);
	  CHECK(y.itv().ub()>1.e8); // upper bound of tan(x) is close to +oo
  }
  { // tan(-pi/2,pi/4)=(-oo,1]
	  AffineTVarVector ax(1);
	  ax[0] = Interval(-(1-1e-10)*piL/2.0,piL/4.0);
	  Interval y= (tan( ax[0] )).itv();
	  CHECK(y.lb()<= -1.e8); // lower bound is close to -oo
	  CHECK(y.ub()>= 1.0);
  }
  {
	  AffineTVarVector ax(1);
	  ax[0] = Interval::pi()/2.0;
	  CHECK_affine_inclu<AA>(tan(ax[0]),Interval());
  }
  {
	  AffineTVarVector ax(1);
	  ax[0] = Interval::pi();
	  AffineT y=tan(-ax[0]);
	  CHECK(y==Approx<AffineT>(Interval(0), ERROR));
	  CHECK(y.contains(0.));
	  CHECK(y.diam() < 1e-8);
  }
  {
	  AffineTVarVector ax(1);
	  ax[0] = Interval(3*piL/4.0 , 5*piU/4.0);
	  AffineT  y = tan(ax[0]);
	  CHECK_affine_inclu<AA>(y, Interval(-1,1));
	  CHECK(y.lb() == Approx(-1.0,1.e-5));
	  CHECK(y.ub() == Approx(1.0,1.e-5));

	  ax[0] = Interval(-oo,oo);
	  y= tan (ax[0]);
	  CHECK(y==Approx<AffineT>(Interval(), ERROR));

	  ax[0] = ((-Interval::pi()/4.0)|(Interval::pi()/4.0));
	  y= tan(ax[0]);
	  CHECK_affine_inclu<AA>(y, Interval(-1,1));
	  CHECK(y.lb() == Approx(-1.0,1.e-5));
	  CHECK(y.ub() == Approx(1.0,1.e-5));
  }
  {
	  // tan(pi/4,pi/2)=[1,+oo)
	  AffineTVarVector ax(1);
	  ax[0] = Interval(piL/4.0,(1-1e-10)*piL/2.0);
	  AffineT  y = tan(ax[0]);
	  CHECK(y.lb() <=1.0);
	  CHECK(y.ub() >= 1.e8); // upper bound of tan(x) is close to +oo
  }
  {
	  // tan(-pi/2,pi/4)=(-oo,1]
	  AffineTVarVector ax(1);
	  ax[0] = Interval(-(1-1e-10)*piL/2.0,piL/4.0);
	  AffineT  y = tan(ax[0]);
	  CHECK_affine_inclu<AA>(y, Interval(-1.e8,1));
  }

  {
    AffineTVarVector ax(1);
    ax[0] = (3*Interval::pi()/4.0 | 5*Interval::pi()/4.0);
    CHECK_affine_inclu<AA>(tan(ax[0]), Interval(-1,1));
  }

  CHECK_pow<AA>(Interval(), 4, Interval(0,oo));
  CHECK_pow<AA>(Interval::empty(), 4, Interval::empty());
  CHECK_pow<AA>(Interval(2,3), 4, Interval(16,81));
  CHECK_pow<AA>(Interval(-2,3), 4, Interval(0,81));
  CHECK_pow<AA>(Interval(-3,2), 4, Interval(0,81));
  CHECK_pow<AA>(Interval(2,oo), 4, Interval(16,oo));
  CHECK_pow<AA>(Interval(), 3, Interval());
  CHECK_pow<AA>(Interval::empty(), 3, Interval::empty());
  CHECK_pow<AA>(Interval(2,3), 3, Interval(8,27));
  CHECK_pow<AA>(Interval(-2,3), 3, Interval(-8,27));
  CHECK_pow<AA>(Interval(-3,2), 3, Interval(-27,8));
  CHECK_pow<AA>(Interval(2,oo), 3, Interval(8,oo));
  CHECK_pow<AA>(Interval(-10,10), -2, Interval(1.0/100,oo));

  CHECK_root<AA>(Interval(0,1), -1, Interval(1.0,oo));
  CHECK_root<AA>(Interval(-27,-8), 3, Interval(-3,-2));
  CHECK_root<AA>(Interval(-4,1), 2, Interval(0,1));
  CHECK_root<AA>(Interval(-8,1), 3, Interval(-2,1));

  /* test: sqrt(const Interval& x); */
  CHECK_sqrt<AA>(Interval(), Interval(0,oo));
  CHECK_sqrt<AA>(Interval(-oo,0), Interval::zero());
  CHECK_sqrt<AA>(Interval(-9,4), Interval(0,2));
  CHECK_sqrt<AA>(Interval(4,9), Interval(2,3));
  CHECK_sqrt<AA>(Interval(-9,-4), Interval::empty());
  CHECK_sqrt<AA>(Interval(-9,oo), Interval(0,oo));

  CHECK_sinh<AA>(Interval());
  CHECK_sinh<AA>(Interval(0,oo));
  CHECK_sinh<AA>(Interval(0,1));
  CHECK_sinh<AA>(Interval(1,oo));
  CHECK_sinh<AA>(Interval(1,1));
  CHECK_sinh<AA>(Interval(2,3));
  CHECK_sinh<AA>(Interval(4,5));
  CHECK_sinh<AA>(Interval(0,0.5));
  CHECK_sinh<AA>(Interval(-0.2,0.5));
  CHECK_sinh<AA>(Interval(0.6,0.9));
  CHECK_sinh<AA>(Interval(-0.5,-0.3));

  CHECK_cosh<AA>(Interval());
  CHECK_cosh<AA>(Interval(0,oo));
  CHECK_cosh<AA>(Interval(0,1));
  CHECK_cosh<AA>(Interval(1,oo));
  CHECK_cosh<AA>(Interval(1,1));
  CHECK_cosh<AA>(Interval(2,3));
  CHECK_cosh<AA>(Interval(4,5));
  CHECK_cosh<AA>(Interval(0,0.5));
  CHECK_cosh<AA>(Interval(-0.2,0.5));
  CHECK_cosh<AA>(Interval(0.6,0.9));
  CHECK_cosh<AA>(Interval(-0.5,-0.3));

  /* test: tanh(const AffineMain<T>& x); defined everywhere, increasing, range (-1,1) */
  CHECK_tanh<AA>(Interval::empty());
  CHECK_tanh<AA>(Interval());
  CHECK_tanh<AA>(Interval(0,oo));
  CHECK_tanh<AA>(Interval(-oo,0));
  CHECK_tanh<AA>(Interval(0,0.5));
  CHECK_tanh<AA>(Interval(-0.2,0.5));
  CHECK_tanh<AA>(Interval(0.6,0.9));
  CHECK_tanh<AA>(Interval(-0.5,-0.3));

  { // issue 248
    AffineTVarVector aff(1);
    aff[0] = Interval(-1.57079632679489678, 1.1780972450961728626);
    CHECK_FALSE((tan(aff[0]).is_empty()));
  }
}


TEST_CASE("AffineForm additional operations")
{

  /* test: sqr(const AffineMain<T>& x); */
  CHECK_sqr<AA>(Interval::empty(), Interval::empty());
  CHECK_sqr<AA>(Interval(), Interval(0,oo));
  CHECK_sqr<AA>(Interval(2,3), Interval(4,9));
  CHECK_sqr<AA>(Interval(-2,3), Interval(0,9));
  CHECK_sqr<AA>(Interval(-3,2), Interval(0,9));
  CHECK_sqr<AA>(Interval(2,oo), Interval(4,oo));

  /* test: asin(const AffineMain<T>& x); domain [-1,1] */
  CHECK_asin<AA>(Interval::empty(), Interval::empty());
  CHECK_asin<AA>(Interval(2,3), Interval::empty()); // entirely out of domain
  CHECK_asin<AA>(Interval(0,1), Interval(0,Interval::half_pi().ub()));
  CHECK_asin<AA>(Interval(-1,0), Interval(-Interval::half_pi().ub(),0));
  CHECK_asin<AA>(Interval(-2,0.5), Interval(-Interval::half_pi().ub(),std::asin(0.5))); // truncated on the left

  {
    for (Interval x_itv :
         {Interval(-0.95, 0.20),
          Interval(-0.75, 0.85),
          Interval(0.10, 0.99)}) {

      AffineTVarVector variables(IntervalVector({x_itv}));
      const AffineT y = asin(variables[0]);

      for (int k = 0; k < 10000; ++k) {
        const double x = x_itv.lb() + (x_itv.ub() - x_itv.lb()) * k / 10000.0;
        CAPTURE(x_itv, x, y.itv());
        CHECK(y.contains(std::asin(x)));
      }
    }
  }

  /* test: acos(const AffineMain<T>& x); domain [-1,1], decreasing */
  CHECK_acos<AA>(Interval::empty(), Interval::empty());
  CHECK_acos<AA>(Interval(2,3), Interval::empty()); // entirely out of domain
  CHECK_acos<AA>(Interval(-1,1), Interval(0,Interval::pi().ub()));
  CHECK_acos<AA>(Interval(0,1), Interval(0,Interval::half_pi().ub()));

  /* test: atan(const AffineMain<T>& x); defined everywhere, increasing */
  CHECK_atan<AA>(Interval::empty(), Interval::empty());
  CHECK_atan<AA>(Interval(), Interval(-Interval::half_pi().ub(),Interval::half_pi().ub()));
  CHECK_atan<AA>(Interval(0,1), Interval(0,std::atan(1.)));
  CHECK_atan<AA>(Interval(-1,0), Interval(-std::atan(1.),0));
  {
    AffineTVarVector ax(1);
    ax[0] = Interval(2., 3.);
    CHECK(acosh(ax[0]).itv().is_superset(Interval(std::acosh(2.), std::acosh(3.))));
  }
  {
    AffineTVarVector ax(1);
    ax[0] = Interval(-1., 1.);
    CHECK(asinh(ax[0]).itv().is_superset(Interval(std::asinh(-1.), std::asinh(1.))));
  }
  {
    AffineTVarVector ax(1);
    ax[0] = Interval(-0.5, 0.5);
    CHECK(atanh(ax[0]).itv().is_superset(Interval(std::atanh(-0.5), std::atanh(0.5))));
  }

  /* test: min, max, operator&, operator| between two affine forms
   * (and between an affine form and a plain Interval) */
  {
    AffineTVarVector ax(2);
    ax[0] = Interval(1,3);
    ax[1] = Interval(2,5);

    CHECK(min(ax[0],ax[1]) == Interval(1,3));
    CHECK(max(ax[0],ax[1]) == Interval(2,5));
    CHECK((ax[0] & ax[1]) == Interval(2,3));
    CHECK((ax[0] | ax[1]) == Interval(1,5));
  }

  {
    // disjoint intervals: intersection is empty, union is the hull
    AffineTVarVector ax(2);
    ax[0] = Interval(0,1);
    ax[1] = Interval(2,3);

    CHECK((ax[0] & ax[1]).is_empty());
    CHECK((ax[0] | ax[1]) == Interval(0,3));
  }

  {
    // mixed operand: AffineMain and plain Interval, both operand orders
    AffineTVarVector ax(1);
    ax[0] = Interval(1,4);
    const Interval b(2,6);

    CHECK((ax[0] & b) == Interval(2,4));
    CHECK((b & ax[0]) == Interval(2,4));
    CHECK((ax[0] | b) == Interval(1,6));
    CHECK((b | ax[0]) == Interval(1,6));
    CHECK(min(ax[0],b) == Interval(1,4));
    CHECK(max(ax[0],b) == Interval(2,6));
    CHECK(min(b,ax[0]) == Interval(1,4));
    CHECK(max(b,ax[0]) == Interval(2,6));
  }
}



TEST_CASE("Critical nonlinear regressions")
{



  // Racine impaire d'un singleton négatif
  {
    AffineTVarVector x(IntervalVector({{-8.}}));
    CHECK(root(x[0], 3).itv() == Approx(Interval(-2.),ERROR));
  }

  // Indice zéro
  {
    AffineTVarVector x(IntervalVector({{1., 2.}}));
    CHECK(root(x[0], 0).is_empty());
  }

  // floor ne doit pas conserver une fausse dépendance
  {
    AffineTVarVector variables(
      IntervalVector({{1., 2.}})
    );

    const AffineT x = variables[0];
    const AffineT residual = floor(x) - x;

    CHECK_FALSE(residual.itv() == Interval(0.));
    CHECK(Interval(-1., 0.).is_subset(residual.itv()));
  }
}





TEST_CASE("Large and exceptional arguments")
{

	const double int_max =
			static_cast<double>(std::numeric_limits<int>::max());

	const std::vector<Interval> inputs = {
			Interval(1e100, 1e100 + 1e90),
			Interval(-1e100, -1e100 + 1e90),
			Interval(int_max, int_max + 2.0),
			Interval(-oo, oo),
			Interval::empty()
	};

	for (const Interval& input : inputs) {
		AffineTVarVector x(IntervalVector({input}));

		CHECK_NOTHROW(sin(x[0]));
		CHECK_NOTHROW(cos(x[0]));
		CHECK_NOTHROW(floor(x[0]));
		CHECK_NOTHROW(ceil(x[0]));
	}
}




TEST_CASE("floor and ceil of non-integer singleton")
{

	AffineTVarVector variables(
			IntervalVector({{1.5}, {-1.5}, {2.0}})
	);

	CHECK(floor(variables[0]) == Interval(1.0));
	CHECK(ceil(variables[0]) == Interval(2.0));

	CHECK(floor(variables[1]) == Interval(-2.0));
	CHECK(ceil(variables[1]) == Interval(-1.0));

	CHECK(floor(variables[2]) == Interval(2.0));
	CHECK(ceil(variables[2]) == Interval(2.0));
}

namespace
{

double odd_root_value(double x, int n)
{
  return std::copysign(
    std::pow(std::fabs(x), 1.0/static_cast<double>(n)), x
  );
}

template<class Function>
void CHECK_SAMPLED_ENCLOSURE(
  const Interval& input,
  const AffineT& output,
  Function function,
  int sample_count = SAMPLE_SIZE)
{
  REQUIRE_FALSE(output.is_empty());
  REQUIRE(sample_count > 100);

  for (int k = 0; k <= sample_count; ++k) {
    double x;

    if (k == 0) {
      x = input.lb();
    }
    else if (k == sample_count) {
      x = input.ub();
    }
    else {
      const long double ratio =
        static_cast<long double>(k) /
        static_cast<long double>(sample_count);

      const long double lower =
        static_cast<long double>(input.lb());
      const long double upper =
        static_cast<long double>(input.ub());

      const long double interpolated =
        (1.0L-ratio)*lower + ratio*upper;

      x = static_cast<double>(interpolated);
      x = std::max(input.lb(), std::min(input.ub(), x));
    }

    CAPTURE(input, output.itv(), k, x);
    REQUIRE(input.contains(x));

    const double expected = function(x);
    CAPTURE(expected);
    CHECK(output.contains(expected));
  }
}

template<class IntervalFunction>
void CHECK_SAMPLED_INTERVAL_ENCLOSURE(
  const Interval& input,
  const AffineT& output,
  const Interval& global_reference,
  IntervalFunction interval_function,
  int sample_count = SAMPLE_SIZE)
{
  REQUIRE_FALSE(output.is_empty());
  REQUIRE(sample_count > 100);

  const Interval output_interval = output.itv();
  int tested_point_count = 0;

  for (int k = 0; k <= sample_count; ++k) {
    double x;

    if (k == 0) {
      x = input.lb();
    }
    else if (k == sample_count) {
      x = input.ub();
    }
    else {
      const long double ratio =
        static_cast<long double>(k) /
        static_cast<long double>(sample_count);
      const long double lower =
        static_cast<long double>(input.lb());
      const long double upper =
        static_cast<long double>(input.ub());
      const long double interpolated =
        (1.0L-ratio)*lower + ratio*upper;

      x = static_cast<double>(interpolated);
      x = std::max(input.lb(), std::min(input.ub(), x));
    }

    CAPTURE(input, output_interval, k, x);
    REQUIRE(input.contains(x));

    // Both operands are certified enclosures of the same point image:
    // interval_function([x]) is the local singleton evaluation, while
    // global_reference encloses the image of the complete valid domain.
    // Their intersection remains a certified enclosure and removes harmless
    // overshoots such as cos([0]) extending by a few ULP beyond 1.
    const Interval local_reference = interval_function(Interval(x));

    // An empty singleton image means that x itself is outside the real
    // pointwise domain. This occurs, for example, at x=0 for log: the
    // interval image log([0,b]) contains -infinity as a one-sided limit,
    // while log([0,0]) is empty because log(0) is not a real value.
    // Such a point imposes no pointwise inclusion condition.
    if (local_reference.is_empty()) {
      continue;
    }

    const Interval expected = local_reference & global_reference;
    CAPTURE(local_reference, global_reference, expected);
    REQUIRE_FALSE(expected.is_empty());
    CHECK(expected.is_subset(output_interval));
    ++tested_point_count;
  }

  // Prevent a wrongly specified domain or interval function from making the
  // whole sampling test vacuous.
  REQUIRE(tested_point_count > 100);
}

} // namespace


TEST_CASE("Odd affine roots crossing zero enclose sampled values")
{
  const std::vector<Interval> inputs = {
    Interval(-8.0, 27.0),
    Interval(-1.e-12, 1.e6),
    Interval(-1.e6, 1.e-12),
    Interval(-1.e-300, 1.e-300),
    Interval(-1.0, 1024.0),
    Interval(-125.0, 32.0)
  };
  const std::vector<int> orders = {3, 5, 7, 9, 101};

  for (const Interval& input : inputs) {
    for (const int order : orders) {
      CAPTURE(input, order);
      AffineTVarVector variables(IntervalVector({input}));
      const AffineT result = root(variables[0], order);
      const Interval reference = root(input, order);

      CHECK_affine_inclu<AA>(result, reference);
      CHECK_SAMPLED_ENCLOSURE(
        input,
        result,
        [order](double x) { return odd_root_value(x, order); }
      );
    }
  }
}


TEST_CASE("Odd affine roots crossing zero preserve symmetry")
{
  const std::vector<int> orders = {3, 5, 7, 9};

  for (const int order : orders) {
    CAPTURE(order);
    const Interval input(-64.0, 64.0);
    AffineTVarVector variables(IntervalVector({input}));
    const AffineT result = root(variables[0], order);

    CHECK(result.contains(0.0));
    CHECK(result.lb() == Approx(-result.ub(), 1.e-10));
    CHECK_SAMPLED_ENCLOSURE(
      input,
      result,
      [order](double x) { return odd_root_value(x, order); }
    );
  }
}


TEST_CASE("Partially valid domains keep a sound affine enclosure")
{
  SECTION("sqrt with a negative invalid part") {
    const Interval input(-9.0, 4.0);
    const Interval valid_domain = input & Interval(0.0, oo);
    AffineTVarVector variables(IntervalVector({input}));
    const AffineT result = sqrt(variables[0]);

    CHECK_affine_inclu<AA>(result, sqrt(valid_domain));
    CHECK_SAMPLED_ENCLOSURE(
      valid_domain, result, [](double x) { return std::sqrt(x); }
    );
  }

  SECTION("asin with values below its domain") {
    const Interval input(-2.0, 0.75);
    const Interval valid_domain = input & Interval(-1.0, 1.0);
    AffineTVarVector variables(IntervalVector({input}));
    const AffineT result = asin(variables[0]);

    CHECK_affine_inclu<AA>(result, asin(valid_domain));
    CHECK_SAMPLED_ENCLOSURE(
      valid_domain, result, [](double x) { return std::asin(x); }
    );
  }

  SECTION("acos with values above its domain") {
    const Interval input(-0.75, 2.0);
    const Interval valid_domain = input & Interval(-1.0, 1.0);
    AffineTVarVector variables(IntervalVector({input}));
    const AffineT result = acos(variables[0]);

    CHECK_affine_inclu<AA>(result, acos(valid_domain));
    CHECK_SAMPLED_ENCLOSURE(
      valid_domain, result, [](double x) { return std::acos(x); }
    );
  }

  SECTION("acosh with values below its domain") {
    const Interval input(0.0, 4.0);
    const Interval valid_domain = input & Interval(1.0, oo);
    AffineTVarVector variables(IntervalVector({input}));
    const AffineT result = acosh(variables[0]);

    CHECK_affine_inclu<AA>(result, acosh(valid_domain));
    CHECK_SAMPLED_ENCLOSURE(
      valid_domain, result, [](double x) { return std::acosh(x); }
    );
  }
}


TEST_CASE("Entirely invalid nonlinear domains remain empty")
{
  AffineTVarVector variables(
    IntervalVector({{-9.0, -4.0}, {2.0, 3.0}, {-3.0, 0.5}})
  );

  CHECK(sqrt(variables[0]).is_empty());
  CHECK(asin(variables[1]).is_empty());
  CHECK(acos(variables[1]).is_empty());
  CHECK(acosh(variables[2]).is_empty());
}


template<class AffineFunction, class IntervalFunction, class ScalarFunction>
void CHECK_NONLINEAR_OPERATION(
  const char* operation_name,
  const Interval& input,
  const Interval& valid_domain,
  AffineFunction affine_function,
  IntervalFunction interval_function,
  ScalarFunction /* scalar_function */,
  int sample_count = SAMPLE_SIZE)
{
  CAPTURE(operation_name, input, valid_domain);

  AffineTVarVector variables(IntervalVector({input}));
  const AffineT result = affine_function(variables[0]);
  const Interval reference = interval_function(valid_domain);

  CHECK_affine_inclu<AA>(result, reference);
  if (!valid_domain.is_empty() && !valid_domain.is_degenerated() &&
      !valid_domain.is_unbounded()) {
    // Use the certified interval image of each singleton as the oracle.
    // A scalar std::* result may differ by a few ULP from the outward-rounded
    // CODAC interval implementation and is therefore not a rigorous oracle.
    CHECK_SAMPLED_INTERVAL_ENCLOSURE(
      valid_domain, result, reference, interval_function, sample_count
    );
  }
}


template<class AffineFunction, class IntervalFunction>
void CHECK_NONLINEAR_INTERVAL_IMAGE(
  const char* operation_name,
  const Interval& input,
  AffineFunction affine_function,
  IntervalFunction interval_function)
{
  CAPTURE(operation_name, input);
  AffineTVarVector variables(IntervalVector({input}));
  const AffineT result = affine_function(variables[0]);
  CHECK_affine_inclu<AA>(result, interval_function(input));
}


TEST_CASE("All nonlinear affine functions enclose reference values")
{
  using Mode = AffineT::Affine_Mode;
  const std::vector<Mode> modes = {
    AffineT::AF_Lin_Chebyshev,
    AffineT::AF_Lin_MinRange
  };

  for (const Mode mode : modes) {
    AffineT::change_mode(mode);
    CAPTURE(static_cast<int>(mode));

    // Reciprocal: both connected components of R\\{0}, including points
    // close to the singularity.
    for (const Interval& input : {
      Interval(1.e-8, 3.0), Interval(-4.0, -1.e-8),
      Interval(1.e-150, 1.e-140), Interval(-1.e8, -1.e-4)}) {
      CHECK_NONLINEAR_OPERATION(
        "inv", input, input,
        [](const AffineT& x) { return inv(x); },
        [](const Interval& x) { return 1.0/x; },
        [](double x) { return 1.0/x; }
      );
    }

    // Square and absolute value, including intervals crossing zero.
    for (const Interval& input : {
      Interval(-4.0, 3.0), Interval(-1.e-12, 1.e12),
      Interval(-3.0, -0.25), Interval(0.25, 5.0)}) {
      CHECK_NONLINEAR_OPERATION(
        "sqr", input, input,
        [](const AffineT& x) { return sqr(x); },
        [](const Interval& x) { return sqr(x); },
        [](double x) { return x*x; }
      );
      CHECK_NONLINEAR_OPERATION(
        "abs", input, input,
        [](const AffineT& x) { return abs(x); },
        [](const Interval& x) { return abs(x); },
        [](double x) { return std::fabs(x); }
      );
    }

    // Square root, including a large invalid negative part.
    for (const Interval& input : {
      Interval(0.0, 4.0), Interval(1.e-24, 1.e6),
      Interval(-1.e12, 1.e-12), Interval(-1.e300, 4.0)}) {
      const Interval domain = input & Interval(0.0, oo);
      CHECK_NONLINEAR_OPERATION(
        "sqrt", input, domain,
        [](const AffineT& x) { return sqrt(x); },
        [](const Interval& x) { return sqrt(x); },
        [](double x) { return std::sqrt(x); }, 8192
      );
    }

    // Exponential and logarithm. The chosen exponential intervals avoid a
    // floating overflow in the scalar reference while still spanning very
    // different scales.
    for (const Interval& input : {
      Interval(-20.0, 20.0), Interval(-700.0, -600.0),
      Interval(-1.e-12, 1.e-12), Interval(1.0, 10.0)}) {
      CHECK_NONLINEAR_OPERATION(
        "exp", input, input,
        [](const AffineT& x) { return exp(x); },
        [](const Interval& x) { return exp(x); },
        [](double x) { return std::exp(x); }
      );
    }
    for (const Interval& input : {
      Interval(1.e-300, 1.e-12), Interval(1.e-12, 1.e12),
      Interval(0.25, 4.0), Interval(-1.e6, 3.0)}) {
      const Interval domain = input & Interval(0.0, oo);
      CHECK_NONLINEAR_OPERATION(
        "log", input, domain,
        [](const AffineT& x) { return log(x); },
        [](const Interval& x) { return log(x); },
        [](double x) { return std::log(x); }, 8192
      );
    }

    // Integer and real powers, plus odd roots crossing zero.
    for (const Interval& input : {
      Interval(-8.0, 27.0), Interval(-1.e6, 1.e-12),
      Interval(-3.0, 2.0)}) {
      CHECK_NONLINEAR_OPERATION(
        "pow(x,3)", input, input,
        [](const AffineT& x) { return pow(x, 3); },
        [](const Interval& x) { return pow(x, 3); },
        [](double x) { return x*x*x; }
      );
      CHECK_NONLINEAR_OPERATION(
        "root(x,3)", input, input,
        [](const AffineT& x) { return root(x, 3); },
        [](const Interval& x) { return root(x, 3); },
        [](double x) { return odd_root_value(x, 3); }, 8192
      );
    }
    for (const Interval& input : {
      Interval(1.e-12, 4.0), Interval(0.25, 1.e6)}) {
      CHECK_NONLINEAR_OPERATION(
        "pow(x,1.5)", input, input,
        [](const AffineT& x) { return pow(x, 1.5); },
        [](const Interval& x) { return pow(x, 1.5); },
        [](double x) { return std::pow(x, 1.5); }
      );
    }

    // Interval and affine exponents deliberately lose the exponent's affine
    // dependency, but their interval image must still be enclosed.
    for (const Interval& input : {
      Interval(0.25, 4.0), Interval(1.e-6, 1.e3)}) {
      const Interval exponent(0.5, 2.0);
      CHECK_NONLINEAR_INTERVAL_IMAGE(
        "pow(x,Interval)", input,
        [exponent](const AffineT& x) { return pow(x, exponent); },
        [exponent](const Interval& x) { return pow(x, exponent); }
      );

      AffineTVarVector variables(
        IntervalVector({input, exponent})
      );
      const AffineT result = pow(variables[0], variables[1]);
      CHECK_affine_inclu<AA>(result, pow(input, exponent));
    }

    // Trigonometric functions. Tangent intervals stay inside one continuity
    // branch, while sine and cosine additionally cover several extrema.
    for (const Interval& input : {
      Interval(-20.0, 20.0), Interval(-1.e-8, 1.e-8),
      Interval(0.25, 5.75)}) {
      CHECK_NONLINEAR_OPERATION(
        "sin", input, input,
        [](const AffineT& x) { return sin(x); },
        [](const Interval& x) { return sin(x); },
        [](double x) { return std::sin(x); }, 8192
      );
      CHECK_NONLINEAR_OPERATION(
        "cos", input, input,
        [](const AffineT& x) { return cos(x); },
        [](const Interval& x) { return cos(x); },
        [](double x) { return std::cos(x); }, 8192
      );
    }
    for (const Interval& input : {
      Interval(-1.4, 1.4), Interval(1.7, 4.5),
      Interval(-1.e-8, 1.e-8)}) {
      CHECK_NONLINEAR_OPERATION(
        "tan", input, input,
        [](const AffineT& x) { return tan(x); },
        [](const Interval& x) { return tan(x); },
        [](double x) { return std::tan(x); }, 8192
      );
    }

    // Inverse trigonometric functions, including largely invalid input
    // intervals and narrow valid parts near the singular endpoints.
    for (const Interval& input : {
      Interval(-1.0, 1.0), Interval(-1.e12, -1.0 + 1.e-12),
      Interval(1.0 - 1.e-12, 1.e12),
      Interval(-1.0 - 1.e-12, 1.0 + 1.e-12)}) {
      const Interval domain = input & Interval(-1.0, 1.0);
      CHECK_NONLINEAR_OPERATION(
        "asin", input, domain,
        [](const AffineT& x) { return asin(x); },
        [](const Interval& x) { return asin(x); },
        [](double x) { return std::asin(x); }, 8192
      );
      CHECK_NONLINEAR_OPERATION(
        "acos", input, domain,
        [](const AffineT& x) { return acos(x); },
        [](const Interval& x) { return acos(x); },
        [](double x) { return std::acos(x); }, 8192
      );
    }
    for (const Interval& input : {
      Interval(-1.e12, 1.e12), Interval(-1.e-12, 1.e-12),
      Interval(-10.0, 3.0)}) {
      CHECK_NONLINEAR_OPERATION(
        "atan", input, input,
        [](const AffineT& x) { return atan(x); },
        [](const Interval& x) { return atan(x); },
        [](double x) { return std::atan(x); }, 8192
      );
    }

    // Hyperbolic functions. Inputs stay below scalar overflow for sinh/cosh.
    for (const Interval& input : {
      Interval(-20.0, 20.0), Interval(-1.e-8, 1.e-8),
      Interval(-5.0, 2.0)}) {
      CHECK_NONLINEAR_OPERATION(
        "sinh", input, input,
        [](const AffineT& x) { return sinh(x); },
        [](const Interval& x) { return sinh(x); },
        [](double x) { return std::sinh(x); }
      );
      CHECK_NONLINEAR_OPERATION(
        "cosh", input, input,
        [](const AffineT& x) { return cosh(x); },
        [](const Interval& x) { return cosh(x); },
        [](double x) { return std::cosh(x); }
      );
      CHECK_NONLINEAR_OPERATION(
        "tanh", input, input,
        [](const AffineT& x) { return tanh(x); },
        [](const Interval& x) { return tanh(x); },
        [](double x) { return std::tanh(x); }
      );
      CHECK_NONLINEAR_OPERATION(
        "asinh", input, input,
        [](const AffineT& x) { return asinh(x); },
        [](const Interval& x) { return asinh(x); },
        [](double x) { return std::asinh(x); }
      );
    }

    // acosh and atanh with their restricted real domains.
    for (const Interval& input : {
      Interval(1.0, 10.0), Interval(-1.e12, 1.0 + 1.e-12),
      Interval(1.0 - 1.e-12, 1.0 + 1.e-8),
      Interval(-1.e12, 2.0)}) {
      const Interval domain = input & Interval(1.0, oo);
      CHECK_NONLINEAR_OPERATION(
        "acosh", input, domain,
        [](const AffineT& x) { return acosh(x); },
        [](const Interval& x) { return acosh(x); },
        [](double x) { return std::acosh(x); }, 8192
      );
    }
    for (const Interval& input : {
      Interval(-0.99, 0.99), Interval(-1.e-12, 1.e-12),
      Interval(-0.999999, 0.5), Interval(-0.5, 0.999999)}) {
      CHECK_NONLINEAR_OPERATION(
        "atanh", input, input,
        [](const AffineT& x) { return atanh(x); },
        [](const Interval& x) { return atanh(x); },
        [](double x) { return std::atanh(x); }, 8192
      );
    }
  }

  AffineT::change_mode(AffineT::AF_Lin_Chebyshev);
}


TEST_CASE("All nonlinear functions accept correlated affine expressions")
{
  using Mode = AffineT::Affine_Mode;
  const std::vector<Mode> modes = {
    AffineT::AF_Lin_Chebyshev,
    AffineT::AF_Lin_MinRange
  };

  for (const Mode mode : modes) {
    AffineT::change_mode(mode);
    CAPTURE(static_cast<int>(mode));

    AffineTVarVector variables(IntervalVector({{-1.0, 1.0}}));
    const AffineT& t = variables[0];

    const AffineT general = 3.0*t - 0.5;       // [-3.5,2.5]
    const AffineT positive = 1.5*t + 2.0;      // [0.5,3.5]
    const AffineT nonzero = 0.25*t + 1.0;      // [0.75,1.25]
    const AffineT unit = 0.75*t;               // [-0.75,0.75]
    const AffineT sqrt_arg = 5.0*t - 4.0;      // [-9,1]
    const AffineT acosh_arg = 5.0*t - 1.0;     // [-6,4]

    CHECK_affine_inclu<AA>(inv(nonzero), 1.0/nonzero.itv());
    CHECK_affine_inclu<AA>(sqr(general), sqr(general.itv()));
    CHECK_affine_inclu<AA>(abs(general), abs(general.itv()));
    CHECK_affine_inclu<AA>(sqrt(sqrt_arg), sqrt(sqrt_arg.itv()));
    CHECK_affine_inclu<AA>(exp(general), exp(general.itv()));
    CHECK_affine_inclu<AA>(log(positive), log(positive.itv()));
    CHECK_affine_inclu<AA>(pow(general, 3), pow(general.itv(), 3));
    CHECK_affine_inclu<AA>(pow(positive, 1.5), pow(positive.itv(), 1.5));
    CHECK_affine_inclu<AA>(root(general, 3), root(general.itv(), 3));

    CHECK_affine_inclu<AA>(sin(general), sin(general.itv()));
    CHECK_affine_inclu<AA>(cos(general), cos(general.itv()));
    CHECK_affine_inclu<AA>(tan(unit), tan(unit.itv()));
    CHECK_affine_inclu<AA>(asin(general), asin(general.itv()));
    CHECK_affine_inclu<AA>(acos(general), acos(general.itv()));
    CHECK_affine_inclu<AA>(atan(general), atan(general.itv()));

    CHECK_affine_inclu<AA>(sinh(general), sinh(general.itv()));
    CHECK_affine_inclu<AA>(cosh(general), cosh(general.itv()));
    CHECK_affine_inclu<AA>(tanh(general), tanh(general.itv()));
    CHECK_affine_inclu<AA>(asinh(general), asinh(general.itv()));
    CHECK_affine_inclu<AA>(acosh(acosh_arg), acosh(acosh_arg.itv()));
    CHECK_affine_inclu<AA>(atanh(unit), atanh(unit.itv()));

    const Interval exponent(0.5, 2.0);
    CHECK_affine_inclu<AA>(pow(positive, exponent),
                           pow(positive.itv(), exponent));

    AffineTVarVector exponent_variable(
      IntervalVector({exponent})
    );
    CHECK_affine_inclu<AA>(pow(positive, exponent_variable[0]),
                           pow(positive.itv(), exponent));
  }

  AffineT::change_mode(AffineT::AF_Lin_Chebyshev);
}



TEST_CASE("Non-smooth nonlinear affine functions enclose reference values")
{
  using Mode = AffineT::Affine_Mode;
  const std::vector<Mode> modes = {
    AffineT::AF_Lin_Chebyshev,
    AffineT::AF_Lin_MinRange
  };

  for (const Mode mode : modes) {
    AffineT::change_mode(mode);
    CAPTURE(static_cast<int>(mode));

    for (const Interval& input : {
      Interval(-4.75, 3.25), Interval(-1.0, 1.0),
      Interval(1.25, 1.75), Interval(-2.0, -2.0)}) {
      CHECK_NONLINEAR_OPERATION(
        "floor", input, input,
        [](const AffineT& x) { return floor(x); },
        [](const Interval& x) { return floor(x); },
        [](double x) { return std::floor(x); }
      );
      CHECK_NONLINEAR_OPERATION(
        "ceil", input, input,
        [](const AffineT& x) { return ceil(x); },
        [](const Interval& x) { return ceil(x); },
        [](double x) { return std::ceil(x); }
      );

      CHECK_NONLINEAR_INTERVAL_IMAGE(
        "integer", input,
        [](const AffineT& x) { return integer(x); },
        [](const Interval& x) { return integer(x); }
      );
    }

    for (const Interval& input : {
      Interval(-4.0, 3.0), Interval(-3.0, -0.25),
      Interval(0.25, 5.0), Interval(0.0, 0.0)}) {
      CHECK_NONLINEAR_INTERVAL_IMAGE(
        "sign", input,
        [](const AffineT& x) { return sign(x); },
        [](const Interval& x) { return sign(x); }
      );
    }
  }

  AffineT::change_mode(AffineT::AF_Lin_Chebyshev);
}


TEST_CASE("atan2 encloses the reference two-argument arctangent")
{
  const std::vector<std::pair<Interval,Interval>> cases = {
    {Interval(1.0, 2.0), Interval(1.0, 2.0)},     // first quadrant
    {Interval(1.0, 2.0), Interval(-2.0, -1.0)},   // second quadrant
    {Interval(-2.0, -1.0), Interval(-2.0, -1.0)}, // third quadrant
    {Interval(-2.0, -1.0), Interval(1.0, 2.0)},   // fourth quadrant
    {Interval(-0.1, 0.1), Interval(1.0, 2.0)}     // straddles y=0, x>0
  };

  for (const auto& input : cases) {
    const Interval& y_itv = input.first;
    const Interval& x_itv = input.second;
    CAPTURE(y_itv, x_itv);

    AffineTVarVector variables(IntervalVector({y_itv, x_itv}));
    const AffineT result = atan2(variables[0], variables[1]);
    const Interval reference = atan2(y_itv, x_itv);

    CHECK_affine_inclu<AA>(result, reference);

    for (int k = 0; k <= 50; ++k) {
      const double y = y_itv.lb() + (y_itv.ub()-y_itv.lb())*k/50.0;
      const double x = x_itv.lb() + (x_itv.ub()-x_itv.lb())*k/50.0;
      CAPTURE(y, x, result.itv());
      CHECK(result.contains(std::atan2(y, x)));
    }
  }
}


TEST_CASE("chi selects a branch according to the sign of the condition")
{
  AffineTVarVector bx(2);
  bx[0] = Interval(10.0, 11.0);  // branch b
  bx[1] = Interval(20.0, 21.0);  // branch c

  const AffineT negative(Interval(-2.0, -1.0));
  const AffineT positive(Interval(1.0, 2.0));

  // Condition strictly negative: selects b, for every Interval/AffineMain
  // combination of the operands that the API exposes.
  CHECK(chi(negative, bx[0], bx[1]).itv() == bx[0].itv());
  CHECK(chi(negative.itv(), bx[0], bx[1]).itv() == bx[0].itv());
  CHECK(chi(negative.itv(), bx[0].itv(), bx[1]).itv() == bx[0].itv());
  CHECK(chi(negative.itv(), bx[0], bx[1].itv()).itv() == bx[0].itv());
  CHECK(chi(negative, Interval(30.0, 31.0), bx[1]).itv() == Interval(30.0, 31.0));

  // Condition strictly positive: selects c.
  CHECK(chi(positive, bx[0], bx[1]).itv() == bx[1].itv());
  CHECK(chi(positive.itv(), bx[0], bx[1]).itv() == bx[1].itv());
  CHECK(chi(positive, bx[0], Interval(30.0, 31.0)).itv() == Interval(30.0, 31.0));

  // Condition straddling zero: both branches remain possible, so the
  // result must enclose their union.
  const AffineT straddling = chi(Interval(-1.0, 1.0), bx[0], bx[1]);
  CHECK(straddling.itv().is_superset(bx[0].itv() | bx[1].itv()));

  // An empty condition yields an empty result.
  CHECK(chi(Interval::empty(), bx[0], bx[1]).is_empty());
}


// NOTE: distance(x1,x2) is declared and documented on AffineMain<T>
// (codac2_AffineMain.h) for 3 overloads (AffineMain-AffineMain,
// Interval-AffineMain, AffineMain-Interval), and each one forwards to a
// free function codac2::distance(const Interval&, const Interval&) that
// does not exist anywhere in the codebase (verified by a full search of
// src/core). Instantiating any of the 3 overloads fails to compile, so
// distance() is left untested here: this is a library gap, not a gap in
// this test. See the reply that added this test suite for the full report.


TEST_CASE("inflate widens the enclosure by the given radius on each side")
{
  AffineTVarVector ax(1);
  ax[0] = Interval(1.0, 2.0);
  const Interval before = ax[0].itv();

  AffineT inflated = ax[0];
  inflated.inflate(0.5);

  CHECK_affine_inclu<AA>(inflated, before + Interval(-0.5, 0.5));
  CHECK(inflated.noise_count() == ax[0].noise_count());

  // Inflating by zero must not change the enclosure.
  AffineT unchanged = ax[0];
  unchanged.inflate(0.0);
  CHECK(unchanged.itv() == before);
}


TEST_CASE("mig, mag, smag, smig and volume delegate to the interval enclosure")
{
  AffineTVarVector ax(3);
  ax[0] = Interval(-3.0, 5.0);   // straddles zero
  ax[1] = Interval(2.0, 4.0);    // strictly positive
  ax[2] = Interval(-5.0, -1.0);  // strictly negative

  for (Index i = 0; i < ax.size(); ++i) {
    CAPTURE(i, ax[i].itv());
    CHECK(ax[i].mig() == ax[i].itv().mig());
    CHECK(ax[i].mag() == ax[i].itv().mag());
    CHECK(ax[i].smag() == ax[i].itv().smag());
    CHECK(ax[i].smig() == ax[i].itv().smig());
    CHECK(ax[i].volume() == ax[i].itv().volume());
    CHECK(ax[i].volume() == ax[i].diam());
  }

  // The documented formula: mig = 0 when the interval straddles zero,
  // +lb when strictly positive, -ub when strictly negative.
  CHECK(ax[0].mig() == 0.0);
  CHECK(ax[1].mig() == 2.0);
  CHECK(ax[2].mig() == 1.0);
  CHECK(ax[0].mag() == 5.0);
}

