/* ============================================================================
 * I B E X - Test of Affine operations
 * ============================================================================
 * Copyright   : ENSTA Bretagne (FRANCE)
 * License     : This program can be distributed under the terms of the GNU LGPL.
 *               See the file LICENSE
 *
 * Author(s)   : Jordan NININ
 * Created     : Juin 20, 2021
 * ---------------------------------------------------------------------------- */

#include <catch2/catch_test_macros.hpp>
#include "codac2_Affine.h"
#include "codac2_Approx.h"
#include <limits>
#include "codac2_Matrix.h"
#include <iostream>
using namespace codac2;

const double MAX_DOUBLE = std::numeric_limits<double>::max();
const double ERROR = std::numeric_limits<double>::epsilon()*10; //__builtin_powi(2.0, -50);

#define piL Interval::pi().lb()
#define piU Interval::pi().ub()

template<class T>
void CHECK_change_mode_MinRange()
{
  AffineMain<T>::change_mode(AffineMain<T>::AF_Lin_MinRange);
  CHECK(AffineMain<T>::AF_Lin_MinRange==AffineMain<T>::get_mode());

  AffineMain<T>::change_mode(AffineMain<T>::AF_Lin_Chebyshev);
  CHECK(AffineMain<T>::AF_Lin_Chebyshev==AffineMain<T>::get_mode());

  AffineMain<T>::change_mode(AffineMain<T>::AF_Lin_Default);
  CHECK(AffineMain<T>::AF_Lin_Chebyshev==AffineMain<T>::get_mode());
}

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
	if (x[0].is_actif() && std::fabs(z)<1.e100) {
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
	if (x.is_actif() && z.is_actif()) {
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
	if (x.is_actif() && !z.is_unbounded() && !z.is_empty()) {
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
void CHECK_sinh(const Interval& x)
{
  double xl=x.lb();
  double xu=x.ub();
  double yl = (xl==-oo)? -oo : 0.5*(std::exp(xl)-std::exp(-xl));
  double yu = (xu==oo)? oo : 0.5*(std::exp(xu)-std::exp(-xu));

  AffineVarMainVector<T> ax(1);
  ax[0] = x;
  CHECK_affine_inclu<T>(sinh(ax[0]), Interval(yl,yu));
  CHECK_affine_inclu<T>(sinh(-ax[0]), Interval(-yu,-yl));
}

template<class T>
void CHECK_cosh(const Interval& x)
{
  Interval y;
  if (x.ub()==oo) {
    if (x.lb()<=0) y=Interval(1,oo);
    else y=Interval((std::cosh(x.lb())),oo);
  }
  else if (x.lb()==-oo) {
    if (x.ub()>=0) y=Interval(1,oo);
    else y=Interval((std::cosh(x.ub())),oo);
  }
  else if (x.lb()>=0)
    y=Interval((std::cosh(x.lb())),(std::cosh(x.ub())));
  else if (x.ub()<=0)
    y=Interval((std::cosh(x.ub())),(std::cosh(x.lb())));
  else
    y=((std::fabs(x.lb())> std::fabs(x.ub())) ? Interval(1,(std::cosh(x.lb()))) :Interval(1,(std::cosh(x.ub()))));

  AffineVarMainVector<T> ax(1);
  ax[0] = x;
  CHECK_affine_inclu<T>(cosh(ax[0]), y);
  CHECK_affine_inclu<T>(cosh(-ax[0]), y);
}





TEST_CASE("Approx<AffineMain<T>>")
{
  typedef AF_Default AA;

  // empty
  {
    AffineVarMainVector<AA> ax(1);
    ax[0] = Interval::empty();
    CHECK(ax[0] == Approx<AffineMain<AA>>(Interval::empty()));
  }

  // unbounded
  {
    AffineVarMainVector<AA> ax(1);
    CHECK(ax[0] == Approx<AffineMain<AA>>(Interval()));
  }
  {

	  AffineVarMainVector<AA> ax(IntervalVector({{0,1},{0,1},{0,1}}));
	  CHECK(ax[0] == Approx<AffineMain<AA>>(Interval(0,1)));
	  CHECK(ax[1] == Approx<AffineMain<AA>>(Interval(0,1)));
	  CHECK(ax[2] == Approx<AffineMain<AA>>(Interval(0,1)));
  }

  // simple linear cases: x in [0,1], 1 noise var
  {
    AffineVarMainVector<AA> ax(IntervalVector({{0,1}}));
    CHECK(ax[0] == Approx<AffineMain<AA>>(Interval(0,1)));
    CHECK((ax[0]+1.0) == Approx<AffineMain<AA>>(Interval(1,2)));
    CHECK((-ax[0]) == Approx<AffineMain<AA>>(Interval(-1,0)));
    CHECK((2.0*ax[0]) == Approx<AffineMain<AA>>(Interval(0,2)));
  }

  // two affine variables summed together (2 noise variables)
  {
	  AffineVarMainVector<AA> ax(IntervalVector({{0,1},{0,1}}));
	  AffineMain<AA> y = ax[0] + ax[1];
	  CHECK(y == Approx<AffineMain<AA>>(Interval(0,2)));
  }

  // must not match a wrong expected interval
  {
    AffineVarMainVector<AA> bx(2);
    bx[0] =Interval(0,1);
    CHECK_FALSE(bx[0] == Approx<AffineMain<AA>>(Interval(5,6)));
    bx[0] =Interval(5,6);
    CHECK(bx[0] == Approx<AffineMain<AA>>(Interval(5,6)));
  }
}





TEST_CASE("AffineForm operations")
{
  typedef AF_Default AA;

  CHECK_change_mode_MinRange<AA>();

  /* test: operator-() */
  {
    AffineVarMainVector<AA> ax(1);
    ax[0] =  Interval(0,1);
    CHECK(-(ax[0]) == Approx<AffineMain<AA>>(Interval(-1,0), ERROR));
  }
  {
    AffineVarMainVector<AA> ax(1);
    CHECK(-(ax[0]) == Approx<AffineMain<AA>>(Interval(), ERROR));
  }
  {
    AffineVarMainVector<AA> ax(1);
    ax[0] = Interval(-oo,0);
    CHECK(-(ax[0]) == Approx<AffineMain<AA>>(Interval(0,oo), ERROR));
  }
  {
    AffineVarMainVector<AA> ax(4);
    ax[0] = Interval(-oo,1);
    CHECK(-(ax[0]) == Approx<AffineMain<AA>>(Interval(-1,oo), ERROR));
    CHECK(ax[2] == Approx<AffineMain<AA>>(Interval(), ERROR));
  }

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
    AffineVarMainVector<AA> ax(1);
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
	  AffineVarMainVector<AA> ax(1);
	  CHECK_affine_inclu<AA>(tan(ax[0]), Interval());
  }
  {
	  AffineVarMainVector<AA> ax(1);
	  ax[0] = (-Interval::pi()/4.0 | Interval::pi()/4.0);
	  CHECK_affine_inclu<AA>(tan(ax[0]), Interval(-1,1));
  }
  { // tan(pi/4,pi/2)=[1,+oo)
	  AffineVarMainVector<AA> x(1); // upper bound of x is close to pi/2
	  x[0] = Interval(piL/4.0,(1-1e-10)*piL/2.0);
	  AffineMain<AA> y=tan(x[0]);
	  CHECK(y.itv().lb()<=1.0);
	  CHECK(y.itv().ub()>1.e8); // upper bound of tan(x) is close to +oo
  }
  { // tan(-pi/2,pi/4)=(-oo,1]
	  AffineVarMainVector<AA> ax(1);
	  ax[0] = Interval(-(1-1e-10)*piL/2.0,piL/4.0);
	  Interval y= (tan( ax[0] )).itv();
	  CHECK(y.lb()<= -1.e8); // lower bound is close to -oo
	  CHECK(y.ub()>= 1.0);
  }
  {
	  AffineVarMainVector<AA> ax(1);
	  ax[0] = Interval::pi()/2.0;
	  CHECK_affine_inclu<AA>(tan(ax[0]),Interval());
  }
  {
	  AffineVarMainVector<AA> ax(1);
	  ax[0] = Interval::pi();
	  AffineMain<AA> y=tan(-ax[0]);
	  CHECK(y==Approx<AffineMain<AA>>(Interval(0), ERROR));
	  CHECK(y.contains(0.));
	  CHECK(y.diam() < 1e-8);
  }
  {
	  AffineVarMainVector<AA> ax(1);
	  ax[0] = Interval(3*piL/4.0 , 5*piU/4.0);
	  AffineMain<AA>  y = tan(ax[0]);
	  CHECK_affine_inclu<AA>(y, Interval(-1,1));
	  CHECK(y.lb() == Approx(-1.0,1.e-5));
	  CHECK(y.ub() == Approx(1.0,1.e-5));

	  ax[0] = Interval(-oo,oo);
	  y= tan (ax[0]);
	  CHECK(y==Approx<AffineMain<AA>>(Interval(), ERROR));

	  ax[0] = ((-Interval::pi()/4.0)|(Interval::pi()/4.0));
	  y= tan(ax[0]);
	  CHECK_affine_inclu<AA>(y, Interval(-1,1));
	  CHECK(y.lb() == Approx(-1.0,1.e-5));
	  CHECK(y.ub() == Approx(1.0,1.e-5));
  }
  {
	  // tan(pi/4,pi/2)=[1,+oo)
	  AffineVarMainVector<AA> ax(1);
	  ax[0] = Interval(piL/4.0,(1-1e-10)*piL/2.0);
	  AffineMain<AA>  y = tan(ax[0]);
	  CHECK(y.lb() <=1.0);
	  CHECK(y.ub() >= 1.e8); // upper bound of tan(x) is close to +oo
  }
  {
	  // tan(-pi/2,pi/4)=(-oo,1]
	  AffineVarMainVector<AA> ax(1);
	  ax[0] = Interval(-(1-1e-10)*piL/2.0,piL/4.0);
	  AffineMain<AA>  y = tan(ax[0]);
	  CHECK_affine_inclu<AA>(y, Interval(-1.e8,1));
  }

  {
    AffineVarMainVector<AA> ax(1);
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

  CHECK_cosh<AA>(Interval());
  CHECK_cosh<AA>(Interval(0,oo));
  CHECK_cosh<AA>(Interval(0,1));
  CHECK_cosh<AA>(Interval(1,oo));
  CHECK_cosh<AA>(Interval(1,1));
  CHECK_cosh<AA>(Interval(2,3));
  CHECK_cosh<AA>(Interval(4,5));

  { // issue 248
    AffineVarMainVector<AA> aff(1);
    aff[0] = Interval(-1.57079632679489678, 1.1780972450961728626);
    CHECK_FALSE((tan(aff[0]).is_empty()));
  }
}
