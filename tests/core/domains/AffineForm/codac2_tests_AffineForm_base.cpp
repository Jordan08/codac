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

using namespace codac2;

const double MAX_DOUBLE = std::numeric_limits<double>::max();

const double ERROR = __builtin_powi(2.0, -55);

template<class T>
void  CHECK_change_mode_MinRange() {
	AffineMain<T>::change_mode(AffineMain<T>::AF_Lin_MinRange);
	CHECK(AffineMain<T>::AF_Lin_MinRange==AffineMain<T>::get_mode());

	AffineMain<T>::change_mode(AffineMain<T>::AF_Lin_Chebyshev);
	CHECK(AffineMain<T>::AF_Lin_Chebyshev==AffineMain<T>::get_mode());

	AffineMain<T>::change_mode(AffineMain<T>::AF_Lin_Default);
	CHECK(AffineMain<T>::AF_Lin_Chebyshev==AffineMain<T>::get_mode());
}


template<class T>
void  CHECK_affine_eq(const AffineMain<T>& y_actual, const Interval& y_expected, double err = ERROR) {
	//std::cout << "CHECK:    " << y_expected << " (expected)        " << y_actual << " (actual)"<< std::endl;
	if (y_expected.is_empty()) { CHECK(y_actual.is_empty()); return; }
	if (y_expected.is_unbounded()) {
		CHECK(y_actual.is_unbounded());
		CHECK(y_actual.itv()==y_expected);
	} else {
		CHECK_FALSE(y_actual.is_empty());
		CHECK(Approx(y_actual.itv().lb(),err)==y_expected.lb());
		CHECK(Approx(y_actual.itv().ub(),err)==y_expected.ub());
		CHECK(Approx(y_actual.itv().mid(),err)==y_expected.mid());
		CHECK(Approx(std::fabs(y_actual.val(0)),err)==y_expected.rad());
	}
}


template<class T>
void  CHECK_affine_eq2(const AffineMain<T>& y_actual, const Interval& y_expected, double err = ERROR) {
	//std::cout << "CHECK:    " << y_expected << " (expected)        " << y_actual << " (actual)"<< std::endl;
	if (y_expected.is_empty()) { CHECK(y_actual.is_empty()); return; }
	if (y_expected.is_unbounded()) {
		CHECK(y_actual.is_unbounded());
		CHECK(y_actual.itv()==y_expected);
	} else {
		CHECK_FALSE(y_actual.is_empty());
		CHECK(Approx(y_actual.itv().lb(),err)==y_expected.lb());
		CHECK(Approx(y_actual.itv().ub(),err)==y_expected.ub());
		CHECK(Approx(y_actual.itv().mid(),err)==y_expected.mid());
		CHECK(Approx(std::fabs(y_actual.val(0))+fabs(y_actual.val(1)),err)==y_expected.rad());
	}
}

template<class T>
void  CHECK_affine_inclu(const AffineMain<T>& y_actual, const Interval& y_expected) {
	//std::cout << "CHECK:    " << y_expected << " (expected)        " << y_actual << " (actual)"<< std::endl;
	if (y_expected.is_empty()) { CHECK(y_actual.is_empty()); return; }

	CHECK_FALSE(y_actual.is_empty());
	CHECK(y_expected.lb()>=y_actual.itv().lb());
	CHECK(y_expected.ub()<=y_actual.itv().ub());
}

template<class T>
void  CHECK_add(const Interval& xi, double z, const Interval& y_expected) {
	AffineVarMainVector<T> x(1,xi);
	AffineMain<T>  y_actual=x[0]+z;
	CHECK_affine_eq<T>(y_actual, y_expected);

	// test the symmetrical case
	y_actual = z+x[0];
	CHECK_affine_eq<T>(y_actual, y_expected);

	// test the +=operator
	y_actual = x[0];
	y_actual += z;
	CHECK_affine_eq<T>(y_actual, y_expected);

	// test the +=operator in the other direction
	y_actual = z;
	y_actual += x[0];
	CHECK_affine_eq<T>(y_actual, y_expected);

	// test subtraction
	y_actual=(-x[0])-z;
	CHECK_affine_eq<T>(y_actual, -y_expected);

	// test the symmetrical case
	y_actual = (-z)-x[0];
	CHECK_affine_eq<T>(y_actual, -y_expected);

	// test the -=operator
	y_actual = -x[0];
	y_actual -= z;
	CHECK_affine_eq<T>(y_actual, -y_expected);

	// test the -=operator in the other direction
	y_actual = -z;
	y_actual -= x[0];
	CHECK_affine_eq<T>(y_actual, -y_expected);

	// test the linear simplification
	if (x[0].is_actif() && std::fabs(z)<1.e100) {
		y_actual = z ;
		y_actual += x[0] ;
		y_actual -= x[0] ;
		y_actual += x[0] ;
		CHECK_affine_eq<T>(y_actual, y_expected);

		y_actual = x[0] ;
		y_actual += z ;
		y_actual -= z ;
		y_actual += z ;
		CHECK_affine_eq<T>(y_actual, y_expected);
	}
}

template<class T>
void  CHECK_add(const AffineMain<T>& x, const AffineMain<T>& z, const Interval& y_expected) {
	AffineMain<T>  y_actual=x+z;
	CHECK_affine_eq2<T>(y_actual, y_expected);

	// test the symmetrical case
	y_actual = z+x;
	CHECK_affine_eq2<T>(y_actual, y_expected);

	// test the +=operator
	y_actual = x;
	y_actual += z;
	CHECK_affine_eq2<T>(y_actual, y_expected);

	// test the +=operator in the other direction
	y_actual = z;
	y_actual += x;
	CHECK_affine_eq2<T>(y_actual, y_expected);

	// test subtraction
	y_actual=(-x)-z;
	CHECK_affine_eq2<T>(y_actual, -y_expected);

	// test the symmetrical case
	y_actual = (-z)-x;
	CHECK_affine_eq2<T>(y_actual, -y_expected);

	// test the -=operator
	y_actual = -x;
	y_actual -= z;
	CHECK_affine_eq2<T>(y_actual, -y_expected);

	// test the -=operator in the other direction
	y_actual = -z;
	y_actual -= x;
	CHECK_affine_eq2<T>(y_actual, -y_expected);

	// test the linear simplification
	if (x.is_actif() && z.is_actif()) {
		y_actual = z ;
		y_actual += x ;
		y_actual -= x ;
		y_actual += x ;
		CHECK_affine_eq2<T>(y_actual, y_expected);

		y_actual = x ;
		y_actual += z ;
		y_actual -= z ;
		y_actual += z ;
		CHECK_affine_eq2<T>(y_actual, y_expected);
	}
}

template<class T>
void  CHECK_add(const AffineMain<T>& x, const Interval& z, const Interval& y_expected) {
	AffineMain<T>  y_actual=x+z;
	CHECK_affine_eq<T>(y_actual, y_expected);

	// test the symmetrical case
	y_actual = z+x;
	CHECK_affine_eq<T>(y_actual, y_expected);

	// test the +=operator
	y_actual = x;
	y_actual += z;
	CHECK_affine_eq<T>(y_actual, y_expected);

	// test the +=operator in the other direction
	y_actual = z;
	y_actual += x;
	CHECK_affine_eq<T>(y_actual, y_expected);

	// test subtraction
	y_actual=(-x)-z;
	CHECK_affine_eq<T>(y_actual, -y_expected);

	// test the symmetrical case
	y_actual = (-z)-x;
	CHECK_affine_eq<T>(y_actual, -y_expected);

	// test the -=operator
	y_actual = -x;
	y_actual -= z;
	CHECK_affine_eq<T>(y_actual, -y_expected);

	// test the -=operator in the other direction
	y_actual = -z;
	y_actual -= x;
	CHECK_affine_eq<T>(y_actual, -y_expected);

	// test the linear simplification
	if (x.is_actif() && !z.is_unbounded() && !z.is_empty()) {
		y_actual = z ;
		y_actual += x ;
		y_actual -= x ;
		y_actual += x ;
		CHECK_affine_eq<T>(y_actual, y_expected);

		y_actual = x ;
		y_actual += z ;
		y_actual -= z ;
		y_actual += z ;
		CHECK_affine_eq<T>(y_actual, y_expected);
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
void  CHECK_add_scal(const Interval& x, double z, const Interval& y_expected) {
	CHECK_add<T>(x,z, y_expected);
	CHECK_add<T>((-x),-z, -y_expected);
}

template<class T>
void  CHECK_add_scal(const Interval& x, const Interval& z, const Interval& y_expected) {
	AffineVarMainVector<T> xa(1,x);
	CHECK_add<T>(xa[0],z,y_expected);
	CHECK_add<T>((-xa[0]),-z, -y_expected);
}


template<class T>
void  CHECK_mul(const AffineMain<T>& x, const AffineMain<T>& z, const Interval& y_expected) {
	AffineMain<T> y_actual=x*z;
	//cout << "CHECK:    " << x << " * " << z << ", " << y_expected << endl;
	//cout << "      out:" << y_actual << endl;
	CHECK_affine_inclu<T>(y_actual, y_expected);

	// test the symmetrical case
	y_actual = z*x;
	CHECK_affine_inclu<T>(y_actual, y_expected);

	// test the *=operator
	y_actual = x;
	y_actual *= z;
	CHECK_affine_inclu<T>(y_actual, y_expected);
}

template<class T>
void  CHECK_mul(const Interval& x, const Interval& z, const Interval& y_expected) {
	AffineVarMainVector<T> ax(2);
	ax[0]=x;
	ax[1]=z;
	CHECK_mul<T>(ax[0],ax[1],y_expected);
}

template<class T>
void  CHECK_mul_scal(const Interval& x, double z, const Interval& y_expected) {
	AffineVarMainVector<T> xa(1,x);
	CHECK_affine_eq<T>(xa[0]*z, y_expected);
	CHECK_affine_eq<T>(z*xa[0], y_expected);
}

template<class T>
void  CHECK_div(const AffineMain<T>& x, const AffineMain<T>& z, const Interval& y_expected) {
	AffineMain<T> y_actual=x/z;
	//cout << "CHECK:    " << x << " / " << z << ", " << y_expected << " (expected)        " << y_actual << " (actual)"<< endl;
	CHECK_affine_inclu<T>(y_actual, y_expected);

	// test the /=operator
	y_actual = x;
	y_actual /= z;
	CHECK_affine_inclu<T>(y_actual, y_expected);
}

template<class T>
void  CHECK_div(const Interval& x, const Interval& z, const Interval& y_expected) {
	AffineVarMainVector<T> ax(2);
	ax[0]=x;
	ax[1]=z;
	CHECK_div<T>(ax[0],ax[1],y_expected);
}

template<class T>
void  CHECK_div_scal(const Interval& x, double z, const Interval& y_expected) {
	AffineVarMainVector<T> xa(1,x);
	AffineMain<T> y_actual=xa[0]/z;
	CHECK_affine_inclu<T>(y_actual, y_expected);

	// test the /=operator
	y_actual = xa[0];
	y_actual /= z;
	CHECK_affine_inclu<T>(y_actual, y_expected);
}


template<class T>
void  minus01() {
	AffineVarMainVector<T> ax(1, Interval(0,1));
	CHECK_affine_eq<T>(-(ax[0]), Interval(-1,0)); }
template<class T>
void  minus02() {
	AffineVarMainVector<T> ax(1);
	CHECK_affine_eq<T>(-(ax[0]), Interval()); }
template<class T>
void  minus03() {
	AffineVarMainVector<T> ax(1, Interval(-oo,0));
	CHECK_affine_eq<T>(-(ax[0]), Interval(0,oo)); }
template<class T>
void  minus04() {
	AffineVarMainVector<T> ax(1, Interval(-oo,1));
	CHECK_affine_eq<T>(-(ax[0]), Interval(-1,oo)); }
template<class T>
void  add01() { CHECK_add<T>(Interval::empty(),    Interval(0,1),      Interval::empty()); }
template<class T>
void  add02() { CHECK_add<T>(Interval(0,1),		Interval::empty(), Interval::empty()); }
template<class T>
void  add03() {CHECK_add<T>(Interval(-oo,1),		Interval(0,1),      Interval(-oo, 2)); }
template<class T>
void  add04() {CHECK_add<T>(Interval(1,oo), Interval(0,1),      Interval(1,oo)); }
template<class T>
void  add05() { CHECK_add<T>(Interval(),      Interval(0,1),      Interval()); }
template<class T>
void  add06() { CHECK_add_scal<T>(Interval::empty(),      oo,       Interval::empty()); }
template<class T>
void  add07() { CHECK_add_scal<T>(Interval::empty(),      0,                  Interval::empty()); }
template<class T>
void  add08() { CHECK_add_scal<T>(Interval(0,1),            1,                  Interval(1,2)); }
template<class T>
void  add09() { CHECK_add_scal<T>(Interval(0,1),            -oo,       Interval::empty()); }
template<class T>
void  add10() { CHECK_add_scal<T>(Interval(0,1),            oo,       Interval::empty()); }
template<class T>
void  add11() { CHECK_add_scal<T>(Interval(-oo,1), 1,                  Interval(-oo,2)); }
template<class T>
void  add12() { CHECK_add<T>(Interval(MAX_DOUBLE,oo), 1, Interval(MAX_DOUBLE,oo)); }
template<class T>
void  add13() { CHECK_add<T>(Interval(MAX_DOUBLE,oo), -1, Interval(codac2::prev_float(MAX_DOUBLE),oo)); }
template<class T>
void  add14() { CHECK_add<T>(Interval(MAX_DOUBLE,oo), Interval(MAX_DOUBLE,oo), Interval(MAX_DOUBLE,oo)); }
template<class T>
void  add15() { CHECK_add<T>(Interval(MAX_DOUBLE,oo), -oo, Interval::empty()); }
template<class T>
void  add16() { CHECK_add<T>(Interval(MAX_DOUBLE,oo), oo, Interval::empty()); }
template<class T>
void  add17() { CHECK_add<T>(Interval(-oo,-MAX_DOUBLE),  1, Interval(-oo,codac2::next_float(-MAX_DOUBLE))); }
template<class T>
void  add18() { CHECK_add<T>(Interval(-oo,-MAX_DOUBLE), -1, Interval(-oo,-MAX_DOUBLE)); }
template<class T>
void  add19() { CHECK_add<T>(Interval(-oo,-MAX_DOUBLE),  Interval(-oo,-MAX_DOUBLE), Interval(-oo,-MAX_DOUBLE)); }

template<class T>
void  mul01() { CHECK_mul<T>(Interval::empty(),         Interval(0,1),               Interval::empty()); }
template<class T>
void  mul02() { CHECK_mul<T>(Interval::zero(),              Interval(),         Interval::zero()); }
template<class T>
void  mul03() { CHECK_mul<T>(Interval(-1,1),              Interval(-oo,0), 	     Interval()); }
template<class T>
void  mul04() { CHECK_mul<T>(Interval(-oo,-1),   Interval(-1,0),              Interval(0,oo)); }
template<class T>
void  mul05() { CHECK_mul<T>(Interval(-oo, 1),   Interval(-1,0),              Interval(-1,oo)); }
template<class T>
void  mul06() { CHECK_mul<T>(Interval(0, 1),              Interval(1,oo), 	 Interval(0,oo)); }
template<class T>
void  mul07() { CHECK_mul<T>(Interval(0, 1),              Interval(-1,oo),   Interval(-1,oo)); }
template<class T>
void  mul08() { CHECK_mul<T>(Interval(-oo,-1),   Interval(0,1),               Interval(-oo,0)); }
template<class T>
void  mul09() { CHECK_mul<T>(Interval(-oo, 1),   Interval(0,1),               Interval(-oo,1)); }
template<class T>
void  mul10() { CHECK_mul<T>(Interval(0, 1),              Interval(-oo,-1),   Interval(-oo,0)); }
template<class T>
void  mul11() { CHECK_mul<T>(Interval(0, 1),              Interval(-oo,1),    Interval(-oo,1)); }
template<class T>
void  mul12() { CHECK_mul<T>(Interval(1,oo),    Interval(0,1),               Interval(0,oo)); }
template<class T>
void  mul13() { CHECK_mul<T>(Interval(-1,oo),   Interval(0,1),               Interval(-1,oo)); }
template<class T>
void  mul14() { CHECK_mul<T>(Interval(1,2),               Interval(1,2), 	        	 Interval(1,4)); }
template<class T>
void  mul15() { CHECK_mul<T>(Interval(1,2),               Interval(-2,3), 	         Interval(-4,6)); }
template<class T>
void  mul16() { CHECK_mul<T>(Interval(-1,1),              Interval(-2,3), 	         Interval(-3,3)); }
template<class T>
void  mul17() { CHECK_mul_scal<T>(Interval(1,2),          -oo, 	        	 Interval::empty()); }
template<class T>
void  mul18() { CHECK_mul_scal<T>(Interval(1,2),          oo, 	        	 Interval::empty()); }
template<class T>
void  mul19() { CHECK_mul_scal<T>(Interval(1,2),          -1, 	        	         Interval(-2,-1)); }

template<class T>
void  div01() { CHECK_div<T>(Interval::empty(),         Interval(0,1),               Interval::empty()); }
template<class T>
void  div02() { CHECK_div<T>(Interval::zero(),              Interval::zero(),              Interval::empty()); }
template<class T>
void  div03() { CHECK_div<T>(Interval(1,2),               Interval::zero(),              Interval::empty()); }
template<class T>
void  div04() { CHECK_div<T>(Interval(),         Interval::zero(),              Interval::empty()); }
template<class T>
void  div05() { CHECK_div<T>(Interval::zero(),              Interval(0,1),               Interval::zero()); }
template<class T>
void  div06() { CHECK_div<T>(Interval::zero(),              Interval(),         Interval::zero()); }
template<class T>
void  div07() { CHECK_div<T>(Interval(6,12),              Interval(2,3),               Interval(2,6)); }
template<class T>
void  div08() { CHECK_div<T>(Interval(6,12),              Interval(-3,-2),             Interval(-6,-2)); }
template<class T>
void  div09() { CHECK_div<T>(Interval(6,12),              Interval(-2,3),              Interval()); }
template<class T>
void  div10() { CHECK_div<T>(Interval(-oo,-1),   Interval(-1,0),              Interval(1,oo)); }
template<class T>
void  div11() { CHECK_div<T>(Interval(-oo,-1),   Interval(0,1),               Interval(-oo,-1)); }
template<class T>
void  div12() { CHECK_div<T>(Interval(1,oo),    Interval(-1,0),              Interval(-oo,-1)); }
template<class T>
void  div13() { CHECK_div<T>(Interval(1,oo),    Interval(0,1),               Interval(1,oo)); }
template<class T>
void  div14() { CHECK_div<T>(Interval(-1,1),              Interval(-1,1), 	         Interval()); }
template<class T>
void  div15() { CHECK_div_scal<T>(Interval(1,2),          -oo, 	        	 Interval::empty()); }
template<class T>
void  div16() { CHECK_div_scal<T>(Interval(1,2),          oo, 	        	 Interval::empty()); }
template<class T>
void  div17() { CHECK_div_scal<T>(Interval(1,2),          -1, 	        	         Interval(-2,-1)); }





template<class T>
void  sqrt01() {
	AffineVarMainVector<T> ax(1);
	CHECK_affine_inclu<T>(sqrt(ax[0]), 	Interval(0,oo)); }
template<class T>
void  sqrt02() {
	AffineVarMainVector<T> ax(1,Interval(-oo,0));
	CHECK_affine_inclu<T>(sqrt(ax[0]), 	Interval::zero()); }
template<class T>
void  sqrt03() {
	AffineVarMainVector<T> ax(1,Interval(-9,4));
	CHECK_affine_inclu<T>(sqrt(ax[0]),     Interval(0,2)); }
template<class T>
void  sqrt04() {
	AffineVarMainVector<T> ax(1,Interval(4,9));
	CHECK_affine_inclu<T>(sqrt(ax[0]),     Interval(2,3)); }
template<class T>
void  sqrt05() {
	AffineVarMainVector<T> ax(1,Interval(-9,-4));
	CHECK_affine_inclu<T>(sqrt(ax[0]),     Interval::empty()); }
template<class T>
void  sqrt06() {
	AffineVarMainVector<T> ax(1,Interval(-9,oo));
	CHECK_affine_inclu<T>(sqrt(ax[0]),     Interval(0,oo)); }

#define piL Interval::pi().lb()
#define piU Interval::pi().ub()

template<class T>
void  CHECK_sinh(const Interval& x) {
	double xl=x.lb();
	double xu=x.ub();
	double yl = (xl==-oo)? -oo : 0.5*(std::exp(xl)-std::exp(-xl));
	double yu = (xu==oo)? oo : 0.5*(std::exp(xu)-std::exp(-xu));

	AffineVarMainVector<T> ax(1,x);
	CHECK_affine_inclu<T>(sinh(ax[0]), Interval(yl,yu));
	CHECK_affine_inclu<T>(sinh(-ax[0]), Interval(-yu,-yl));
}
template<class T>
void  sinh01() { CHECK_sinh<T>(Interval()); }
template<class T>
void  sinh02() { CHECK_sinh<T>(Interval(0,oo)); }
template<class T>
void  sinh03() { CHECK_sinh<T>(Interval(0,1)); }
template<class T>
void  sinh04() { CHECK_sinh<T>(Interval(1,oo)); }
template<class T>
void  sinh05() { CHECK_sinh<T>(Interval(1,1)); }
template<class T>
void  sinh06() { CHECK_sinh<T>(Interval(2,3)); }
template<class T>
void  sinh07() { CHECK_sinh<T>(Interval(4,5)); }

template<class T>
void  CHECK_cosh(const Interval& x) {

	Interval y ;
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


	//	std::cout << x<<" : " << cosh(x)<< "|||  "<<y<<std::endl;
	AffineVarMainVector<T> ax(1,x);
	CHECK_affine_inclu<T>(cosh(ax[0]), y);
	CHECK_affine_inclu<T>(cosh(-ax[0]), y);
	//	CHECK(x,acosh(cosh(x)));
	//	CHECK(x,acosh(cosh(-x)));
}
template<class T>
void  cosh01() { CHECK_cosh<T>(Interval()); }
template<class T>
void  cosh02() { CHECK_cosh<T>(Interval(0,oo)); }
template<class T>
void  cosh03() { CHECK_cosh<T>(Interval(0,1)); }
template<class T>
void  cosh04() { CHECK_cosh<T>(Interval(1,oo)); }
template<class T>
void  cosh05() { CHECK_cosh<T>(Interval(1,1)); }
template<class T>
void  cosh06() { CHECK_cosh<T>(Interval(2,3)); }
template<class T>
void  cosh07() { CHECK_cosh<T>(Interval(4,5)); }

template<class T>
void  CHECK_trigo(const Interval& x, const Interval& sin_x_expected) {

	AffineVarMainVector<T> ax(1,x);
	CHECK_affine_inclu<T>(sin(ax[0]), sin_x_expected);
	CHECK_affine_inclu<T>(sin(Interval::pi()-ax[0]), sin_x_expected);
	CHECK_affine_inclu<T>(sin(ax[0]+Interval::two_pi()), sin_x_expected);
	CHECK_affine_inclu<T>(sin(-ax[0]), -sin_x_expected);
	CHECK_affine_inclu<T>(cos(ax[0]-Interval::half_pi()), sin_x_expected);
	CHECK_affine_inclu<T>(cos(Interval::half_pi()-ax[0]), sin_x_expected);
	CHECK_affine_inclu<T>(cos(ax[0]+Interval::half_pi()), -sin_x_expected);
	CHECK_affine_inclu<T>(cos(ax[0]+Interval::two_pi()-Interval::half_pi()), sin_x_expected);
}
template<class T>
void  log01() {
	AffineVarMainVector<T> ax(1,Interval::empty());
	CHECK_affine_inclu<T>(log(ax[0]), Interval::empty()); }
template<class T>
void  log02() {
	AffineVarMainVector<T> ax(1,Interval());
	CHECK_affine_inclu<T>(log(ax[0]), Interval()); }
template<class T>
void  log03() {
	AffineVarMainVector<T> ax(1,Interval(0,oo));
	CHECK_affine_inclu<T>(log(ax[0]), Interval()); }
template<class T>
void  log04() {
	AffineVarMainVector<T> ax(1,Interval(-oo,0));
	CHECK_affine_inclu<T>(log(ax[0]),  Interval::empty()); /*Interval(-oo,-MAX_DOUBLE));*/ }
template<class T>
void  log05() {
	AffineVarMainVector<T> ax(1,Interval(1,2));
	CHECK_affine_inclu<T>(log(ax[0]),       Interval(0,std::log(2))); }
template<class T>
void  log06() {
	AffineVarMainVector<T> ax(1,Interval(-1,1));
	CHECK_affine_inclu<T>(log(ax[0]),     Interval(-oo,0)); }
template<class T>
void  log07() {
	AffineVarMainVector<T> ax(1,Interval(0,codac2::next_float(0)));
	CHECK(log(ax[0]).itv().ub()> -744.5); }
template<class T>
void  log08() {
	AffineVarMainVector<T> ax(1,Interval(0,1));
	CHECK_affine_inclu<T>(log(ax[0]),       Interval(-oo,0)); }
template<class T>
void  log09() {
	AffineVarMainVector<T> ax(1,Interval(1,oo));
	CHECK_affine_inclu<T>(log(ax[0]),  Interval(0,oo)); }
template<class T>
void  log10() {
	AffineVarMainVector<T> ax(1,Interval(0));
	CHECK_affine_inclu<T>(log(ax[0]), Interval::empty()); /* Interval(-oo,-MAX_DOUBLE)); */ }
template<class T>
void  log11() {
	AffineVarMainVector<T> ax(1,Interval(-2,-1));
	CHECK_affine_inclu<T>(log(ax[0]), Interval::empty()); }
template<class T>
void  exp01() {
	AffineVarMainVector<T> ax(1,Interval::empty());
	CHECK_affine_inclu<T>(exp(ax[0]), Interval::empty()); }
template<class T>
void  exp02() {
	AffineVarMainVector<T> ax(1,Interval());
	CHECK_affine_inclu<T>(exp(ax[0]), Interval(0,oo)); }
template<class T>
void  exp03() {
	AffineVarMainVector<T> ax(1,Interval(0,oo));
	CHECK_affine_inclu<T>(exp(ax[0]), Interval(1,oo)); }
template<class T>
void  exp04() {
	AffineVarMainVector<T> ax(1,Interval(-oo,0));
	CHECK_affine_inclu<T>(exp(ax[0]),  Interval(0,1)); }
template<class T>
void  exp05() {
	AffineVarMainVector<T> ax(1,Interval(0,2));
	CHECK_affine_inclu<T>(exp(ax[0]),    Interval(1,std::exp(2))); }
template<class T>
void  exp06() {
	AffineVarMainVector<T> ax(1,Interval(-1,1));
	CHECK_affine_inclu<T>(exp(ax[0]),     Interval(std::exp(-1),std::exp(1))); }
template<class T>
void  exp07() {
	AffineVarMainVector<T> ax(1,Interval(1.e100,1.e111));
	CHECK_affine_inclu<T>(exp(ax[0]), Interval(MAX_DOUBLE,oo)); }
template<class T>
void  exp08() {
	AffineVarMainVector<T> ax(1,Interval(MAX_DOUBLE,oo));
	CHECK_affine_inclu<T>(exp(ax[0]), Interval(MAX_DOUBLE,oo)); }
template<class T>
void  exp09() {
	AffineVarMainVector<T> ax(1,Interval(0, MAX_DOUBLE));
	CHECK_affine_inclu<T>(exp(ax[0]), Interval(1,oo)); }
template<class T>
void  sin01() { CHECK_trigo<T>(Interval(), Interval(-1,1)); }
template<class T>
void  sin02() { CHECK_trigo<T>(Interval::empty(), Interval::empty()); }
template<class T>
void  sin03() { CHECK_trigo<T>(Interval(0,piU/2.0), Interval(0,1)); }
template<class T>
void  sin04() { CHECK_trigo<T>(Interval(0,piU), Interval(0,1)); }
template<class T>
void  sin05() { CHECK_trigo<T>(Interval(0,3*piU/2.0), Interval(-1,1)); }
template<class T>
void  sin06() { CHECK_trigo<T>(Interval(piL,3*piU/2.0), Interval(-1,0)); }
template<class T>
void  sin07() { CHECK_trigo<T>(Interval(0.5,1.5), Interval(std::sin(0.5),std::sin(1.5))); }
template<class T>
void  sin08() { CHECK_trigo<T>(Interval(1.5,3), Interval(std::sin(3.0),1)); }
template<class T>
void  sin09() { CHECK_trigo<T>(Interval(3,4), Interval(std::sin(4.0),std::sin(3.0))); }
template<class T>
void  sin10() { CHECK_trigo<T>(Interval(3,5), Interval(-1,std::sin(3.0))); }
template<class T>
void  sin11() { CHECK_trigo<T>(Interval(3,2*piU+1.5), Interval(-1,std::sin(1.5))); }
template<class T>
void  sin12() { CHECK_trigo<T>(Interval(5,2*piU+1.5), Interval(std::sin(5.0),std::sin(1.5))); }
template<class T>
void  sin13() { CHECK_trigo<T>(Interval(5,2*piU+3), Interval(std::sin(5.0),1)); }
template<class T>
void  tan01() {
	AffineVarMainVector<T> ax(1,Interval());
	CHECK_affine_inclu<T>(tan(ax[0]), Interval()); }
template<class T>
void  tan02() {
	AffineVarMainVector<T> ax(1,(-Interval::pi()/4.0 | Interval::pi()/4.0));
	CHECK_affine_inclu<T>(tan(ax[0]), Interval(-1,1)); }
template<class T>
void  tan03() {	// tan(pi/4,pi/2)=[1,+oo)
	AffineVarMainVector<T> x(1,Interval(piL/4.0,(1-1e-10)*piL/2.0)); // upper bound of x is close to pi/2
	AffineMain<T>  y=tan(x[0]);
	CHECK(y.itv().lb()<=1.0);
	CHECK(y.itv().ub()>1000); // upper bound of tan(x) is close to +oo
}
template<class T>
void  tan04() {	// tan(-pi/2,pi/4)=(-oo,1]
	AffineVarMainVector<T> ax(1,Interval(-(1-1e-10)*piL/2.0,piL/4.0));
	Interval y= (tan( ax[0] )).itv();
	CHECK(y.lb()<-1000); // lower bound is close to -oo
	CHECK(y.ub()>=1.0);
}
template<class T>
void  tan05() {
	AffineVarMainVector<T> ax(1,Interval::pi()/2.0);
	CHECK_affine_inclu<T>(tan(ax[0]),Interval());
}
template<class T>
void  tan06() {
	AffineVarMainVector<T> ax(1,Interval::pi());
	CHECK_affine_inclu<T>(tan(-ax[0]),Interval());
}
template<class T>
void  tan07() {
	AffineVarMainVector<T> ax(1,(3*Interval::pi()/4.0 | 5*Interval::pi()/4.0));
	CHECK_affine_inclu<T>(tan(ax[0]), Interval(-1,1));
}

template<class T>
void  CHECK_pow(const Interval& x, int p, const Interval& y_expected) {
	AffineVarMainVector<T> ax(1,x);
	CHECK_affine_inclu<T>(pow(ax[0],p),y_expected);
	CHECK_affine_inclu<T>(pow(-ax[0],p),(p%2==0)? y_expected : -y_expected);
}
template<class T>
void  pow01() { CHECK_pow<T>(Interval(), 4, Interval(0,oo)); }
template<class T>
void  pow02() { CHECK_pow<T>(Interval::empty(), 4, Interval::empty()); }
template<class T>
void  pow03() { CHECK_pow<T>(Interval(2,3), 4, Interval(16,81)); }
template<class T>
void  pow04() { CHECK_pow<T>(Interval(-2,3), 4, Interval(0,81)); }
template<class T>
void  pow05() { CHECK_pow<T>(Interval(-3,2), 4, Interval(0,81)); }
template<class T>
void  pow06() { CHECK_pow<T>(Interval(2,oo), 4, Interval(16,oo)); }
template<class T>
void  pow07() { CHECK_pow<T>(Interval(), 3, Interval()); }
template<class T>
void  pow08() { CHECK_pow<T>(Interval::empty(), 3, Interval::empty()); }
template<class T>
void  pow09() { CHECK_pow<T>(Interval(2,3), 3, Interval(8,27)); }
template<class T>
void  pow10() { CHECK_pow<T>(Interval(-2,3), 3, Interval(-8,27)); }
template<class T>
void  pow11() { CHECK_pow<T>(Interval(-3,2), 3, Interval(-27,8)); }
template<class T>
void  pow12() { CHECK_pow<T>(Interval(2,oo), 3, Interval(8,oo)); }
template<class T>
void  pow13() { CHECK_pow<T>(Interval(-10,10), -2, Interval(1.0/100,oo)); }
template<class T>
void  root01() {
	AffineVarMainVector<T> ax(1, Interval(0,1));
	CHECK_affine_inclu<T>(root(ax[0],-1), Interval(1.0,oo));
}
template<class T>
void  root02() {
	AffineVarMainVector<T> ax(1, Interval(-27,-8));
	CHECK_affine_inclu<T>(root(ax[0], 3),Interval(-3,-2));
}
template<class T>
void  root03() {
	AffineVarMainVector<T> ax(1, Interval(-4,1));
	CHECK_affine_inclu<T>(root(ax[0],2), Interval(0,1));
}
template<class T>
void  root04() {
	AffineVarMainVector<T> ax(1, Interval(-8,1));
	CHECK_affine_inclu<T>(root(ax[0],3), Interval(-2,1));
}


template<class T>
void  tan_issue248() {
	//Interval itv = Interval(-Interval::pi().lb()/2,3*Interval::pi().ub()/8);
	AffineVarMainVector<T> aff(1,Interval(-1.57079632679489678, 1.1780972450961728626));
	CHECK_FALSE((tan(aff[0]).is_empty()));
}


TEST_CASE("AffineForm operations")
{

typedef AF_Default  AA;

CHECK_change_mode_MinRange<AA>();

/* test: operator&=(const Interval& x) */
 minus01<AA>();
 minus02<AA>();
 minus03<AA>();
 minus04<AA>();

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
 add01<AA>();
 add02<AA>();
 add03<AA>();
 add04<AA>();
 add05<AA>();
 add06<AA>();
 add07<AA>();
 add08<AA>();
 add09<AA>();
 add10<AA>();
 add11<AA>();
 add12<AA>();
 add13<AA>();
 add14<AA>();
 add15<AA>();
 add16<AA>();
 add17<AA>();
 add18<AA>();
 add19<AA>();

/* test:
 * =======
 *   operator*=(double d)
 *   operator*=(const Interval& x)
 *   operator*(const Interval& x, double d)
 *   operator*(double d,const Interval& x)
 *   operator*(const Interval& x1, const Interval& x2)
 */
 mul01<AA>();
 mul02<AA>();
 mul03<AA>();
 mul04<AA>();
 mul05<AA>();
 mul06<AA>();
 mul07<AA>();
 mul08<AA>();
 mul09<AA>();
 mul10<AA>();
 mul11<AA>();
 mul12<AA>();
 mul13<AA>();
 mul14<AA>();
 mul15<AA>();
 mul16<AA>();
 mul17<AA>();
 mul18<AA>();
 mul19<AA>();

/* test:
 * =======
 *   operator/=(double d)
 *   operator/=(const Interval& x)
 *   operator/(const Interval& x, double d)
 *   operator/(double d,const Interval& x)
 *   operator/(const Interval& x1, const Interval& x2)
 */

 div01<AA>();
 div02<AA>();
 div03<AA>();
 div04<AA>();
 div05<AA>();
 div06<AA>();
 div07<AA>();
 div08<AA>();
 div09<AA>();
 div10<AA>();
 div11<AA>();
 div12<AA>();
 div13<AA>();
 div14<AA>();
 div15<AA>();
 div16<AA>();
 div17<AA>();


 log01<AA>();
 log02<AA>();
 log03<AA>();
 log04<AA>();
 log05<AA>();
 log06<AA>();
 log07<AA>();
 log08<AA>();
 log09<AA>();
 log10<AA>();
 log11<AA>();

 exp01<AA>();
 exp02<AA>();
 exp03<AA>();
 exp04<AA>();
 exp05<AA>();
 exp06<AA>();
 exp07<AA>();
 exp08<AA>();
 exp09<AA>();

/* test:
 * sin
 * cos
 */
 sin01<AA>();
 sin02<AA>();
 sin03<AA>();
 sin04<AA>();
 sin05<AA>();
 sin06<AA>();
 sin07<AA>();
 sin08<AA>();
 sin09<AA>();
 sin10<AA>();
 sin11<AA>();
 sin12<AA>();
 sin13<AA>();

 tan01<AA>();
 tan02<AA>();
 tan03<AA>();
 tan04<AA>();
 tan05<AA>();
 tan06<AA>();
 tan07<AA>();

 pow01<AA>();
 pow02<AA>();
 pow03<AA>();
 pow04<AA>();
 pow05<AA>();
 pow06<AA>();
 pow07<AA>();
 pow08<AA>();
 pow09<AA>();
 pow10<AA>();
 pow11<AA>();
 pow12<AA>();
 pow13<AA>();

 root01<AA>();
 root02<AA>();
 root03<AA>();
 root04<AA>();

/* test: sqrt(const Interval& x); */
 sqrt01<AA>();
 sqrt02<AA>();
 sqrt03<AA>();
 sqrt04<AA>();
 sqrt05<AA>();
 sqrt06<AA>();

 sinh01<AA>();
 sinh02<AA>();
 sinh03<AA>();
 sinh04<AA>();
 sinh05<AA>();
 sinh06<AA>();
 sinh07<AA>();

 cosh01<AA>();
 cosh02<AA>();
 cosh03<AA>();
 cosh04<AA>();
 cosh05<AA>();
 cosh06<AA>();
 cosh07<AA>();


 tan_issue248<AA>();



}


