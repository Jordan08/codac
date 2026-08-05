/**
 *  \file codac2_Affine2_fAF2.h
 *
 *  Core storage and exact floating-point helpers for the AF_fAF2 affine model.
 *
 * ----------------------------------------------------------------------------
 *  \date       2020
 *  \author     Jordan Ninin
 *  \copyright  Copyright 2020 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */


#pragma once

#include "codac2_math.h"

#ifndef _MSC_VER
	#ifdef __FP_FAST_FMA
		#define codac2_FMA
	#endif
	#ifdef FP_FAST_FMA
		#define codac2_FMA
	#endif

#else
#if (_MSC_VER >= 1800)
	#define codac2_FMA
#endif
#endif


namespace codac2 {

class AF_fAF2;

//typedef AF_fAF1  AF_Default;
typedef AF_fAF2  AF_Default;
//typedef AF_fAF2_fma  AF_Default;
//typedef AF_iAF  AF_Default;
//typedef AF_sAF  AF_Default;
//typedef AF_No  AF_Default;
//typedef AF_fAFFullI AF_Default;



template<class T>  class AffineMain;
template<class T>  class AffineVarMain;

class AF_fAF2 {

private:
	
	friend class AffineVarMain<AF_fAF2>;
	friend class AffineMain<AF_fAF2>;

	/**
	 * Code for the particular case:
	 * if the affine form is actif, _actif=1  and _n is the size of the affine form
	 * if the set is degenerate, _actif = 0 and itv().diam()< AF_EC
	 * if the set is empty, _actif = -1
	 * if the set is ]-oo,+oo[, _actif = -2 and _err =]-oo,+oo[
	 * if the set is [a, +oo[ , _actif = -3 and _err = [a, +oo[
	 * if the set is ]-oo, a] , _actif = -4 and _err = ]-oo, a]
	 *
	 */

	double * _val; 		// vector of elements of the affine form
	double _err; 	// error of the affine form, corresponded to the last term
	//	bool _actif; // boolean to know if the affine form is actif or not. This is to manage the particular case of EMPTY and an unbounded Interval

	/**
	 * \brief Returns the exact rounding error of the addition of two floating-point values.
	 */
	static double twoSum(double a, double b, double *res);

	/**
	 * \brief Returns the exact rounding error of the multiplication of two floating-point values.
	 */
	static double twoProd(double a, double b, double *res);
	static void Split(double x, int sp, double *x_high, double *x_low);

	/** \brief Creates an affine core from coefficients and remainder error. */
	AF_fAF2(double * val, double err);

public:
	/** \brief Destroys the affine core. */
	virtual ~AF_fAF2();

};


inline AF_fAF2::AF_fAF2(double * val, double err) :
	_val	(val ),
	_err	(err) {

}



inline AF_fAF2::~AF_fAF2() {
	if (_val!=nullptr) delete[] _val;
}



/////////////////////
// CODE extract from "Handbook of Floating-Point Arithmetic" p.132-139
inline void AF_fAF2::Split(double x, int sp, double *x_high, double *x_low)
{
	unsigned long C = (1UL << sp) + 1;
	double gamma = (C * x);
	double delta = (x - gamma);
	*x_high= (gamma + delta);
	*x_low= (x - *x_high);
}

inline double AF_fAF2::twoProd(double x, double y, double *r_1)
{
#ifdef codac2_FMA
	*r_1 = (x * y);
	return std::fma(x,y,-(*r_1));
#else

	int SHIFT_POW = 27; //  53 / 2 for double precision.
	double x_high, x_low;
	double y_high, y_low;
	double t_1;
	double t_2;
	double t_3;
	Split(x, SHIFT_POW, &x_high, &x_low);
	Split(y, SHIFT_POW, &y_high, &y_low);
	*r_1 = (x * y);
	t_1 = (-*r_1 + x_high * y_high);
	t_2 =   (t_1 + x_high * y_low );
	t_3 =	(t_2 + x_low  * y_high);
	return  (t_3 + x_low  * y_low );

#endif

}



// CODE extract from "Handbook of Floating-Point Arithmetic" p.130
inline double AF_fAF2::twoSum(double a, double b, double *res) {
	*res = (a+b);
	double a2 = (*res - b);
	double b2 = (*res - a2);
	double delta_a = (a - a2);
	double delta_b = (b - b2);
	return (delta_a + delta_b);
}

//////////////////////

}