/**
 *  \file codac2_Affine2_fAF2.h
 *
 *  Core storage and exact floating-point helpers for the AF_fAF2 affine model.
 *
 * ----------------------------------------------------------------------------
 *  \date       2026
 *  \author     Jordan Ninin
 *  \copyright  Copyright 2026 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */


#pragma once

#include "codac2_math.h"
#include <memory>
#include <utility>


namespace codac2 {


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

	std::unique_ptr<double[]> _val; 		// vector of elements of the affine form
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

	AF_fAF2(const AF_fAF2&) = delete;
	AF_fAF2& operator=(const AF_fAF2&) = delete;
	AF_fAF2(AF_fAF2&& other) noexcept;
	AF_fAF2& operator=(AF_fAF2&& other) noexcept;

	/** \brief Destroys the affine core. */
	~AF_fAF2(){
	// if (_val!=nullptr) delete[] _val;
	// unique_ptr<double[]> _val cleans itself up; nothing left to do here.
	};
	

};

inline AF_fAF2::AF_fAF2(AF_fAF2&& other) noexcept
  : _val(std::exchange(other._val, nullptr)),
    _err(other._err)
{
}

inline AF_fAF2& AF_fAF2::operator=(AF_fAF2&& other) noexcept
{
  if(this != &other)
  {
    _val = std::exchange(other._val, nullptr);
    _err = other._err;
  }

  return *this;
}


inline AF_fAF2::AF_fAF2(double * val, double err) :
	_val	(val ),
	_err	(err) {

}





// check one time if  fma is available at runtime, and store the result in a static variable.
// This avoids repeated checks and ensures that the FMA support is determined only once during the program's execution.
inline bool codac_fma_runtime() noexcept
{
    static const bool has_fma = []() noexcept {
#if defined(_MSC_VER) && (defined(_M_X64) || defined(_M_IX86))

        int cpu_info[4] = {};
        __cpuid(cpu_info, 0);
        if (cpu_info[0] < 1) {
            return false;
        }

        __cpuidex(cpu_info, 1, 0);
        constexpr int FMA_BIT     = 1 << 12;
        constexpr int OSXSAVE_BIT = 1 << 27;
        constexpr int AVX_BIT     = 1 << 28;

        if ((cpu_info[2] & FMA_BIT) == 0 ||
            (cpu_info[2] & OSXSAVE_BIT) == 0 ||
            (cpu_info[2] & AVX_BIT) == 0) {
            return false;
        }

        constexpr unsigned __int64 AVX_STATE =(1ULL << 1) | (1ULL << 2);

        return (_xgetbv(0) & AVX_STATE) == AVX_STATE;

#elif (defined(__GNUC__) || defined(__clang__)) && (defined(__x86_64__) || defined(__i386__))
        __builtin_cpu_init();
        return __builtin_cpu_supports("fma");

#elif defined(__aarch64__) || defined(_M_ARM64)
        return true;

#else
        return false;

#endif
    }();
    return has_fma;
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
	if (codac_fma_runtime() ) {
		*r_1 = (x * y);
		return std::fma(x,y,-(*r_1));

	} else {

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

	}

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
