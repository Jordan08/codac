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

/**
 * \class AF_fAF2
 * \brief Core coefficient storage of the "fAF2" affine-arithmetic model.
 *
 * Stores the vector of noise-symbol coefficients and the accumulated
 * rounding/linearization error of an affine form, together with the
 * error-free floating-point transformations (\c twoSum, \c twoProd) used to
 * keep the model rigorous. The buffer is exclusively owned (non-copyable,
 * movable only), which is why \c AffineMain and \c AffineVarMain always
 * manipulate it through move semantics.
 */
class AF_fAF2 {

private:
	
	friend class AffineVarMain<AF_fAF2>;
	friend class AffineMain<AF_fAF2>;

	std::unique_ptr<double[]> _val; 		// vector of noise symbols of the affine form
	double _err; 	// error of the affine form, corresponded to the last term

	/**
	 * \brief Returns the exact rounding error of the addition of two floating-point values.
	 *
	 * \param a first operand
	 * \param b second operand
	 * \param res rounded floating-point sum \f$a+b\f$
	 * \return the exact rounding error \f$\varepsilon\f$ such that \f$a+b=\mathrm{res}+\varepsilon\f$ in infinite precision
	 */
	static double twoSum(double a, double b, double *res);

	/**
	 * \brief Returns the exact rounding error of the multiplication of two floating-point values.
	 *
	 * Uses a hardware FMA instruction when available at runtime (see
	 * \c codac_fma_runtime), and falls back to a Split-based Dekker product
	 * otherwise.
	 *
	 * \param x first operand
	 * \param y second operand
	 * \param r_1 rounded floating-point product \f$x \cdot y\f$
	 * \return the exact rounding error \f$\varepsilon\f$ such that \f$x \cdot y=r_1+\varepsilon\f$ in infinite precision
	 */
	static double twoProd(double x, double y, double *r_1);

	/**
	 * \brief Splits a double-precision value into a high and a low part.
	 *
	 * Used internally by the non-FMA fallback implementation of \c twoProd
	 * (Dekker's error-free splitting, see "Handbook of Floating-Point
	 * Arithmetic" p.132-139).
	 *
	 * \param x value to split
	 * \param sp number of significant bits kept in \p x_high (splitting exponent)
	 * \param x_high high-order part of \p x
	 * \param x_low low-order part of \p x, such that \f$x=x_{high}+x_{low}\f$
	 */
	static void Split(double x, int sp, double *x_high, double *x_low);

	/**
	 * \brief Creates an affine core from coefficients and remainder error.
	 *
	 * \param val vector of noise-symbol coefficients, whose ownership is transferred
	 * \param err accumulated remainder error of the affine form
	 */
	AF_fAF2(std::unique_ptr<double[]> val, double err);

public:

	/** \brief Copy constructor, explicitly disabled: an \c AF_fAF2 core exclusively owns its coefficient buffer. */
	AF_fAF2(const AF_fAF2&) = delete;

	/** \brief Copy assignment, explicitly disabled: an \c AF_fAF2 core exclusively owns its coefficient buffer. */
	AF_fAF2& operator=(const AF_fAF2&) = delete;

	/**
	 * \brief Move constructor.
	 *
	 * Transfers ownership of the coefficient buffer from \p other, which is
	 * left with a null buffer afterwards.
	 *
	 * \param other affine core to move from
	 */
	AF_fAF2(AF_fAF2&& other) noexcept;

	/**
	 * \brief Move assignment.
	 *
	 * Transfers ownership of the coefficient buffer from \p other, which is
	 * left with a null buffer afterwards.
	 *
	 * \param other affine core to move from
	 * \return a reference to this
	 */
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

inline AF_fAF2& AF_fAF2::operator=(AF_fAF2&& other) noexcept {
  if(this != &other)  {
    _val = std::exchange(other._val, nullptr);
    _err = other._err;
  }
  return *this;
}


inline AF_fAF2::AF_fAF2(std::unique_ptr<double[]> val, double err) :
	_val	(std::move(val)),
	_err	(err) {
}




/**
 * \brief Checks, once per program execution, whether the CPU supports the FMA instruction set.
 *
 * The runtime CPU-feature check is expensive, so it is performed only once
 * (on first call) and its result is cached in a function-local static
 * variable for all subsequent calls. Note that this only reflects hardware
 * support: the FMA code path is only ever emitted by the compiler for
 * translation units also built with the matching flag (e.g. ``-mfma``).
 *
 * \return true if FMA is available at runtime, false otherwise
 */
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
