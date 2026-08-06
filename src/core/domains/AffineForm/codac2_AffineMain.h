/**
 *  \file codac2_AffineMain.h
 *
 *  Main affine-form class and arithmetic operators.
 *
 * ----------------------------------------------------------------------------
 *  \date       2026
 *  \author     Jordan Ninin
 *  \copyright  Copyright 2026 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include "codac2_Interval.h"
#include <cmath>
#include <ostream>
#include <iostream>


#include "codac2_Affine2_fAF2.h"
#include "codac2_matrices.h"


namespace codac2 {


template<class T> class AffineVarMain;
template<class T> class AffineMain;


template<class T>
struct is_interval_based< AffineMain<T> > : std::true_type {};

template<class T>
struct is_ctc< AffineMain<T> > : std::false_type {};

template<class T>
struct is_sep< AffineMain<T> > : std::false_type {};

} // namespace codac2

/*
 * The following specializations make ``AffineMain<T>`` usable as an Eigen
 * Scalar type (as done for ``codac2::Interval`` in codac2_matrices.h), which
 * is required for ``AffineMainVector<T>``/``AffineMainMatrix<T>`` (Eigen
 * matrices of AffineMain elements) to compile: without them, generic Eigen
 * algorithms (e.g. matrix printing, precision handling) instantiate things
 * like ``log2(NumTraits<AffineMain<T>>::epsilon())`` which do not compile,
 * since AffineMain<T> has no arithmetic conversion to double.
 * Every member below is inherited unchanged from ``NumTraits<double>``
 * (only ``Real``/``Nested``/``Scalar`` are aliased to AffineMain<T>), exactly
 * like it is done for ``codac2::Interval``.
 */
namespace Eigen
{
  template<class T>
  struct NumTraits<codac2::AffineMain<T>>
   : NumTraits<double> // permits to get the epsilon, dummy_precision, lowest, highest functions
  {
    typedef codac2::AffineMain<T> Real;
    typedef codac2::AffineMain<T> NonInteger;
    typedef codac2::AffineMain<T> Nested;
    typedef codac2::AffineMain<T> Scalar;
    typedef double RealScalar;

    enum {
      IsComplex = 0,
      IsInteger = 0,
      IsSigned = 1,
      RequireInitialization = 1,
      ReadCost = 1,
      AddCost = 3,
      MulCost = 3
    };
  };

  template<class T,typename BinOp>
  struct ScalarBinaryOpTraits<codac2::AffineMain<T>,double,BinOp>
  { typedef codac2::AffineMain<T> ReturnType; };

  template<class T,typename BinOp>
  struct ScalarBinaryOpTraits<double,codac2::AffineMain<T>,BinOp>
  { typedef codac2::AffineMain<T> ReturnType; };
}

namespace codac2 {

/**
 * \ingroup arithmetic
 *
 * \brief Main affine-form domain based on \c AF_Default.
 *
 * Affine forms rely on a compact status code:
 * - \c _actif=1: active affine form, \c _n is the number of variables
 * - \c _actif=0: degenerate interval
 * - \c _actif=-1: empty set
 * - \c _actif=-2: \f$[-\infty,+\infty]\f$
 * - \c _actif=-3: \f$[a,+\infty]\f$ encoded with \c _err=a
 * - \c _actif=-4: \f$[-\infty,a]\f$ encoded with \c _err=a
 */


//=================================================================================================================
//=================================================================================================================
//=================================================================================================================

template<class T>
class AffineMain  : public DomainInterface<Interval,double>
{

protected:


	/** \brief Tolerance for the default compacting procedure. */
	static const double AF_COMPAC_Tol;
	static const double AF_EM;
	static const double AF_EC;
	static const double AF_EE;

	static bool mode;

	/**
	 * \brief Encodes the current affine-state mode.
	 *
	 * See class description above for the exact code-to-state mapping.
	 */
	int _actif;		// boolean to know if the affine form is actif or not.
	// This is to manage the particular case of EMPTY and an unbounded Interval

	int _n; 		// dimension (size of _elt._val)-1  , ie number of variable

	T _elt;			// core of the affine form

	/**
	 * \brief Creates an affine form with \p size variables.
	 *
	 * The \p var-th variable is initialized from \p itv.
	 */
	AffineMain(int size, int var, const Interval& itv);

	/**
	 * \brief Changes the number of affine variables.
	 *
	 * \param n new number of variables
	 */
	void resize(int n);

public:

	typedef enum {
		AF_Lin_Default=0, AF_Lin_Chebyshev=1, AF_Lin_MinRange=2
	} Affine_Mode; // ...etc...

	/**
	 * \brief Changes the linearization mode globally.
	 *
	 * \param tt approximation mode (Chebyshev by default)
	 */
	static void change_mode(Affine_Mode tt=AF_Lin_Default);
	/**
	 * \brief Get the linearization mode globally.
	 *
	 * \return the approximation mode
	 */
	static Affine_Mode get_mode();

	/** \brief Creates an unbounded affine form (like \c Interval()). */
	AffineMain();

	/**
	 * \brief Creates an affine form by copy.
	 *
	 * \param x affine form to copy
	 */
	AffineMain(const AffineMain& x);

	/** \brief Destroys the affine form. */
	virtual ~AffineMain() { };

	/** \brief Returns \f$-*\mathrm{this}\f$. */
	virtual AffineMain operator-() const;

	/**
	 * \brief Tests exact equality with another affine form.
	 *
	 * \param x affine form to compare with
	 * \return true iff both interval enclosures are exactly equal
	 */
	bool operator==(const AffineMain& x) const;

	/**
	 * \brief Tests exact equality with an interval.
	 *
	 * \param x interval to compare with
	 * \return true iff both intervals are exactly equal
	 */
	bool operator==(const Interval& x) const;

	/**
	 * \brief Tests non-equality with another affine form.
	 *
	 * \param x affine form to compare with
	 * \return true iff interval enclosures are different
	 */
	bool operator!=(const AffineMain& x) const;

	/**
	 * \brief Tests non-equality with an interval.
	 *
	 * \param x interval to compare with
	 * \return true iff intervals are different
	 */
	bool operator!=(const Interval& x) const ;

	/** \brief Sets this affine form to the empty set. */
	void set_empty();

	/**
	 * \brief Assigns from another affine form.
	 *
	 * \param x affine form to copy
	 * \return a reference to this
	 */
	virtual AffineMain& operator=(const AffineMain& x);

	/**
	 * \brief Assigns from a real value.
	 *
	 * \param x real value
	 * \return a reference to this
	 */
	AffineMain& operator=(double x);

	/**
	 * \brief Assigns from an interval.
	 *
	 * \param itv interval value
	 * \return a reference to this
	 */
	virtual AffineMain& operator=(const Interval& itv);

	/* Union and Intersection of two Affine form must not be implemented
	 * That could produce to much confusion.
	 */
	/** \brief Intersection of *this and x.
	 * \param x - the interval to compute the intersection with.*/
	//	AffineMain& operator&=(const AffineMain& x);
	/** \brief Union of *this and I.
	 * \param x - the interval to compute the hull with.*/
	//	AffineMain& operator|=(const AffineMain& x);

	/**
	 * \brief Adds \f$[-\mathrm{rad},+\mathrm{rad}]\f$ to this.
	 *
	 * \param radd inflation radius
	 * \return a reference to this
	 */
	AffineMain& inflate(double radd);

	/**
	 * \brief Returns the lower bound of this.
	 *
	 * \return lower bound
	 */
	double lb() const;

	/**
	 * \brief Returns the upper bound of this.
	 *
	 * \return upper bound
	 */
	double ub() const;

	/**
	 * \brief Returns the midpoint of this.
	 *
	 * \return midpoint
	 */
	double mid() const;

	/**
	 * \brief Returns the radius of this.
	 *
	 * \return radius, or 0 if this is empty
	 */
	double rad() const;

	/**
	 * \brief Returns the diameter of this.
	 *
	 * \return diameter, or 0 if this is empty
	 */
	double diam() const;

	/**
	 * \brief Returns the mignitude of this.
	 *
	 * The value is:
	 * - +lb(this) if this > 0
	 * - -ub(this) if this < 0
	 * - 0 otherwise
	 *
	 * \return mignitude
	 */
	double mig() const;

	/**
	 * \brief Returns the magnitude of this.
	 *
	 * \return \f$\max(|\mathrm{lb}(x)|,|\mathrm{ub}(x)|)\f$
	 */
	double mag() const;

	/**
	 * \brief Returns the signed magnitude of this.
	 *
	 * \return signed magnitude
	 */
	double smag() const;

	/**
	 * \brief Returns the signed mignitude of this.
	 *
	 * \return signed mignitude
	 */
	double smig() const;

	/**
	 * \brief Volume (alias of diam for 1D domains).
	 */
	double volume() const;

	/**
	 * \brief Returns the number of represented noise variables.
	 *
	 * \return number of variables
	 */
	int size() const;

	/**
	 * \brief Returns the interval enclosure of this affine form.
	 *
	 * \return interval enclosure
	 */
	const Interval itv() const ;

	/**
	 * \brief Returns the i-th noise coefficient.
	 *
	 * \param i coefficient index
	 * \return coefficient value
	 */
	double val(int i) const;

	/**
	 * \brief Returns the remainder error term.
	 *
	 * \return remainder error
	 */
	double err() const;
	/**
	 * \brief Tests whether this affine form is active and valid.
	 *
	 * \return true iff this affine form is active
	 */
	bool is_actif() const;

	/**
	 * \brief Tests whether this affine form is empty.
	 *
	 * \return true iff this is empty
	 */
	bool is_empty() const;

	/**
	 * \brief Tests whether this affine form is degenerate.
	 *
	 * \return true iff this represents a singleton interval
	 */
	bool is_degenerated() const;

	/**
	 * \brief True if one bound of *this is infinite.
	 *
	 * \note An empty affine form is always bounded.
	 */
	bool is_unbounded() const;

	/**
	 * \brief Compacts small coefficients into the remainder term.
	 *
	 * \param tol threshold below which noise coefficients are absorbed
	 */
	void compact(double tol);
	/** \brief Compacts with the default tolerance. */
	void compact();

	/**
	 * \brief Tests whether this interval enclosure is a subset of \p x.
	 *
	 * \note Always returns true if this is empty.
	 */
	bool is_subset(const Interval& x) const;

	/**
	 * \brief Tests whether this interval enclosure is a strict subset of \p x.
	 *
	 * \note In particular, (-oo,oo) is not a strict subset of (-oo,oo)
	 * and the empty set is not a strict subset of the empty set although
	 * in both cases, the first is inside the interior of the second.
	 */
	bool is_strict_subset(const Interval& x) const;

	/**
	 * \brief Tests whether this interval enclosure is inside the interior of \p x.
	 *
	 * \note In particular, (-oo,oo) is in the interior of (-oo,oo)
	 * and the empty set is in the interior of the empty set.
	 * \note Always returns true if this is empty.
	 */
	bool is_interior_subset(const Interval& x) const;

	/**
	 * \brief Tests whether this interval enclosure is in the relative interior of \p x.
	 *
	 * When x is degenerated, the relative interior of x is x itself
	 * (it is not an open set). Otherwise, the relative interior of x
	 * is the interior (in the usual meaning).
	 */
	bool is_relative_interior_subset(const Interval& x) const;

	/**
	 * \brief Tests whether this interval enclosure is a strict interior subset of \p x.
	 *
	 * \note In particular, (-oo,oo) is not "strictly" in the interior of (-oo,oo)
	 * and the empty set is not "strictly" in the interior of the empty set.
	 */
	bool is_strict_interior_subset(const Interval& x) const;

	/**
	 * \brief Tests whether this interval enclosure is a superset of \p x.
	 *
	 * \note Always returns true if x is empty.
	 */
	bool is_superset(const Interval& x) const;

	/**
	 * \brief Tests whether this interval enclosure is a strict superset of \p x.
	 *
	 * \see #is_strict_subset(const Interval&) const.
	 */
	bool is_strict_superset(const Interval& x) const;

	/**
	 * \brief Tests whether this interval enclosure contains \p d.
	 *
	 * \note d can also be an "open bound", i.e., infinity.
	 * So this function is not restricted to a set-membership
	 * interpretation.
	 */
	bool contains(const double& d) const;

	/**
	 * \brief Tests whether the interior of this contains \p d.
	 *
	 * \param d real value
	 * \return true iff the interior contains \p d
	 */
	bool interior_contains(const double& d) const;

	/**
	 * \brief Tests whether this interval enclosure intersects \p x.
	 */
	bool intersects(const Interval &x) const;

	/**
	 * \brief Tests whether this and \p x overlap with non-null volume.
	 *
	 * Equivalently, some interior points (of this or x) must belong to the intersection.
	 */
	bool overlaps(const Interval &x) const;

	/**
	 * \brief Tests whether this and \p x are disjoint.
	 *
	 * \param x interval value
	 * \return true iff this and \p x do not intersect
	 */
	bool is_disjoint(const Interval &x) const;


	/**
	 * \brief Adds \p d to this.
	 *
	 * \param d real value
	 * \return a reference to this
	 */
	virtual AffineMain& operator+=(double d);

	/**
	 * \brief Subtracts \p d from this.
	 *
	 * \param d real value
	 * \return a reference to this
	 */
	virtual AffineMain& operator-=(double d);

	/**
	 * \brief Multiplies this by \p d.
	 *
	 * \param d real value
	 * \return a reference to this
	 */
	virtual AffineMain& operator*=(double d);

	/**
	 * \brief Divides this by \p d.
	 *
	 * \param d real value
	 * \return a reference to this
	 */
	virtual AffineMain& operator/=(double d) ;

	/**
	 * \brief Adds interval \p x to this.
	 *
	 * \param x interval value
	 * \return a reference to this
	 */
	virtual AffineMain& operator+=(const Interval& x);

	/**
	 * \brief Subtracts interval \p x from this.
	 *
	 * \param x interval value
	 * \return a reference to this
	 */
	virtual AffineMain& operator-=(const Interval& x);

	/**
	 * \brief Multiplies this by interval \p x.
	 *
	 * \param x interval value
	 * \return a reference to this
	 */
	virtual AffineMain& operator*=(const Interval& x);

	/**
	 * \brief Divides this by interval \p x.
	 *
	 * \param x interval value
	 * \return a reference to this
	 */
	virtual AffineMain& operator/=(const Interval& x);

	/**
	 * \brief Adds affine form \p x to this.
	 *
	 * \param x affine form value
	 * \return a reference to this
	 */
	virtual AffineMain& operator+=(const AffineMain& x);

	/**
	 * \brief Subtracts affine form \p x from this.
	 *
	 * \param x affine form value
	 * \return a reference to this
	 */
	virtual AffineMain& operator-=(const AffineMain& x);

	/**
	 * \brief Multiplies this by affine form \p x.
	 *
	 * \param x affine form value
	 * \return a reference to this
	 */
	virtual AffineMain& operator*=(const AffineMain& x);

	/**
	 * \brief Divides this by affine form \p x.
	 *
	 * \param x affine form value
	 * \return a reference to this
	 */
	virtual AffineMain& operator/=(const AffineMain& x);





private:
	/** \brief Internal helper: compute square. */
	virtual AffineMain& Asqr(const Interval& itv);
	/** \brief Internal helper: negate. */
	virtual AffineMain&  Aneg();
	/** \brief Internal helper: compute inverse (Chebyshev/MinRange). */
	virtual AffineMain&  Ainv(const Interval& itv);
	/** \brief Internal helper: compute square root. */
	virtual AffineMain&  Asqrt(const Interval& itv);
	/** \brief Internal helper: compute exponential. */
	virtual AffineMain&  Aexp(const Interval& itv);
	/** \brief Internal helper: compute logarithm. */
	virtual AffineMain&  Alog(const Interval& itv);
	/** \brief Internal helper: compute power (integer exponent). */
	virtual AffineMain&  Apow(int n, const Interval& itv);
	/** \brief Internal helper: compute power (real exponent). */
	virtual AffineMain&  Apow(double d, const Interval& itv);
	/** \brief Internal helper: compute power (interval exponent). */
	virtual AffineMain&  Apow(const Interval &y, const Interval& itvx);
	/** \brief Internal helper: compute n-th root. */
	virtual AffineMain&  Aroot(int n, const Interval& itv);
	/** \brief Internal helper: compute cosine. */
	virtual AffineMain&  Acos(const Interval& itv);
	/** \brief Internal helper: compute sine. */
	virtual AffineMain&  Asin(const Interval& itv);
	/** \brief Internal helper: compute tangent. */
	virtual AffineMain&  Atan(const Interval& itv);
	/** \brief Internal helper: compute arc-cosine. */
	virtual AffineMain&  Aacos(const Interval& itv);
	/** \brief Internal helper: compute arc-sine. */
	virtual AffineMain&  Aasin(const Interval& itv);
	/** \brief Internal helper: compute arc-tangent. */
	virtual AffineMain&  Aatan(const Interval& itv);
	/** \brief Internal helper: compute hyperbolic cosine. */
	virtual AffineMain&  Acosh(const Interval& itv);
	/** \brief Internal helper: compute hyperbolic sine. */
	virtual AffineMain&  Asinh(const Interval& itv);
	/** \brief Internal helper: compute hyperbolic tangent. */
	virtual AffineMain&  Atanh(const Interval& itv);
	/** \brief Internal helper: compute absolute value. */
	virtual AffineMain&  Aabs(const Interval& itv);
	/** \brief Internal helper: Chebyshev inverse. */
	virtual AffineMain&  Ainv_CH(const Interval& itv);
	/** \brief Internal helper: Chebyshev square root. */
	virtual AffineMain&  Asqrt_CH(const Interval& itv);
	/** \brief Internal helper: Chebyshev exponential. */
	virtual AffineMain&  Aexp_CH(const Interval& itv);
	/** \brief Internal helper: Chebyshev logarithm. */
	virtual AffineMain&  Alog_CH(const Interval& itv);
	/** \brief Internal helper: MinRange inverse. */
	virtual AffineMain&  Ainv_MR(const Interval& itv);
	/** \brief Internal helper: MinRange square root. */
	virtual AffineMain&  Asqrt_MR(const Interval& itv);
	/** \brief Internal helper: MinRange exponential. */
	virtual AffineMain&  Aexp_MR(const Interval& itv);
	/** \brief Internal helper: MinRange logarithm. */
	virtual AffineMain&  Alog_MR(const Interval& itv);



//	template<class A>   friend AffineMain<A> AffineVarMain<A>::operator-() const ;
	template<class A>   friend AffineMain<A> operator/(double d, const AffineMain<A>& x);
	template<class A>	friend AffineMain<A> inv(const AffineMain<A>&  x);
	template<class A>	friend AffineMain<A> sqr(const AffineMain<A>&  x);
	template<class A>	friend AffineMain<A> sqrt(const AffineMain<A>&  x);
	template<class A>	friend AffineMain<A> exp(const AffineMain<A>&  x);
	template<class A>	friend AffineMain<A> log(const AffineMain<A>&  x);
	template<class A>	friend AffineMain<A> pow(const AffineMain<A>&  x, int n);
	template<class A>	friend AffineMain<A> pow(const AffineMain<A>&  x, double d);
	template<class A>	friend AffineMain<A> pow(const AffineMain<A>& x, const Interval &y);
	template<class A>	friend AffineMain<A> pow(const AffineMain<A>& x, const AffineMain<A>& y);
	template<class A>	friend AffineMain<A> root(const AffineMain<A>&  x, int n);
	template<class A>	friend AffineMain<A> cos(const AffineMain<A>&  x);
	template<class A>	friend AffineMain<A> sin(const AffineMain<A>&  x);
	template<class A>	friend AffineMain<A> tan(const AffineMain<A>&  x);
	template<class A>	friend AffineMain<A> acos(const AffineMain<A>&  x);
	template<class A>	friend AffineMain<A> asin(const AffineMain<A>&  x);
	template<class A>	friend AffineMain<A> atan(const AffineMain<A>&  x);
	template<class A>	friend AffineMain<A> cosh(const AffineMain<A>&  x);
	template<class A>	friend AffineMain<A> sinh(const AffineMain<A>&  x);
	template<class A>	friend AffineMain<A> tanh(const AffineMain<A>&  x);
	template<class A>	friend AffineMain<A> abs(const AffineMain<A>& x);
	template<class A>	friend AffineMain<A> floor(const AffineMain<A>&  x);
	template<class A>	friend AffineMain<A> ceil(const AffineMain<A>&  x);

	/** \note Not yet implemented. */
	template<class A>	friend AffineMain<A> acosh(const AffineMain<A>&  x);
	/** \note Not yet implemented. */
	template<class A>	friend AffineMain<A> asinh(const AffineMain<A>&  x);
	/** \note Not yet implemented. */
	template<class A>	friend AffineMain<A> atanh(const AffineMain<A>&  x);


	template<class A> friend std::ostream& operator<<(std::ostream& os, const AffineMain<A>& x);
};


/** \ingroup arithmetic */
/*@{*/

/**
 * \brief Streams out an affine form.
 *
 * \param os output stream
 * \param x affine form value
 * \return a reference to the updated stream
 */
template<class T>
std::ostream& operator<<(std::ostream& os, const AffineMain<T>&  x);

template<class T>
inline double AffineMain<T>::lb() const { return this->itv().lb(); }

template<class T>
inline double AffineMain<T>::ub() const { return this->itv().ub(); }

template<class T>
inline double AffineMain<T>::rad() const { return this->itv().rad(); }

template<class T>
inline double AffineMain<T>::diam() const { return this->itv().diam(); }

template<class T>
inline double AffineMain<T>::mig() const { return this->itv().mig(); }

template<class T>
inline double AffineMain<T>::mag() const { return this->itv().mag(); }

template<class T>
inline double AffineMain<T>::smag() const { return this->itv().smag(); }

template<class T>
inline double AffineMain<T>::smig() const { return this->itv().smig(); }

template<class T>
inline double AffineMain<T>::volume() const { return this->itv().volume(); }

template<class T>
inline bool AffineMain<T>::is_subset(const Interval& x) const { return this->itv().is_subset(x); }

template<class T>
inline bool AffineMain<T>::is_strict_subset(const Interval& x) const { return this->itv().is_strict_subset(x); }

template<class T>
inline bool AffineMain<T>::is_interior_subset(const Interval& x) const { return this->itv().is_interior_subset(x); }

template<class T>
inline bool AffineMain<T>::is_relative_interior_subset(const Interval& x) const { return this->itv().is_relative_interior_subset(x); }

template<class T>
inline bool AffineMain<T>::is_strict_interior_subset(const Interval& x) const { return this->itv().is_strict_interior_subset(x); }

template<class T>
inline bool AffineMain<T>::is_superset(const Interval& x) const { return this->itv().is_superset(x); }

template<class T>
inline bool AffineMain<T>::is_strict_superset(const Interval& x) const { return this->itv().is_strict_superset(x); }

template<class T>
inline bool AffineMain<T>::contains(const double& d) const { return this->itv().contains(d); }

template<class T>
inline bool AffineMain<T>::interior_contains(const double& d) const { return this->itv().interior_contains(d); }

template<class T>
inline bool AffineMain<T>::intersects(const Interval &x) const { return this->itv().intersects(x); }

template<class T>
inline bool AffineMain<T>::overlaps(const Interval &x) const { return this->itv().overlaps(x); }

template<class T>
inline bool AffineMain<T>::is_disjoint(const Interval &x) const { return this->itv().is_disjoint(x); }


/** \brief Returns \f$-x\f$. */
template<class T>
inline AffineMain<T> AffineMain<T>::operator-() const {
	AffineMain<T> res(*this);
	res.Aneg();
	return res;
}

/** \brief Returns \f$\mathrm{AF}[x]_1+\mathrm{AF}[x]_2\f$. */
template<class T>
AffineMain<T> operator+(const AffineMain<T>&  x1, const AffineMain<T>&  x2);

/** \brief Returns \f$\mathrm{AF}[x]+d\f$. */
template<class T>
AffineMain<T> operator+(const AffineMain<T>&  x, double d);

/** \brief Returns \f$d+\mathrm{AF}[x]\f$. */
template<class T>
AffineMain<T> operator+(double d, const AffineMain<T>&  x);

/** \brief Returns \f$\mathrm{AF}[x]_1+[x]_2\f$. */
template<class T>
AffineMain<T> operator+(const AffineMain<T>&  x1, const Interval& x2);

/** \brief Returns \f$[x]_1+\mathrm{AF}[x]_2\f$. */
template<class T>
AffineMain<T> operator+(const Interval& x1, const AffineMain<T>&  x2);

/** \brief AF[x]_1-AF[x]_2. */
template<class T>
AffineMain<T> operator-(const AffineMain<T>&  x1, const AffineMain<T>&  x2);

/** \brief Returns \f$\mathrm{AF}[x]-d\f$. */
template<class T>
AffineMain<T> operator-(const AffineMain<T>&  x, double d);

/** \brief Returns \f$d-\mathrm{AF}[x]\f$. */
template<class T>
AffineMain<T> operator-(double d, const AffineMain<T>&  x);

/** \brief Returns \f$\mathrm{AF}[x]_1-[x]_2\f$. */
template<class T>
AffineMain<T> operator-(const AffineMain<T>&  x1, const Interval& x2);

/** \brief Returns \f$[x]_1-\mathrm{AF}[x]_2\f$. */
template<class T>
AffineMain<T> operator-(const Interval& x1, const AffineMain<T>&  x2);

/** \brief Returns \f$\mathrm{AF}[x]_1\times\mathrm{AF}[x]_2\f$. */
template<class T>
AffineMain<T> operator*(const AffineMain<T>&  x1, const AffineMain<T>&  x2);

/** \brief Returns \f$\mathrm{AF}[x]\times d\f$. */
template<class T>
AffineMain<T> operator*(const AffineMain<T>&  x, double d);

/** \brief Returns \f$d\times \mathrm{AF}[x]\f$. */
template<class T>
AffineMain<T> operator*(double d, const AffineMain<T>&  x);

/** \brief Returns \f$\mathrm{AF}[x]_1\times[x]_2\f$. */
template<class T>
AffineMain<T> operator*(const AffineMain<T>&  x1, const Interval& x2);

/** \brief Returns \f$[x]_1\times\mathrm{AF}[x]_2\f$. */
template<class T>
AffineMain<T> operator*(const Interval& x1, const AffineMain<T>&  x2);

/** \brief Returns \f$\mathrm{AF}[x]_1/\mathrm{AF}[x]_2\f$. */
template<class T>
AffineMain<T> operator/(const AffineMain<T>&  x1, const AffineMain<T>&  x2);

/** \brief Returns \f$\mathrm{AF}[x]/d\f$. */
template<class T>
AffineMain<T> operator/(const AffineMain<T>&  x, double d);

/** \brief Returns \f$d/\mathrm{AF}[x]\f$. */
template<class T>
AffineMain<T> operator/(double d, const AffineMain<T>&  x);

/** \brief Returns \f$\mathrm{AF}[x]_1/[x]_2\f$. */
template<class T>
AffineMain<T> operator/(const AffineMain<T>&  x1, const Interval& x2);

/** \brief Returns \f$[x]_1/\mathrm{AF}[x]_2\f$. */
template<class T>
AffineMain<T> operator/(const Interval& x1, const AffineMain<T>&  x2);


/** \brief Returns the Hausdorff distance between two affine forms. */
template<class T>
double distance(const AffineMain<T>& x1, const AffineMain<T>& x2);

/** \brief Returns the Hausdorff distance between an interval and an affine form. */
template<class T>
double distance(const Interval &x1, const AffineMain<T>& x2);

/** \brief Returns the Hausdorff distance between an affine form and an interval. */
template<class T>
double distance(const AffineMain<T>& x1, const Interval &x2);


/** \brief Returns \f$1/\mathrm{AF}[x]\f$. */
template<class T>
AffineMain<T> inv(const AffineMain<T>&  x);

/** \brief Returns \f$\mathrm{AF}[x]^2\f$. */
template<class T>
AffineMain<T> sqr(const AffineMain<T>&  x);

/** \brief Returns \f$\sqrt{\mathrm{AF}[x]}\f$. */
template<class T>
AffineMain<T> sqrt(const AffineMain<T>&  x);

/** \brief Returns \f$\exp(\mathrm{AF}[x])\f$. */
template<class T>
AffineMain<T> exp(const AffineMain<T>&  x);

/** \brief Returns \f$\log(\mathrm{AF}[x])\f$. */
template<class T>
AffineMain<T> log(const AffineMain<T>&  x);

/** \brief Returns \f$\mathrm{AF}[x]^n\f$. */
template<class T>
AffineMain<T> pow(const AffineMain<T>&  x, int n);

/** \brief Returns \f$\mathrm{AF}[x]^d\f$. */
template<class T>
AffineMain<T> pow(const AffineMain<T>&  x, double d);

/** \brief Returns \f$\mathrm{AF}[x]^{[y]}\f$. */
template<class T>
AffineMain<T> pow(const AffineMain<T>& x, const Interval &y);

/** \brief Returns \f$\mathrm{AF}[x]^{\mathrm{AF}[y]}\f$. */
template<class T>
AffineMain<T> pow(const AffineMain<T>& x, const AffineMain<T>& y);

/** \brief Returns the n-th root of \f$\mathrm{AF}[x]\f$. */
template<class T>
AffineMain<T> root(const AffineMain<T>&  x, int n);

/** \brief Returns \f$\cos(\mathrm{AF}[x])\f$. */
template<class T>
AffineMain<T> cos(const AffineMain<T>&  x);

/** \brief Returns \f$\sin(\mathrm{AF}[x])\f$. */
template<class T>
AffineMain<T> sin(const AffineMain<T>&  x);

/** \brief Returns \f$\tan(\mathrm{AF}[x])\f$. */
template<class T>
AffineMain<T> tan(const AffineMain<T>&  x);

/** \brief Returns \f$\arccos(\mathrm{AF}[x])\f$. */
template<class T>
AffineMain<T> acos(const AffineMain<T>&  x);

/** \brief Returns \f$\arcsin(\mathrm{AF}[x])\f$. */
template<class T>
AffineMain<T> asin(const AffineMain<T>&  x);

/** \brief Returns \f$\arctan(\mathrm{AF}[x])\f$. */
template<class T>
AffineMain<T> atan(const AffineMain<T>&  x);

/** \brief Returns \f$\cosh(\mathrm{AF}[x])\f$. */
template<class T>
AffineMain<T> cosh(const AffineMain<T>&  x);

/** \brief Returns \f$\sinh(\mathrm{AF}[x])\f$. */
template<class T>
AffineMain<T> sinh(const AffineMain<T>&  x);

/** \brief Returns \f$\tanh(\mathrm{AF}[x])\f$. */
template<class T>
AffineMain<T> tanh(const AffineMain<T>&  x);

/** \brief Returns \f$|\mathrm{AF}[x]|\f$. */
template<class T>
AffineMain<T> abs(const AffineMain<T>& x);

template<class T>
Interval max(const AffineMain<T>&  x, const AffineMain<T>&  y);
template<class T>
Interval max(const Interval& x, const AffineMain<T>&  y);
template<class T>
Interval max(const AffineMain<T>&  x, const Interval& y);

template<class T>
Interval min(const AffineMain<T>&  x, const AffineMain<T>&  y);
template<class T>
Interval min(const Interval& x, const AffineMain<T>&  y) ;
template<class T>
Interval min(const AffineMain<T>&  x, const Interval& y);


/**
 * \brief Returns \f$[x]_1\cap[x]_2\f$.
 *
 * \return \c Interval::EMPTY if the intersection is empty
 */
template<class T>
Interval operator&(const AffineMain<T>&  x1, const AffineMain<T>&  x2);
template<class T>
Interval operator&(const Interval& x1, const AffineMain<T>&  x2);
template<class T>
Interval operator&(const AffineMain<T>&  x1, const Interval& x2);


/** \brief Returns \f$\square([x]_1\cup[x]_2)\f$. */
template<class T>
Interval operator|(const AffineMain<T>&  x1, const AffineMain<T>&  x2);
template<class T>
Interval operator|(const Interval& x1, const AffineMain<T>&  x2);
template<class T>
Interval operator|(const AffineMain<T>&  x1, const Interval& x2);


/** \brief Returns the largest integer interval included in \p x. */
template<class T>
AffineMain<T> integer(const AffineMain<T>&  x);
/**
 * \brief Returns \f$\lfloor[x]\rfloor\f$.
 *
 * \param x affine form value
 * \return operation result
 */
template<class T>
AffineMain<T> floor(const AffineMain<T>& x);

/**
 * \brief Returns \f$\lceil[x]\rceil\f$.
 *
 * \param x affine form value
 * \return operation result
 */
template<class T>
AffineMain<T> ceil(const AffineMain<T>& x);

/** \brief Sign of AF[x].
 *
 *  Return \f$sign(AF[x]) = hull \{ sign{x}, x\inAF[x] \}\f$.
 * \remark By convention, \f$ 0\inAF[x] \Longrightarrow sign(AF[x])=[-1,1]\f$. */
template<class T>
AffineMain<T> sign(const AffineMain<T>&  x);



/** \brief Chi of AF[a], AF[b] and AF[c].
 *
 *  Return \f$chi(AF[a],AF[b],AF[c]) = 0.5*(1-\sign(AF[a]))*AF[b] + 0.5*(\sign(AF[a])+1)*AF[c]. \f$
 * \remark  chi(AF[a],AF[b],AF[c]) =AF[b] if [a]<=0, AF[c] if AF[a]>0, hull \{AF[b], AF[c]\} else.  */
template<class T>
AffineMain<T> chi(const AffineMain<T>&  a,const AffineMain<T>&  b,const AffineMain<T>&  c);
template<class T>
AffineMain<T> chi(const Interval&  a,const AffineMain<T>&  b,const AffineMain<T>&  c);
template<class T>
AffineMain<T> chi(const Interval&  a,const Interval&  b,const AffineMain<T>&  c);
template<class T>
AffineMain<T> chi(const Interval&  a,const AffineMain<T>&  b,const Interval&  c);
template<class T>
AffineMain<T> chi(const AffineMain<T>&  a,const Interval&  b,const AffineMain<T>&  c);
template<class T>
AffineMain<T> chi(const AffineMain<T>&  a,const AffineMain<T>&  b,const Interval&  c);








/*@}*/

/*============================================ inline implementation ============================================ */




template<class T> const double AffineMain<T>::AF_COMPAC_Tol = 1.e-6;
template<class T> const double AffineMain<T>::AF_EM = __builtin_powi(2.0, -51);
template<class T> const double AffineMain<T>::AF_EC = __builtin_powi(2.0, -55);
template<class T> const double AffineMain<T>::AF_EE = 2.0;
template<class T> bool AffineMain<T>::mode=true;



template<class T>
inline void AffineMain<T>::change_mode(Affine_Mode tt) {
	switch(tt) {
	case AF_Lin_Default:
	case AF_Lin_Chebyshev: {
		mode =true;
		break;
	}
	case AF_Lin_MinRange:
		mode =false;
		break;
	}
}
template<class T>
inline AffineMain<T>::Affine_Mode AffineMain<T>::get_mode() {
	if (mode==true) {
		return AF_Lin_Chebyshev;
	} else {
		return AF_Lin_MinRange;
	}
}

template<class T>
inline void AffineMain<T>::compact(){	compact(AF_COMPAC_Tol); }

template<class T>
inline AffineMain<T>& AffineMain<T>::operator=(double d) {
	*this = Interval(d);
	return *this;
}



template<class T>
inline Interval operator&(const AffineMain<T>& x1, const AffineMain<T>& x2) {
	return ((x1.itv())&(x2.itv()));
}

template<class T>
inline Interval operator&(const Interval& x1, const AffineMain<T>& x2) {
	return ((x1)&(x2.itv()));
}

template<class T>
inline Interval operator&(const AffineMain<T>& x1, const Interval& x2) {
	return ((x1.itv())&(x2));
}

template<class T>
inline Interval operator|(const AffineMain<T>& x1, const AffineMain<T>& x2) {
	return ((x1.itv())|(x2.itv()));
}

template<class T>
inline Interval operator|(const Interval& x1, const AffineMain<T>& x2) {
	return ((x1)|(x2.itv()));
}

template<class T>
inline Interval operator|(const AffineMain<T>& x1, const Interval& x2) {
	return ((x1.itv())|(x2));
}

template<class T>
inline bool AffineMain<T>::operator==(const AffineMain<T>& x) const{
	return (this->itv() == x.itv());
}

template<class T>
inline bool AffineMain<T>::operator==(const Interval& x) const{
	return (this->itv() == x);
}

template<class T>
inline bool AffineMain<T>::operator!=(const AffineMain<T>& x) const{
	return (this->itv() != x.itv());
}

template<class T>
inline bool AffineMain<T>::operator!=(const Interval& x) const {
	return (this->itv() != x);
}

template<class T>
inline void AffineMain<T>::set_empty(){
	*this = Interval::empty();
}

//template<class T>
//inline AffineMain<T>& AffineMain<T>::inflate(double radd){
//	if (fabs(radd)>= POS_INFINITY) {
//		*this = Interval::all_reals();
//	} else {
//		saxpy(1.0, AffineMain<T>(),0.0, radd, false, false, false, true);
//	}
//	return *this;
//}

template<class T>
inline int AffineMain<T>::size() const{
	return _n;
}


template<class T>
inline bool AffineMain<T>::is_actif() const{
	return (_actif>-1);
}

template<class T>
inline bool AffineMain<T>::is_empty() const{
	return (_actif==-1);
}

template<class T>
inline bool AffineMain<T>::is_degenerated() const {
	//return (itv().diam() <	AF_EC);
	return (_actif==0);
}

template<class T>
inline bool AffineMain<T>::is_unbounded() const{
	return ((-1>_actif)&&(_actif>-5));
}

template<class T>
inline AffineMain<T>& AffineMain<T>::operator-=(double d){
	*this += (-d);
	return *this;
}

template<class T>
inline 	AffineMain<T>& AffineMain<T>::operator/=(double d) {
	return *this *= (1.0 / Interval(d)) ;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::operator+=(const Interval& x){
	if ( x.is_empty()) {
		*this = Interval::empty();
	} else if (x.is_unbounded() || (!is_actif())) {
		*this = this->itv() + x;
	} else {
		*this += x.mid();
		this->inflate(x.rad());
	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::operator-=(const Interval& x){
	return *this += (-x);
}

template<class T>
inline AffineMain<T>& AffineMain<T>::operator/=(const Interval& x){
	return *this *= (1.0/x) ;
}

//template<class T>
//inline AffineMain<T>& AffineMain<T>::operator+=(const AffineMain<T>& x){
//	return saxpy(1.0, x, 0.0, 0.0, false, true, false, false);
//}

template<class T>
inline AffineMain<T>& AffineMain<T>::operator-=(const AffineMain<T>& x){
	return *this += (-AffineMain<T>(x));
}

template<class T>
inline AffineMain<T>& AffineMain<T>::operator/=(const AffineMain<T>& x){
	return *this *= (AffineMain<T>(x).Ainv(x.itv()));
}

template<class T>
inline AffineMain<T> operator+(const AffineMain<T>& x1, const AffineMain<T>& x2){
	return AffineMain<T>(x1) += x2;
}

template<class T>
inline AffineMain<T> operator+(const AffineMain<T>& x, double d){
	return AffineMain<T>(x) += d;
}

template<class T>
inline AffineMain<T> operator+(double d, const AffineMain<T>& x){
	return AffineMain<T>(x) += d;
}

template<class T>
inline AffineMain<T> operator+(const AffineMain<T>& x1, const Interval& x2){
	return AffineMain<T>(x1) += x2;
}

template<class T>
inline AffineMain<T> operator+(const Interval& x1, const AffineMain<T>& x2){
	return AffineMain<T>(x2) += x1;
}

template<class T>
inline AffineMain<T> operator-(const AffineMain<T>& x1, const AffineMain<T>& x2){
	return AffineMain<T>(x1) += (-x2);
}

template<class T>
inline AffineMain<T> operator-(const AffineMain<T>& x, double d){
	return AffineMain<T>(x) -= d;
}

template<class T>
inline AffineMain<T> operator-(double d, const AffineMain<T>& x){
	AffineMain<T> res = (-x);
	return res += d;
}

template<class T>
inline AffineMain<T> operator-(const AffineMain<T>& x1, const Interval& x2) {
	return AffineMain<T>(x1) -= x2;
}

template<class T>
inline AffineMain<T> operator-(const Interval& x1, const AffineMain<T>& x2) {
	AffineMain<T> res = (- x2);
	return res += x1;
}

template<class T>
inline AffineMain<T> operator*(const AffineMain<T>& x1, const AffineMain<T>& x2) {
	return AffineMain<T>(x1) *= x2;
}

template<class T>
inline AffineMain<T> operator*(const AffineMain<T>& x, double d){
	return AffineMain<T>(x) *= d;
}

template<class T>
inline AffineMain<T> operator*(double d, const AffineMain<T>& x){
	return AffineMain<T>(x) *= d;
}

template<class T>
inline AffineMain<T> operator*(const AffineMain<T>& x1, const Interval& x2){
	return AffineMain<T>(x1) *= x2;
}

template<class T>
inline AffineMain<T> operator*(const Interval& x1, const AffineMain<T>& x2){
	return AffineMain<T>(x2) *= x1;
}

template<class T>
inline AffineMain<T> operator/(const AffineMain<T>& x1, const AffineMain<T>& x2){
	return AffineMain<T>(x1) /= x2;
}

template<class T>
inline AffineMain<T> operator/(const AffineMain<T>& x, double d){
	return AffineMain<T>(x) /= d;
}

template<class T>
inline AffineMain<T> operator/(double d, const AffineMain<T>& x){
	AffineMain<T> out;
	out = d;
	return out *= (AffineMain<T>(x).Ainv(x.itv()));
}

template<class T>
inline AffineMain<T> operator/(const AffineMain<T>& x1, const Interval& x2){
	return AffineMain<T>(x1) /= x2;
}

template<class T>
inline AffineMain<T> operator/(const Interval& x1, const AffineMain<T>& x2){
	return AffineMain<T>(x1) *= (AffineMain<T>(x2).Ainv(x2.itv()));
}

template<class T>
inline double distance(const AffineMain<T> &x1, const AffineMain<T> &x2){
	return distance(x1.itv(), x2.itv());
}

template<class T>
inline double distance(const Interval &x1, const AffineMain<T> &x2){
	return distance(x1, x2.itv());
}

template<class T>
inline double distance(const AffineMain<T> &x1, const Interval &x2){
	return distance(x1.itv(), x2);
}

template<class T>
inline AffineMain<T> inv(const AffineMain<T>& x){
	return AffineMain<T>(x).Ainv(x.itv());
}

template<class T>
inline AffineMain<T> sqr(const AffineMain<T>& x){
	return AffineMain<T>(x).Asqr(x.itv());
}

template<class T>
inline AffineMain<T> sqrt(const AffineMain<T>& x){
	return AffineMain<T>(x).Asqrt(x.itv());
}

template<class T>
inline AffineMain<T> exp(const AffineMain<T>& x){
	return AffineMain<T>(x).Aexp(x.itv());
}

template<class T>
inline AffineMain<T> log(const AffineMain<T>& x){
	return AffineMain<T>(x).Alog(x.itv());
}

template<class T>
inline AffineMain<T> pow(const AffineMain<T>& x, int n) {
	return AffineMain<T>(x).Apow(n,x.itv());
}
template<class T>
inline AffineMain<T> pow(const AffineMain<T>& x, double d){
	return AffineMain<T>(x).Apow(d,x.itv());
}
template<class T>
inline AffineMain<T> pow(const AffineMain<T> &x, const Interval &y){
	// return exp(y * log(x));
	return AffineMain<T>(x).Apow(y,x.itv());
}
template<class T>
inline AffineMain<T> pow(const AffineMain<T> &x, const AffineMain<T> &y){
	return AffineMain<T>(x).Apow(y.itv(),x.itv());
}

template<class T>
inline AffineMain<T> root(const AffineMain<T>& x, int n) {
	return AffineMain<T>(x).Aroot(n,x.itv());
}

template<class T>
inline AffineMain<T> cos(const AffineMain<T>& x){
	return AffineMain<T>(x).Acos(x.itv());
}

template<class T>
inline AffineMain<T> sin(const AffineMain<T>& x){
	return AffineMain<T>(x).Asin(x.itv());
}

template<class T>
inline AffineMain<T> tan(const AffineMain<T>& x){
	return AffineMain<T>(x).Atan(x.itv());
}

template<class T>
inline AffineMain<T> acos(const AffineMain<T>& x){
	return AffineMain<T>(x).Aacos(x.itv());
}

template<class T>
inline AffineMain<T> asin(const AffineMain<T>& x){
	return AffineMain<T>(x).Aasin(x.itv());
}

template<class T>
inline AffineMain<T> atan(const AffineMain<T>& x){
	return AffineMain<T>(x).Aatan(x.itv());
}

template<class T>
inline AffineMain<T> cosh(const AffineMain<T>& x){
	return AffineMain<T>(x).Acosh(x.itv());
}

template<class T>
inline AffineMain<T> sinh(const AffineMain<T>& x){
	return AffineMain<T>(x).Asinh(x.itv());
}

template<class T>
inline AffineMain<T> tanh(const AffineMain<T>& x){
	return AffineMain<T>(x).Atanh(x.itv());
}

template<class T>
inline AffineMain<T> abs(const AffineMain<T> &x){
	return AffineMain<T>(x).Aabs(x.itv());
}


template<class T>
inline Interval max(const AffineMain<T>& x, const AffineMain<T>& y) {
	return max(x.itv(), y.itv());
}
template<class T>
inline Interval max(const Interval& x, const AffineMain<T>& y) {
	return max(x, y.itv());
}
template<class T>
inline Interval max(const AffineMain<T>& x, const Interval& y) {
	return max(x.itv(), y);
}


template<class T>
inline Interval min(const AffineMain<T>& x, const AffineMain<T>& y) {
	return min(x.itv(), y.itv());
}
template<class T>
inline Interval min(const Interval& x, const AffineMain<T>& y) {
	return min(x, y.itv());
}
template<class T>
inline Interval min(const AffineMain<T>& x, const Interval& y) {
	return min(x.itv(), y);
}


template<class T>
inline AffineMain<T> integer(const AffineMain<T>& x){
	Interval tmp= x.itv();
	if ((tmp.lb()==(int) tmp.lb())&&(tmp.ub()==(int) tmp.ub())) {
		return AffineMain<T>(x);
	} else  {
		AffineMain<T> out;
		out =integer(x.itv());
		return out;
	}
}
template<class T>
inline AffineMain<T> floor(const AffineMain<T>& x){
	Interval tmp= x.itv();
	if ((tmp.lb()==(int) tmp.lb())&&(tmp.ub()==(int) tmp.ub())) {
		return AffineMain<T>(x);
	} else  {
		AffineMain<T> out;
		out =floor(x.itv());
		return out;
	}
}
template<class T>
inline AffineMain<T> ceil(const AffineMain<T>& x){
	Interval tmp= x.itv();
	if ((tmp.lb()==(int) tmp.lb())&&(tmp.ub()==(int) tmp.ub())) {
		return AffineMain<T>(x);
	} else  {
		AffineMain<T> out;
		out =ceil(x.itv());
		return out;
	}
}

template<class T>
inline AffineMain<T> sign(const AffineMain<T>& x) {
	AffineMain<T> out;
	out =sign(x.itv());
	return out;
}


template<class T>
inline AffineMain<T> chi(const AffineMain<T>&  a,const AffineMain<T>&  b,const AffineMain<T>&  c) {
	return chi(a.itv(),b,c);
}
template<class T>
inline AffineMain<T> chi(const Interval&  a,const Interval&  b,const AffineMain<T>&  c) {
	AffineMain<T> out(c);
	out = b;
	return chi(a,out,c);
}
template<class T>
inline AffineMain<T> chi(const Interval&  a,const AffineMain<T>&  b,const Interval&  c) {
	AffineMain<T> out(b);
	out = c;
	return chi(a,b,out);
}
template<class T>
inline AffineMain<T> chi(const AffineMain<T>&  a,const Interval&  b,const AffineMain<T>&  c) {
	AffineMain<T> out(c);
	out = b;
	return chi(a.itv(),out,c);
}
template<class T>
inline AffineMain<T> chi(const AffineMain<T>&  a,const AffineMain<T>&  b,const Interval&  c){
	AffineMain<T> out(b);
	out = c;
	return chi(a.itv(),b,out);
}
template<class T>
inline AffineMain<T> chi(const Interval&  a,const AffineMain<T>&  b,const AffineMain<T>&  c){
	if (a.ub()<=0) {
		return AffineMain<T>(b);
	} else if (a.lb()>0) {
		return AffineMain<T>(c);
	} else {
		Interval b_itv = b.itv();
		Interval c_itv = c.itv();
		if (b_itv.is_subset(c_itv)) {
			return AffineMain<T>(c);
		} else if (c_itv.is_subset(b_itv)) {
			return AffineMain<T>(b);
		} else {
			AffineMain<T> out;
			out = b|c;
			return  out;
		}
	}
}

template<class T>
inline std::ostream& operator<<(std::ostream& os, const AffineMain<T>& x) {
	{
		os << x.itv() << " : ";
		if (x.is_actif()) {
			os << x.mid();
			for (int i = 0; i < x.size(); i++) {
				os << " + " << x.val(i) << " eps_" << i;
			}
			os << " + " << x.err() << " [-1,1] ";
		} else {
			os << "Affine Form is not enabled. ";
		}
		return os;

	}
}


template<class T>
inline AffineMain<T>& AffineMain<T>::Ainv(const Interval& itv){
	if (mode) {
		this->Ainv_CH(itv);
	} else {
		this->Ainv_MR(itv);
	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Aexp(const Interval& itv){
	if (mode) {
		this->Aexp_CH(itv);
	} else {
		this->Aexp_MR(itv);
	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Alog(const Interval& itv){
	if (mode) {
		this->Alog_CH(itv);
	} else {
		this->Alog_MR(itv);
	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Asqrt(const Interval& itv){
	if (mode) {
		this->Asqrt_CH(itv);
	} else {
		this->Asqrt_MR(itv);
	}
	return *this;
}


// Chebyshev formula:
// alpha = (f(ub(x)) - f(lb(x))) / diam(x)
// u     = (f')^{-1}(alpha)
// d_a   = f(lb(x)) - alpha*lb(x)
// d_b   = f(ub(x)) - alpha*ub(x)
// d_min = min(d_a, d_b)
// d_max = f(u) - alpha*u
// beta  = mid([d_min, d_max])
// zeta  = rad([d_min, d_max])

template<class T>
inline AffineMain<T>& AffineMain<T>::Ainv_CH(const Interval& itv){
	Interval res_itv = 1.0/(itv);;

	// Particular case
	if ((itv.is_unbounded()) || res_itv.is_empty() || res_itv.is_unbounded() || (!is_actif()) || (itv.diam()<AF_EC)) {
		*this = res_itv;
	}  else  {
		// General case
		double alpha, beta, ddelta, t1, t2;
		Interval dmm(0.0), TEMP1(0.0), TEMP2(0.0), band(0.0);
		Interval itv2 =abs(itv);
		if (itv2.diam()<AF_EC) {
			*this = res_itv;
			return *this;
		}
		dmm = (1.0/(itv2));
		alpha = -(dmm.diam()/itv2.diam());

		//u = 1/sqrt(-alpha);
		TEMP1 = (1.0/Interval((itv2).lb()))-alpha*Interval((itv2).lb());
		TEMP2 = (1.0/Interval((itv2).ub()))-alpha*Interval((itv2).ub());
		if (TEMP1.ub()>TEMP2.ub()) {
			band = Interval((2*sqrt(-Interval(alpha))).lb(),TEMP1.ub());
		}
		else {
			band = Interval((2*sqrt(-Interval(alpha))).lb(),TEMP2.ub());
		}

		beta = band.mid();
		t1 = (beta -band).ub();
		t2 = (band -beta).ub();
		ddelta = (t1>t2)? t1 : t2;

		if (itv.lb()<0.0) beta = -beta;

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);

	}


	return *this;
}


template<class T>
inline AffineMain<T>& AffineMain<T>::Asqrt_CH(const Interval& itv){
	Interval itv2 = itv & Interval(0,oo);
	Interval res_itv = sqrt(itv2);;

	// Particular case
	if (res_itv.is_empty() || res_itv.is_unbounded() || (!is_actif()) || (itv2.diam()<AF_EC)) {
		*this = res_itv;
	}  else  {
		// General case
		double alpha, beta, ddelta, t1, t2;
		Interval TEMP1(0.0), TEMP2(0.0), band(0.0);

		alpha = res_itv.diam()/itv2.diam();

		//u = 1/(4*alpha^2);
		TEMP1 = res_itv.lb()-alpha*Interval(itv2.lb());
		TEMP2 = res_itv.ub()-alpha*Interval(itv2.ub());
		if (TEMP1.lb()>TEMP2.lb()) {
			band = Interval(TEMP2.lb(),(1.0/(4*Interval(alpha))).ub());
		} else {
			band = Interval(TEMP1.lb(),(1.0/(4*Interval(alpha))).ub());
		}

		beta = band.mid();
		t1 = (beta -band).ub();
		t2 = (band-beta).ub();
		ddelta = (t1>t2)? t1 : t2;

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);


	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Aexp_CH(const Interval& itv){
	Interval res_itv = exp(itv);;

	// Particular case
	if (res_itv.is_empty() || res_itv.is_unbounded() || (!is_actif()) || (itv.diam()<AF_EC)) {
		*this = res_itv;
	}  else  {
		// General case
		double alpha, beta, ddelta, t1, t2;
		Interval TEMP1(0.0), TEMP2(0.0), band(0.0);
		alpha = res_itv.diam()/itv.diam();

		//u = log(alpha);
		TEMP1 = res_itv.lb()-alpha*Interval(itv.lb());
		TEMP2 = res_itv.ub()-alpha*Interval(itv.ub());
		if (TEMP1.ub()>TEMP2.ub()) {
			band = Interval((alpha*(1-log(Interval(alpha)))).lb(),TEMP1.ub());
		} else {
			band = Interval((alpha*(1-log(Interval(alpha)))).lb(),TEMP2.ub());
		}

		beta = band.mid();
		t1 = (beta -band).ub();
		t2 = (band -beta).ub();
		ddelta = (t1>t2)? t1 : t2;

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);

	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Alog_CH(const Interval& itv){
	Interval res_itv = log(itv);;

	// Particular case
	if (res_itv.is_empty() || res_itv.is_unbounded() || (!is_actif()) || (itv.diam()<AF_EC)) {
		*this = res_itv;
	}  else  {
		// General case
		double alpha, beta, ddelta, t1, t2;
		Interval TEMP1(0.0), TEMP2(0.0), band(0.0);
		alpha = res_itv.diam()/itv.diam();

		//u = 1/alpha;
		TEMP1 = res_itv.lb()-alpha*Interval(itv.lb());
		TEMP2 = res_itv.ub()-alpha*Interval(itv.ub());
		if (TEMP1.lb()>TEMP2.lb()) {
			band = Interval(TEMP2.lb(),(-log(Interval(alpha))-1).ub());
		}
		else {
			band = Interval(TEMP1.lb(),(-log(Interval(alpha))-1).ub());
		}

		beta = band.mid();
		t1 = (beta -band).ub();
		t2 = (band -beta).ub();
		ddelta = (t1>t2)? t1 : t2;

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);


	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Apow(const Interval &y, const Interval& itvx){
	// return exp(y * log(x));
	this->Alog(itvx);
	*this *= y;
	this->Aexp((y*log(itvx)));
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Acos(const Interval& itv){
	Interval res_itv = cos(itv);;

	// Particular case
	if (res_itv.is_empty() || res_itv.is_unbounded() || (!is_actif()) || (itv.diam()<AF_EC)) {
		*this = res_itv;
	}  else  {
		// General case
		double alpha, beta, ddelta, t1 ;//, t2;
		Interval dmm(0.0), TEMP1(0.0), TEMP2(0.0), band(0.0);
		if (itv.diam()>=Interval::two_pi().lb()) {
			*this = Interval(-1,1);
			return *this;
		}
		//  pour _itv = [a,b]
		// x0 = 1/sqrt(2)
		// x1= - x0
		// xb0 = 0.5*((b-a)*x0 +(a+b))
		// xb1 = 0.5*((b-a)*x1 +(a+b))
		// c0 = 0.5 (f(xb0)+f(xb1))
		// c1 = x0*f(xb0)+x1*f(xb1)
		// alpha = 2*c1/(b-a)
		// beta = c0-c1*(a+b)/(b-a)
		//  old : ddelta = (b-a)^2 * f''(_itv)/16
		//  new : ddelta = evaluate the error at the bound and the points when f'(x)=alpha

		double x0,xb0,xb1,fxb0,fxb1,c0,c1;

		x0 = 1.0/std::sqrt(2.);
		xb0 = (0.5)*(itv.diam()*x0 +itv.lb()+itv.ub());
		xb1 = (0.5)*(itv.diam()*(-x0) +itv.lb()+itv.ub());


		fxb0 = std::cos(xb0);
		fxb1 = std::cos(xb1);

		c0 = (0.5)*(fxb0+fxb1);
		c1 = x0*fxb0-x0*fxb1;

		alpha  = 2*c1/(itv.diam());
		beta   = c0-c1*((itv.lb()+itv.ub())/(itv.diam()));
		//ddelta = ((size()*Interval(TEMP1.rad())) + Interval(TEMP2.rad())).ub();

		// compute the maximal error
		ddelta= 0.0;
		Interval u,nb_period;

		// compute the error at _itv.lb() and _itv.ub() and compute the first point such as f'(u) = alpha

		ddelta = (abs(cos(Interval(itv.lb()))-(alpha*Interval(itv.lb())+beta))).ub();
		t1     = (abs(cos(Interval(itv.ub()))-(alpha*Interval(itv.ub())+beta))).ub();
		if (t1>ddelta)  ddelta= t1;
		u = asin(-Interval(alpha));
		nb_period = (itv+Interval::half_pi()) / Interval::pi();

		// evaluate the error at the point such as f'(u) = alpha
		int p1 = ((int) nb_period.lb())-2;
		int p2 = ((int) nb_period.ub())+2;

		int i = p1;

		while (i<=p2) { // looking for a point
			TEMP1 = (itv & (i%2==0? (u + i*Interval::pi()) : (i*Interval::pi() - u)));
			if (!(TEMP1.is_empty())) { // check if maximize the error
				t1 = (abs(cos(TEMP1)-(alpha*TEMP1+beta))).ub();
				if (t1>ddelta)  ddelta= t1;
			}
			i++;
		}

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);


	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Asin(const Interval& itv){
	Interval res_itv = sin(itv);;

	// Particular case
	if (res_itv.is_empty() || res_itv.is_unbounded() || (!is_actif()) || (itv.diam()<AF_EC)) {
		*this = res_itv;
	}  else  {
		// General case
		double alpha, beta, ddelta, t1;//, t2;
		Interval dmm(0.0), TEMP1(0.0), TEMP2(0.0), band(0.0);
		if (itv.diam()>=Interval::two_pi().lb()) {
			*this = Interval(-1,1);
			return *this;
		}
		//  pour _itv = [a,b]
		// x0 = 1/sqrt(2)
		// x1= - x0
		// xb0 = 0.5*((b-a)*x0 +(a+b))
		// xb1 = 0.5*((b-a)*x1 +(a+b))
		// c0 = 0.5 (f(xb0)+f(xb1))
		// c1 = x0*f(xb0)+x1*f(xb1)
		// alpha = 2*c1/(b-a)
		// beta = c0-c1*(a+b)/(b-a)
		//  old : ddelta = (b-a)^2 * f''(_itv)/16
		//  new : ddelta = evaluate the error at the bound and the points when f'(x)=alpha

		double x0,xb0,xb1,fxb0,fxb1,c0,c1;

		x0 = 1.0/std::sqrt(2.);
		xb0 = (0.5)*(itv.diam()*x0 +itv.lb()+itv.ub());
		xb1 = (0.5)*(itv.diam()*(-x0) +itv.lb()+itv.ub());


		fxb0 = std::sin(xb0);
		fxb1 = std::sin(xb1);

		c0 = (0.5)*(fxb0+fxb1);
		c1 = x0*fxb0-x0*fxb1;

		alpha  = 2*c1/(itv.diam());
		beta   = c0-c1*((itv.lb()+itv.ub())/(itv.diam()));
		//ddelta = ((size()*Interval(TEMP1.rad())) + Interval(TEMP2.rad())).ub();

		// compute the maximal error
		ddelta= 0.0;
		Interval u,nb_period;

		// compute the error at _itv.lb() and _itv.ub() and compute the first point such as f'(u) = alpha

		ddelta = (abs(sin(Interval(itv.lb()))-(alpha*Interval(itv.lb())+beta))).ub();
		t1     = (abs(sin(Interval(itv.ub()))-(alpha*Interval(itv.ub())+beta))).ub();
		if (t1>ddelta)  ddelta= t1;
		u = acos(Interval(alpha));
		nb_period = (itv) / Interval::pi();

		// evaluate the error at the point such as f'(u) = alpha
		int p1 = ((int) nb_period.lb())-2;
		int p2 = ((int) nb_period.ub())+2;

		int i = p1;

		while (i<=p2) { // looking for a point
			TEMP1 = (itv & (i%2==0? (u + i*Interval::pi()) : ((i+1)*Interval::pi() - u)));
			if (!(TEMP1.is_empty())) {
				t1 = (abs(sin(TEMP1)-(alpha*TEMP1+beta))).ub();
				if (t1>ddelta)  ddelta= t1;
			}
			i++;
		}

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);


	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Atan(const Interval& itv){
	Interval res_itv = tan(itv);;

	// Particular case
	if (res_itv.is_empty() || res_itv.is_unbounded() || (!is_actif()) || (itv.diam()<AF_EC)) {
		*this = res_itv;
	}  else  {
		// General case
		double alpha, beta, ddelta, t1;//, t2;
		Interval dmm(0.0), TEMP1(0.0), TEMP2(0.0), band(0.0);
		if (itv.diam()>=Interval::two_pi().lb()) {
			*this = Interval(-1,1);
			return *this;
		}
		//  pour _itv = [a,b]
		// x0 = 1/sqrt(2)
		// x1= - x0
		// xb0 = 0.5*((b-a)*x0 +(a+b))
		// xb1 = 0.5*((b-a)*x1 +(a+b))
		// c0 = 0.5 (f(xb0)+f(xb1))
		// c1 = x0*f(xb0)+x1*f(xb1)
		// alpha = 2*c1/(b-a)
		// beta = c0-c1*(a+b)/(b-a)
		//  old : ddelta = (b-a)^2 * f''(_itv)/16
		//  new : ddelta = evaluate the error at the bound and the points when f'(x)=alpha

		double x0,xb0,xb1,fxb0,fxb1,c0,c1;

		x0 = 1.0/std::sqrt(2.);
		xb0 = (0.5)*(itv.diam()*  x0  +itv.lb()+itv.ub());
		xb1 = (0.5)*(itv.diam()*(-x0) +itv.lb()+itv.ub());

		fxb0 = std::tan(xb0);
		fxb1 = std::tan(xb1);

		c0 = (0.5)*(fxb0+fxb1);
		c1 = x0*fxb0-x0*fxb1;

		alpha  = 2*c1/(itv.diam());
		beta   = c0-c1*((itv.lb()+itv.ub())/(itv.diam()));
		//ddelta = ((size()*Interval(TEMP1.rad())) + Interval(TEMP2.rad())).ub();

		// compute the maximal error
		ddelta= 0.0;
		Interval u,nb_period;

		// compute the error at _itv.lb() and _itv.ub() and compute the first point such as f'(u) = alpha

		ddelta = (abs(tan(Interval(itv.lb()))-(alpha*Interval(itv.lb())+beta))).ub();
		t1     = (abs(tan(Interval(itv.ub()))-(alpha*Interval(itv.ub())+beta))).ub();
		if (t1>ddelta)  ddelta= t1;
		u = acos(1/sqrt(Interval(alpha)));
		nb_period = (itv) / Interval::pi();

		// evaluate the error at the point such as f'(u) = alpha
		int p1 = ((int) nb_period.lb())-2;
		int p2 = ((int) nb_period.ub())+2;

		int i = p1;

		while (i<=p2) { // looking for a point
			TEMP1 = (itv & ( i*Interval::pi() + u));
			if ((!(TEMP1.is_empty()))) {
				t1 = (abs(tan(TEMP1)-(alpha*TEMP1+beta))).ub();
				if (t1>ddelta)  ddelta= t1;
			}
			TEMP1 = (itv & ( i*Interval::pi() - u ));
			if ((!(TEMP1.is_empty()))) {
				t1 = (abs(tan(TEMP1)-(alpha*TEMP1+beta))).ub();
				if (t1>ddelta)  ddelta= t1;
			}
			i++;
		}

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);

	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Aacos(const Interval& itv){
	Interval res_itv = acos(itv);;
	Interval itv2 = itv & Interval(-1,1);
	// Particular case
	if ( res_itv.is_empty() || res_itv.is_unbounded() || (!is_actif()) || (itv2.diam()<AF_EC)) {
		*this = res_itv;
	}  else  {
		// General case
		double alpha, beta, ddelta, t1;//, t2;

		//  pour _itv = [a,b]
		// x0 = 1/sqrt(2)
		// x1= - x0
		// xb0 = 0.5*((b-a)*x0 +(a+b))
		// xb1 = 0.5*((b-a)*x1 +(a+b))
		// c0 = 0.5 (f(xb0)+f(xb1))
		// c1 = x0*f(xb0)+x1*f(xb1)
		// alpha = 2*c1/(b-a)
		// beta = c0-c1*(a+b)/(b-a)
		//  old : ddelta = (b-a)^2 * f''(_itv)/16
		//  new : ddelta = evaluate the error at the bound and the points when f'(x)=alpha

		double x0,xb0,xb1,fxb0,fxb1,c0,c1;

		x0 = 1.0/std::sqrt(2.);
		xb0 = (0.5)*(itv2.diam()*  x0  +itv2.lb()+itv2.ub());
		xb1 = (0.5)*(itv2.diam()*(-x0) +itv2.lb()+itv2.ub());

		fxb0 = std::acos(xb0);
		fxb1 = std::acos(xb1);

		c0 = (0.5)*(fxb0+fxb1);
		c1 = x0*fxb0-x0*fxb1;

		alpha  = 2*c1/(itv2.diam());
		beta   = c0-c1*((itv2.lb()+itv2.ub())/(itv2.diam()));
		//ddelta = ((_n*Interval(TEMP1.rad())) + Interval(TEMP2.rad())).ub();

		// compute the maximal error
		ddelta= 0.0;

		// compute the error at _itv.lb(), _itv.ub() and u such as f'(u) =alpha

		ddelta = (abs(acos(Interval(itv2.lb()))-(alpha*Interval(itv2.lb())+beta))).ub();
		t1     = (abs(acos(Interval(itv2.ub()))-(alpha*Interval(itv2.ub())+beta))).ub();
		if (t1>ddelta)  ddelta= t1;
		// acos'(u)=-1/sqrt(1-u^2) = alpha
		// u = +-sqrt(1-1/(alpha^2))
		Interval TEMP2 = sqrt(1-1/(pow(Interval(alpha),2)));
		if (!((TEMP2 & itv2).is_empty())) {
			t1 = (abs(acos(TEMP2)-(alpha*TEMP2+beta))).ub();
			if (t1>ddelta)  ddelta= t1;
		}
		if (!(((-TEMP2) & itv2).is_empty())) {
			t1 = (abs(acos(-TEMP2)-(alpha*(-TEMP2)+beta))).ub();
			if (t1>ddelta)  ddelta= t1;
		}

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);

	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Aasin(const Interval& itv){
	Interval res_itv = asin(itv);;
	Interval itv2 =itv & Interval(-1,1);
	// Particular case
	if (res_itv.is_empty() || res_itv.is_unbounded() || (!is_actif()) || (itv2.diam()<AF_EC)) {
		*this = res_itv;
	}  else  {
		// General case
		double alpha, beta, ddelta, t1;//, t2;

		//  pour _itv = [a,b]
		// x0 = 1/sqrt(2)
		// x1= - x0
		// xb0 = 0.5*((b-a)*x0 +(a+b))
		// xb1 = 0.5*((b-a)*x1 +(a+b))
		// c0 = 0.5 (f(xb0)+f(xb1))
		// c1 = x0*f(xb0)+x1*f(xb1)
		// alpha = 2*c1/(b-a)
		// beta = c0-c1*(a+b)/(b-a)
		//  old : ddelta = (b-a)^2 * f''(_itv)/16
		//  new : ddelta = evaluate the error at the bound and the points when f'(x)=alpha

		double x0,xb0,xb1,fxb0,fxb1,c0,c1;

		x0 = 1.0/std::sqrt(2.);
		xb0 = (0.5)*(itv2.diam()*  x0  +itv2.lb()+itv2.ub());
		xb1 = (0.5)*(itv2.diam()*(-x0) +itv2.lb()+itv2.ub());

		fxb0 = std::asin(xb0);
		fxb1 = std::asin(xb1);

		c0 = (0.5)*(fxb0+fxb1);
		c1 = x0*fxb0-x0*fxb1;

		alpha  = 2*c1/(itv2.diam());
		beta   = c0-c1*((itv2.lb()+itv2.ub())/(itv2.diam()));
		//ddelta = ((_n*Interval(TEMP1.rad())) + Interval(TEMP2.rad())).ub();

		// compute the maximal error
		ddelta= 0.0;

		// compute the error at _itv.lb(), _itv.ub() and u such as f'(u) =alpha

		ddelta = (abs(asin(Interval(itv2.lb()))-(alpha*Interval(itv2.lb())+beta))).ub();
		t1     = (abs(asin(Interval(itv2.ub()))-(alpha*Interval(itv2.ub())+beta))).ub();
		if (t1>ddelta)  ddelta= t1;
		// asin'(u)=1/sqrt(1-u^2) = alpha
		// u = sqrt(1-1/(alpha^2))
		Interval TEMP2 = sqrt(1/Interval(alpha)-1);
		if (!((TEMP2 & itv2).is_empty())) {
			// tanh(acosh(x)) = sqrt(sqr(x)-1) :x
			t1 = (abs(asin(TEMP2)-(alpha*TEMP2+beta))).ub();
			if (t1>ddelta)  ddelta= t1;
		}
		if (!(((-TEMP2) & itv2).is_empty())) {
			// sinh(acosh(-x)) = -sqrt(sqr(x)-1)
			t1 = (abs(asin(-TEMP2)-(alpha*(-TEMP2)+beta))).ub();
			if (t1>ddelta)  ddelta= t1;
		}

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);

	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Aatan(const Interval& itv){
	Interval res_itv = atan(itv);;

	// Particular case
	if (itv.is_unbounded() || res_itv.is_empty() || res_itv.is_unbounded() || (!is_actif()) || (itv.diam()<AF_EC)) {
		*this = res_itv;
	}  else  {
		// General case
		double alpha, beta, ddelta, t1;//, t2;
		Interval  TEMP2(0.0);
		//  pour _itv = [a,b]
		// x0 = 1/sqrt(2)
		// x1= - x0
		// xb0 = 0.5*((b-a)*x0 +(a+b))
		// xb1 = 0.5*((b-a)*x1 +(a+b))
		// c0 = 0.5 (f(xb0)+f(xb1))
		// c1 = x0*f(xb0)+x1*f(xb1)
		// alpha = 2*c1/(b-a)
		// beta = c0-c1*(a+b)/(b-a)
		//  old : ddelta = (b-a)^2 * f''(_itv)/16
		//  new : ddelta = evaluate the error at the bound and the points when f'(x)=alpha

		double x0,xb0,xb1,fxb0,fxb1,c0,c1;

		x0 = 1.0/std::sqrt(2.);
		xb0 = (0.5)*(itv.diam()*  x0  +itv.lb()+itv.ub());
		xb1 = (0.5)*(itv.diam()*(-x0) +itv.lb()+itv.ub());

		fxb0 = std::atan(xb0);
		fxb1 = std::atan(xb1);

		c0 = (0.5)*(fxb0+fxb1);
		c1 = x0*fxb0-x0*fxb1;

		alpha  = 2*c1/(itv.diam());
		beta   = c0-c1*((itv.lb()+itv.ub())/(itv.diam()));
		//ddelta = ((_n*Interval(TEMP1.rad())) + Interval(TEMP2.rad())).ub();

		// compute the maximal error
		ddelta= 0.0;

		// compute the error at _itv.lb(), _itv.ub() and u such as f'(u) =alpha

		ddelta = (abs(atan(Interval(itv.lb()))-(alpha*Interval(itv.lb())+beta))).ub();
		t1     = (abs(atan(Interval(itv.ub()))-(alpha*Interval(itv.ub())+beta))).ub();
		if (t1>ddelta)  ddelta= t1;
		// atan'(u)=1/(u^2+1) = alpha
		// u = +-sqrt(1/alpha -1)
		TEMP2 = sqrt(1/Interval(alpha)-1);
		if (!((TEMP2 & itv).is_empty())) {
			t1 = (abs(atan(TEMP2)-(alpha*TEMP2+beta))).ub();
			if (t1>ddelta)  ddelta= t1;
		}
		if (!(((-TEMP2) & itv).is_empty())) {
			t1 = (abs(atan(-TEMP2)-(alpha*(-TEMP2)+beta))).ub();
			if (t1>ddelta)  ddelta= t1;
		}

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);

	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Acosh(const Interval& itv){
	Interval res_itv = cosh(itv);;

	// Particular case
	if (res_itv.is_empty() || res_itv.is_unbounded() || (!is_actif()) || (itv.diam()<AF_EC)) {
		*this = res_itv;
	}  else  {
		// General case
		double alpha, beta, ddelta, t1, t2;
		Interval TEMP1(0.0), TEMP2(0.0), band(0.0);
		alpha = ((cosh(Interval(itv.ub()))-cosh(Interval(itv.lb())))/itv.diam()).lb();

		//u = asinh(alpha);
		TEMP1 = cosh(Interval(itv.lb()))-alpha*Interval(itv.lb());
		TEMP2 = cosh(Interval(itv.ub()))-alpha*Interval(itv.ub());
		if (TEMP1.ub()>TEMP2.ub()) {
			// cosh(asinh(alpha)) = sqrt(sqr(alpha)+1)
			band = Interval((sqrt(pow(Interval(alpha),2)+1)-alpha*asinh(Interval(alpha))).lb(),TEMP1.ub());
		}
		else {
			band = Interval((sqrt(pow(Interval(alpha),2)+1)-alpha*asinh(Interval(alpha))).lb(),TEMP2.ub());
		}

		beta = band.mid();
		t1 = (beta -band).ub();
		t2 = (band -beta).ub();
		ddelta = (t1>t2)? t1 : t2;

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);


	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Asinh(const Interval& itv){
	Interval res_itv = sinh(itv);;

	// Particular case
	if (res_itv.is_empty() || res_itv.is_unbounded() || (!is_actif()) || (itv.diam()<AF_EC)) {
		*this = res_itv;
	}  else  {
		// General case
		double alpha, beta, ddelta, t1;//, t2;
		Interval  TEMP2(0.0);
		//  pour _itv = [a,b]
		// x0 = 1/sqrt(2)
		// x1= - x0
		// xb0 = 0.5*((b-a)*x0 +(a+b))
		// xb1 = 0.5*((b-a)*x1 +(a+b))
		// c0 = 0.5 (f(xb0)+f(xb1))
		// c1 = x0*f(xb0)+x1*f(xb1)
		// alpha = 2*c1/(b-a)
		// beta = c0-c1*(a+b)/(b-a)
		//  old : ddelta = (b-a)^2 * f''(_itv)/16
		//  new : ddelta = evaluate the error at the bound and the points when f'(x)=alpha

		double x0,xb0,xb1,fxb0,fxb1,c0,c1;

		x0 = 1.0/std::sqrt(2.);
		xb0 = (0.5)*(itv.diam()*  x0  +itv.lb()+itv.ub());
		xb1 = (0.5)*(itv.diam()*(-x0) +itv.lb()+itv.ub());

		fxb0 = std::sinh(xb0);
		fxb1 = std::sinh(xb1);

		c0 = (0.5)*(fxb0+fxb1);
		c1 = x0*fxb0-x0*fxb1;

		alpha  = 2*c1/(itv.diam());
		beta   = c0-c1*((itv.lb()+itv.ub())/(itv.diam()));
		//ddelta = ((_n*Interval(TEMP1.rad())) + Interval(TEMP2.rad())).ub();

		// compute the maximal error
		ddelta= 0.0;

		// compute the error at _itv.lb(), _itv.ub() and u such as f'(u) =alpha
		ddelta = (abs(sinh(Interval(itv.lb()))-(alpha*Interval(itv.lb())+beta))).ub();
		t1     = (abs(sinh(Interval(itv.ub()))-(alpha*Interval(itv.ub())+beta))).ub();
		if (t1>ddelta)  ddelta= t1;
		// u = acosh(alpha)
		TEMP2 = acosh(Interval(alpha));
		if (!((TEMP2 & itv).is_empty())) {
			// sinh(acosh(x)) = sqrt(sqr(x)-1)
			t1 = (abs(sqrt(pow(Interval(alpha),2)-1)-(alpha*TEMP2+beta))).ub();
			if (t1>ddelta)  ddelta= t1;
		}
		if (!(((-TEMP2) & itv).is_empty())) {
			// sinh(acosh(-x)) = -sqrt(sqr(x)-1)
			t1 = (abs((-sqrt(pow(Interval(alpha),2)-1))-(alpha*(-TEMP2)+beta))).ub();
			if (t1>ddelta)  ddelta= t1;
		}

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);

	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Atanh(const Interval& itv){
	Interval res_itv = tanh(itv);;

	// Particular case
	if (itv.is_unbounded() || res_itv.is_empty() || res_itv.is_unbounded() || (!is_actif()) || (itv.diam()<AF_EC)) {
		*this = res_itv;
	}  else  {
		// General case
		double alpha, beta, ddelta, t1;//, t2;
		Interval  TEMP2(0.0);
		//  pour _itv = [a,b]
		// x0 = 1/sqrt(2)
		// x1= - x0
		// xb0 = 0.5*((b-a)*x0 +(a+b))
		// xb1 = 0.5*((b-a)*x1 +(a+b))
		// c0 = 0.5 (f(xb0)+f(xb1))
		// c1 = x0*f(xb0)+x1*f(xb1)
		// alpha = 2*c1/(b-a)
		// beta = c0-c1*(a+b)/(b-a)
		//  old : ddelta = (b-a)^2 * f''(_itv)/16
		//  new : ddelta = evaluate the error at the bound and the points when f'(x)=alpha

		double x0,xb0,xb1,fxb0,fxb1,c0,c1;

		x0 = 1.0/std::sqrt(2.);
		xb0 = (0.5)*(itv.diam()*  x0  +itv.lb()+itv.ub());
		xb1 = (0.5)*(itv.diam()*(-x0) +itv.lb()+itv.ub());

		fxb0 = std::tanh(xb0);
		fxb1 = std::tanh(xb1);

		c0 = (0.5)*(fxb0+fxb1);
		c1 = x0*fxb0-x0*fxb1;

		alpha  = 2*c1/(itv.diam());
		beta   = c0-c1*((itv.lb()+itv.ub())/(itv.diam()));
		//ddelta = ((_n*Interval(TEMP1.rad())) + Interval(TEMP2.rad())).ub();

		// compute the maximal error
		ddelta= 0.0;

		// compute the error at _itv.lb(), _itv.ub() and u such as f'(u) =alpha
		ddelta = (abs(tanh(Interval(itv.lb()))-(alpha*Interval(itv.lb())+beta))).ub();
		t1     = (abs(tanh(Interval(itv.ub()))-(alpha*Interval(itv.ub())+beta))).ub();
		if (t1>ddelta)  ddelta= t1;
		// tanh'(u)=alpha
		// cosh(u)= -2/alpha -1
		// u = +-acosh(-2/alpha -1)
		TEMP2 = acosh(-2/Interval(alpha) -1);
		if (!((TEMP2 & itv).is_empty())) {
			t1 = (abs(tanh(TEMP2)-(alpha*TEMP2+beta))).ub();
			if (t1>ddelta)  ddelta= t1;
		}
		if (!(((-TEMP2) & itv).is_empty())) {
			t1 = (abs(tanh(-TEMP2)-(alpha*(-TEMP2)+beta))).ub();
			if (t1>ddelta)  ddelta= t1;
		}

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);

	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Aabs(const Interval& itv){
	Interval res_itv = abs(itv);;

	// Particular case
	if (0<=itv.lb()) return *this;
	if (itv.ub()<=0) {
		this->Aneg();
		return *this;
	}
	if (res_itv.is_empty() || res_itv.is_unbounded() || (!is_actif()) || (itv.diam()<AF_EC)) {
		*this = res_itv;
	}  else  {
		// General case
		double alpha, beta, ddelta, t1, t2;
		Interval TEMP1(0.0), TEMP2(0.0), band(0.0);

		alpha = ((abs(Interval(itv.ub()))-abs(Interval(itv.lb())))/itv.diam()).ub();

		TEMP1 = res_itv.lb()-alpha*Interval(itv.lb());
		TEMP2 = res_itv.ub()-alpha*Interval(itv.ub());
		if (TEMP1.ub()>TEMP2.ub()) {
			// u = 0
			band = Interval(0.0,TEMP1.ub());
		}
		else {
			band = Interval(0.0,TEMP2.ub());
		}

		beta = band.mid();
		t1 = (beta -band).ub();
		t2 = (band -beta).ub();
		ddelta = (t1>t2)? t1 : t2;

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);


	}
	return *this;
}




template<class T>
inline AffineMain<T>& AffineMain<T>::Apow(int n, const Interval& itv) {
	//	std::cout << "in power "<<std::endl;

	// Particular case
	if (itv.is_empty()||itv.is_unbounded()) {
		*this = pow(itv,n);
	} else if (!is_actif()) {
		*this = pow(itv,n);
	} else if (itv.diam()< AF_EC) {
		*this = pow(itv,n);
	} else {
		// General Case
		if (n == 0) {
			*this = Interval::one();
			return *this;
		} else if (n == 1)
			return *this;
		else if (n == 2)
			return this->Asqr(itv);
		else if (n<0) {
			this->Apow(-n,itv);
			return this->Ainv(pow(itv,-n));

		} else if (n % 2 == 0) {
			// alpha = (F(sup(x)) - F(inf(x)))/diam(X)
			// u = (f')^{-1}(alpha)
			// d_a = f(inf(x)) -alpha*inf(X)
			// d_b = f(sup(x)) -alpha*sup(x)
			// d_min = min(d_a,d_b)
			// d_max = f(u) - alpha*u
			// beta = Interval(d_min,d_max).mid()
			// zeta = Interval(d_min,d_max).rad()
			double alpha, beta, ddelta, t1, t2;
			Interval dmm(0.0), TEMP1(0.0), TEMP2(0.0), band(0.0);

			dmm = pow(itv, n);
			alpha = ((__builtin_powi(itv.ub(),n)-__builtin_powi(itv.lb(),n))/itv.diam());

			TEMP1 = (dmm.lb()) - alpha * Interval(itv.lb());
			TEMP2 = (dmm.ub()) - alpha * Interval(itv.ub());
			// u = (alpha/n)^(1/(n-1))
			if (TEMP1.ub() > TEMP2.ub()) {
				TEMP2 = Interval(alpha) / n;
				band = Interval(
						((1 - n) * TEMP2 * (root(TEMP2, n - 1))).lb(),
						TEMP1.ub());
			} else {
				TEMP1 = Interval(alpha) / n;
				band = Interval(
						((1 - n) * TEMP1 * (root(TEMP1, n - 1))).lb(),
						TEMP2.ub());
			}

			beta = band.mid();
			t1 = (beta - band).ub();
			t2 = (band - beta).ub();
			ddelta = (t1 > t2) ? t1 : t2;

			*this *= alpha;
			*this += beta;
			this->inflate(ddelta);
			//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);

		} else {
			// for _itv = [a,b]
			// x0 = 1/sqrt(2)
			// x1= - x0
			// xb0 = 0.5*((b-a)*x0 +(a+b))
			// xb1 = 0.5*((b-a)*x1 +(a+b))
			// c0 = 0.5 (f(xb0)+f(xb1))
			// c1 = x0*f(xb0)+x1*f(xb1)
			// alpha = 2*c1/(b-a)
			// beta = c0-c1*(a+b)/(b-a)
			//  old : ddelta = (b-a)^2 * f''(_itv)/16
			//  new : ddelta = evaluate the error at the bound and the points when f'(x)=alpha

			Interval  TEMP1, TEMP2;
			double alpha, beta, ddelta, t1, x0, xb0, xb1, fxb0, fxb1, c0, c1;

			x0  = 1.0/std::sqrt(2.);
			xb0 = (0.5) * (itv.diam() * ( x0) + itv.lb() + itv.ub());
			xb1 = (0.5) * (itv.diam() * (-x0) + itv.lb() + itv.ub());
			fxb0 = __builtin_powi(xb0, n);
			fxb1 = __builtin_powi(xb1, n);
			c0 = (0.5) * (fxb0 + fxb1);
			c1 = x0 * fxb0 - x0 * fxb1;

			alpha = 2 * c1 / (itv.diam());
			beta  = c0 - c1 * ((itv.lb() + itv.ub()) / (itv.diam()));
			//ddelta = ((_n * Interval(TEMP1.rad())) + Interval(TEMP2.rad())).ub(); //

			// compute the maximal error

			// compute the error at _itv.lb() and _itv.ub()
			ddelta = (abs(
					pow(Interval(itv.lb()), n)
					- (alpha * Interval(itv.lb()) + beta))).ub();
			t1 = (abs(
					pow(Interval(itv.ub()), n)
					- (alpha * Interval(itv.ub()) + beta))).ub();
			if (t1 > ddelta) ddelta= t1 ;
			// u = (alpha/n)^(1/(n-1))
			TEMP2 = pow(Interval(alpha) / n, 1.0 / Interval(n - 1));
			if (!((TEMP2 & itv).is_empty())) {
				t1 = (abs(pow(TEMP2, n) - (alpha * TEMP2 + beta))).ub();
				if (t1 > ddelta) ddelta= t1 ;
			}
			if (!(((-TEMP2) & itv).is_empty())) {
				t1 = (abs(pow(-TEMP2, n) - (alpha * (-TEMP2) + beta))).ub();
				if (t1 > ddelta) ddelta= t1 ;
			}

			*this *= alpha;
			*this += beta;
			this->inflate(ddelta);
			//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);
		}

	}
	//	std::cout << "out power "<<std::endl;
	return *this;
}


template<class T>
inline AffineMain<T>& AffineMain<T>::Apow(double d, const Interval& itv) {
	if ( ((int) (d)) == d) {
		this->Apow((int) (d),itv);
		return *this;
	} else if (d<0) {
		this->Apow(Interval(-d), itv);
		this->Ainv(pow(itv,-d));
		return *this;
	} else {
		this->Apow(Interval(d),itv);
		return *this;
	}
}


template<class T>
inline AffineMain<T>& AffineMain<T>::Aroot(int n, const Interval& itv) {

	if (is_empty()) return *this;
	else if (n==0)  return *this = Interval::one();
	else if (n==1)  return *this;
	else if (is_degenerated()) {
		return *this = pow(Interval(mid()),1.0/n);
	}
	else if (n<0) {
		this->Aroot(-n,itv);
		this->Ainv(root(itv,-n));
		return *this;
	}
	else if (n % 2 == 0) {
		this->Apow(Interval::one()/n,itv);
		return *this; // the negative part of x should be removed
	}
	else if (0 <= itv.lb()) {
		this->Apow(Interval::one()/n,itv);
		return  *this;
	}
	else if (itv.ub() <= 0) {
		this->Aneg();
		this->Apow(Interval::one()/n,-itv);
		this->Aneg();
		return  *this;
	}
	else {
		// TODO do the root when x contains ZERO more properly
		//codac2_error("warning: Affine2Main ROOT non completely well implemented yet");
		//		y=pow(x,e) |  // the negative part of x should be removed
		//	    (-pow(-x,e)); // the positive part of x should be removed
		// BE CAREFULL the result of this union is an INTERVAL, so y lost all its affine form
		return *this = ((pow(itv & Interval(0,oo), Interval::one()/n)) | (-pow(-(itv & Interval({-oo,0})),Interval::one()/n)));  // BE CAREFULL the result of this union is an INTERVAL, so y lost all its affine form
		// BE CAREFULL the result of this union is an INTERVAL, so y lost all its affine form
	/*
	 *  TODO check the proposition of GPT-5.3-Codex
	 * 		// Odd root across zero: use a secant affine model and bound the nonlinearity
		// with an interval remainder instead of collapsing directly to an interval hull.
		Interval root_itv = root(itv, n);
		double alpha = root_itv.diam()/itv.diam();
		double beta = root_itv.mid() - alpha*itv.mid();
		Interval rem = root_itv - (alpha*itv + beta);
		double t1 = std::fabs(rem.lb());
		double t2 = std::fabs(rem.ub());
		double ddelta = (t1>t2)? t1 : t2;

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		return *this;
	 */
	
	}

}

//////////////////////////////////////////////////////////////////////////////////
// MinRange formula


// For x = [a,b]:
// if f is increasing:
//   alpha = f'(a)
//   band  = [f(a)-alpha*a, f(b)-alpha*b]
// if f is decreasing:
//   alpha = f'(b)
//   band  = [f(a)-alpha*b, f(b)-alpha*a]
// then:
//   beta   = mid(band)
//   ddelta = rad(band)

template<class T>
inline AffineMain<T>& AffineMain<T>::Asqrt_MR(const Interval& itv) {

	Interval res_itv = sqrt(itv);
	if (itv.is_unbounded() || res_itv.is_empty() || res_itv.is_unbounded() || (!is_actif()) || (itv.diam()<AF_EC)) {
		*this = res_itv;
	}  else  {  // _actif && b

		double alpha, beta, ddelta, t1, t2;
		Interval dmm(0.0), TEMP1(0.0), TEMP2(0.0), band(0.0), itv2;
		//if ((itv.ub() == POS_INFINITY) || (itv.ub()<0) ) {
		//	_err = sqrt(itv);
		//	_actif = false;
		//	delete[] _val;
		//	_val = NULL;
		//}
		//else {
		if (itv.lb()<0) {
			itv2 = Interval(0.0,itv.ub());
		} else {
			itv2 =itv;
		}

		dmm = sqrt(itv2);
		if (itv2.diam()< AF_EC) {
			alpha = 0.0;
			band =dmm;
		}
		else {
			alpha = (1.0/(2.0*dmm)).lb();  // compute the derivative
			if (alpha<=0) {
				alpha = 0.0;
				band = dmm;
			}
			else {
				TEMP1 = Interval(dmm.lb()) -alpha*itv2.lb();
				TEMP2 = Interval(dmm.ub()) -alpha*itv2.ub();
				if (TEMP1.lb()>TEMP2.ub()) {
					band = Interval(TEMP2.lb(),TEMP1.ub());
					// normally this case never happen
				}
				else {
					band = Interval(TEMP1.lb(),TEMP2.ub());
				}
			}
		}

		beta = band.mid();
		t1 = (beta -band).ub();
		t2 = (band -beta).ub();
		ddelta = (t1>t2)? t1 : t2;

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);

	}
	return *this;
}

template<class T>
AffineMain<T>& AffineMain<T>::Aexp_MR(const Interval& itv) {

	Interval res_itv = exp(itv);
	if (itv.is_unbounded() || res_itv.is_empty() || res_itv.is_unbounded() || (!is_actif()) || (itv.diam()<AF_EC)) {
		*this = res_itv;
	}  else  {  // _actif && b

		double alpha, beta, ddelta, t1, t2;
		Interval dmm(0.0), TEMP1(0.0), TEMP2(0.0), band(0.0), itv2;


		dmm = res_itv;
		if (itv.diam()< AF_EC) {
			alpha = 0.0;
			band =dmm;
		}
		else {
			alpha = dmm.lb();  // compute the derivative
			if (alpha<=0) {
				alpha = 0.0;
				band = dmm;
			}
			else {
				TEMP1 = Interval(dmm.lb()) -alpha*itv.lb();
				TEMP2 = Interval(dmm.ub()) -alpha*itv.ub();
				if (TEMP1.lb()>TEMP2.ub()) {
					band = Interval(TEMP2.lb(),TEMP1.ub());
					// normally this case never happen
				}
				else {
					band = Interval(TEMP1.lb(),TEMP2.ub());
				}
			}
		}

		beta = band.mid();
		t1 = (beta -band).ub();
		t2 = (band -beta).ub();
		ddelta = (t1>t2)? t1 : t2;

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);

	}


	return *this;
}

template<class T>
AffineMain<T>& AffineMain<T>::Alog_MR(const Interval& itv) {

	Interval res_itv = log(itv);
	if (itv.is_unbounded() || res_itv.is_empty() || res_itv.is_unbounded() || (!is_actif()) || (itv.diam()<AF_EC)) {
		*this = res_itv;
	}  else  {  // _actif && b

		double alpha, beta, ddelta, t1, t2;
		Interval dmm(0.0), TEMP1(0.0), TEMP2(0.0), band(0.0), itv2;

		//if ((itv.lb()<= 0) ||(itv.ub()==POS_INFINITY) ) {
		//	_err = log(itv);
		//	_actif = false;
		//	delete[] _val;
		//	_val = NULL;
		//}
		//else {

		dmm = res_itv;
		if (itv.diam()< AF_EC) {
			alpha = 0.0;
			band =dmm;
		}
		else {
			alpha = (1.0/itv).lb();  // compute the derivative
			if (alpha<=0) {
				alpha = 0.0;
				band = dmm;
			}
			else {
				TEMP1 = Interval(dmm.lb()) -alpha*itv.lb();
				TEMP2 = Interval(dmm.ub()) -alpha*itv.ub();
				if (TEMP1.lb()>TEMP2.ub()) {
					band = Interval(TEMP2.lb(),TEMP1.ub());
					// normally this case never happen
				}
				else {
					band = Interval(TEMP1.lb(),TEMP2.ub());
				}
			}
		}

		beta = band.mid();
		t1 = (beta -band).ub();
		t2 = (band -beta).ub();
		ddelta = (t1>t2)? t1 : t2;

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);

		//}

	}
	return *this;
}

template<class T>
AffineMain<T>& AffineMain<T>::Ainv_MR(const Interval& itv) {

	Interval res_itv = (1.0/itv);

	if (itv.is_unbounded() || res_itv.is_empty() || res_itv.is_unbounded() || (!is_actif()) || (itv.diam()<AF_EC)) {
		*this = res_itv;
	}  else  {  // _actif && b

		double alpha, beta, ddelta, t1, t2;
		Interval dmm(0.0), TEMP1(0.0), TEMP2(0.0), band(0.0), itv2;

		dmm = (1.0/abs(itv));
		if (itv.diam()< AF_EC) {
			alpha = 0.0;
			band = dmm;
		}
		else {
			alpha = ((-1.0)/pow(itv,2)).ub(); // compute the derivative
			if (alpha<=0) {
				alpha = 0.0;
				band = dmm;
			}
			else {
				TEMP1 = Interval(dmm.lb()) -alpha*itv.ub();
				TEMP2 = Interval(dmm.ub()) -alpha*itv.lb();
				if (TEMP1.lb()>TEMP2.ub()) {
					band = Interval(TEMP2.lb(),TEMP1.ub());
					// normally this case never happen
				}
				else {
					band = Interval(TEMP1.lb(),TEMP2.ub());
				}
			}
		}

		beta = band.mid();
		t1 = (beta -band).ub();
		t2 = (band -beta).ub();
		ddelta = (t1>t2)? t1 : t2;

		if (itv.lb()<0.0) beta = -beta;

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);

	}
	return *this;
}





} // end namespace codac2


/** \brief atan2(AF[y],AF[x]). */
//Affine2 atan2(const Affine2& y, const Affine2& x);
/** \brief atan2([y],AF[x]). */
//Affine2 atan2(const Interval& y, const Affine2& x);
/** \brief atan2(AF[y],[x]). */
//Affine2 atan2(const Affine2& y, const Interval& x);
/** \brief cosh(AF[x]). */
