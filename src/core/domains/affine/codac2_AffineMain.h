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

#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <list>
#include <ostream>
#include <type_traits>
#include <utility>

#include "codac2_Interval.h"
#include "codac2_Affine2_fAF2.h"
#include "codac2_assert.h"
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

template<class T>
struct is_affine_based< AffineMain<T> > : std::true_type {};

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
      ReadCost = 4,
      AddCost = 20,
      MulCost = 50
    };
  };

  template<class T,typename BinOp>
  struct ScalarBinaryOpTraits<codac2::AffineMain<T>,double,BinOp>
  { typedef codac2::AffineMain<T> ReturnType; };

  template<class T,typename BinOp>
  struct ScalarBinaryOpTraits<double,codac2::AffineMain<T>,BinOp>
  { typedef codac2::AffineMain<T> ReturnType; };

  template<class T,typename BinOp>
  struct ScalarBinaryOpTraits<codac2::AffineMain<T>,codac2::Interval,BinOp>
  { typedef codac2::AffineMain<T> ReturnType; };

  template<class T,typename BinOp>
  struct ScalarBinaryOpTraits<codac2::Interval,codac2::AffineMain<T>,BinOp>
  { typedef codac2::AffineMain<T> ReturnType; };

  template<class T,typename BinOp>
  struct ScalarBinaryOpTraits<codac2::AffineMain<T>,codac2::AffineMain<T>,BinOp>
  { typedef codac2::AffineMain<T> ReturnType; };

}

namespace codac2 {


/**
 *
 * \brief Main affine-form domain based on \c AF_Default.
 *
 * The affine form stores an interval enclosure together with an optional
 * affine dependency representation. Its state is encoded by \c AffineStatus:
 * - \c Active: active affine form with noise variables
 * - \c Degenerate: singleton interval
 * - \c Empty: empty set
 * - \c AllReals: \f$[-\infty,+\infty]\f$
 * - \c UpperUnbounded: \f$[a,+\infty]\f$
 * - \c LowerUnbounded: \f$[-\infty,a]\f$
 */
enum class AffineStatus : std::int8_t
{
  LowerUnbounded = -4,
  UpperUnbounded = -3,
  AllReals       = -2,
  Empty          = -1,
  Degenerate     = 0,
  Active         = 1
};

//=================================================================================================================
//=================================================================================================================
//=================================================================================================================

template<class T>
class AffineMain  : public DomainInterface<Interval,double>, public DomainInterface<AffineMain<T>,double>
{


public:

	typedef enum {
		AF_Lin_Chebyshev=1, 
		AF_Lin_MinRange=0
	} Affine_Mode; // ...etc...

	/**
	 * \brief Changes the linearization mode for the current thread.
	 *
	 * \param tt approximation mode; Chebyshev is used by default
	 */
	static void change_mode(Affine_Mode tt=AF_Lin_Chebyshev);
	/**
	 * \brief Returns the linearization mode for the current thread.
	 *
	 * \return the approximation mode
	 */
	static Affine_Mode get_mode();

	/** \brief Creates an unbounded affine form (like \c Interval()). */
	AffineMain();

	/**
	 * \brief Creates a degenerate affine form from a real value.
	 *
	 * \param value real value
	 */
	AffineMain(double value);
	/**
	 * \brief Creates an affine form from an interval.
	 *
	 * \param itv interval enclosure
	 */
	AffineMain(const Interval& itv);
	/**
	 * \brief Converts this affine form to its interval enclosure.
	 *
	 * \return interval enclosure of this affine form
	 */
	operator Interval() const;

	/**
	 * \brief Creates an affine form by copy.
	 *
	 * \param x affine form to copy
	 */
	AffineMain(const AffineMain<T>& x);

	/** \brief Destroys the affine form. */
	virtual ~AffineMain() { };

	/** \brief Returns \f$-*\mathrm{this}\f$. */
	AffineMain<T> operator-() const;

	/**
	 * \brief Tests exact equality with another affine form.
	 * Compares the interval enclosures of two affine forms.
	 * Affine coefficients and dependency structures are not compared.
	 *
	 * \param x affine form to compare with
	 * \return true iff both interval enclosures are exactly equal
	 */
	bool operator==(const AffineMain<T>& x) const;

	/**
	 * \brief Tests exact equality with an interval.
	 *
	 * \param x interval to compare with
	 * \return true iff both intervals are exactly equal
	 */
	bool operator==(const Interval& x) const;

	/**
	 * \brief Tests exact equality with a real value.
	 *
	 * \param x interval to compare with
	 * \return true iff both intervals are exactly equal
	 */
	bool operator==(const double x) const;

	/**
	 * \brief Tests non-equality with another affine form.
	 *
	 * \param x affine form to compare with
	 * \return true iff interval enclosures are different
	 */
	bool operator!=(const AffineMain<T>& x) const;

	/**
	 * \brief Tests non-equality with an interval.
	 *
	 * \param x interval to compare with
	 * \return true iff intervals are different
	 */
	bool operator!=(const Interval& x) const ;

     /**
     * \brief Comparison (strict less-than) between this and x
     * 
     * The returned ``BoolInterval`` encloses the truth value of
     * \f$t<s\f$ for \f$t\in[\mathrm{this}]\f$ and \f$s\in[x]\f$.
     * 
     * \note Returns:
     *       - ``BoolInterval::EMPTY`` if this or x is empty
     *       - ``BoolInterval::TRUE`` iff \f$\mathrm{ub}([\mathrm{this}])<\mathrm{lb}([x])\f$
     *       - ``BoolInterval::FALSE`` iff \f$\mathrm{ub}([x])\leq\mathrm{lb}([\mathrm{this}])\f$
     *       - ``BoolInterval::UNKNOWN`` otherwise
     * 
     * \param x interval to be compared with
     * \return interval Boolean result
     */
	BoolInterval operator<(const Interval& x) const;

     /**
     * \brief Comparison (strict greater-than) between this and x
     * 
     * The returned BoolInterval encloses the truth value of
     * \f$t>s\f$ for \f$t\in[\mathrm{this}]\f$ and \f$s\in[x]\f$.
     * 
     * \note Returns:
     *       - ``BoolInterval::EMPTY`` if this or x is empty
     *       - ``BoolInterval::TRUE`` iff \f$\mathrm{lb}([\mathrm{this}])>\mathrm{ub}([x])\f$
     *       - ``BoolInterval::FALSE`` iff \f$\mathrm{lb}([x])\geq\mathrm{ub}([\mathrm{this}])\f$
     *       - ``BoolInterval::UNKNOWN`` otherwise
     * 
     * \param x interval to be compared with
     * \return interval Boolean result
     */
	BoolInterval operator>(const Interval& x) const;

	/**
	 * \brief Comparison (strict less-than) between two affine forms.
	 *
	 * The comparison is performed on their interval enclosures.
	 *
	 * \param x affine form to be compared with
	 * \return interval Boolean result
	 */
	BoolInterval operator<(const AffineMain<T>& x) const;

	/**
	 * \brief Comparison (strict greater-than) between two affine forms.
	 *
	 * The comparison is performed on their interval enclosures.
	 *
	 * \param x affine form to be compared with
	 * \return interval Boolean result
	 */
	BoolInterval operator>(const AffineMain<T>& x) const;


	/**
	 * \brief Tests non-equality with a real value.
	 *
	 * \param x interval to compare with
	 * \return true iff intervals are different
	 */
	bool operator!=(const double x) const;

	/** \brief Sets this affine form to the empty set. */
	void set_empty();

	/**
	 * \brief Assigns from another affine form.
	 *
	 * \param x affine form to copy
	 * \return a reference to this
	 */
	AffineMain<T>& operator=(const AffineMain<T>& x);

	/**
	 * \brief Assigns from a real value.
	 *
	 * \param x real value
	 * \return a reference to this
	 */
	AffineMain<T>& operator=(double x);

	/**
	 * \brief Assigns from an interval.
	 *
	 * \param itv interval value
	 * \return a reference to this
	 */
	virtual AffineMain<T>& operator=(const Interval& itv);

	/**
	 * \brief Sets this affine form to [-oo,oo]
	 *
	 * \note This function is used for template purposes.
	 *
	 * \return a reference to this
	 */
	AffineMain<T>& init();

	/**
	 * \brief Sets this affine form to x
	 *
	 * \note This function is used for template purposes.
	 *
	 * \param x the value for re-initialization
	 * \return a reference to this
	 */
	AffineMain<T>& init(const Interval& x);

	/**
	 * \brief Sets the bounds of this affine form as the hull of a list of values.
	 *
	 * \note This function is separated from the constructor for py binding purposes.
	 *
	 * \param l list of one or two values contained in the resulting affine form
	 * \return a reference to this
	 */
	AffineMain<T>& init_from_list(const std::list<double>& l);


	// AffineMain(AffineMain<T>&&) noexcept = default;

	// AffineMain<T>& operator=(AffineMain<T>&&) noexcept = default;

	/* Union and Intersection of two Affine form must not be implemented
	 * That could produce to much confusion.
	 */
	/**
	 * \brief Intersection assignment is disabled for affine forms.
	 *
	 * This operation is intentionally deleted because the intersection of two
	 * affine forms cannot be represented without discarding affine information.
	 *
	 * \param x affine form
	 */
	AffineMain<T>& operator&=(const AffineMain<T>& x) = delete;

	/**
	 * \brief Hull assignment is disabled for affine forms.
	 *
	 * This operation is intentionally deleted because the hull of two affine
	 * forms would discard affine dependency information.
	 *
	 * \param x affine form
	 */
	AffineMain<T>& operator|=(const AffineMain<T>& x) = delete;

	/**
	 * \brief Adds \f$[-\mathrm{rad},+\mathrm{rad}]\f$ to this.
	 *
	 * \param radd inflation radius
	 * \return a reference to this
	 */
	AffineMain<T>& inflate(double radd);

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
	 * \brief Returns the dimension of this (which is always 1).
	 *
	 * \note This function is used for template purposes.
	 *
	 * \return 1
	 */
	Index size() const;

	/**
	 * \brief Returns the number of represented noise variables.
	 *
	 * \return number of variables
	 */
	Index noise_count() const;

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
	double noise(Index i) const;

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
	bool is_active() const;

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
	 * \brief Tests whether the interval enclosure of this affine form is a subset of \p x.
	 *
	 * \param x interval or affine form to compare with
	 * \note Returns true when this affine form is empty.
	 * \return true iff the interval enclosure is a subset of \p x
	 */
	bool is_subset(const AffineMain<T>& x) const;

	/**
	 * \brief Tests whether this interval enclosure is a strict subset of \p x.
	 *
	 * \note In particular, (-oo,oo) is not a strict subset of (-oo,oo)
	 * and the empty set is not a strict subset of the empty set although
	 * in both cases, the first is inside the interior of the second.
	 */
	bool is_strict_subset(const Interval& x) const;
	/**
	 * \brief Tests whether the interval enclosure of this affine form is a strict subset of \p x.
	 *
	 * \param x interval or affine form to compare with
	 * \return true iff the enclosure is a strict subset of \p x
	 */
	bool is_strict_subset(const AffineMain<T>& x) const;

	/**
	 * \brief Tests whether this interval enclosure is inside the interior of \p x.
	 *
	 * \note In particular, (-oo,oo) is in the interior of (-oo,oo)
	 * and the empty set is in the interior of the empty set.
	 * \note Always returns true if this is empty.
	 */
	bool is_interior_subset(const Interval& x) const;
	/**
	 * \brief Tests whether the interval enclosure of this affine form is contained in the interior of \p x.
	 *
	 * \param x interval or affine form to compare with
	 * \note Returns true when this affine form is empty.
	 * \return true iff the enclosure is contained in the interior of \p x
	 */
	bool is_interior_subset(const AffineMain<T>& x) const;

	// NOTE: is_relative_interior_subset() was removed: it forwarded to
	// Interval::is_relative_interior_subset(), which codac2::Interval does
	// not provide, so it could never be instantiated. See the Python
	// bindings (codac2_py_Affine.cpp) for the same finding.

	/**
	 * \brief Tests whether this interval enclosure is a strict interior subset of \p x.
	 *
	 * \note In particular, (-oo,oo) is not "strictly" in the interior of (-oo,oo)
	 * and the empty set is not "strictly" in the interior of the empty set.
	 */
	bool is_strict_interior_subset(const Interval& x) const;
	/**
	 * \brief Tests whether the interval enclosure is strictly contained in the interior of \p x.
	 *
	 * \param x interval or affine form to compare with
	 * \return true iff the enclosure is strictly contained in the interior of \p x
	 */
	bool is_strict_interior_subset(const AffineMain<T>& x) const;

	/**
	 * \brief Tests whether this interval enclosure is a superset of \p x.
	 *
	 * \param x interval or affine form to compare with
	 * \return true iff the enclosure is a superset of \p x
	 */
	bool is_superset(const Interval& x) const;
	/**
	 * \brief Tests whether the interval enclosure of this affine form contains \p x.
	 *
	 * \param x interval or affine form to compare with
	 * \return true iff the enclosure is a superset of \p x
	 */
	bool is_superset(const AffineMain<T>& x) const;

	/**
	 * \brief Tests whether this interval enclosure is a strict superset of \p x.
	 *
	 * \param x interval or affine form to compare with
	 * \return true iff the enclosure is a strict superset of \p x
	 */
	bool is_strict_superset(const Interval& x) const;
	/**
	 * \brief Tests whether the interval enclosure of this affine form is a strict superset of \p x.
	 *
	 * \param x interval or affine form to compare with
	 * \return true iff the enclosure is a strict superset of \p x
	 */
	bool is_strict_superset(const AffineMain<T>& x) const;

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
	 * \brief Tests whether this affine form intersects \p x.
	 *
	 * \param x interval or affine form to compare with
	 * \return true iff the interval enclosures have a non-empty intersection
	 */
	bool intersects(const AffineMain<T>& x) const;

	/**
	 * \brief Tests whether this and \p x overlap with non-null volume.
	 *
	 * Equivalently, some interior points (of this or x) must belong to the intersection.
	 */
	bool overlaps(const Interval &x) const;

	/**
	 * \brief Tests whether this affine form and \p x overlap with non-zero volume.
	 *
	 * \param x interval or affine form to compare with
	 * \return true iff their intersection has non-zero volume
	 */
	bool overlaps(const AffineMain<T> &x) const;

	/**
	 * \brief Tests whether this and \p x are disjoint.
	 *
	 * \param x interval value
	 * \return true iff this and \p x do not intersect
	 */
	bool is_disjoint(const Interval &x) const;

	/**
	 * \brief Tests whether this affine form and \p x are disjoint.
	 *
	 * \param x interval or affine form to compare with
	 * \return true iff their interval enclosures do not intersect
	 */
	bool is_disjoint(const AffineMain<T> &x) const;

	/**
	 * \brief Adds \p d to this.
	 *
	 * \param d real value
	 * \return a reference to this
	 */
	AffineMain<T>& operator+=(double d);

	/**
	 * \brief Subtracts \p d from this.
	 *
	 * \param d real value
	 * \return a reference to this
	 */
	AffineMain<T>& operator-=(double d);

	/**
	 * \brief Multiplies this by \p d.
	 *
	 * \param d real value
	 * \return a reference to this
	 */
	AffineMain<T>& operator*=(double d);

	/**
	 * \brief Divides this by \p d.
	 *
	 * \param d real value
	 * \return a reference to this
	 */
	AffineMain<T>& operator/=(double d) ;

	/**
	 * \brief Adds interval \p x to this.
	 *
	 * \param x interval value
	 * \return a reference to this
	 */
	AffineMain<T>& operator+=(const Interval& x);

	/**
	 * \brief Subtracts interval \p x from this.
	 *
	 * \param x interval value
	 * \return a reference to this
	 */
	AffineMain<T>& operator-=(const Interval& x);

	/**
	 * \brief Multiplies this by interval \p x.
	 *
	 * \param x interval value
	 * \return a reference to this
	 */
	AffineMain<T>& operator*=(const Interval& x);

	/**
	 * \brief Divides this by interval \p x.
	 *
	 * \param x interval value
	 * \return a reference to this
	 */
	AffineMain<T>& operator/=(const Interval& x);

	/**
	 * \brief Adds affine form \p x to this.
	 *
	 * \param x affine form value
	 * \return a reference to this
	 */
	AffineMain<T>& operator+=(const AffineMain<T>& x);

	/**
	 * \brief Subtracts affine form \p x from this.
	 *
	 * \param x affine form value
	 * \return a reference to this
	 */
	AffineMain<T>& operator-=(const AffineMain<T>& x);

	/**
	 * \brief Multiplies this by affine form \p x.
	 *
	 * \param x affine form value
	 * \return a reference to this
	 */
	AffineMain<T>& operator*=(const AffineMain<T>& x);

	/**
	 * \brief Divides this by affine form \p x.
	 *
	 * \param x affine form value
	 * \return a reference to this
	 */
	AffineMain<T>& operator/=(const AffineMain<T>& x);

    /**
     * \brief Provides an empty Affine Form
     *
     * \return an empty set
     */
    static AffineMain<T> empty();


	private:
	/**
	 * \brief Internal helper that computes the square of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>& Asqr(const Interval& itv);
	/**
	 * \brief Internal helper that computes the negation of the current affine form.
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Aneg();
	/**
	 * \brief Internal helper that computes the inverse of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Ainv(const Interval& itv);
	/**
	 * \brief Internal helper that computes the square root of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Asqrt(const Interval& itv);
	/**
	 * \brief Internal helper that computes the exponential of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Aexp(const Interval& itv);
	/**
	 * \brief Internal helper that computes the natural logarithm of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Alog(const Interval& itv);
	/**
	 * \brief Internal helper that computes the power of the current affine form. *
	 * \param n integer exponent
	 * \param itv interval enclosure of the base
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Apow(int n, const Interval& itv);
	/**
	 * \brief Internal helper that computes the power of the current affine form. *
	 * \param d real exponent
	 * \param itv interval enclosure of the base
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Apow(double d, const Interval& itv);
	/**
	 * \brief Internal helper that computes the power of the current affine form. *
	 * \param y interval exponent
	 * \param itvx interval enclosure of the base
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Apow(const Interval &y, const Interval& itvx);
	/**
	 * \brief Internal helper that computes the n-th root of the current affine form. *
	 * \param n root degree
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Aroot(int n, const Interval& itv);
	/**
	 * \brief Internal helper that computes the cosine of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Acos(const Interval& itv);
	/**
	 * \brief Internal helper that computes the sine of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Asin(const Interval& itv);
	/**
	 * \brief Internal helper that computes the tangent of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Atan(const Interval& itv);
	/**
	 * \brief Internal helper that computes the arc-cosine of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Aacos(const Interval& itv);
	/**
	 * \brief Internal helper that computes the arc-sine of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Aasin(const Interval& itv);
	/**
	 * \brief Internal helper that computes the arc-tangent of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Aatan(const Interval& itv);
	/**
	 * \brief Internal helper that computes the hyperbolic cosine of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Acosh(const Interval& itv);
	/**
	 * \brief Internal helper that computes the hyperbolic sine of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Asinh(const Interval& itv);
	/**
	 * \brief Internal helper that computes the hyperbolic tangent of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Atanh(const Interval& itv);
	/**
	 * \brief Internal helper that computes the absolute value of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Aabs(const Interval& itv);
	/**
	 * \brief Internal helper that computes the inverse using Chebyshev linearization of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Ainv_CH(const Interval& itv);
	/**
	 * \brief Internal helper that computes the square root using Chebyshev linearization of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Asqrt_CH(const Interval& itv);
	/**
	 * \brief Internal helper that computes the exponential using Chebyshev linearization of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Aexp_CH(const Interval& itv);
	/**
	 * \brief Internal helper that computes the natural logarithm using Chebyshev linearization of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Alog_CH(const Interval& itv);
	/**
	 * \brief Internal helper that computes the inverse using MinRange linearization of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Ainv_MR(const Interval& itv);
	/**
	 * \brief Internal helper that computes the square root using MinRange linearization of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Asqrt_MR(const Interval& itv);
	/**
	 * \brief Internal helper that computes the exponential using MinRange linearization of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Aexp_MR(const Interval& itv);
	/**
	 * \brief Internal helper that computes the natural logarithm using MinRange linearization of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Alog_MR(const Interval& itv);
	/**
	 * \brief Internal helper that computes the inverse hyperbolic cosine of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>& Aacosh(const Interval& itv);
	/**
	 * \brief Internal helper that computes the inverse hyperbolic sine of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>& Aasinh(const Interval& itv);
	/**
	 * \brief Internal helper that computes the inverse hyperbolic tangent of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>& Aatanh(const Interval& itv);
	/**
	 * \brief Internal helper that computes the two-argument arctangent of the current affine form. *
	 * \param x affine form representing the x-coordinate
	 * \param itvY interval enclosure of the y-coordinate
	 * \param itvX interval enclosure of the x-coordinate
	 *
	 * \return reference to this affine form
	 */
	AffineMain<T>& Aatan2(const AffineMain<T>& x, const Interval& itvY, const Interval& itvX);

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
	template<class A>	friend AffineMain<A> pow(double  d, const AffineMain<A>&  x);
	template<class A>	friend AffineMain<A> pow(const Interval& itv, const AffineMain<A>&  x);
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

	template<class A>	friend AffineMain<A> acosh(const AffineMain<A>&  x);
	template<class A>	friend AffineMain<A> asinh(const AffineMain<A>&  x);
	template<class A>	friend AffineMain<A> atanh(const AffineMain<A>&  x);

	template<class A> 	friend AffineMain<A> atan2(const AffineMain<A>& y, const AffineMain<A>& x);

	template<class A> friend std::ostream& operator<<(std::ostream& os, const AffineMain<A>& x);



protected:

	static thread_local Affine_Mode mode;	// one mode per thread, Chebyshev by default


	/** \brief Tolerance for the default compacting procedure. */
	inline static constexpr double AF_COMPAC_Tol = 1.e-6;
	inline static constexpr double AF_EM = 0x1p-51;
	inline static constexpr double AF_EC = 0x1p-55;
	inline static constexpr double AF_EE = 2.0;


	/**
	 * \brief Encodes the current affine-state mode.
	 *
	 * See class description above for the exact code-to-state mapping.
	 */
	AffineStatus _status;		// Current affine-form status.

	Index _n_noise; 		// dimension (size of _elt._val)-1  , ie number of noise symbols

	T _elt;			// core of the affine form

	/**
	 * \brief Creates an affine form with \p size variables.
	 *
	 * The \p var-th variable is initialized from \p itv.
	 */
	AffineMain(Index noise_count, Index var, const Interval& itv);

	/**
	 * \brief Changes the number of affine variables.
	 *
	 * \param n new number of variables
	 */
	void resize_noise(Index n);

};




template<class T>
inline AffineMain<T> AffineMain<T>::empty()
{
	AffineMain<T> result;
	result.set_empty();
	return result;
}


template<class T>
thread_local typename AffineMain<T>::Affine_Mode AffineMain<T>::mode = AffineMain<T>::AF_Lin_Chebyshev;// default mode is Chebyshev


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
inline Index AffineMain<T>::size() const { return 1; }

template<class T>
inline AffineMain<T>& AffineMain<T>::init()
{
	*this = Interval(-oo,oo);
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::init(const Interval& x)
{
	*this = x;
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::init_from_list(const std::list<double>& l)
{
	if (l.size() == 1)
		*this = Interval(*l.begin());

	else if (l.size() == 2)
		*this = Interval(*l.begin(), *std::prev(l.end()));

	else
	{
		assert_release("'AffineMain' can only be defined by one or two 'double' values.");
	}

	return *this;
}

template<class T>
inline bool AffineMain<T>::is_subset(const Interval& x) const { return this->itv().is_subset(x); }

template<class T>
inline bool AffineMain<T>::is_strict_subset(const Interval& x) const { return this->itv().is_strict_subset(x); }

template<class T>
inline bool AffineMain<T>::is_interior_subset(const Interval& x) const { return this->itv().is_interior_subset(x); }

template<class T>
inline bool AffineMain<T>::is_strict_interior_subset(const Interval& x) const { return this->itv().is_strict_interior_subset(x); }

template<class T>
inline bool AffineMain<T>::is_superset(const Interval& x) const { return this->itv().is_superset(x); }

template<class T>
inline bool AffineMain<T>::is_strict_superset(const Interval& x) const { return this->itv().is_strict_superset(x); }



template<class T>
inline bool AffineMain<T>::is_subset(const AffineMain<T>& x) const { return this->itv().is_subset(x.itv()); }

template<class T>
inline bool AffineMain<T>::is_strict_subset(const AffineMain<T>& x) const { return this->itv().is_strict_subset(x.itv()); }

template<class T>
inline bool AffineMain<T>::is_interior_subset(const AffineMain<T>& x) const { return this->itv().is_interior_subset(x.itv()); }

template<class T>
inline bool AffineMain<T>::is_strict_interior_subset(const AffineMain<T>& x) const { return this->itv().is_strict_interior_subset(x.itv()); }

template<class T>
inline bool AffineMain<T>::is_superset(const AffineMain<T>& x) const { return this->itv().is_superset(x.itv()); }

template<class T>
inline bool AffineMain<T>::is_strict_superset(const AffineMain<T>& x) const { return this->itv().is_strict_superset(x.itv()); }


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

template<class T>
inline bool AffineMain<T>::intersects(const AffineMain<T> &x) const { return this->itv().intersects(x.itv()); }

template<class T>
inline bool AffineMain<T>::overlaps(const AffineMain<T> &x) const { return this->itv().overlaps(x.itv()); }

template<class T>
inline bool AffineMain<T>::is_disjoint(const AffineMain<T> &x) const { return this->itv().is_disjoint(x.itv()); }


/**
 * \brief Compares an interval with an affine form using strict less-than.
 *
 * The comparison is performed on their interval enclosures.
 *
 * \param x interval
 * \param y affine form
 * \return interval Boolean enclosing the truth value of \f$x < [y]\f$
 */
template<class T>
BoolInterval operator<(const Interval& x, const AffineMain<T>& y);

/**
 * \brief Compares an interval with an affine form using strict greater-than.
 *
 * The comparison is performed on their interval enclosures.
 *
 * \param x interval
 * \param y affine form
 * \return interval Boolean enclosing the truth value of \f$x > [y]\f$
 */
template<class T>
BoolInterval operator>(const Interval& x, const AffineMain<T>& y);



/**
 * \brief Returns the additive inverse of an affine form.
 *
 * \return \f$-x\f$
 */
template<class T>
inline AffineMain<T> AffineMain<T>::operator-() const {
	AffineMain<T> res(*this);
	res.Aneg();
	return res;
}

/**
 * \brief Adds two affine forms.
 *
 * \param x1 first affine form
 * \param x2 second affine form
 * \return affine sum
 */
template<class T>
AffineMain<T> operator+(const AffineMain<T>&  x1, const AffineMain<T>&  x2);

/**
 * \brief Adds a real value to an affine form.
 *
 * \param x affine form
 * \param d real value
 * \return \f$x+d\f$
 */
template<class T>
AffineMain<T> operator+(const AffineMain<T>&  x, double d);

/**
 * \brief Adds a real value to an affine form.
 *
 * \param d real value
 * \param x affine form
 * \return \f$d+x\f$
 */
template<class T>
AffineMain<T> operator+(double d, const AffineMain<T>&  x);

/**
 * \brief Adds an interval to an affine form.
 *
 * \param x1 affine form
 * \param x2 interval
 * \return \f$x_1+x_2\f$
 */
template<class T>
AffineMain<T> operator+(const AffineMain<T>&  x1, const Interval& x2);

/**
 * \brief Adds an interval to an affine form.
 *
 * \param x1 interval
 * \param x2 affine form
 * \return affine sum
 */

template<class T>
AffineMain<T> operator+(const Interval& x1, const AffineMain<T>&  x2);

/**
 * \brief Subtracts two affine forms.
 *
 * \param x1 first affine form
 * \param x2 second affine form
 * \return affine difference
 */
template<class T>
AffineMain<T> operator-(const AffineMain<T>&  x1, const AffineMain<T>&  x2);

/**
 * \brief Returns the difference between an affine form and a real value.
 *
 * \param x affine form
 * \param d real value
 * \return \f$x-d\f$
 */
template<class T>
AffineMain<T> operator-(const AffineMain<T>&  x, double d);

/**
 * \brief Subtracts an affine form from a real value.
 *
 * \param d real value
 * \param x affine form
 * \return \f$d-x\f$
 */

template<class T>
AffineMain<T> operator-(double d, const AffineMain<T>&  x);

/**
 * \brief Subtracts an interval from an affine form.
 *
 * \param x1 affine form
 * \param x2 interval
 * \return \f$x_1-x_2\f$
 */
template<class T>
AffineMain<T> operator-(const AffineMain<T>&  x1, const Interval& x2);

/**
 * \brief Subtracts an affine form from an interval.
 *
 * \param x1 interval
 * \param x2 affine form
 * \return \f$x_1-x_2\f$
 */

template<class T>
AffineMain<T> operator-(const Interval& x1, const AffineMain<T>&  x2);

/**
 * \brief Multiplies two affine forms.
 *
 * \param x1 first affine form
 * \param x2 second affine form
 * \return affine product
 */
template<class T>
AffineMain<T> operator*(const AffineMain<T>&  x1, const AffineMain<T>&  x2);

/**
 * \brief Multiplies an affine form by a real value.
 *
 * \param x affine form
 * \param d real value
 * \return \f$x d\f$
 */
template<class T>
AffineMain<T> operator*(const AffineMain<T>&  x, double d);

/**
 * \brief Multiplies an affine form by a real value.
 *
 * \param d real value
 * \param x affine form
 * \return \f$d x\f$
 */
template<class T>
AffineMain<T> operator*(double d, const AffineMain<T>&  x);

/**
 * \brief Multiplies an affine form by an interval.
 *
 * \param x1 affine form
 * \param x2 interval
 * \return affine product
 */
template<class T>
AffineMain<T> operator*(const AffineMain<T>&  x1, const Interval& x2);

/**
 * \brief Multiplies an interval by an affine form.
 *
 * \param x1 interval
 * \param x2 affine form
 * \return affine product
 */

template<class T>
AffineMain<T> operator*(const Interval& x1, const AffineMain<T>&  x2);

/**
 * \brief Divides two affine forms.
 *
 * \param x1 numerator affine form
 * \param x2 denominator affine form
 * \return affine quotient
 */
template<class T>
AffineMain<T> operator/(const AffineMain<T>&  x1, const AffineMain<T>&  x2);

/**
 * \brief Divides an affine form by a real value.
 *
 * \param x affine form
 * \param d real divisor
 * \return affine quotient
 */
template<class T>
AffineMain<T> operator/(const AffineMain<T>&  x, double d);

/**
 * \brief Divides a real value by an affine form.
 *
 * \param d real value
 * \param x affine form
 * \return affine quotient
 */
template<class T>
AffineMain<T> operator/(double d, const AffineMain<T>&  x);

/**
 * \brief Divides an affine form by an interval.
 *
 * \param x1 affine form
 * \param x2 interval
 * \return affine quotient
 */
template<class T>
AffineMain<T> operator/(const AffineMain<T>&  x1, const Interval& x2);

/**
 * \brief Divides an interval by an affine form.
 *
 * \param x1 interval numerator
 * \param x2 affine form denominator
 * \return affine quotient
 */

template<class T>
AffineMain<T> operator/(const Interval& x1, const AffineMain<T>&  x2);


// NOTE: distance(x1,x2) was removed: every overload forwarded to a free
// function codac2::distance(const Interval&, const Interval&) that
// codac2 does not provide, so none of them could ever be instantiated.
// See the Python bindings (codac2_py_Affine_operations.cpp) for the same
// finding.


/**
 * \brief Returns the multiplicative inverse of an affine form.
 *
 * \param x affine form
 * \return affine result
 */
template<class T>
AffineMain<T> inv(const AffineMain<T>&  x);

/**
 * \brief Returns the square of an affine form.
 *
 * \param x affine form
 * \return affine result
 */
template<class T>
AffineMain<T> sqr(const AffineMain<T>&  x);

/**
 * \brief Returns the square root of an affine form.
 *
 * \param x affine form
 * \return affine result
 */
template<class T>
AffineMain<T> sqrt(const AffineMain<T>&  x);

/**
 * \brief Returns the exponential of an affine form.
 *
 * \param x affine form
 * \return affine result
 */
template<class T>
AffineMain<T> exp(const AffineMain<T>&  x);

/**
 * \brief Returns the natural logarithm of an affine form.
 *
 * \param x affine form
 * \return affine result
 */
template<class T>
AffineMain<T> log(const AffineMain<T>&  x);

/**
 * \brief Returns an integer power of an affine form.
 *
 * \param x affine form
 * \param n integer exponent
 * \return affine result
 */
template<class T>
AffineMain<T> pow(const AffineMain<T>&  x, int n);

/**
 * \brief Returns a real power of an affine form.
 *
 * \param x affine form
 * \param d real exponent
 * \return affine result
 */
template<class T>
AffineMain<T> pow(const AffineMain<T>&  x, double d);

/**
 * \brief Returns an affine form raised to an interval exponent.
 *
 * \param x affine form base
 * \param y interval exponent
 * \return affine result
 *
 * \note The dependency structure of the exponent is not represented when
 * the exponent is supplied as an interval.
 */
template<class T>
AffineMain<T> pow(const AffineMain<T>& x, const Interval &y);

/**
 * \brief Returns an affine form raised to an affine-form exponent.
 *
 * \param x affine form base
 * \param y affine form exponent
 * \return affine result
 *
 * \note The affine dependency structure of the exponent is not preserved;
 * the exponent is evaluated through its interval enclosure.
 */
template<class T>
AffineMain<T> pow(const AffineMain<T>& x, const AffineMain<T>& y);

/**
 * \brief Returns a real value raised to an affine-form exponent.
 *
 * \param d real base
 * \param x affine form exponent
 * \return affine result
 *
 * \note The affine dependency structure of the exponent is evaluated through
 * its interval enclosure.
 */
template<class T>
AffineMain<T> pow(double d, const AffineMain<T>& x);

/**
 * \brief Returns an interval raised to an affine-form exponent.
 *
 * \param x interval base
 * \param y affine form exponent
 * \return affine result
 *
 * \note The affine dependency structure of the exponent is evaluated through
 * its interval enclosure.
 */
template<class T>
AffineMain<T> pow(const Interval& x, const AffineMain<T>& y);

/**
 * \brief Returns the n-th root of an affine form.
 *
 * \param x affine form
 * \param n root degree
 * \return affine result
 */
template<class T>
AffineMain<T> root(const AffineMain<T>&  x, int n);

/**
 * \brief Returns the cosine of an affine form.
 *
 * \param x affine form
 * \return affine result
 */
template<class T>
AffineMain<T> cos(const AffineMain<T>&  x);

/**
 * \brief Returns the sine of an affine form.
 *
 * \param x affine form
 * \return affine result
 */
template<class T>
AffineMain<T> sin(const AffineMain<T>&  x);

/**
 * \brief Returns the tangent of an affine form.
 *
 * \param x affine form
 * \return affine result
 */
template<class T>
AffineMain<T> tan(const AffineMain<T>&  x);

/**
 * \brief Returns the arc-cosine of an affine form.
 *
 * \param x affine form
 * \return affine result
 */
template<class T>
AffineMain<T> acos(const AffineMain<T>&  x);

/**
 * \brief Returns the arc-sine of an affine form.
 *
 * \param x affine form
 * \return affine result
 */
template<class T>
AffineMain<T> asin(const AffineMain<T>&  x);

/**
 * \brief Returns the arc-tangent of an affine form.
 *
 * \param x affine form
 * \return affine result
 */
template<class T>
AffineMain<T> atan(const AffineMain<T>&  x);

/**
 * \brief Returns the hyperbolic cosine of an affine form.
 *
 * \param x affine form
 * \return affine result
 */
template<class T>
AffineMain<T> cosh(const AffineMain<T>&  x);

/**
 * \brief Returns the hyperbolic sine of an affine form.
 *
 * \param x affine form
 * \return affine result
 */
template<class T>
AffineMain<T> sinh(const AffineMain<T>&  x);

/**
 * \brief Returns the hyperbolic tangent of an affine form.
 *
 * \param x affine form
 * \return affine result
 */
template<class T>
AffineMain<T> tanh(const AffineMain<T>&  x);

/**
 * \brief Returns the inverse hyperbolic cosine of an affine form.
 *
 * \param x affine form
 * \return affine result
 */
template<class T>
AffineMain<T> acosh(const AffineMain<T>& x);

/**
 * \brief Returns the inverse hyperbolic sine of an affine form.
 *
 * \param x affine form
 * \return affine result
 */
template<class T>
AffineMain<T> asinh(const AffineMain<T>& x);

/**
 * \brief Returns the inverse hyperbolic tangent of an affine form.
 *
 * \param x affine form
 * \return affine result
 */
template<class T>
AffineMain<T> atanh(const AffineMain<T>& x);


/**
 * \brief Returns the two-argument arctangent.
 *
 * Computes \f$\operatorname{atan2}(y,x)\f$ from two affine forms.
 *
 * \param y affine form representing the y-coordinate
 * \param x affine form representing the x-coordinate
 * \return affine result
 */
template<class T>
AffineMain<T> atan2(const AffineMain<T>& y, const AffineMain<T>& x);


/**
 * \brief Returns the absolute value of an affine form.
 *
 * \param x affine form
 * \return affine result
 */
template<class T>
AffineMain<T> abs(const AffineMain<T>& x);

/**
 * \brief Returns the interval maximum of two arguments.
 *
 * \param x first interval or affine form
 * \param y second interval or affine form
 * \return interval containing the maximum
 *
 * \note The result is an \c Interval and does not preserve affine
 * dependency information.
 */
template<class T>
Interval max(const AffineMain<T>&  x, const AffineMain<T>&  y);
/**
 * \brief Returns the interval maximum of two arguments.
 *
 * \param x first interval or affine form
 * \param y second interval or affine form
 * \return interval containing the maximum
 *
 * \note The result is an \c Interval and does not preserve affine
 * dependency information.
 */
template<class T>
Interval max(const Interval& x, const AffineMain<T>&  y);
/**
 * \brief Returns the interval maximum of two arguments.
 *
 * \param x first interval or affine form
 * \param y second interval or affine form
 * \return interval containing the maximum
 *
 * \note The result is an \c Interval and does not preserve affine
 * dependency information.
 */
template<class T>
Interval max(const AffineMain<T>&  x, const Interval& y);

/**
 * \brief Returns the interval minimum of two arguments.
 *
 * \param x first interval or affine form
 * \param y second interval or affine form
 * \return interval containing the minimum
 *
 * \note The result is an \c Interval and does not preserve affine
 * dependency information.
 */
template<class T>
Interval min(const AffineMain<T>&  x, const AffineMain<T>&  y);
/**
 * \brief Returns the interval minimum of two arguments.
 *
 * \param x first interval or affine form
 * \param y second interval or affine form
 * \return interval containing the minimum
 *
 * \note The result is an \c Interval and does not preserve affine
 * dependency information.
 */
template<class T>
Interval min(const Interval& x, const AffineMain<T>&  y) ;
/**
 * \brief Returns the interval minimum of two arguments.
 *
 * \param x first interval or affine form
 * \param y second interval or affine form
 * \return interval containing the minimum
 *
 * \note The result is an \c Interval and does not preserve affine
 * dependency information.
 */
template<class T>
Interval min(const AffineMain<T>&  x, const Interval& y);


/**
 * \brief Returns the intersection of two interval enclosures.
 *
 * \param x1 first affine form or interval
 * \param x2 second affine form or interval
 * \return intersection as an \c Interval; \c Interval::EMPTY if empty
 */
template<class T>
Interval operator&(const AffineMain<T>&  x1, const AffineMain<T>&  x2);
/**
 * \brief Returns the intersection of two interval enclosures.
 *
 * \param x1 first affine form or interval
 * \param x2 second affine form or interval
 * \return intersection as an \c Interval; \c Interval::EMPTY if empty
 */
template<class T>
Interval operator&(const Interval& x1, const AffineMain<T>&  x2);
/**
 * \brief Returns the intersection of two interval enclosures.
 *
 * \param x1 first affine form or interval
 * \param x2 second affine form or interval
 * \return intersection as an \c Interval; \c Interval::EMPTY if empty
 */
template<class T>
Interval operator&(const AffineMain<T>&  x1, const Interval& x2);


/**
 * \brief Returns the interval hull of the union of two interval enclosures.
 *
 * \param x1 first affine form or interval
 * \param x2 second affine form or interval
 * \return hull of the union as an \c Interval
 */
template<class T>
Interval operator|(const AffineMain<T>&  x1, const AffineMain<T>&  x2);
/**
 * \brief Returns the interval hull of the union of two interval enclosures.
 *
 * \param x1 first affine form or interval
 * \param x2 second affine form or interval
 * \return hull of the union as an \c Interval
 */
template<class T>
Interval operator|(const Interval& x1, const AffineMain<T>&  x2);
/**
 * \brief Returns the interval hull of the union of two interval enclosures.
 *
 * \param x1 first affine form or interval
 * \param x2 second affine form or interval
 * \return hull of the union as an \c Interval
 */
template<class T>
Interval operator|(const AffineMain<T>&  x1, const Interval& x2);


/**
 * \brief Returns the largest integer interval contained in an affine form.
 *
 * \param x affine form
 * \return integer interval
 *
 * \note The result is an \c Interval and does not preserve affine
 * dependency information.
 */
template<class T>
Interval integer(const AffineMain<T>&  x);
/**
 * \brief Returns the interval floor of an affine form.
 *
 * \param x affine form
 * \return \f$\lfloor[x]\rfloor\f$ as an \c Interval
 *
 * \note The result does not preserve affine dependency information.
 */
template<class T>
Interval floor(const AffineMain<T>& x);

/**
 * \brief Returns the interval ceiling of an affine form.
 *
 * \param x affine form
 * \return \f$\lceil[x]\rceil\f$ as an \c Interval
 *
 * \note The result does not preserve affine dependency information.
 */
template<class T>
Interval ceil(const AffineMain<T>& x);

/**
 * \brief Returns the interval sign enclosure of an affine form.
 *
 * \param x affine form
 * \return affine result
 */
template<class T>
AffineMain<T> sign(const AffineMain<T>&  x);



/**
 * \brief Selects between two affine forms according to the sign of a condition.
 *
 * For a condition \p a that is strictly non-positive, the result is \p b;
 * for a condition that is strictly positive, the result is \p c. If both
 * signs are possible, the result encloses both branches.
 *
 * \param a condition, as an interval or affine form
 * \param b value selected when \p a is non-positive
 * \param c value selected when \p a is positive
 * \return affine result enclosing the selected branch or both branches
 */
template<class T>
AffineMain<T> chi(const AffineMain<T>&  a,const AffineMain<T>&  b,const AffineMain<T>&  c);
/**
 * \brief Selects between two affine forms according to the sign of a condition.
 *
 * For a condition \p a that is strictly non-positive, the result is \p b;
 * for a condition that is strictly positive, the result is \p c. If both
 * signs are possible, the result encloses both branches.
 *
 * \param a condition, as an interval or affine form
 * \param b value selected when \p a is non-positive
 * \param c value selected when \p a is positive
 * \return affine result enclosing the selected branch or both branches
 */
template<class T>
AffineMain<T> chi(const Interval&  a,const AffineMain<T>&  b,const AffineMain<T>&  c);
/**
 * \brief Selects between two affine forms according to the sign of a condition.
 *
 * For a condition \p a that is strictly non-positive, the result is \p b;
 * for a condition that is strictly positive, the result is \p c. If both
 * signs are possible, the result encloses both branches.
 *
 * \param a condition, as an interval or affine form
 * \param b value selected when \p a is non-positive
 * \param c value selected when \p a is positive
 * \return affine result enclosing the selected branch or both branches
 */
template<class T>
AffineMain<T> chi(const Interval&  a,const Interval&  b,const AffineMain<T>&  c);
/**
 * \brief Selects between two affine forms according to the sign of a condition.
 *
 * For a condition \p a that is strictly non-positive, the result is \p b;
 * for a condition that is strictly positive, the result is \p c. If both
 * signs are possible, the result encloses both branches.
 *
 * \param a condition, as an interval or affine form
 * \param b value selected when \p a is non-positive
 * \param c value selected when \p a is positive
 * \return affine result enclosing the selected branch or both branches
 */
template<class T>
AffineMain<T> chi(const Interval&  a,const AffineMain<T>&  b,const Interval&  c);
/**
 * \brief Selects between two affine forms according to the sign of a condition.
 *
 * For a condition \p a that is strictly non-positive, the result is \p b;
 * for a condition that is strictly positive, the result is \p c. If both
 * signs are possible, the result encloses both branches.
 *
 * \param a condition, as an interval or affine form
 * \param b value selected when \p a is non-positive
 * \param c value selected when \p a is positive
 * \return affine result enclosing the selected branch or both branches
 */
template<class T>
AffineMain<T> chi(const AffineMain<T>&  a,const Interval&  b,const AffineMain<T>&  c);
/**
 * \brief Selects between two affine forms according to the sign of a condition.
 *
 * For a condition \p a that is strictly non-positive, the result is \p b;
 * for a condition that is strictly positive, the result is \p c. If both
 * signs are possible, the result encloses both branches.
 *
 * \param a condition, as an interval or affine form
 * \param b value selected when \p a is non-positive
 * \param c value selected when \p a is positive
 * \return affine result enclosing the selected branch or both branches
 */
template<class T>
AffineMain<T> chi(const AffineMain<T>&  a,const AffineMain<T>&  b,const Interval&  c);







/*@}*/

/*============================================ inline implementation ============================================ */




template<class T>
inline void AffineMain<T>::change_mode(Affine_Mode tt) {
	switch (tt) {
	case Affine_Mode::AF_Lin_Chebyshev:
	case Affine_Mode::AF_Lin_MinRange:
		mode = tt;
		break;
	default:
		mode = Affine_Mode::AF_Lin_Chebyshev;
		break;
	}
}
template<class T>
inline AffineMain<T>::Affine_Mode AffineMain<T>::get_mode() {
	return mode;
}

template<class T>
inline void AffineMain<T>::compact(){	compact(AF_COMPAC_Tol); }



template<class T>
AffineMain<T>::AffineMain(double d) : AffineMain() {
	*this = Interval(d);
 }

template<class T>
AffineMain<T>::AffineMain(const Interval& x) : AffineMain() {
	*this = x;
 }

template<typename AF>
AffineMain<AF>::operator Interval() const {
    return this->itv();
}


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
bool AffineMain<T>::operator==(double x) const{
    return (this->itv() == x);
}

template<class T>
bool AffineMain<T>::operator!=(double x) const{
    return (this->itv() != x);
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
inline BoolInterval AffineMain<T>::operator<(const Interval& x) const
{
    return this->itv() < x;
}

template<class T>
inline BoolInterval AffineMain<T>::operator>(const Interval& x) const
{
    return this->itv() > x;
}

template<class T>
inline BoolInterval AffineMain<T>::operator<(const AffineMain& x) const
{
    return this->itv() < x.itv();
}

template<class T>
inline BoolInterval AffineMain<T>::operator>(const AffineMain& x) const
{
    return this->itv() > x.itv();
}

template<class T>
inline void AffineMain<T>::set_empty(){
	*this = Interval::empty();
}


template<class T>
inline Index AffineMain<T>::noise_count() const{
	return _n_noise;
}

template<class T>
inline bool AffineMain<T>::is_active() const{
	return (_status==AffineStatus::Active||_status==AffineStatus::Degenerate);
}

template<class T>
inline bool AffineMain<T>::is_empty() const{
	return (_status==AffineStatus::Empty);
}

template<class T>
inline bool AffineMain<T>::is_degenerated() const {
	//return (itv().diam() <	AF_EC);
	return (_status==AffineStatus::Degenerate);
}

template<class T>
inline bool AffineMain<T>::is_unbounded() const{
	return ((_status==AffineStatus::AllReals)||(_status==AffineStatus::UpperUnbounded)||(_status==AffineStatus::LowerUnbounded));
}

template<class T>
inline BoolInterval operator<(const Interval& x, const AffineMain<T>& y) {
    return x < y.itv();
}

template<class T>
inline BoolInterval operator>(const Interval& x, const AffineMain<T>& y){
    return x > y.itv();
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
	} else if (x.is_unbounded() || (!is_active())) {
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

template<class T>
inline AffineMain<T>& AffineMain<T>::operator-=(const AffineMain<T>& x){
	return *this += (AffineMain<T>(x).Aneg());
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
inline AffineMain<T> inv(const AffineMain<T>& x){
	AffineMain<T> out(x);
	out.Ainv(x.itv());
	return out;
}

template<class T>
inline AffineMain<T> sqr(const AffineMain<T>& x){
	AffineMain<T> out(x);
	out.Asqr(x.itv());
	return out;
}

template<class T>
inline AffineMain<T> sqrt(const AffineMain<T>& x){
	AffineMain<T> out(x);
	out.Asqrt(x.itv());
	return out;
}

template<class T>
inline AffineMain<T> exp(const AffineMain<T>& x){
	AffineMain<T> out(x);
	out.Aexp(x.itv());
	return out;
}

template<class T>
inline AffineMain<T> log(const AffineMain<T>& x){
	AffineMain<T> out(x);
	out.Alog(x.itv());
	return out;
}

template<class T>
inline AffineMain<T> pow(const AffineMain<T>& x, int n) {
	AffineMain<T> out(x);
	out.Apow(n,x.itv());
	return out;
}
template<class T>
inline AffineMain<T> pow(const AffineMain<T>& x, double d){
	AffineMain<T> out(x);
	out.Apow(d,x.itv());
	return out;
}
template<class T>
inline AffineMain<T> pow(const AffineMain<T> &x, const Interval &y){
	// return exp(y * log(x));
	AffineMain<T> out(x);
	out.Apow(y,x.itv());
	return out;
}
template<class T>
inline AffineMain<T> pow(const AffineMain<T>& x, const AffineMain<T>& y)
{
  // Deliberate design choice: the affine dependency structure of the
  // exponent y is not preserved. The exponent is evaluated through y.itv()
  // as an independent interval. The affine structure of the base x is kept
  // whenever Apow(const Interval&, const Interval&) can safely preserve it.
  AffineMain<T> out(x);
  out.Apow(y.itv(), x.itv());
  return out;
}
template<class T>
inline AffineMain<T> pow(const Interval& x, const AffineMain<T>& y)
{
  // Deliberate design choice: the affine dependency structure of the
  // exponent y is not preserved. The exponent is evaluated through y.itv()
  // as an independent interval. The affine structure of the base x is kept
  // whenever Apow(const Interval&, const Interval&) can safely preserve it.
  AffineMain<T> out(x);
  out.Apow(y.itv(), x);
  return out;
}

template<class T>
inline AffineMain<T> pow(double x, const AffineMain<T>& y)
{
  // Deliberate design choice: the affine dependency structure of the
  // exponent y is not preserved. The exponent is evaluated through y.itv()
  // as an independent interval. The affine structure of the base x is kept
  // whenever Apow(const Interval&, const Interval&) can safely preserve it.
  AffineMain<T> out(x);
  out.Apow(y.itv(), Interval(x));
  return out;
}

template<class T>
inline AffineMain<T> root(const AffineMain<T>& x, int n) {
	AffineMain<T> out(x);
	out.Aroot(n,x.itv());
	return out;
}

template<class T>
inline AffineMain<T> cos(const AffineMain<T>& x){
	AffineMain<T> out(x);
	out.Acos(x.itv());
	return out;
}

template<class T>
inline AffineMain<T> sin(const AffineMain<T>& x){
	AffineMain<T> out(x);
	out.Asin(x.itv());
	return out;
}

template<class T>
inline AffineMain<T> tan(const AffineMain<T>& x){
	AffineMain<T> out(x);
	out.Atan(x.itv());
	return out;
}

template<class T>
inline AffineMain<T> acos(const AffineMain<T>& x){
	AffineMain<T> out(x);
	out.Aacos(x.itv());
	return out;
}

template<class T>
inline AffineMain<T> asin(const AffineMain<T>& x){
	AffineMain<T> out(x);
	out.Aasin(x.itv());
	return out;
}

template<class T>
inline AffineMain<T> atan(const AffineMain<T>& x){
	AffineMain<T> out(x);
	out.Aatan(x.itv());
	return out;
}

template<class T>
inline AffineMain<T> cosh(const AffineMain<T>& x){
	AffineMain<T> out(x);
	out.Acosh(x.itv());
	return out;
}

template<class T>
inline AffineMain<T> sinh(const AffineMain<T>& x){
	AffineMain<T> out(x);
	out.Asinh(x.itv());
	return out;
}

template<class T>
inline AffineMain<T> tanh(const AffineMain<T>& x){
	AffineMain<T> out(x);
	out.Atanh(x.itv());
	return out;
}

template<class T>
inline AffineMain<T> acosh(const AffineMain<T>& x){
	AffineMain<T> out(x);
	out.Aacosh(x.itv());
	return out;
}

template<class T>
inline AffineMain<T> asinh(const AffineMain<T>& x){
	AffineMain<T> out(x);
	out.Aasinh(x.itv());
	return out;
}

template<class T>
inline AffineMain<T> atanh(const AffineMain<T>& x){
	AffineMain<T> out(x);
	out.Aatanh(x.itv());
	return out;
}

template<class T>
inline AffineMain<T> atan2(const AffineMain<T>& y, const AffineMain<T>& x)
{
    AffineMain<T> res(y);
    res.Aatan2(x, y.itv(), x.itv());
    return res;
}


template<class T>
inline AffineMain<T> abs(const AffineMain<T> &x){
	AffineMain<T> out(x);
	out.Aabs(x.itv());
	return out;
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
inline Interval integer(const AffineMain<T>& x)
{
  return integer(x.itv());
}

template<class T>
inline Interval floor(const AffineMain<T>& x)
{
  return floor(x.itv());
}

template<class T>
inline Interval ceil(const AffineMain<T>& x)
{
  return ceil(x.itv());
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
inline AffineMain<T> chi(const Interval& a,
                         const AffineMain<T>& b,
                         const AffineMain<T>& c)
{
  if (a.is_empty()) {
    AffineMain<T> out;
    out = Interval::empty();
    return out;
  }
  if (a.ub() <= 0.0) {
    return AffineMain<T>(b);
  }
  if (a.lb() > 0.0) {
    return AffineMain<T>(c);
  }

  // Both branches are possible. Returning the interval hull deliberately
  // drops either branch's dependency rather than selecting an arbitrary one.
		Interval b_itv = b.itv();
		Interval c_itv = c.itv();
		if (b_itv.is_subset(c_itv)) {
			return AffineMain<T>(c);
		} else if (c_itv.is_subset(b_itv)) {
			return AffineMain<T>(b);
		} else {
			AffineMain<T> out;
			out = b_itv|c_itv;
			return  out;
		}
}

template<class T>
inline std::ostream& operator<<(std::ostream& os, const AffineMain<T>& x) {
	{
		os << x.itv() << " : ";
		if (x.is_active()) {
			os << x.mid();
			for (int i = 0; i < x.noise_count(); i++) {
				os << " + " << x.noise(i) << " eps_" << i;
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
	switch (mode) {
	case Affine_Mode::AF_Lin_Chebyshev:
		this->Ainv_CH(itv);
		break;
	case Affine_Mode::AF_Lin_MinRange:
		this->Ainv_MR(itv);
		break;
	default:
		*this = 1.0/itv;
		break;
	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Aexp(const Interval& itv){
	switch (mode) {
	case Affine_Mode::AF_Lin_Chebyshev:
		this->Aexp_CH(itv);
		break;
	case Affine_Mode::AF_Lin_MinRange:
		this->Aexp_MR(itv);
		break;
	default:
		*this = exp(itv);
		break;
	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Alog(const Interval& itv){
	switch (mode) {
	case Affine_Mode::AF_Lin_Chebyshev:
		this->Alog_CH(itv);
		break;
	case Affine_Mode::AF_Lin_MinRange:
		this->Alog_MR(itv);
		break;
	default:
		*this = log(itv);
		break;
	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Asqrt(const Interval& itv){
	switch (mode) {
	case Affine_Mode::AF_Lin_Chebyshev:
		this->Asqrt_CH(itv);
		break;
	case Affine_Mode::AF_Lin_MinRange:
		this->Asqrt_MR(itv);
		break;
	default:
		*this = sqrt(itv);
		break;
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
	Interval res_itv = 1.0/(itv);

	// Particular case
	if ((itv.is_unbounded()) || res_itv.is_empty() || res_itv.is_unbounded() || (!is_active()) || (itv.diam()<AF_EC)) {
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
		if (!std::isfinite(alpha) || alpha >= 0.0) {
			*this = res_itv;
			return *this;
		}

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
	Interval res_itv = sqrt(itv2);

	// Particular case
	if (res_itv.is_empty() || res_itv.is_unbounded() || (!is_active()) || (itv2.diam()<AF_EC)) {
		*this = res_itv;
	}  else  {
		// General case
		double alpha, beta, ddelta, t1, t2;
		Interval TEMP1(0.0), TEMP2(0.0), band(0.0);

		alpha = res_itv.diam()/itv2.diam();
		if (!std::isfinite(alpha) || alpha <= 0.0) {
			*this = res_itv;
			return *this;
		}

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
	Interval res_itv = exp(itv);

	// Particular case
	if (res_itv.is_empty() || res_itv.is_unbounded() || (!is_active()) || (itv.diam()<AF_EC)) {
		*this = res_itv;
	}  else  {
		// General case
		double alpha, beta, ddelta, t1, t2;
		Interval TEMP1(0.0), TEMP2(0.0), band(0.0);
		alpha = res_itv.diam()/itv.diam();
		if (!std::isfinite(alpha) || alpha <= 0.0) {
			*this = res_itv;
			return *this;
		}

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
	Interval res_itv = log(itv);

	// Particular case
	if (res_itv.is_empty() || res_itv.is_unbounded() || (!is_active()) || (itv.diam()<AF_EC)) {
		*this = res_itv;
	}  else  {
		// General case
		double alpha, beta, ddelta, t1, t2;
		Interval TEMP1(0.0), TEMP2(0.0), band(0.0);
		alpha = res_itv.diam()/itv.diam();
		if (!std::isfinite(alpha) || alpha <= 0.0) {
			*this = res_itv;
			return *this;
		}

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
inline AffineMain<T>&
AffineMain<T>::Apow(const Interval& y, const Interval& itvx)
{
  const Interval reference = pow(itvx, y);
  if (y.is_empty() || y.is_unbounded() || itvx.is_empty() ||
      itvx.is_unbounded()) {
    *this = reference;
    return *this;
  }

  // An interval singleton that is exactly representable as an int must use
  // the integer-power algorithm. This also preserves valid negative bases.
  if (y.is_degenerated() && std::isfinite(y.lb()) &&
      y.lb() >= static_cast<double>(std::numeric_limits<int>::min()) &&
      y.lb() <= static_cast<double>(std::numeric_limits<int>::max()) &&
      std::trunc(y.lb()) == y.lb()) {
    return this->Apow(static_cast<int>(y.lb()), itvx);
  }

  // exp(y*log(x)) is a real-valued identity only for a strictly positive
  // base. For zero or negative bases, use the Interval implementation,
  // whose domain conventions define the public result safely.
  if (itvx.lb() <= 0.0) {
    *this = reference;
    return *this;
  }

  this->Alog(itvx);
  *this *= y;
  this->Aexp(y*log(itvx));
  return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Acos(const Interval& itv){
	Interval res_itv = cos(itv);

	// Particular case
	if (res_itv.is_empty() || res_itv.is_unbounded() || (!is_active()) || (itv.diam()<AF_EC)) {
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
		//ddelta = ((noise_count()*Interval(TEMP1.rad())) + Interval(TEMP2.rad())).ub();

		// compute the maximal error
		ddelta= 0.0;
		Interval u,nb_period;

		// compute the error at _itv.lb() and _itv.ub() and compute the first point such as f'(u) = alpha

		ddelta = (abs(cos(Interval(itv.lb()))-(alpha*Interval(itv.lb())+beta))).ub();
		t1     = (abs(cos(Interval(itv.ub()))-(alpha*Interval(itv.ub())+beta))).ub();
		if (t1>ddelta)  ddelta= t1;
		
		Interval alpha_domain = Interval(alpha) & Interval(-1.0, 1.0);
		if (alpha_domain.is_empty()) {
			*this = res_itv;
			return *this;
		}
		u = asin(-alpha_domain);
		nb_period = (itv+Interval::half_pi()) / Interval::pi();

		// Large arguments cannot be converted safely to an integer period.
		const double period_limit =
				static_cast<double>(std::numeric_limits<int>::max() - 2);
		if (!std::isfinite(nb_period.lb()) || !std::isfinite(nb_period.ub()) ||
			nb_period.lb() < -period_limit || nb_period.ub() > period_limit) {
			*this = res_itv;
			return *this;
		}

		// evaluate the error at the points such that f'(u) = alpha
		int p1 = static_cast<int>(nb_period.lb()) - 2;
		int p2 = static_cast<int>(nb_period.ub()) + 2;

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
	Interval res_itv = sin(itv);

	// Particular case
	if (res_itv.is_empty() || res_itv.is_unbounded() || (!is_active()) || (itv.diam()<AF_EC)) {
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
		//ddelta = ((noise_count()*Interval(TEMP1.rad())) + Interval(TEMP2.rad())).ub();

		// compute the maximal error
		ddelta= 0.0;
		Interval u,nb_period;

		// compute the error at _itv.lb() and _itv.ub() and compute the first point such as f'(u) = alpha

		ddelta = (abs(sin(Interval(itv.lb()))-(alpha*Interval(itv.lb())+beta))).ub();
		t1     = (abs(sin(Interval(itv.ub()))-(alpha*Interval(itv.ub())+beta))).ub();
		if (t1>ddelta)  ddelta= t1;
		
		Interval alpha_domain =Interval(alpha) & Interval(-1.0, 1.0);
		if (alpha_domain.is_empty()) {
			*this = res_itv;
			return *this;
		}
		u = acos(alpha_domain);
		nb_period = itv / Interval::pi();

		const double period_limit =
				static_cast<double>(std::numeric_limits<int>::max() - 2);
		if (!std::isfinite(nb_period.lb()) || !std::isfinite(nb_period.ub()) ||
			nb_period.lb() < -period_limit || nb_period.ub() > period_limit) {
			*this = res_itv;
			return *this;
		}

		// evaluate the error at the points such that f'(u) = alpha
		int p1 = static_cast<int>(nb_period.lb()) - 2;
		int p2 = static_cast<int>(nb_period.ub()) + 2;

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
	Interval res_itv = tan(itv);

	// Particular case
	if (res_itv.is_empty() || res_itv.is_unbounded() || (!is_active()) || (itv.diam()<AF_EC)) {
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
		//ddelta = ((noise_count()*Interval(TEMP1.rad())) + Interval(TEMP2.rad())).ub();

		// compute the maximal error
		ddelta= 0.0;
		Interval u,nb_period;

		// compute the error at _itv.lb() and _itv.ub() and compute the first point such as f'(u) = alpha

		ddelta = (abs(tan(Interval(itv.lb()))-(alpha*Interval(itv.lb())+beta))).ub();
		t1     = (abs(tan(Interval(itv.ub()))-(alpha*Interval(itv.ub())+beta))).ub();
		if (t1>ddelta)  ddelta= t1;
		if (!std::isfinite(alpha) || alpha < 1.0) {
			*this = res_itv;
			return *this;
		}
		const Interval cosine_domain =(Interval::one()/sqrt(Interval(alpha))) & Interval(0.0, 1.0);
		if (cosine_domain.is_empty()) {
			*this = res_itv;
			return *this;
		}
		u = acos(cosine_domain);
		nb_period = itv / Interval::pi();

		const double period_limit =
				static_cast<double>(std::numeric_limits<int>::max() - 2);
		if (!std::isfinite(nb_period.lb()) || !std::isfinite(nb_period.ub()) ||
			nb_period.lb() < -period_limit || nb_period.ub() > period_limit) {
			*this = res_itv;
			return *this;
		}

		// evaluate the error at the points such that f'(u) = alpha
		int p1 = static_cast<int>(nb_period.lb()) - 2;
		int p2 = static_cast<int>(nb_period.ub()) + 2;

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
	Interval res_itv = acos(itv);
	Interval itv2 = itv & Interval(-1,1);
	// Particular case
	if ( res_itv.is_empty() || res_itv.is_unbounded() || (!is_active()) || (itv2.diam()<AF_EC)) {
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
		//ddelta = ((_n_noise*Interval(TEMP1.rad())) + Interval(TEMP2.rad())).ub();

		// compute the maximal error
		ddelta= 0.0;

		// compute the error at _itv.lb(), _itv.ub() and u such as f'(u) =alpha

		ddelta = (abs(acos(Interval(itv2.lb()))-(alpha*Interval(itv2.lb())+beta))).ub();
		t1     = (abs(acos(Interval(itv2.ub()))-(alpha*Interval(itv2.ub())+beta))).ub();
		if (t1>ddelta)  ddelta= t1;
		// acos'(u)=-1/sqrt(1-u^2) = alpha
		// u = +-sqrt(1-1/(alpha^2))
		Interval TEMP2 = sqrt(1-1/(pow(Interval(alpha),2)));
		Interval critical = TEMP2 & itv2;
		if (!critical.is_empty()) {
			t1 = (abs(acos(critical)-(alpha*critical+beta))).ub();
			if (t1>ddelta)  ddelta= t1;
		}
		critical = (-TEMP2) & itv2;
		if (!critical.is_empty()) {
			t1 = (abs(acos(critical)-(alpha*critical+beta))).ub();
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
	Interval res_itv = asin(itv);
	Interval itv2 =itv & Interval(-1,1);
	// Particular case
	if (res_itv.is_empty() || res_itv.is_unbounded() || (!is_active()) || (itv2.diam()<AF_EC)) {
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
		//ddelta = ((_n_noise*Interval(TEMP1.rad())) + Interval(TEMP2.rad())).ub();

		// compute the maximal error
		ddelta= 0.0;

		// compute the error at _itv.lb(), _itv.ub() and u such as f'(u) =alpha

		ddelta = (abs(asin(Interval(itv2.lb()))-(alpha*Interval(itv2.lb())+beta))).ub();
		t1     = (abs(asin(Interval(itv2.ub()))-(alpha*Interval(itv2.ub())+beta))).ub();
		if (t1>ddelta)  ddelta= t1;
		// asin'(u)=1/sqrt(1-u^2) = alpha
		// u = sqrt(1-1/(alpha^2))
		Interval TEMP2 = sqrt(1-1/sqr(Interval(alpha)));
		if (!((TEMP2 & itv2).is_empty())) {
			t1 = (abs(asin(TEMP2)-(alpha*TEMP2+beta))).ub();
			if (t1>ddelta)  ddelta= t1;
		}
		if (!(((-TEMP2) & itv2).is_empty())) {
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
	Interval res_itv = atan(itv);

	// Particular case
	if (itv.is_unbounded() || res_itv.is_empty() || res_itv.is_unbounded() || (!is_active()) || (itv.diam()<AF_EC)) {
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
		//ddelta = ((_n_noise*Interval(TEMP1.rad())) + Interval(TEMP2.rad())).ub();

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
	Interval res_itv = cosh(itv);

	// Particular case
	if (res_itv.is_empty() || res_itv.is_unbounded() || (!is_active()) || (itv.diam()<AF_EC)) {
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
	Interval res_itv = sinh(itv);

	// Particular case
	if (res_itv.is_empty() || res_itv.is_unbounded() || (!is_active()) || (itv.diam()<AF_EC)) {
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
		//ddelta = ((_n_noise*Interval(TEMP1.rad())) + Interval(TEMP2.rad())).ub();

		// compute the maximal error
		ddelta= 0.0;

		// compute the error at _itv.lb(), _itv.ub() and u such as f'(u) =alpha
		ddelta = (abs(sinh(Interval(itv.lb()))-(alpha*Interval(itv.lb())+beta))).ub();
		t1     = (abs(sinh(Interval(itv.ub()))-(alpha*Interval(itv.ub())+beta))).ub();
		if (t1>ddelta)  ddelta= t1;
		// u = acosh(alpha)
		if (!std::isfinite(alpha) || alpha < 1.0) {
			*this = res_itv;
			return *this;
		}
		TEMP2 = acosh(Interval(alpha));
		Interval critical = TEMP2 & itv;
		if (!critical.is_empty()) {
			t1 = (abs(sinh(critical)-(alpha*critical+beta))).ub();
			if (t1>ddelta)  ddelta= t1;
		}
		critical = (-TEMP2) & itv;
		if (!critical.is_empty()) {
			t1 = (abs(sinh(critical)-(alpha*critical+beta))).ub();
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
	Interval res_itv = tanh(itv);

	// Particular case
	if (itv.is_unbounded() || res_itv.is_empty() || res_itv.is_unbounded() || (!is_active()) || (itv.diam()<AF_EC)) {
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
		//ddelta = ((_n_noise*Interval(TEMP1.rad())) + Interval(TEMP2.rad())).ub();

		// compute the maximal error
		ddelta= 0.0;

		// compute the error at _itv.lb(), _itv.ub() and u such as f'(u) =alpha
		ddelta = (abs(tanh(Interval(itv.lb()))-(alpha*Interval(itv.lb())+beta))).ub();
		t1     = (abs(tanh(Interval(itv.ub()))-(alpha*Interval(itv.ub())+beta))).ub();
		if (t1>ddelta)  ddelta= t1;
		// tanh'(u)=alpha
		// 1/cosh(u)^2= alpha
		// u = +-acosh(1/sqrt(alpha))
		TEMP2 = acosh(1/(sqrt(Interval(alpha))));
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


// acosh is concave on [1,+oo) (acosh''(x) = -x/(x^2-1)^{3/2} < 0), so it
// follows the same Chebyshev pattern as Asqrt_CH: band = [chord intercept
// (min), tangent intercept (max)]. Since acosh is strictly increasing,
// alpha = diam(res_itv)/diam(domain) is always > 0.
template<class T>
inline AffineMain<T>&
AffineMain<T>::Aacosh(const Interval& itv)
{
  const Interval domain = itv & Interval(1.0, oo);
  const Interval res_itv = acosh(domain);

  // The model is built on the real domain intersection. Tests below check
  // that applying it to the original affine form still encloses all values
  // belonging to that valid domain when the input is only partly admissible.
  if (res_itv.is_empty() || res_itv.is_unbounded() ||
      !is_active() || domain.diam() < AF_EC) {
    *this = res_itv;
    return *this;
  }

  const double alpha = res_itv.diam()/domain.diam();
  if (!std::isfinite(alpha) || alpha <= 0.0) {
    *this = res_itv;
    return *this;
  }

  const Interval left = acosh(Interval(domain.lb())) -
                        alpha*Interval(domain.lb());
  const Interval right = acosh(Interval(domain.ub())) -
                         alpha*Interval(domain.ub());
  Interval band = left | right;

  // acosh'(u)=alpha gives u=sqrt(1+1/alpha^2).
  const Interval candidate =
      sqrt(1.0 + 1.0/sqr(Interval(alpha))) & domain;
  if (!candidate.is_empty()) {
    band = band | (acosh(candidate) - alpha*candidate);
  }

  const double beta = band.mid();
  const double ddelta = band.rad();
  *this *= alpha;
  *this += beta;
  this->inflate(ddelta);
  return *this;
}


template<class T>
inline AffineMain<T>&
AffineMain<T>::Aasinh(const Interval& itv)
{
  const Interval res_itv = asinh(itv);
  if (res_itv.is_empty() || res_itv.is_unbounded() || !is_active() ||
      itv.diam() < AF_EC) {
    *this = res_itv;
    return *this;
  }

  const double x0 = 1.0/std::sqrt(2.0);
  const double xb0 = 0.5*(itv.diam()*x0 + itv.lb() + itv.ub());
  const double xb1 = 0.5*(-itv.diam()*x0 + itv.lb() + itv.ub());
  const double fxb0 = std::asinh(xb0);
  const double fxb1 = std::asinh(xb1);
  const double c0 = 0.5*(fxb0 + fxb1);
  const double c1 = x0*(fxb0 - fxb1);
  const double alpha = 2.0*c1/itv.diam();
  const double beta0 = c0-c1*(itv.lb()+itv.ub())/itv.diam();

  if (!std::isfinite(alpha) || alpha <= 0.0 || alpha > 1.0) {
    *this = res_itv;
    return *this;
  }

  Interval residual =
      (asinh(Interval(itv.lb())) -
       (alpha*Interval(itv.lb()) + beta0)) |
      (asinh(Interval(itv.ub())) -
       (alpha*Interval(itv.ub()) + beta0));

  // asinh'(u)=alpha gives u=+-sqrt(1/alpha^2-1).
  const Interval critical =
      sqrt(1.0/sqr(Interval(alpha)) - 1.0);
  for (const Interval& candidate : {critical, -critical}) {
    const Interval point = candidate & itv;
    if (!point.is_empty()) {
      residual = residual | (asinh(point) - (alpha*point + beta0));
    }
  }

  const double beta = beta0 + residual.mid();
  const double ddelta = residual.rad();
  *this *= alpha;
  *this += beta;
  this->inflate(ddelta);
  return *this;
}


// atanh also has an inflection point at x=0 (same S-shape as tanh/asinh),
// but unlike them its derivative 1/(1-x^2) is minimal (=1) at x=0 and grows
// without bound towards the domain edges +-1, so alpha is always >= 1
// (mean value theorem), guaranteeing 1-1/alpha >= 0 below.
template<class T>
inline AffineMain<T>&
AffineMain<T>::Aatanh(const Interval& itv)
{
  const Interval domain = itv & Interval(-1.0, 1.0);
  const Interval res_itv = atanh(itv);
  if (domain != itv || res_itv.is_empty() || res_itv.is_unbounded() ||
      !is_active() || itv.diam() < AF_EC) {
    *this = res_itv;
    return *this;
  }

  const double x0 = 1.0/std::sqrt(2.0);
  const double xb0 = 0.5*(itv.diam()*x0 + itv.lb() + itv.ub());
  const double xb1 = 0.5*(-itv.diam()*x0 + itv.lb() + itv.ub());
  const double fxb0 = std::atanh(xb0);
  const double fxb1 = std::atanh(xb1);
  const double c0 = 0.5*(fxb0 + fxb1);
  const double c1 = x0*(fxb0 - fxb1);
  const double alpha = 2.0*c1/itv.diam();
  const double beta0 = c0-c1*(itv.lb()+itv.ub())/itv.diam();

  if (!std::isfinite(alpha) || alpha < 1.0) {
    *this = res_itv;
    return *this;
  }

  Interval residual =
      (atanh(Interval(itv.lb())) -
       (alpha*Interval(itv.lb()) + beta0)) |
      (atanh(Interval(itv.ub())) -
       (alpha*Interval(itv.ub()) + beta0));

  // atanh'(u)=alpha gives u=+-sqrt(1-1/alpha).
  const Interval critical = sqrt(1.0 - 1.0/Interval(alpha));
  for (const Interval& candidate : {critical, -critical}) {
    const Interval point = candidate & itv;
    if (!point.is_empty()) {
      residual = residual | (atanh(point) - (alpha*point + beta0));
    }
  }

  const double beta = beta0 + residual.mid();
  const double ddelta = residual.rad();
  *this *= alpha;
  *this += beta;
  this->inflate(ddelta);
  return *this;
}


template<class T>
inline AffineMain<T>& AffineMain<T>::Aatan2(const AffineMain<T>& x, const Interval& itvY, const Interval& itvX){
    if(itvY.is_empty() || itvX.is_empty()){
		*this = Interval::empty();
		return *this;
	}
    // We handle the special case x=[0,0] separately
    else if(itvX == Interval::zero())  {
      if(itvY.lb() >= 0) {
        if(itvY.ub() == 0)
          *this = Interval::empty(); // atan2(0,0) is undefined
        else 
          *this = Interval::half_pi();
      }
      else if(itvY.ub() <= 0) 
        *this = (-Interval::half_pi());
      else
        *this = Interval(-1,1)*Interval::half_pi();
	  return *this;
    }
    else if(itvX.lb() >= 0) {
		// Zone sûre (x>0 sur toute la boîte): atan2(y,x) == atan(y/x)
    	// exactement, donc on réutilise la division affine et l'atan affine
    	// déjà validés, sans nouvelle linéarisation.
		*this /= x;                 // *this contient AF[y] -> devient AF[y]/AF[x]
    	this->Aatan(this->itv());
      	return *this; // now, x.ub()>0 -> atan does not give an empty set
	}
    else if(itvX.ub() <= 0) {
      if(itvY.lb() >= 0) {
		*this /= x;                 // *this contient AF[y] -> devient AF[y]/AF[x]
    	this->Aatan(this->itv());
		*this += Interval::pi(); // x.ub()<0
      	return *this;
		}
      else if(itvY.ub() < 0) {
		*this /= x;         
    	this->Aatan(this->itv());
		*this -= Interval::pi(); 
      	return *this;
		}
      else {
        *this = Interval(-1,1)*Interval::pi();
		return *this;
	  }
    }
    else  {     
	// Cas singulier ou coupure de branche potentielle: on ne tente pas de
    // préserver la corrélation affine, on redescend sur l'Interval.
		Interval res_itv = atan2(itvY, itvX);
        *this = res_itv;
        return *this;
    }
  }





template<class T>
inline AffineMain<T>& AffineMain<T>::Aabs(const Interval& itv){
	Interval res_itv = abs(itv);

	// Particular case
	if (0<=itv.lb()) return *this;
	if (itv.ub()<=0) {
		this->Aneg();
		return *this;
	}
	if (res_itv.is_empty() || res_itv.is_unbounded() || (!is_active()) || (itv.diam()<AF_EC)) {
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



namespace detail {

/** Portable integer power for finite scalar interpolation points. */
inline double affine_powi(double base, unsigned int exponent)
{
  double result = 1.0;
  while (exponent != 0U) {
    if ((exponent & 1U) != 0U) {
      result *= base;
    }
    exponent >>= 1U;
    if (exponent != 0U) {
      base *= base;
    }
  }
  return result;
}

} // namespace detail

template<class T>
inline AffineMain<T>& AffineMain<T>::Apow(int n, const Interval& itv) {
	//	std::cout << "in power "<<std::endl;

	// Particular case
	if (itv.is_empty()||itv.is_unbounded()) {
		*this = pow(itv,n);
	} else if (!is_active()) {
		*this = pow(itv,n);
	} else if (itv.diam()< AF_EC) {
		*this = pow(itv,n);
	} else {
		// General Case
		if (n == 0) {
			*this = Interval::one();
		} else if (n == 1)
			return *this;
		else if (n == 2)
			this->Asqr(itv);
		else if (n<0) {
			if (n == std::numeric_limits<int>::min()) {
				*this = pow(itv,n);
			} else {
				this->Apow(-n,itv);
				this->Ainv(pow(itv,-n));
			}

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
			const double f_lb = detail::affine_powi(itv.lb(), static_cast<unsigned int>(n));
			const double f_ub = detail::affine_powi(itv.ub(), static_cast<unsigned int>(n));
			alpha = (f_ub - f_lb)/itv.diam();
			if (!std::isfinite(f_lb) || !std::isfinite(f_ub) ||
					!std::isfinite(alpha)) {
				*this = pow(itv, n);
				return *this;
			}

			TEMP1 = Interval(f_lb) - alpha * Interval(itv.lb());
			TEMP2 = Interval(f_ub) - alpha * Interval(itv.ub());
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
			fxb0 = detail::affine_powi(xb0, static_cast<unsigned int>(n));
			fxb1 = detail::affine_powi(xb1, static_cast<unsigned int>(n));
			c0 = (0.5) * (fxb0 + fxb1);
			c1 = x0 * fxb0 - x0 * fxb1;

			alpha = 2 * c1 / (itv.diam());
			beta  = c0 - c1 * ((itv.lb() + itv.ub()) / (itv.diam()));
			if (!std::isfinite(fxb0) || !std::isfinite(fxb1) ||
					!std::isfinite(alpha) || !std::isfinite(beta)) {
				*this = pow(itv, n);
				return *this;
			}
			//ddelta = ((_n_noise * Interval(TEMP1.rad())) + Interval(TEMP2.rad())).ub(); //

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
	if (!std::isfinite(d)) {
		*this = pow(itv, d);
		return *this;
	}
	if (d >= static_cast<double>(std::numeric_limits<int>::min()) &&
	    d <= static_cast<double>(std::numeric_limits<int>::max()) &&
	    std::trunc(d) == d) {
	  this->Apow(static_cast<int>(d), itv);
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
	else if (n==0)  return *this = Interval::empty();
	else if (n==1)  return *this;
	else if (is_degenerated()) {
		return *this = root(itv,n);
	}
	else if (n<0) {
		if (n == std::numeric_limits<int>::min()) {
			return *this = root(itv, n);
		}
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
		// Odd root across zero. The derivative is singular at zero, but the
		// extrema of r(x)=root(x,n)-alpha*x-beta are the endpoints, zero,
		// and the two points where root'(x)=alpha.
		const Interval root_itv = root(itv, n);
		const double alpha = root_itv.diam()/itv.diam();
		if (!std::isfinite(alpha) || alpha <= 0.0) {
			return *this = root_itv;
		}

		const double beta0 = root_itv.mid() - alpha*itv.mid();
		Interval residual =
				(root(Interval(itv.lb()), n) -
				 (alpha*Interval(itv.lb()) + beta0)) |
				(root(Interval(itv.ub()), n) -
				 (alpha*Interval(itv.ub()) + beta0));
		residual = residual | Interval(-beta0); // residual at x=0

		const double exponent =
				static_cast<double>(n)/static_cast<double>(n-1);
		const Interval critical =
				pow(Interval::one()/(static_cast<double>(n)*alpha), exponent);

		for (const Interval& candidate : {critical, -critical}) {
			const Interval point = candidate & itv;
			if (!point.is_empty()) {
				residual = residual | (root(point, n) - (alpha*point + beta0));
			}
		}

		const double beta = beta0 + residual.mid();
		const double ddelta = residual.rad();
		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		return *this;
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
	if (itv.is_unbounded() || res_itv.is_empty() || res_itv.is_unbounded() || (!is_active()) || (itv.diam()<AF_EC)) {
		*this = res_itv;
	}  else  {  // _actif && b

		double alpha, beta, ddelta, t1, t2;
		Interval dmm(0.0), TEMP1(0.0), TEMP2(0.0), band(0.0), itv2;

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
	if (itv.is_unbounded() || res_itv.is_empty() || res_itv.is_unbounded() || (!is_active()) || (itv.diam()<AF_EC)) {
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
	if (itv.is_unbounded() || res_itv.is_empty() || res_itv.is_unbounded() || (!is_active()) || (itv.diam()<AF_EC)) {
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
  const Interval res_itv = 1.0/itv;
  if (itv.is_unbounded() || res_itv.is_empty() || res_itv.is_unbounded() ||
      !is_active() || itv.diam() < AF_EC) {
    *this = res_itv;
    return *this;
  }

  // On either side of zero, 1/x is decreasing. The endpoint farthest from
  // zero has the derivative of smallest magnitude and gives the MinRange
  // slope. Keep the negative slope instead of collapsing it to zero.
  const double endpoint = (itv.lb() > 0.0) ? itv.ub() : itv.lb();
  const Interval derivative = -1.0/sqr(Interval(endpoint));
  const double alpha = derivative.mid();
  if (!std::isfinite(alpha) || alpha >= 0.0) {
    *this = res_itv;
    return *this;
  }

  Interval band =
      (1.0/Interval(itv.lb()) - alpha*Interval(itv.lb())) |
      (1.0/Interval(itv.ub()) - alpha*Interval(itv.ub()));

  // If rounding moved the stationary point slightly inside the interval,
  // include it explicitly in the residual range.
  const Interval stationary_abs = sqrt(-1.0/Interval(alpha));
  for (const Interval& candidate : {stationary_abs, -stationary_abs}) {
    const Interval point = candidate & itv;
    if (!point.is_empty()) {
      band = band | (1.0/point - alpha*point);
    }
  }

  const double beta = band.mid();
  const double ddelta = band.rad();

  *this *= alpha;
  *this += beta;
  this->inflate(ddelta);
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
