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

	/**
	 * \brief Creates an unbounded affine form (like \c Interval()).
	 *
	 * \note Implemented in codac2_Affine2_fAF2.cpp (for T = AF_fAF2).
	 */
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
	 * \note Implemented in codac2_Affine2_fAF2.cpp (for T = AF_fAF2).
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
	 * \note Implemented in codac2_Affine2_fAF2.cpp (for T = AF_fAF2).
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
	 * \note Implemented in codac2_Affine2_fAF2.cpp (for T = AF_fAF2).
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
	 * \note Implemented in codac2_Affine2_fAF2.cpp (for T = AF_fAF2).
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
	 * \note Implemented in codac2_Affine2_fAF2.cpp (for T = AF_fAF2).
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
	 * \note Implemented in codac2_Affine2_fAF2.cpp (for T = AF_fAF2).
	 */
	const Interval itv() const ;

	/**
	 * \brief Returns the i-th noise coefficient.
	 *
	 * \param i coefficient index
	 * \return coefficient value
	 * \note Implemented in codac2_Affine2_fAF2.cpp (for T = AF_fAF2).
	 */
	double noise(Index i) const;

	/**
	 * \brief Returns the remainder error term.
	 *
	 * \return remainder error
	 * \note Implemented in codac2_Affine2_fAF2.cpp (for T = AF_fAF2).
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
	 * \note Implemented in codac2_Affine2_fAF2.cpp (for T = AF_fAF2).
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
	 * \note Implemented in codac2_Affine2_fAF2.cpp (for T = AF_fAF2).
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
	 * \note Implemented in codac2_Affine2_fAF2.cpp (for T = AF_fAF2).
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
	 * \note Implemented in codac2_Affine2_fAF2.cpp (for T = AF_fAF2).
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
	 * \note Implemented in codac2_Affine2_fAF2.cpp (for T = AF_fAF2).
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
	 * \note Implemented in codac2_Affine2_fAF2.cpp (for T = AF_fAF2).
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
	 * \note Implemented in codac2_Affine2_fAF2.cpp (for T = AF_fAF2).
	 */
		AffineMain<T>& Asqr(const Interval& itv);
	/**
	 * \brief Internal helper that computes the negation of the current affine form.
	 *
	 * \return reference to this affine form
	 * \note Implemented in codac2_Affine2_fAF2.cpp (for T = AF_fAF2).
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
	 * \brief Internal helper that computes the hyperbolic cosine using Chebyshev linearization of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Acosh_CH(const Interval& itv);
	/**
	 * \brief Internal helper that computes the hyperbolic tangent using Chebyshev linearization of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Atanh_CH(const Interval& itv);
	/**
	 * \brief Internal helper that computes the arctangent using Chebyshev linearization of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Aatan_CH(const Interval& itv);
	/**
	 * \brief Internal helper that computes the arcsine using Chebyshev linearization of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Aasin_CH(const Interval& itv);
	/**
	 * \brief Internal helper that computes the arccosine using Chebyshev linearization of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Aacos_CH(const Interval& itv);
	/**
	 * \brief Internal helper that computes the hyperbolic sine using Chebyshev linearization of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Asinh_CH(const Interval& itv);
	/**
	 * \brief Internal helper that computes the inverse hyperbolic sine using Chebyshev linearization of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Aasinh_CH(const Interval& itv);
	/**
	 * \brief Internal helper that computes the inverse hyperbolic cosine using Chebyshev linearization of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Aacosh_CH(const Interval& itv);
	/**
	 * \brief Internal helper that computes the inverse hyperbolic tangent using Chebyshev linearization of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Aatanh_CH(const Interval& itv);
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
	 * \brief Internal helper that computes the hyperbolic cosine using MinRange linearization of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Acosh_MR(const Interval& itv);
	/**
	 * \brief Internal helper that computes the hyperbolic tangent using MinRange linearization of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Atanh_MR(const Interval& itv);
	/**
	 * \brief Internal helper that computes the arctangent using MinRange linearization of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Aatan_MR(const Interval& itv);
	/**
	 * \brief Internal helper that computes the arcsine using MinRange linearization of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Aasin_MR(const Interval& itv);
	/**
	 * \brief Internal helper that computes the arccosine using MinRange linearization of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Aacos_MR(const Interval& itv);
	/**
	 * \brief Internal helper that computes the hyperbolic sine using MinRange linearization of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Asinh_MR(const Interval& itv);
	/**
	 * \brief Internal helper that computes the inverse hyperbolic sine using MinRange linearization of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Aasinh_MR(const Interval& itv);
	/**
	 * \brief Internal helper that computes the inverse hyperbolic cosine using MinRange linearization of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Aacosh_MR(const Interval& itv);
	/**
	 * \brief Internal helper that computes the inverse hyperbolic tangent using MinRange linearization of the current affine form. *
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Aatanh_MR(const Interval& itv);
	/**
	 * \brief Internal helper that computes an integer power using Chebyshev linearization of the current affine form. *
	 * \param n integer exponent, at least 3
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Apow_CH(int n, const Interval& itv);
	/**
	 * \brief Internal helper that computes an integer power using MinRange linearization of the current affine form. *
	 * \param n integer exponent, at least 3
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Apow_MR(int n, const Interval& itv);
	/**
	 * \brief Internal helper that computes an odd integer root using MinRange linearization of the current affine form. *
	 * \param n odd integer index, at least 3
	 * \param itv interval enclosure of the input
	 *
	 * \return reference to this affine form
	 */
		AffineMain<T>&  Aroot_MR(int n, const Interval& itv);
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
	 * \note Implemented in codac2_Affine2_fAF2.cpp (for T = AF_fAF2).
	 */
	AffineMain(Index noise_count, Index var, const Interval& itv);

	/**
	 * \brief Changes the number of affine variables.
	 *
	 * \param n new number of variables
	 * \note Implemented in codac2_Affine2_fAF2.cpp (for T = AF_fAF2).
	 */
	void resize_noise(Index n);

};


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

} // namespace codac2


#include "codac2_AffineMain_impl.h"
#include "codac2_AffineMain_operations_impl.h"




/** \brief atan2(AF[y],AF[x]). */
//Affine2 atan2(const Affine2& y, const Affine2& x);
/** \brief atan2([y],AF[x]). */
//Affine2 atan2(const Interval& y, const Affine2& x);
/** \brief atan2(AF[y],[x]). */
//Affine2 atan2(const Affine2& y, const Interval& x);
/** \brief cosh(AF[x]). */
