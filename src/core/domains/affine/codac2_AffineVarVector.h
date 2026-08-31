/**
 *  \file codac2_AffineVarVector.h
 *
 *  Vector container dedicated to affine variables.
 *
 * ----------------------------------------------------------------------------
 *  \date       2026
 *  \author     Jordan Ninin
 *  \copyright  Copyright 2026 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */
#pragma once

#include "codac2_Interval.h"
#include "codac2_IntervalVector.h"
#include <cassert>
#include <limits>
#include <ostream>
#include <concepts>

#include "codac2_AffineVar.h"
#include "codac2_AffineVector.h"

namespace codac2 {


/**
 * \class AffineVarMainVector
 * \brief Vector of affine variables.
 *
 * Each component is an \c AffineVarMain sharing a coherent number of
 * noise symbols with the other components.
 *
 * \note \b Invariant: position \c i in the vector always owns the dedicated
 *       noise symbol of index \c i (see \c AffineVarMain::noise_index()).
 *       Every factory-style method of this class (the constructors,
 *       \c resize(), \c conservativeResize(), \c operator=(const
 *       IntervalVector&), \c init(const Interval&)) establishes or restores
 *       this invariant by (re)building each component as
 *       ``AffineVarMain<T>(size(), i, itv)``.
 *
 *       This invariant is also what the many deleted mutating operators
 *       (``operator+=``, ``-=``, ``*=``, ``/=``) exist to protect: a
 *       declared affine variable must never be updated in place through an
 *       arithmetic combination, since the result would no longer be a pure
 *       function of its own noise symbol alone.
 *
 *       The remaining opening, plain element assignment through
 *       ``operator[]``/``operator()`` (e.g. ``v[i] = w[j]``), is guarded the
 *       same way at the \c AffineVarMain level: \c AffineVarMain::operator=
 *       only adopts the source's noise symbol identity while the
 *       destination is still an unbound Eigen placeholder (freshly
 *       (re)allocated, not yet contextualized); once a component has an
 *       established identity, assigning it another \c AffineVarMain (same
 *       vector or not, same size or not) absorbs only that source's
 *       interval enclosure, rebuilt on the destination's own noise symbol.
 *       Without this guard, ``v[3] = v[7]`` would alias positions 3 and 7
 *       onto the same noise symbol, making them silently -- and
 *       incorrectly -- perfectly correlated in every later computation
 *       (e.g. ``v[3] - v[7]`` would wrongly collapse to a single point
 *       instead of enclosing the actual range of two independent
 *       components). See the regression test
 *       "AffineVarMain : assigning between two components of the same
 *       vector does not alias their noise symbols" in
 *       codac2_tests_AffineVariables.cpp.
 */


template<class T>
class AffineVarMainVector :public Eigen::Matrix<AffineVarMain<T>, -1, 1> {


public:

	using Base = Eigen::Matrix<AffineVarMain<T>, -1, 1>;

	// Expose the selected Eigen arithmetic operators for affine-variable
	// vectors. Constructors and assignment operators are intentionally not
	// inherited from Eigen: some generic Eigen constructors perform scalar
	// initialization that is not compatible with AffineVarMain<T>.
	using Base::operator+;
	using Base::operator-;
	using Base::operator*;
	using Base::operator/;

	/**
	 * \brief Creates a vector of size \p n.
	 *
	 * Each component is set to the default (unbounded) affine form, i.e.
	 * the equivalent of ``Interval()`` (]-oo,+oo[).
	 *
	 * \note The generic ``Eigen::Matrix(int n)`` constructor is deliberately
	 * not inherited. It may perform scalar initialization through
	 * ``init(0.)``, which is not compatible with ``AffineVarMain<T>``.
	 * This constructor therefore provides the intended initialization
	 * semantics explicitly.
	 *
	 * \param n vector size
	 */
	explicit AffineVarMainVector(Index n);

	/**
	 * \brief Creates a vector from an interval vector.
	 *
	 * \param x interval vector used for component-wise initialization
	 */
	explicit AffineVarMainVector(const IntervalVector& x);

	/**
	 * \brief Creates a degenerate affine vector from a real-valued vector.
	 *
	 * \param x real-valued vector
	 */
	explicit AffineVarMainVector(const Vector& x);

	/**
	 * \brief Resizes this vector, rebuilding common components from their current interval enclosure.
	 *
	 * \note Existing affine dependency (noise-symbol correlation) is not preserved:
	 *       - components common to both sizes are recreated as independent affine
	 *         variables from their current interval enclosure (\c itv()), empty
	 *         components remaining empty
	 *       - newly created components (when growing) are initialized to \f$[-\infty,+\infty]\f$
	 *       This operation has the strong exception guarantee.
	 *
	 * \param n2 new non-negative vector size
	 */
	void conservativeResize(Index n2);

	/**
	 * \brief Resizes this vector, resetting every component to an independent unbounded affine variable.
	 *
	 * \note Unlike conservativeResize(), existing interval enclosures are discarded: every
	 *       component, including those common to both sizes, is reset to \f$[-\infty,+\infty]\f$.
	 *
	 * \param n2 new non-negative vector size
	 */
	void resize(Index n2);

	/**
	 * \brief Assigns this vector from an interval vector.
	 *
	 * The destination is resized when necessary. Empty components and all
	 * other component values are copied exactly. The operation has the strong
	 * exception guarantee.
	 * \param x interval vector
	 * \return a reference to this
	 */
	AffineVarMainVector<T>& operator=(const IntervalVector& x);


	/**
	 * \brief Broadcasts a single interval to every component.
	 *
	 * \note This function is used for template purposes.
	 *
	 * This overload accepts an ``Interval`` directly and
	 * assigns it component-wise through ``AffineVarMain::operator=(const Interval&)``,
	 * which preserves each component's dedicated noise symbol (\c _var).
	 * The generic scalar-broadcast ``init`` cannot be used here: it would
	 * require an implicit ``Interval -> AffineVarMain<T>`` conversion, which is
	 * intentionally not provided (see the assignment-asymmetry regression test
	 * in codac2_tests_AffineForm_coverage.cpp).
	 *
	 * \param x interval value assigned to every component
	 * \return a reference to this
	 */
	AffineVarMainVector<T>& init(const Interval& x);

	/** \brief Disabled: broadcasting an ``AffineVarMain<T>`` would require duplicating its noise symbol across components. */
	AffineVarMainVector<T>& init(const AffineVarMain<T>& x)= delete;
	/** \brief Disabled: broadcasting an ``AffineMain<T>`` would require an implicit conversion to ``AffineVarMain<T>`` that is intentionally not provided. */
	AffineVarMainVector<T>& init(const AffineMain<T>& x)= delete;

	/** \brief Returns \f$-*\mathrm{this}\f$. */
	AffineMainVector<T> operator-() const;


	/**
	 * \brief Returns \f$(*\mathrm{this}) \cdot \mathrm{rhs}\f$.
	 *
	 * The ``AffineVarMainVector`` is internally converted to an
	 * ``AffineMainVector``. The resulting Eigen object contains ``AffineMain<T>``
	 * coefficients and preserves the dimensions of the matrix product.
	 *
	 * \param rhs matrix expression to multiply this vector by
	 * \return the resulting affine matrix expression, evaluated eagerly
	 */
	template<typename OtherDerived>
	requires requires(const AffineMainVector<T>& lhs,const Eigen::MatrixBase<OtherDerived>&  rhs)
			{	lhs * rhs.derived();	}
	inline auto operator*( const Eigen::MatrixBase<OtherDerived>& rhs)
	{
	    assert(this->cols() == rhs.rows());
	    const AffineMainVector<T> converted_lhs(*this);
	    return (converted_lhs * rhs.derived()).eval();
	}


	/**
	 * \name Disabled mutating operators
	 *
	 * A vector of affine variables must not be mutated component-wise in a
	 * way that bypasses \c AffineVarMain's noise-symbol bookkeeping. Any
	 * arithmetic result must be captured in a plain ``AffineMainVector<T>``
	 * instead (see the free operators and \c operator- above).
	 * @{
	 */
	AffineVarMainVector& operator=(const AffineMainVector<T>&) = delete;

	AffineVarMainVector& operator+=(const Vector&) = delete;
	AffineVarMainVector& operator+=(const IntervalVector&) = delete;
	AffineVarMainVector& operator+=(const AffineMainVector<T>&) = delete;
	AffineVarMainVector& operator+=(const AffineVarMainVector<T>&) = delete;

	AffineVarMainVector& operator-=(const Vector&) = delete;
	AffineVarMainVector& operator-=(const IntervalVector&) = delete;
	AffineVarMainVector& operator-=(const AffineMainVector<T>&) = delete;
	AffineVarMainVector& operator-=(const AffineVarMainVector<T>&) = delete;

	AffineVarMainVector& operator*=(double) = delete;
	AffineVarMainVector& operator*=(const Interval&) = delete;
	AffineVarMainVector& operator*=(const AffineMain<T>&) = delete;
	AffineVarMainVector& operator*=(const AffineVarMain<T>&) = delete;

	AffineVarMainVector& operator/=(double) = delete;
	AffineVarMainVector& operator/=(const Interval&) = delete;
	AffineVarMainVector& operator/=(const AffineMain<T>&) = delete;
	AffineVarMainVector& operator/=(const AffineVarMain<T>&) = delete;

	template<typename OtherDerived>
	AffineVarMainVector& operator+=(const Eigen::EigenBase<OtherDerived>&) = delete;
	template<typename OtherDerived>
	AffineVarMainVector& operator-=(const Eigen::EigenBase<OtherDerived>&) = delete;
	template<typename OtherDerived>
	AffineVarMainVector& operator*=(const Eigen::EigenBase<OtherDerived>&) = delete;
	template<typename OtherDerived>
	AffineVarMainVector& operator/=(const Eigen::EigenBase<OtherDerived>&) = delete;
	/** @} */

};


/**
 * \brief Returns \f$\mathrm{lhs} \cdot \mathrm{rhs}\f$.
 *
 * The ``AffineVarMainVector`` is internally converted to an
 * ``AffineMainVector``. The resulting Eigen object contains ``AffineMain<T>``
 * coefficients and preserves the dimensions of the matrix product.
 *
 * \param lhs matrix expression on the left-hand side
 * \param rhs vector of affine variables on the right-hand side
 * \return the resulting affine matrix expression, evaluated eagerly
 */
template<typename OtherDerived, class T>
requires  requires(const Eigen::MatrixBase<OtherDerived>& lhs,const AffineMainVector<T>& rhs)
		{	lhs.derived() * rhs;	}
inline auto operator*(
    const Eigen::MatrixBase<OtherDerived>& lhs,
    const AffineVarMainVector<T>& rhs)
{
    assert(lhs.cols() == rhs.rows());
    const AffineMainVector<T> converted_rhs(rhs);
    return (lhs.derived() * converted_rhs).eval();
}




template<class T>
AffineVarMainVector<T>::AffineVarMainVector(Index n) : Base() {
	assert(n >= 0);
	Base::resize(n);
	for (Index i = 0; i < n; ++i) {
		(*this)[i] = AffineVarMain<T>(n, i, Interval());
	}
}


template<class T>
AffineVarMainVector<T>::AffineVarMainVector(const IntervalVector& x) : Base()  {
	assert(x.size() <= static_cast<Eigen::Index>(std::numeric_limits<int>::max()));
	 Base::resize(x.size());
	for (Index i = 0; i < x.size(); i++){
		(*this)[i] = AffineVarMain<T>(x.size(), i, x[i]);
	}
}

template<class T>
AffineVarMainVector<T>::AffineVarMainVector(const Vector& x) : Base() {
	assert(x.size() <= static_cast<Eigen::Index>(std::numeric_limits<int>::max()));
	Base::resize(x.size());
	for (Index i = 0; i < x.size(); i++){
		(*this)[i] = AffineVarMain<T>(x.size(), i, Interval(x[i]));
	}
}

template<class T>
void AffineVarMainVector<T>::conservativeResize(Index n2) {
	assert(n2 >= 0);
	assert(n2 <= static_cast<codac2::Index>(std::numeric_limits<int>::max()));
	const Index n1 = this->size();
    if(n2 == 0)   {
        Base::resize(0);
        return;
    }
	if (n2!=n1) {
		IntervalVector old_values(n1);

		for(Index i = 0; i < n1; ++i)
			old_values[i] = (*this)[i].itv();

		Base::resize(n2);
		Index i=0;
		for (; i<n1 && i<n2; i++){
			(*this)[i]=AffineVarMain<T>(n2, i,old_values[i]);
		}
		for (; i<n2; i++) {
			(*this)[i]=AffineVarMain<T>(n2, i, Interval());
		}
	}
}

template<class T>
void AffineVarMainVector<T>::resize(Index n2) {
	assert(n2 >= 0);
	assert(n2 <= static_cast<codac2::Index>(std::numeric_limits<int>::max()));
	const Index n1 = this->size();
    if(n2 == 0)   {
        Base::resize(0);
        return;
    }
	if (n2!=n1) {

		Base::resize(n2);
		for (Index i=0; i<n2; i++) {
			(*this)[i]=AffineVarMain<T>(n2, i, Interval());
		}
	}
}

template<class T>
AffineVarMainVector<T>& AffineVarMainVector<T>::operator=(const IntervalVector& x)  {
	if (this->size()!=x.size()){
		Base::resize(x.size());
	}
	for (Index i=0; i<x.size(); i++){
		(*this)[i]=AffineVarMain<T>(x.size(), i, x[i]);
	}
	return *this;
}

template<class T>
AffineVarMainVector<T>& AffineVarMainVector<T>::init(const Interval& x)
{
    for (Index i = 0; i < this->size(); ++i)
        (*this)[i] = AffineVarMain<T>(this->size(), i, x);
    return *this;
}

template<class T>
AffineMainVector<T> AffineVarMainVector<T>::operator-() const
{
  AffineMainVector<T> result(this->size());
  for(Index i = 0; i < this->size(); ++i)
    result[i] = (-((*this)[i]));

  return result;
}


/**
 * \brief Streams out x
 *
 * \param os the stream to be updated
 * \param x the vector of affine variables to stream out
 * \return a reference to the updated stream
 */
template<class T>
inline std::ostream& operator<<(std::ostream& os, const AffineVarMainVector<T>& x)
{
	if(x.is_empty())
		return os << "[ empty " << x.size() << "d box ]";

	else
	{
		os << x.format(codac_vector_fmt());
		return os;
	}

}

//===============================================================================================

}


/** \brief atan2(AF[y],AF[x]). */
//Affine2 atan2(const Affine2& y, const Affine2& x);
/** \brief atan2([y],AF[x]). */
//Affine2 atan2(const Interval& y, const Affine2& x);
/** \brief atan2(AF[y],[x]). */
//Affine2 atan2(const Affine2& y, const Interval& x);
/** \brief cosh(AF[x]). */
