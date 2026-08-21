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
 * \ingroup arithmetic
 *
 * \brief Vector of affine variables.
 *
 * Each component is an \c AffineVarMain sharing a coherent number of
 * noise symbols with the other components.
 */


template<class T>
class AffineVarMainVector :public Eigen::Matrix<AffineVarMain<T>, -1, 1> {


public:

	using Base = Eigen::Matrix<AffineVarMain<T>, -1, 1>;

	// Inherit the constructors and the assignment operators of the Eigen base
	// class (dynamic size constructor, expression assignment, etc.).
	// Without these "using" declarations, the derived class only gets the
	// implicitly-generated default/copy/move constructors and assignment
	// operators, and every other Eigen::Matrix constructor/assignment
	// (e.g. AffineMainVector<T>(n) or v = v1+v2) would fail to compile.
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
	 * \note This constructor hides (replaces) the generic
	 * ``Eigen::Matrix(int n)`` constructor inherited above via
	 * ``using Base::Base``. That generic constructor always tries to
	 * zero-initialize non-interval scalar types through ``init(0.)``,
	 * which does not compile for ``AffineVarMain<T>`` since there is no
	 * implicit conversion from ``double`` to ``AffineVarMain<T>``.
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
	 * \brief Resizes the vector and reconstructs all components as independent
	 * affine variables from their current interval enclosures.
	 * Existing affine dependency information is not preserved.
	 *
	 * Existing components are reconstructed from their interval enclosures;
	 * newly created components are initialized to (-inf,+inf). Empty
	 * components remain empty, consistently with IntervalVector semantics.
	 * The operation has the strong exception guarantee.
	 *
	 * \param n2 new non-negative vector size
	 */
	void conservativeResize(Index n2);

	/**
	 * \brief Resizes the vector and reconstructs all components as independent
	 * affine variables from their current interval enclosures.
	 * Existing affine dependency information is not preserved.
	 *
	 * All components are initialized to (-inf,+inf).
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
	 * Unlike the generic Eigen::Matrix::init(const Scalar&) inherited via
	 * ``using Base::Base``, this overload accepts an ``Interval`` directly and
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
	AffineVarMainVector<T>& init(const AffineVarMain<T>& x)= delete;
	AffineVarMainVector<T>& init(const AffineMain<T>& x)= delete;

	/** \brief Returns \f$-*\mathrm{this}\f$. */
	AffineMainVector<T> operator-() const;


	/**
	 * Product of a matrix expression by a vector of affine variables.
	 *
	 * The AffineVarMainVector is internally converted to an
	 * AffineMainVector. The resulting Eigen object contains AffineMain<T>
	 * coefficients and preserves the dimensions of the matrix product.
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

};


/**
 * Product of a matrix expression by a vector of affine variables.
 *
 * The AffineVarMainVector is internally converted to an
 * AffineMainVector. The resulting Eigen object contains AffineMain<T>
 * coefficients and preserves the dimensions of the matrix product.
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
 * \brief Stream output operator for ``IntervalVector`` objects.
 *
 * \param os The output stream to write to.
 * \param x The interval vector whose contents are to be printed.
 * \return A reference to the modified output stream.
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
