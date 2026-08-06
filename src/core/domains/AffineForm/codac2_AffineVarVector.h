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
#include <ostream>

#include "codac2_AffineVar.h"
#include "codac2_AffineVector.h"

namespace codac2 {

template<class T>  class AffineVarMainVector;

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
	using Base::Base;
	using Base::operator=;

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
	explicit AffineVarMainVector(int n);

	/**
	 * \brief Creates a vector of size \p n with all components initialized from \p x.
	 *
	 * \pre n > 0
	 * \param n vector size
	 * \param x interval used to initialize each component
	 */
//	explicit AffineVarMainVector(int n, const Interval& x);

	/**
	 * \brief  Create \a n AffineVarMainVector of dimension \a n with
	 * all the components initialized to \a x.
	 * \pre n>0
	 */
//	explicit AffineVarMainVector(int n, const AffineMain<T>& x);

	/**
	 * \brief Creates a copy from an affine vector.
	 *
	 * \param x affine vector to copy
	 */
	AffineVarMainVector(const AffineMainVector<T>& x);

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
	 * \brief Returns the interval hull of each affine component.
	 *
	 * Build an \c IntervalVector with the same size as \c *this, where
	 * the i-th entry is \c (*this)[i].itv().
	 * \pre (*this) must be nonempty
	 */
	IntervalVector itv() const;


	/**
	 * \brief Sets all elements to zero, even if this is empty.
	 *
	 * \note Emptiness is overridden.
	 */
	void clear();

	/**
	 * \brief Initializes all components with the same interval value.
	 *
	 * \note Emptiness is overridden.
	 * \param x interval value
	 */
	void init(const Interval& x);

	/**
	 * \brief Resizes this affine variable vector.
	 *
	 * If the size is increased, the existing components are not
	 * modified and the new ones are set to (-inf,+inf), even if
	 * (*this) is the empty Interval (however, in this case, the status of
	 * (*this) remains "empty").
	 *
	 * \param n2 new vector size
	 */
	void resize(codac2::Index n2);


	/**
	 * \brief Assigns this vector from an interval vector.
	 *
	 * \pre Dimensions of this and x must match.
	 * \note Emptiness is overridden.
	 * \param x interval vector
	 * \return a reference to this
	 */
	AffineVarMainVector<T>& operator=(const IntervalVector& x);



private:
/*	void put(int start_index, const AffineMainVector<T>& subvec) {codac2_error(" AffineVarMainVector : operator put non valid");};
	AffineVarMainVector& operator=(const AffineMainVector<T>& x) {codac2_error(" AffineVarMainVector : operator= non valid");};
	AffineVarMainVector& operator+=(const Vector& x2) {codac2_error(" AffineVarMainVector : operator+= non valid");};
	AffineVarMainVector& operator+=(const IntervalVector& x2) {codac2_error(" AffineVarMainVector : operator+= non valid");};
	AffineVarMainVector& operator+=(const AffineMainVector<T>& x2) {codac2_error(" AffineVarMainVector : operator+= non valid");};
	AffineVarMainVector& operator+=(const AffineVarMainVector<T>& x2) {codac2_error(" AffineVarMainVector : operator+= non valid");};
	AffineVarMainVector& operator-=(const Vector& x2) {codac2_error(" AffineVarMainVector : operator-= non valid");};
	AffineVarMainVector& operator-=(const IntervalVector& x2) {codac2_error(" AffineVarMainVector : operator-= non valid");};
	AffineVarMainVector& operator-=(const AffineMainVector<T>& x2) {codac2_error(" AffineVarMainVector : operator-= non valid");};
	AffineVarMainVector& operator-=(const AffineVarMainVector<T>& x2) {codac2_error(" AffineVarMainVector : operator-= non valid");};
	AffineVarMainVector& operator*=(double d) {codac2_error(" AffineVarMainVector : operator*= non valid");};
	AffineVarMainVector& operator*=(const Interval& x1) {codac2_error(" AffineVarMainVector : operator*= non valid");};
	AffineVarMainVector& operator*=(const AffineMain<T>& x1) {codac2_error(" AffineVarMainVector : operator*= non valid");};
	AffineVarMainVector& operator*=(const AffineVarMain<T>& x1) {codac2_error(" AffineVarMainVector : operator*= non valid");};
	AffineVarMainVector& operator/=(double d) {codac2_error(" AffineVarMainVector : operator/= non valid");};
	AffineVarMainVector& operator/=(const Interval& x1) {codac2_error(" AffineVarMainVector : operator/= non valid");};
	AffineVarMainVector& operator/=(const AffineMain<T>& x1) {codac2_error(" AffineVarMainVector : operator/= non valid");};
	AffineVarMainVector& operator/=(const AffineVarMain<T>& x1) {codac2_error(" AffineVarMainVector : operator/= non valid");};
*/
};



template<class T>
AffineVarMainVector<T>::AffineVarMainVector(int n)  : Base(n,1)
	 {
	assert(n>=1);
	for (int i = 0; i < n; i++){
		(*this)[i] = AffineVarMain<T>(n, i, Interval());
	}
}

//template<class T>
//AffineVarMainVector<T>::AffineVarMainVector(int n, const Interval& x) : Base(n,1)
//	 {
//	assert(n>=1);
//	for (int i = 0; i < n; i++) {
//		(*this)[i] = AffineVarMain<T>(n, i, x);
//	}
//}


template<class T>
AffineVarMainVector<T>::AffineVarMainVector(const AffineMainVector<T>& x): Base(x.size(),1)
	   {
	for (int i = 0; i < x.size(); i++){
		(*this)[i] = AffineVarMain<T>(x.size(),i,(x[i]).itv());
	}
}

template<class T>
AffineVarMainVector<T>::AffineVarMainVector(const IntervalVector& x) : Base(x.size(),1)
	   {
	for (int i = 0; i < x.size(); i++){
		(*this)[i] = AffineVarMain<T>(x.size(), i, x[i]);
	}
}

template<class T>
AffineVarMainVector<T>::AffineVarMainVector(const Vector& x) : Base(x.size(),1)
	{
	for (int i = 0; i < x.size(); i++){
		double val = x[i];
		(*this)[i] = AffineVarMain<T>(x.size(), i, Interval(val));
	}
}

template<class T>
void AffineVarMainVector<T>::resize(codac2::Index n2) {
	Index n1 =this->size();
	if (n2!=n1) {
		AffineVarMainVector<T>::conservativeResize(n2);
		int i=0;
		for (; i<n1 && i<n2; i++){
			(*this)[i]=AffineVarMain<T>(n2, i, ((*this)[i]).itv());
		}
		for (; i<n2; i++) {
			(*this)[i]=AffineVarMain<T>(n2, i, Interval());
		}
	}
}


template<class T>
AffineVarMainVector<T>& AffineVarMainVector<T>::operator=(const IntervalVector& x)  {
	if (this->size()!=x.size()){
		AffineVarMainVector<T>::resize(x.size());
	}
	for (int i=0; i<x.size(); i++){
		(*this)[i]=AffineVarMain<T>(x.size(), i, x[i]);
	}
	return *this;
}


template<class T>
void AffineVarMainVector<T>::init(const Interval& x) {
	for (int i = 0; i < this->size(); i++) {
		(*this)[i] =  AffineVarMain<T>(this->size(), i, x);
	}
}


template<class T>
inline void AffineVarMainVector<T>::clear() {
	this->init(Interval::zero());
}



template<class T>
IntervalVector AffineVarMainVector<T>::itv() const {
	IntervalVector tmp(this->size());
	for (int i = 0; i < this->size(); i++) {
		tmp[i] = (*this)[i].itv();
	}
	return tmp;
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
