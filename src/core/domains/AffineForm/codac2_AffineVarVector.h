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

typedef AffineVarMainVector<AF_Default> Affine2Variables;
/**
 * \ingroup arithmetic
 *
 * \brief Vector of affine variables.
 *
 * Each component is an \c AffineVarMain sharing a coherent number of
 * noise symbols with the other components.
 */


template<class T=AF_Default>
class AffineVarMainVector : public AffineMainVector<T> {

public:

	/**
	 * \brief Creates a vector of size \p n initialized to unbounded affine variables.
	 *
	 * \pre n > 0
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
	explicit AffineVarMainVector(int n, const Interval& x);

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


	/** \brief Destroys this affine variable vector. */
	virtual ~AffineVarMainVector();


	/**
	 * \brief Sets this vector to the empty set.
	 *
	 * The dimension remains unchanged.
	 */
	void set_empty();

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
	 * \brief Initializes components from an interval vector.
	 *
	 * \param x interval vector of values
	 */
	void init(const IntervalVector& x);


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
	void resize(Index n2);


	/**
	 * \brief Assigns this vector from an interval vector.
	 *
	 * \pre Dimensions of this and x must match.
	 * \note Emptiness is overridden.
	 * \param x interval vector
	 * \return a reference to this
	 */
	AffineVarMainVector<T>& operator=(const IntervalVector& x);

	/**
	 * \brief Returns the interval hull of each affine component.
	 *
	 * \pre (*this) must be nonempty
	 * \return interval vector projection
	 */
	IntervalVector itv() const;


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
AffineVarMainVector<T>::AffineVarMainVector(int n)  :
	AffineMainVector<T>(n)   {
	assert(n>=1);
	for (int i = 0; i < n; i++){
		(*this)[i] = AffineVarMain<T>(n, i, Interval());
	}
}

template<class T>
AffineVarMainVector<T>::AffineVarMainVector(int n, const Interval& x) :
	AffineMainVector<T>(n)   {
	assert(n>=1);
	for (int i = 0; i < n; i++) {
		(*this)[i] = AffineVarMain<T>(n, i, x);
	}
}


template<class T>
AffineVarMainVector<T>::AffineVarMainVector(const AffineMainVector<T>& x):
	AffineMainVector<T>(x.size())   {
	for (int i = 0; i < x.size(); i++){
		(*this)[i] = AffineVarMain<T>(x.size(),i,(x[i]).itv());
	}

}

template<class T>
AffineVarMainVector<T>::AffineVarMainVector(const IntervalVector& x) :
	AffineMainVector<T>(x.size())    {
	for (int i = 0; i < x.size(); i++){
		(*this)[i] = AffineVarMain<T>(x.size(), i, x[i]);
	}
}

template<class T>
AffineVarMainVector<T>::AffineVarMainVector(const Vector& x) :
	AffineMainVector<T>(x.size()) {
	for (int i = 0; i < x.size(); i++){
		(*this)[i] = AffineVarMain<T>(x.size(), i, (Interval(x[i])));
	}
}

template<class T>
void AffineVarMainVector<T>::resize(Index n2) {
	Index n1 =this->size();
	if (n2==n1) return;
	else {
		AffineMainVector<T>::conservativeResize(n2);
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
	AffineMainVector<T>::resize(x.size());
	for (int i=0; i<x.size(); i++){
		(*this)[i]=AffineVarMain<T>(x.size(), i, x[i]);
	}
	return *this;
}


template<class T>
void AffineVarMainVector<T>::init(const Interval& x) {
	for (int i = 0; i < this->size(); i++) {
		(*this)[i] =  x;
	}
}

template<class T>
void AffineVarMainVector<T>::init(const IntervalVector& x) {
	assert(x.size() == this->size());
	for (int i = 0; i < this->size(); i++) {
		(*this)[i] = x[i];
	}
}


template<class T>
inline void AffineVarMainVector<T>::set_empty() {
	(*this)[0] = Interval::empty();
}

template<class T>
inline void AffineVarMainVector<T>::clear() {
	this->init(Interval::zero());
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
