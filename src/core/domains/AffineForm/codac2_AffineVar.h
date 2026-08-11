/**
 *  \file codac2_AffineVar.h
 *
 *  Affine variables: affine forms that keep track of a dedicated noise symbol.
 *
 * ----------------------------------------------------------------------------
 *  \date       2026
 *  \author     Jordan Ninin
 *  \copyright  Copyright 2026 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include <cassert>
#include <cstddef>
#include <type_traits>

#include "codac2_AffineMain.h"

// Déclaration anticipée de la fonction interne d'Eigen qui construit
// par défaut chaque élément lors de l'allocation d'un Eigen::Matrix
// dynamique de type non-POD (voir Eigen/src/Core/util/Memory.h, Eigen 5.0.0).
// Signature vérifiée par compilation directe contre les sources d'Eigen 5.0.0 :
// elle retourne T*, pas void.
namespace Eigen { namespace internal {
  template<typename T> T* default_construct_elements_of_array(T* ptr, std::size_t size);
} }



namespace codac2 {

template<class T> class AffineVarMainVector;


template<class T>
struct is_interval_based< AffineVarMain<T> > : std::true_type {};

template<class T>
struct is_ctc< AffineVarMain<T> > : std::false_type {};

template<class T>
struct is_sep< AffineVarMain<T> > : std::false_type {};

template<class T>
struct is_affine_based< AffineVarMain<T> > : std::true_type {};

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
  struct NumTraits<codac2::AffineVarMain<T>>
   : NumTraits<double> // permits to get the epsilon, dummy_precision, lowest, highest functions
  {
    typedef codac2::AffineVarMain<T> Real;
    typedef codac2::AffineVarMain<T> NonInteger;
    typedef codac2::AffineVarMain<T> Nested;
    typedef codac2::AffineVarMain<T> Scalar;
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
  struct ScalarBinaryOpTraits<codac2::AffineVarMain<T>,double,BinOp>
  { typedef codac2::AffineMain<T> ReturnType; };

  template<class T,typename BinOp>
  struct ScalarBinaryOpTraits<double,codac2::AffineVarMain<T>,BinOp>
  { typedef codac2::AffineMain<T> ReturnType; };

  template<class T,typename BinOp>
  struct ScalarBinaryOpTraits<codac2::AffineVarMain<T>,codac2::AffineVarMain<T>,BinOp>
  { typedef codac2::AffineMain<T> ReturnType; };

  template<class T,typename BinOp>
  struct ScalarBinaryOpTraits<codac2::AffineVarMain<T>,codac2::AffineMain<T>,BinOp>
  { typedef codac2::AffineMain<T> ReturnType; };

  template<class T,typename BinOp>
  struct ScalarBinaryOpTraits<codac2::AffineMain<T>,codac2::AffineVarMain<T>,BinOp>
  { typedef codac2::AffineMain<T> ReturnType; };

  template<class T,typename BinOp>
  struct ScalarBinaryOpTraits<codac2::Interval,codac2::AffineVarMain<T>,BinOp>
  { typedef codac2::AffineMain<T> ReturnType; };

  template<class T,typename BinOp>
  struct ScalarBinaryOpTraits<codac2::AffineVarMain<T>,codac2::Interval,BinOp>
  { typedef codac2::AffineMain<T> ReturnType; };

}

namespace codac2 {
/**
 * \ingroup arithmetic
 *
 * \brief Affine variables built on top of \c AF_Default.
 *
 * This class specializes \c AffineMain by associating one variable index
 * (stored in \c _var) to each component.
 */



template<class T>
class AffineVarMain : public AffineMain<T> {


private:
	friend class Eigen::Matrix< AffineVarMain<T>,-1,1>;
	friend class AffineVarMainVector<T>;
	friend AffineVarMain<T>* Eigen::internal::default_construct_elements_of_array<AffineVarMain<T>>(AffineVarMain<T>*, std::size_t); //  <-- pour garder AffineVarMain() en protected

    int _var;

protected:
	// Eigen first creates unbound placeholders. AffineVarMainVector then
	// contextualizes each placeholder through copy assignment. An unbound
	// object must never be assigned an Interval directly.
    AffineVarMain() : AffineMain<T>(), _var(-1) {}

	/**
	 * \brief Creates an affine variable among a set of variables.
	 *
	 * \param size total number of variables
	 * \param var index of the variable represented by this instance
	 * \param itv initial interval value
	 */
	explicit AffineVarMain(int size, int var, const Interval& itv);

public:
	/**
	 * \brief Creates an affine variable by copy.
	 *
	 * \param x affine variable to copy
	 */
    AffineVarMain(const AffineVarMain<T>& x);

	/**
	 * \brief Assigns from another affine variable.
	 *
	 * \param x affine variable to copy
	 * \return a reference to this
	 */
	AffineVarMain& operator=(const AffineVarMain<T>& x) ;

	/**
	 * \brief Sets this affine variable from an interval.
	 *
	 * \param itv interval value assigned to this variable
	 * \return a reference to this
	 */
	AffineVarMain& operator=(const Interval& itv);
	AffineVarMain& operator=(double d);


	// Une variable déclarée ne doit pas être mutée en place : cela casserait la
	// cohérence entre _var et le contexte (taille totale de variables) dans
	// lequel elle a été créée. Toute opération arithmétique doit repasser par
	// une AffineMain<T> ordinaire (voir operator+, operator-, etc. libres).
	AffineVarMain& operator=(const AffineMain<T>&) = delete;

	AffineVarMain& operator+=(double) = delete;
	AffineVarMain& operator+=(const Interval&) = delete;
	AffineVarMain& operator+=(const AffineMain<T>&) = delete;

	AffineVarMain& operator-=(double) = delete;
	AffineVarMain& operator-=(const Interval&) = delete;
	AffineVarMain& operator-=(const AffineMain<T>&) = delete;

	AffineVarMain& operator*=(double) = delete;
	AffineVarMain& operator*=(const Interval&) = delete;
	AffineVarMain& operator*=(const AffineMain<T>&) = delete;

	AffineVarMain& operator/=(double) = delete;
	AffineVarMain& operator/=(const Interval&) = delete;
	AffineVarMain& operator/=(const AffineMain<T>&) = delete;

	template<typename OtherDerived>
	AffineVarMain& operator+=(const Eigen::EigenBase<OtherDerived>&) = delete;
	template<typename OtherDerived>
	AffineVarMain& operator-=(const Eigen::EigenBase<OtherDerived>&) = delete;
	template<typename OtherDerived>
	AffineVarMain& operator*=(const Eigen::EigenBase<OtherDerived>&) = delete;
	template<typename OtherDerived>
	AffineVarMain& operator/=(const Eigen::EigenBase<OtherDerived>&) = delete;

    /**
     * \brief Provides an empty Affine Form
     *
     * \return an empty set
     */
    static AffineVarMain<T> empty();
};

template<class T>
inline AffineVarMain<T> AffineVarMain<T>::empty()
{
	AffineVarMain<T> result(0,-1,Interval::empty());;
	return result;
}


template<class T>
AffineVarMain<T>::AffineVarMain(int size, int var1, const Interval& itv) :
		AffineMain<T>(size, var1, itv),
		_var		(var1) {
	assert(size > 0|| itv.is_empty());
	assert(_var >= 0 || itv.is_empty());
	assert(_var < size);
}



template<class T>
AffineVarMain<T>::AffineVarMain(const AffineVarMain<T>& x) :
		AffineMain<T>(x),
		_var		(x._var) {
}

template<class T>
AffineVarMain<T>& AffineVarMain<T>::operator=(const AffineVarMain<T>& x) {
	if (this == &x) {
		return *this;
	}

	if (_var < 0) {
		// Eigen placeholder: acquire the complete variable identity once.
		AffineMain<T>::operator=(x);
		_var = x._var;
		return *this;
	}

	if (x._var == _var && x.size() == this->size()) {
		// Same affine context and same dedicated noise symbol.
		AffineMain<T>::operator=(x);
		return *this;
	}

	// Assigning a different affine variable must not overwrite the destination
	// identity. Rebuild the destination from the certified interval enclosure,
	// which gives it its own noise symbol and avoids duplicated or out-of-range
	// symbol indices inside an AffineVarMainVector.
	return (*this = x.itv());
}

template<class T>
inline AffineVarMain<T>& AffineVarMain<T>::operator=(double d) {
	*this = Interval(d);
	return *this;
}


//===============================================================================================

} // end namespace codac2


/** \brief atan2(AF[y],AF[x]). */
//Affine2 atan2(const Affine2& y, const Affine2& x);
/** \brief atan2([y],AF[x]). */
//Affine2 atan2(const Interval& y, const Affine2& x);
/** \brief atan2(AF[y],[x]). */
//Affine2 atan2(const Affine2& y, const Interval& x);
/** \brief cosh(AF[x]). */
