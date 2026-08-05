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

#include <cmath>
#include <ostream>
#include <cassert>

#include "codac2_IntervalVector.h"
#include "codac2_AffineMain.h"
#include "codac2_AffineVector.h"

namespace codac2 {

template<class T> class AffineVarMainVector;

/**
 * \ingroup arithmetic
 *
 * \brief Affine variables built on top of \c AF_Default.
 *
 * This class specializes \c AffineMain by associating one variable index
 * (stored in \c _var) to each component.
 */


typedef AffineVarMainVector<AF_Default> Affine2Variables;



template<class T=AF_Default>
class AffineVarMain : public AffineMain<T> {


public:
	/**
	 * \brief Sets this affine variable from an interval.
	 *
	 * \param itv interval value assigned to this variable
	 * \return a reference to this
	 */
	AffineVarMain& operator=(const Interval& itv);


private:
	friend class AffineVarMainVector<T>;

    int _var;

 //  necesary with eigen inAffineVarMainVector
    AffineVarMain() : _var(-1) {};

	/**
	 * \brief Creates an affine variable by copy.
	 *
	 * \param x affine variable to copy
	 */
    AffineVarMain(const AffineVarMain<T>& x);

	/**
	 * \brief Creates an affine variable among a set of variables.
	 *
	 * \param size total number of variables
	 * \param var index of the variable represented by this instance
	 * \param itv initial interval value
	 */
	explicit AffineVarMain(int size, int var, const Interval& itv);

	/**
	 * \brief Assigns from another affine variable.
	 *
	 * \param x affine variable to copy
	 * \return a reference to this
	 */
	AffineVarMain& operator=(const AffineVarMain<T>& x) ;


/*	AffineVarMain& operator+=(const Vector& x2) 		{codac2_error(" AffineVarMain : operator+= non valid");};
	AffineVarMain& operator+=(const Interval& x2) 		{codac2_error(" AffineVarMain : operator+= non valid");};
	AffineVarMain& operator+=(const AffineMain<T>& x2)	{codac2_error(" AffineVarMain : operator+= non valid");};
	AffineVarMain& operator-=(const Vector& x2) 		{codac2_error(" AffineVarMain : operator-= non valid");};
	AffineVarMain& operator-=(const Interval& x2) 		{codac2_error(" AffineVarMain : operator-= non valid");};
	AffineVarMain& operator-=(const AffineMain<T>& x2) 	{codac2_error(" AffineVarMain : operator-= non valid");};
	AffineVarMain& operator*=(double d) 				{codac2_error(" AffineVarMain : operator*= non valid");};
	AffineVarMain& operator*=(const Interval& x1) 		{codac2_error(" AffineVarMain : operator*= non valid");};
	AffineVarMain& operator*=(const AffineMain<T>& x1) 	{codac2_error(" AffineVarMain : operator*= non valid");};

	AffineVarMain&  Asqr(const Interval& itv)  	{codac2_error(" AffineVarMain : operator Asqr non valid");};
	AffineVarMain&  Aneg()  					{codac2_error(" AffineVarMain : operator Aneg non valid");};
	AffineVarMain&  Ainv(const Interval& itv)  	{codac2_error(" AffineVarMain : operator Ainv non valid");};
	AffineVarMain&  Asqrt(const Interval& itv) 	{codac2_error(" AffineVarMain : operator non Asqrt valid");};
	AffineVarMain&  Aexp(const Interval& itv)  	{codac2_error(" AffineVarMain : operator Aexp non valid");};
	AffineVarMain&  Alog(const Interval& itv)  	{codac2_error(" AffineVarMain : operator Alog non valid");};
	AffineVarMain&  Apow(int n, const Interval& itv)  	{codac2_error(" AffineVarMain : operator Apow non valid");};
	AffineVarMain&  Apow(double d, const Interval& itv)	{codac2_error(" AffineVarMain : operator Apow non valid");};
	AffineVarMain&  Apow(const Interval &y, const Interval& itvx)  {codac2_error(" AffineVarMain : operator Apow non valid");};
	AffineVarMain&  Aroot(int n, const Interval& itv)  	{codac2_error(" AffineVarMain : operator Aroot non valid");};
	AffineVarMain&  Acos(const Interval& itv)  	{codac2_error(" AffineVarMain : operator Acos non valid");};
	AffineVarMain&  Asin(const Interval& itv)  	{codac2_error(" AffineVarMain : operator Asin non valid");};
	AffineVarMain&  Atan(const Interval& itv)  	{codac2_error(" AffineVarMain : operator Atan non valid");};
	AffineVarMain&  Aacos(const Interval& itv)  {codac2_error(" AffineVarMain : operator Aacos non valid");};
	AffineVarMain&  Aasin(const Interval& itv)  {codac2_error(" AffineVarMain : operator Aasin non valid");};
	AffineVarMain&  Aatan(const Interval& itv)  {codac2_error(" AffineVarMain : operator Aatan non valid");};
	AffineVarMain&  Acosh(const Interval& itv)  {codac2_error(" AffineVarMain : operator Acosh non valid");};
	AffineVarMain&  Asinh(const Interval& itv)  {codac2_error(" AffineVarMain : operator Asinh non valid");};
	AffineVarMain&  Atanh(const Interval& itv)  {codac2_error(" AffineVarMain : operator Atanh non valid");};
	AffineVarMain&  Aabs(const Interval& itv)  	{codac2_error(" AffineVarMain : operator Aabs non valid");};
	AffineVarMain&  Ainv_CH(const Interval& itv) {codac2_error(" AffineVarMain : operator Ainv_CH non valid");};
	AffineVarMain&  Asqrt_CH(const Interval& itv){codac2_error(" AffineVarMain : operator Asqrt_CH non valid");};
	AffineVarMain&  Aexp_CH(const Interval& itv) {codac2_error(" AffineVarMain : operator Aexp_CH non valid");};
	AffineVarMain&  Alog_CH(const Interval& itv) {codac2_error(" AffineVarMain : operator Alog_CH non valid");};
	AffineVarMain&  Ainv_MR(const Interval& itv) {codac2_error(" AffineVarMain : operator Ainv_CH non valid");};
	AffineVarMain&  Asqrt_MR(const Interval& itv){codac2_error(" AffineVarMain : operator Asqrt_MR non valid");};
	AffineVarMain&  Aexp_MR(const Interval& itv) {codac2_error(" AffineVarMain : operator Aexp_MR non valid");};
	AffineVarMain&  Alog_MR(const Interval& itv) {codac2_error(" AffineVarMain : operator Alog_MR non valid");};
*/
};



template<class T>
AffineVarMain<T>::AffineVarMain(int size, int var1, const Interval& itv) :
		AffineMain<T>(size, var1, itv),
		_var		(var1) {
	assert((size>var1)&&(var1>=0));
}



template<class T>
AffineVarMain<T>::AffineVarMain(const AffineVarMain<T>& x) :
		AffineMain<T>(x),
		_var		(x._var) {
}

template<class T>
AffineVarMain<T>& AffineVarMain<T>::operator=(const AffineVarMain<T>& x) {
	if (this != &x) {
		AffineMain<T>::operator=(x);
		_var = x._var;
	}
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
