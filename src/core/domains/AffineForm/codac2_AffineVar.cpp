/**
 *  \file codac2_AffineVar.cpp
 *
 *  Implementation of affine-variable specializations.
 *
 * ----------------------------------------------------------------------------
 *  \date       2020
 *  \author     Jordan Ninin
 *  \copyright  Copyright 2020 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */


#include <cassert>

#include "codac2_AffineVar.h"
#include "codac2_Affine2_fAF2.h"

namespace codac2 {


/**
 * Code for the particular case:
 * if the affine form is actif, _actif=1  and _n is the size of the affine form
 * if the set is degenerate, _actif = 0 and itv().diam()< AF_EC
 * if the set is empty, _actif = -1
 * if the set is ]-oo,+oo[, _actif = -2 and _err = ]-oo,+oo[
 * if the set is [a, +oo[ , _actif = -3 and _err = [a, +oo[
 * if the set is ]-oo, a] , _actif = -4 and _err = ]-oo, a]
 *
 */

//===========================================================================================
//===========================================================================================



template<>
AffineVarMain<AF_fAF2>& AffineVarMain<AF_fAF2>::operator=(const Interval& x) {
	assert(this->_n > this->_var);
	if (x.is_empty()) {
		this->_actif = -1;
		this->_elt._err = 0.0;
	} else if (x.ub()>= oo && x.lb()<= -oo ) {
		_actif = -2;
		this->_elt._err = 0.0;
	} else if (x.ub()>= oo ) {
		this->_actif = -3;
		this->_elt._err = x.lb();
	} else if (x.lb()<= -oo ) {
		this->_actif = -4;
		this->_elt._err = x.ub();
	} else  {
		if (this->_elt._val==nullptr) { this->_elt._val = new double[_n+1]; }
		this->_elt._val[0] = x.mid();
		for (int i=1; i<=_n;i++) {
			this->_elt._val[i] =0;
		}
		if ( x.is_degenerated()) {
			_actif=0;
			this->_elt._err	= 0;
		} else {
			_actif = 1;
			this->_elt._val[_var+1] = x.rad();
		}
	}
	return *this;
}





}// end namespace codac2




