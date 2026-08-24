/**
 *  \file codac2_Affine2_fAF2.cpp
 *
 *  Implementation of \c AffineMain specializations for the \c AF_fAF2 model.
 *
 * ----------------------------------------------------------------------------
 *  \date       2026
 *  \author     Jordan Ninin
 *  \copyright  Copyright 2026 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */




#include "codac2_Affine2_fAF2.h"

#include <cassert>

#include "codac2_AffineMain.h"
#include "codac2_AffineVar.h"
#include <climits>


namespace codac2 {

/**
 * Code for the particular case:
 * if the affine form is actif, _status=1  and _n_noise is the number of noise symbols of the affine form
 * if the set is degenerate, _status = 0 and itv().diam()< AF_EC
 * if the set is empty, _status = -1
 * if the set is ]-oo,+oo[, _status = -2 and _err = ]-oo,+oo[
 * if the set is [a, +oo[ , _status = -3 and _err = [a, +oo[
 * if the set is ]-oo, a] , _status = -4 and _err = ]-oo, a]
 *
 */
template<>
AffineMain<AF_fAF2>::AffineMain() :
		 _status  (AffineStatus::AllReals     ),
		 _n_noise (0		),
		 _elt	  (nullptr	,0.0	)	{
 }

template<>
AffineVarMain<AF_fAF2>::AffineVarMain() : 
		AffineMain<AF_fAF2>(), 
		_var(-1) {}

template<>
AffineMain<AF_fAF2>& AffineMain<AF_fAF2>::operator=(const Interval& x) {
	assert(_n_noise >= 0);
	if (x.is_empty()) {
		_status = AffineStatus::Empty;
		_elt._err = 0.0;
	} else if (x.ub()>= oo && x.lb()<= -oo ) {
		_status = AffineStatus::AllReals;
		_elt._err = 0.0;
	} else if (x.ub()>= oo ) {
		_status = AffineStatus::UpperUnbounded;
		if (x.lb()>= std::numeric_limits<double>::max())
			this->_elt._err = std::numeric_limits<double>::max();
		else
			this->_elt._err = x.lb();
	} else if (x.lb()<= -oo ) {
		_status = AffineStatus::LowerUnbounded;
		if (x.ub()<= -std::numeric_limits<double>::max())
			this->_elt._err = -std::numeric_limits<double>::max();
		else
			this->_elt._err = x.ub();
	} else  {
		if (_elt._val==nullptr) { _elt._val = std::make_unique<double[]>(_n_noise+1); }
		_elt._val[0] = x.mid();
		for (Index i=1; i<=_n_noise;i++) {
			_elt._val[i] =0;
		}
		if ( x.is_degenerated()) {
			_status=AffineStatus::Degenerate;
			_elt._err	= 0;
		} else {
			_status = AffineStatus::Active;
			_elt._err	= x.rad();
		}
	}
	return *this;
}



template<>
AffineVarMain<AF_fAF2>& AffineVarMain<AF_fAF2>::operator=(const Interval& x) {
	assert(_n_noise >= 0);
	assert(_var >= 0);
	assert(_var < _n_noise);
	if (x.is_empty()) {
		this->_status = AffineStatus::Empty;
		this->_elt._err = 0.0;
	} else if (x.ub()>= oo && x.lb()<= -oo ) {
		this->_status = AffineStatus::AllReals;
		this->_elt._err = 0.0;
	} else if (x.ub()>= oo ) {
		this->_status = AffineStatus::UpperUnbounded;
		if (x.lb()>= std::numeric_limits<double>::max()) this->_elt._err = std::numeric_limits<double>::max();
		else this->_elt._err = x.lb();
	} else if (x.lb()<= -oo ) {
		this->_status = AffineStatus::LowerUnbounded;
		if (x.ub()<= -std::numeric_limits<double>::max()) this->_elt._err = -std::numeric_limits<double>::max();
		else this->_elt._err = x.ub();
	} else  {
		if (this->_elt._val==nullptr) { this->_elt._val = std::make_unique<double[]>(_n_noise+1); }
		this->_elt._val[0] = x.mid();
		for (Index i=1; i<=_n_noise;i++) {
			this->_elt._val[i] =0;
		}
		if ( x.is_degenerated()) {
			this->_status = AffineStatus::Degenerate;
			this->_elt._err = 0.0;
		} else {
			this->_status = AffineStatus::Active;
			this->_elt._val[_var+1] = x.rad();
			this->_elt._err = 0.0;  // uncertainty fully captured in _val[_var+1]
		}
	} 
	return *this;
}





template<>
AffineMain<AF_fAF2>::AffineMain(Index noise_count, Index var, const Interval& itv) :
			_status	 (AffineStatus::Degenerate),
			_n_noise (noise_count),
			_elt	 (nullptr,0.0)
{
	assert((((noise_count>=0) && (var>=0))||itv.is_empty()) && (var<noise_count));
	if (!(itv.is_unbounded()||itv.is_empty())) {
		_elt._val	=std::make_unique<double[]>(noise_count + 1);
		_elt._val[0] = itv.mid();
		for (Index i = 1; i <= noise_count; i++){
			_elt._val[i] = 0.0;
		}
		if (! itv.is_degenerated()) {
			_status = AffineStatus::Active;
			_elt._val[var+1] = itv.rad();
		}
	} else if (itv.is_empty()) {
		_status =AffineStatus::Empty;
		_elt._err = 0.0;
	} else if (itv.ub()>= oo && itv.lb()<= -oo ) {
		_status = AffineStatus::AllReals;
		_elt._err = 0.0;
	} else if (itv.ub()>= oo ) {
		_status = AffineStatus::UpperUnbounded;
		if (itv.lb()>= std::numeric_limits<double>::max())
			this->_elt._err = std::numeric_limits<double>::max();
		else
			this->_elt._err = itv.lb();
	} else if (itv.lb()<= -oo ) {
		_status = AffineStatus::LowerUnbounded;
		if (itv.ub()<= -std::numeric_limits<double>::max())
			this->_elt._err = -std::numeric_limits<double>::max();
		else
			this->_elt._err = itv.ub();
	}
}



template<>
AffineMain<AF_fAF2>::AffineMain(const AffineMain<AF_fAF2>& x) :
		_status	(x._status),
		_n_noise(x._n_noise),
		_elt	(nullptr	,x._elt._err ) {
	if (x.is_active()) {
		_elt._val =std::make_unique<double[]>(x._n_noise + 1);
		for (Index i = 0; i <= x._n_noise; i++){
			_elt._val[i] = x._elt._val[i];
		}
	}
}



template<>
AffineMain<AF_fAF2>&
AffineMain<AF_fAF2>::operator=(const AffineMain<AF_fAF2>& x)
{
  if (this == &x) {
    return *this;
  }

  _status = x._status;
  _elt._err = x._elt._err;

  if (!x.is_active()) {
    // An inactive affine form has no usable coefficient representation.
    // Releasing the old storage prevents _n_noise from becoming inconsistent
    // with the actual allocation capacity.
    _elt._val = nullptr;
    _n_noise = x._n_noise;
    return *this;
  }

  if (_elt._val == nullptr || _n_noise != x._n_noise) {
    // unique_ptr::operator= frees whatever _elt._val previously owned only
    // after the new array has been successfully allocated, so the
    // destination is left unchanged if allocation throws.
    _elt._val = std::make_unique<double[]>(x._n_noise + 1);
    _n_noise = x._n_noise;
  }

  for (Index i = 0; i <= x._n_noise; ++i) {
    _elt._val[i] = x._elt._val[i];
  }

  return *this;
}


template<>
double AffineMain<AF_fAF2>::noise(Index i) const{
	assert(is_active() &&(0<=i) && (i<noise_count()));
	return _elt._val[i+1];
}

template<>
double AffineMain<AF_fAF2>::err() const{
	assert(is_active() );
	return _elt._err;
}


template<>
const Interval AffineMain<AF_fAF2>::itv() const {

	switch(_status) {
	case AffineStatus::Empty : {
		return Interval::empty();
		break;
	}
	case AffineStatus::AllReals : {
		return Interval();
		break;
	}
	case AffineStatus::UpperUnbounded : {
		return Interval(_elt._err,oo);
		break;
	}
	case AffineStatus::LowerUnbounded: {
		return Interval(-oo,_elt._err);
		break;
	}
	case AffineStatus::Degenerate: {
		return Interval(_elt._val[0]);
		break;
	}
	default: { // _status== AffineStatus::Active
		Interval res(_elt._val[0]);
		Interval pmOne(-1.0, 1.0);
		for (Index i = 1; i <= noise_count(); i++){
			res += (_elt._val[i] * pmOne);
		}
		res += _elt._err * pmOne;
		return res;
		break;
	}
	}

}


template<>
double AffineMain<AF_fAF2>::mid() const{
	return (is_active())? _elt._val[0] : itv().mid();
}




/**
 * Code for the particular case:
 * if the affine form is actif, _status=1  and _n_noise is the number of noise symbols of the affine form
 * if the set is degenerate, _status = 0 and itv().diam()< AF_EC
 * if the set is empty, _status = -1
 * if the set is ]-oo,+oo[, _status = -2 and _err =]-oo,+oo[
 * if the set is [a, +oo[ , _status = -3 and _err = [a, +oo[
 * if the set is ]-oo, a] , _status = -4 and _err = ]-oo, a]
 *
 */

template<>
AffineMain<AF_fAF2>& AffineMain<AF_fAF2>::Aneg() {
	switch(_status) {
	case AffineStatus::UpperUnbounded : {
		_elt._err = -_elt._err;
		_status    = AffineStatus::LowerUnbounded;
		break;
	}
	case AffineStatus::LowerUnbounded : {
		_elt._err = -_elt._err;
		_status    = AffineStatus::UpperUnbounded;
		break;
	}
	case AffineStatus::Degenerate :{
		_elt._val[0] = (-_elt._val[0]);
		break;
	}
	case AffineStatus::Active : {
		for (Index i = 0; i <= noise_count(); i++) {
			_elt._val[i] = (-_elt._val[i]);
		}
		break;
	}
	default :
		break;
	}

	return *this;
}




template<>
AffineMain<AF_fAF2>& AffineMain<AF_fAF2>::operator*=(double alpha) {
	if (_status==AffineStatus::Active) {  // multiply by a scalar alpha
		if (alpha==0.0) {
			_status = AffineStatus::Degenerate;
			for (Index i=0; i<=_n_noise;i++) {
				_elt._val[i]=0;
			}
			_elt._err = 0;
		} else if ( std::fabs(alpha) < oo) {
			double temp, ttt, sss, eee;
			ttt= 0.0;
			sss= 0.0;
			for (Index i=0; i<=noise_count();i++) {
				eee = _elt.twoProd(_elt._val[i], alpha, &temp);
				_elt._val[i] = temp;
				ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
				if (std::fabs(_elt._val[i])<AF_EC) {
					sss = (1+2*AF_EM)*(sss+ std::fabs(_elt._val[i]));
					_elt._val[i] = 0.0;
				}
			}
			//	_elt._err = (1+2*AF_EM)*((1+2*AF_EM)*std::fabs(alpha)*_elt._err+AF_EE*AF_EM*ttt + AF_EE*sss);
			_elt._err = (1+2*AF_EM)*( ((1+2*AF_EM)*std::fabs(alpha)*_elt._err) +	((AF_EE*ttt) +	(AF_EE*sss)) );

			bool b = (_elt._err<oo);
			for (Index i=0;i<=noise_count();i++) {
				b &= (std::fabs(_elt._val[i])<oo);
			}
			if (!b) { *this = Interval(); }

		} else {
			*this = itv()*alpha;
		}
	} else {  //scalar alpha
		*this = this->itv()* alpha;
	}
	return *this;
}



template<>
AffineMain<AF_fAF2>& AffineMain<AF_fAF2>::operator+=(double beta) {
	if (beta==0) return *this;

	if (std::fabs(beta)>=oo) {
		// Adding a raw +/-oo scalar.
		// If this affine form currently holds a proper (bounded or degenerate)
		// value, it becomes the canonical "reaches infinity in that direction"
		// state, mirroring what operator=(double) already does: a bounded
		// quantity plus infinity is dominated by that infinity.
		// Otherwise (already empty, or already an extended/half-infinite
		// state), we keep the historical convention -- obtained through the
		// Interval+double operator -- that combining two raw infinities (of
		// possibly conflicting signs) is undefined and yields the empty set.
		if (is_active()) {
			*this = beta;
			return *this;
		}
		*this = this->itv() + beta;
		return *this;
	}

	if (_status==AffineStatus::Active) {
		double temp, ttt, sss, eee;
		ttt=0.0;
		sss=0.0;
		eee = _elt.twoSum(_elt._val[0],beta,&temp);
		ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
		if (std::fabs(temp)<AF_EC) {
			sss = (1+2*AF_EM)*(sss+std::fabs(temp));
			_elt._val[0] = 0.0;
		} else {
			_elt._val[0]=temp;
		}
		//				_elt._err = (1+2*AF_EM)*(_elt._err+ (AF_EE*(AF_EM*ttt)+AF_EE*sss));
		_elt._err = (1+2*AF_EM)*(_elt._err +	(AF_EE*(ttt)+ AF_EE*sss) );

		if (!(_elt._err<oo && (std::fabs(_elt._val[0])<oo))) { *this = Interval(); }

	} else {
		*this = this->itv()+ beta;
	}
	return *this;

}




template<>
AffineMain<AF_fAF2>& AffineMain<AF_fAF2>::inflate(double ddelta) {
	assert(ddelta>=0);
	if (ddelta>0) {
		if (is_active()) {
			if ((ddelta)<oo) {
				_status=AffineStatus::Active;
				double temp, ttt, sss, eee;
				ttt=0.0;
				sss=0.0;
				eee = _elt.twoSum(_elt._err,std::fabs(ddelta), &temp);
				ttt = (1+2*AF_EM)*(std::fabs(eee));
				if (std::fabs(temp)<AF_EC) {
					sss = (1+2*AF_EM)*(std::fabs(temp));
					temp =0;
				}
				//				_elt._err = (1+2*AF_EM)*(temp+ (AF_EE*(AF_EM*ttt)));;
				_elt._err = (1+2*AF_EM)*( temp + (AF_EE*(ttt) + AF_EE*sss) );

				if (!(_elt._err<oo)) { *this = Interval(); }
			}
			else {
				*this = Interval();
			}
		} else {
			*this = this->itv()+Interval(-1,1)*ddelta;
		}
	}
	return *this;
}



template<>
void AffineMain<AF_fAF2>::resize_noise(Index n) {
	assert(n>=1);
	if(n == _n_noise)	return;

	if (_elt._val==nullptr) {
		_n_noise = n;
		_elt._val = std::make_unique<double[]>(n+1);
		for (Index i=0;i<=n;i++) {
			_elt._val[i] = 0;
		}
	} else {
		auto tmp = std::make_unique<double[]>(n+1);
		Index i=0;
		for (;i<=_n_noise;i++) {
			tmp[i] = _elt._val[i];
		}
		for (;i<=n;i++) {
			tmp[i] = 0;
		}
		_n_noise = n;
		_elt._val = std::move(tmp);
	}
	return;
}





template<>
AffineMain<AF_fAF2>& AffineMain<AF_fAF2>::operator+=(const AffineMain<AF_fAF2>& y) {

	if (is_active() && y.is_active()) {
		if (y.is_degenerated()) {
			*this += y._elt._val[0];
		} else if (is_degenerated()) {
			double tmp = _elt._val[0];
			*this = y;
			*this += tmp;
		} else {
			if (_n_noise < y.noise_count()) {
				this->resize_noise(y.noise_count());
			}
			double temp, ttt, sss, eee;
			ttt=0.0;
			sss=0.0;
			for(Index i=0;i<=y.noise_count();i++) {
				eee = _elt.twoSum(_elt._val[i], y._elt._val[i], &temp);
				ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
				if (std::fabs(temp)<AF_EC) {
					sss = (1+2*AF_EM)*(sss+ std::fabs(temp));
					_elt._val[i] = 0.0;
				}
				else {
					_elt._val[i]=temp;
				}
			}
			// _elt._err = (1+2*AF_EM)*((_elt._err+y._elt._err+ (AF_EE*(AF_EM*ttt)+AF_EE*sss));
			_elt._err = (1+2*AF_EM)*( (_elt._err+y._elt._err) + ((AF_EE*(ttt)) + (AF_EE*sss)) );

			bool b = (_elt._err<oo);
			for (Index i=0;i<=_n_noise;i++) {
				b &= (std::fabs(_elt._val[i])<oo);
			}
			if (!b) {*this = Interval(); }

		}
	} else {
		*this = itv() + y.itv();
	}
	//std::cout << "OUT += "<<std::endl<< *this  << std::endl;
	return *this;
}


/**
 * see  Equation(17)  of
 * X.-H. Vu, D. Sam-Haroud, and B. Faltings. Combining multiple inclusion representa-
tions in numerical constraint propagation. In Tools with Artificial Intelligence, IEEE
International Conference on, pages 458–467, Los Alamitos, CA, USA, 2004. IEEE Com-
puter Society.
 */
template<>
AffineMain<AF_fAF2>& AffineMain<AF_fAF2>::operator*=(const AffineMain<AF_fAF2>& y) {

	if (is_active() && (y.is_active())) {
		if (y.is_degenerated()) {
			*this *= y._elt._val[0];
		}	else if (is_degenerated()) {
			double tmp = _elt._val[0];
			*this = y;
			*this *= tmp;
		} else 	 {
			if (_n_noise < y.noise_count()) {
				this->resize_noise(y.noise_count());
			}
			double Sx, Sy, Sxy, Sz, ttt, sss, ppp, tmp, xVal0, eee;
			std::unique_ptr<double[]> xTmp;

			xTmp = std::make_unique<double[]>(_n_noise + 1);
			Sx=0.0; Sy=0.0; Sxy=0.0; Sz=0.0; ttt=0.0; sss=0.0; ppp=0.0; tmp=0.0; xVal0=0.0; eee=0.0;

			// These accumulators may later be multiplied by quantities at the
			// opposite scale. Do not discard them using an absolute threshold.
			// For example, 1e-300 * 1e300 contributes at order one.
			for (Index i = 1; i <= _n_noise; i++) {
				ppp = 0.0;

				if (i <= y.noise_count()) {
					eee = _elt.twoProd(_elt._val[i], y._elt._val[i], &ppp);
					ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));

					eee = _elt.twoSum(Sz, ppp, &tmp);
					ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
					Sz = tmp;

					eee = _elt.twoSum(Sxy, std::fabs(ppp), &tmp);
					ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
					Sxy = tmp;
				}

				eee = _elt.twoSum(Sx, std::fabs(_elt._val[i]), &tmp);
				ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
				Sx = tmp;

				if (i <= y.noise_count()) {
					eee = _elt.twoSum(Sy, std::fabs(y._elt._val[i]), &tmp);
					ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
					Sy = tmp;
				}
			}

			xVal0 = _elt._val[0];
			// RES = X%T(0) * res
			for (Index i = 0; i <= _n_noise; i++) {
				eee = _elt.twoProd(_elt._val[i],y._elt._val[0], &ppp);
				ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
				_elt._val[i] = ppp;

				if (std::fabs(_elt._val[i]) < AF_EC) {
					sss = (1+2*AF_EM)*(sss+ std::fabs(_elt._val[i]));
					_elt._val[i] = 0.0;
				}
			}

			// Xtmp = X%T(0) * Y
			xTmp[0] = 0.0;
			for (Index i = 1; i <= y.noise_count(); i++) {
				eee = _elt.twoProd(xVal0,y._elt._val[i], &ppp);
				ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
				xTmp[i] = ppp;

				if (std::fabs(xTmp[i]) < AF_EC) {
					sss = (1+2*AF_EM)*(sss+ std::fabs(xTmp[i]));
					xTmp[i] = 0.0;
				}

			}

			//RES =  RES + Xtmp = ( Y%(0) * X ) + ( X%T(0) * Y - X%T(0)*Y%(0) )
			for (Index i = 0; i <= y.noise_count(); i++) {

				eee = _elt.twoSum(_elt._val[i],xTmp[i], &tmp);
				ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
				_elt._val[i] = tmp;

				if (std::fabs(_elt._val[i]) < AF_EC) {
					sss = (1+2*AF_EM)*(sss+ std::fabs(_elt._val[i]));
					_elt._val[i] = 0.0;
				}

			}

			eee = _elt.twoProd(0.5,Sz, &ppp);
			ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));

			eee = _elt.twoSum(_elt._val[0],ppp, &tmp);
			ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
			_elt._val[0] = tmp;

			if (std::fabs(_elt._val[0]) < AF_EC) {
				sss = (1+2*AF_EM)*(sss+ std::fabs(_elt._val[0]));
				_elt._val[0] = 0.0;
			}

			eee = _elt.twoSum(_elt._err,Sx, &tmp);
			ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));

			eee = _elt.twoSum(y._elt._err,Sy, &ppp);
			ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));


			_elt._err = (1+ 2*AF_EM) * (
					((1+ 2*AF_EM) *std::fabs(y._elt._val[0]) * _elt._err)  +
					((1+ 2*AF_EM) *std::fabs(xVal0) * y._elt._err)  +
					((1+ 2*AF_EM) *(tmp * ppp)) +
					((1- 2*AF_EM) *(-0.5) *  Sxy)  +
					//					(AF_EE * (AF_EM * ttt))  +
					(AF_EE * (ttt))  +
					(AF_EE * sss)
			);


			bool b = (_elt._err<oo);
			for (Index i=0;i<=_n_noise;i++) {
				b &= (std::fabs(_elt._val[i])<oo);
			}
			if (!b) *this = Interval();
		}

	} else { // y or x is not a valid affine form. So we multiply y.itv() such as an interval
		*this = (itv() * y.itv());
	}

	return *this;
}



template<>
AffineMain<AF_fAF2>& AffineMain<AF_fAF2>::operator*=(const Interval& y) {
	if (	(!is_active())||
			y.is_empty()||
			y.is_unbounded() ) {
		*this = itv()*y;
	} else {
		if (y.is_degenerated()) {
			*this *= y.mid();
		} else {
			AffineMain<AF_fAF2> tmp(y);
			*this *= tmp;
		}
	}
	return *this;
}




template<>
AffineMain<AF_fAF2>& AffineMain<AF_fAF2>::Asqr(const Interval& itv) {

	if (	(!is_active())||
			itv.is_empty()||
			itv.is_unbounded()||
			(itv.diam() < AF_EC)  ) {
		*this = pow(itv,2);

	} else  {

		double Sx, Sx2, ttt, sss, ppp, x0, eee,tmp;
		Sx = 0; Sx2 = 0; ttt = 0; sss = 0; ppp = 0; x0 = 0; eee =0.0; tmp =0.0;

		// compute the error
		for (Index i = 1; i <= _n_noise; i++) {

			eee = _elt.twoProd(_elt._val[i],_elt._val[i], &ppp);
			ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));

			eee = _elt.twoSum(Sx2,ppp, &tmp);
			ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
			Sx2 = tmp;

			if (std::fabs(Sx2) < AF_EC) {
				sss = (1+2*AF_EM)*(sss+ std::fabs(Sx2));
				Sx2 = 0.0;
			}

			eee = _elt.twoSum(Sx,std::fabs(_elt._val[i]), &tmp);
			ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
			Sx = tmp;

			if (std::fabs(Sx) < AF_EC) {
				sss = (1+2*AF_EM)*(sss+ std::fabs(Sx));
				Sx = 0.0;
			}

		}
		// compute 2*_elt._val[0]*(*this)
		x0 = _elt._val[0];

		eee = _elt.twoProd(x0,x0, &ppp);
		ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
		_elt._val[0] = ppp;

		if (std::fabs(_elt._val[0]) < AF_EC) {
			sss = (1+2*AF_EM)*(sss+ std::fabs(_elt._val[0]));
			_elt._val[0] = 0.0;
		}

		// compute 2*_elt._val[0]*(*this)
		for (Index i = 1; i <= _n_noise; i++) {

			eee = _elt.twoProd((2*x0),_elt._val[i], &ppp);
			ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
			_elt._val[i] = ppp;

			if (std::fabs(_elt._val[i]) < AF_EC) {
				sss = (1+2*AF_EM)*(sss+ std::fabs(_elt._val[i]));
				_elt._val[i] = 0.0;
			}

		}

		eee = _elt.twoProd(0.5,Sx2, &ppp);
		ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));

		eee = _elt.twoSum(_elt._val[0],ppp, &tmp);
		ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
		_elt._val[0] = tmp;

		if (std::fabs(_elt._val[0]) < AF_EC) {
			sss = (1+2*AF_EM)*(sss+ std::fabs(_elt._val[0]));
			_elt._val[0] = 0.0;
		}

		eee = _elt.twoSum(_elt._err,Sx, &tmp);
		ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));

		_elt._err = (1+ 2*AF_EM) * (
				((1+ 2*AF_EM) *2*std::fabs(x0) * _elt._err)  +
				((1+ 2*AF_EM) *(tmp * tmp)) +
				((1- 2*AF_EM) *(-0.5) *  Sx2)  +
//					(AF_EE * (AF_EM * ttt))  +
				(AF_EE * (ttt))  +
				(AF_EE * sss)
				);

		{
			bool b = (_elt._err<oo);
			for (Index i=0;i<=_n_noise;i++) {
				b &= (std::fabs(_elt._val[i])<oo);
			}
			if (!b) {
				*this = Interval();
			}
		}

	}

	return *this;
}


template<>
void AffineMain<AF_fAF2>::compact(double tol)
{
  // Empty and unbounded affine forms have no coefficient representation.
  // Their logical dimension may nevertheless be non-zero, so iterating up
  // to _n_noise would dereference a null _elt._val pointer.
  if (!is_active() || _elt._val == nullptr) {
    return;
  }

  // A negative or non-finite threshold has no meaningful compacting
  // semantics. Keep the form unchanged rather than modifying its remainder.
  if (!std::isfinite(tol) || tol < 0.0) {
    return;
  }
  	for (Index i=1;i<=_n_noise;i++) {
  		if (std::fabs(_elt._val[i])<tol) {
  			double temp=0.0;
  			double sss=0.0;
  			double eee = _elt.twoSum(_elt._err,std::fabs(_elt._val[i]), &temp);
  			double ttt = (1+2*AF_EM)*(std::fabs(eee));
  			if (std::fabs(temp)<AF_EC) {
  				sss = (1+2*AF_EM)*(std::fabs(temp));
  				temp =0;
  			}
  //			_elt._err = (1+2*AF_EM)*(temp+ (AF_EE*(AF_EM*ttt)));;
  			_elt._err = (1+2*AF_EM)*( temp + (AF_EE*(ttt) + AF_EE*sss) );
  
  			_elt._val[i] =0;
  		}
  	}
}






//===========================================================================================
//===========================================================================================



}// end namespace codac2





