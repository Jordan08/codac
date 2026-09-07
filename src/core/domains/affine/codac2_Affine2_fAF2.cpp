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

#include <algorithm>
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
void AffineMain<AF_fAF2>::assign_interval_core(const Interval& x) {
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
}

template<>
AffineMain<AF_fAF2>& AffineMain<AF_fAF2>::operator=(const Interval& x) {
	assert(_n_noise >= 0);
	assign_interval_core(x);
	return *this;
}



template<>
AffineVarMain<AF_fAF2>& AffineVarMain<AF_fAF2>::operator=(const Interval& x) {
	assert(_n_noise >= 0);
	assert(_var >= 0);
	assert(_var < _n_noise);
	this->assign_interval_core(x);
	if (this->is_active()) {
		// The radius the shared core left in _elt._err belongs on this
		// variable's own dedicated noise symbol instead.
		this->_elt._val[_var+1] = this->_elt._err;
		this->_elt._err = 0.0;  // uncertainty fully captured in _val[_var+1]
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
		double temp;
		const double eee = _elt.twoSum(_elt._val[0],beta,&temp);
		const double ttt = (1+2*AF_EM)*std::fabs(eee);
		double sss = 0.0;
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
				double temp;
				const double eee = _elt.twoSum(_elt._err,std::fabs(ddelta), &temp);
				const double ttt = (1+2*AF_EM)*std::fabs(eee);
				double sss = 0.0;
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
			const double tmp = _elt._val[0];
			*this = y;
			*this += tmp;
		} else {
			const Index yn = y.noise_count();   // y is const throughout: read it once, not once per iteration
			if (_n_noise < yn) {
				this->resize_noise(yn);
			}
			double temp, ttt, sss, eee;
			ttt=0.0;
			sss=0.0;
			for(Index i=0;i<=yn;i++) {
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
International Conference on, pages 458–467, Los Alamitos, CA, USA, 2004. IEEE Computer Society.
Adding from Jordan NININ:
	The quadratic remainder u*v, with u=sum(a_i.e_i) and
	v=sum(b_i.e_i), admits two valid bands. Neither one dominates,
	so the tightest is their intersection.
		B1 = Sz/2 +- (Sx.Sy - Sxy/2)            the AF2 band
		B2 = [-Sm^2/4 , Sp^2/4]                 by polarisation,
			u*v = ((u+v)^2-(u-v)^2)/4, each square lying in [0,S^2]
	qbeta is the mid point of the intersection, qdelta its radius.
	On a square, u=v gives Sm=0 and Sp=2.Sx, hence [0,Sx^2], the
	Chebyshev band, which always beats B1.
	The accumulated errors behave like two extra independent noise
	symbols of coefficients _err, so they enter every radius below.
 */
template<>
AffineMain<AF_fAF2>& AffineMain<AF_fAF2>::operator*=(const AffineMain<AF_fAF2>& y) {

	if (is_active() && (y.is_active())) {
		if (y.is_degenerated()) {
			*this *= y._elt._val[0];
		}	else if (is_degenerated()) {
			const double tmp = _elt._val[0];
			*this = y;
			*this *= tmp;
		} else 	 {
			const Index yn = y.noise_count();   // y is const throughout: read it once, not once per iteration
			if (_n_noise < yn) {
				this->resize_noise(yn);
			}
			double Sx, Sy, Sxy, Sz, Sp, Sm, ttt, sss, ppp, tmp, eee;
			std::unique_ptr<double[]> xTmp;

			xTmp = std::make_unique<double[]>(_n_noise + 1);
			Sx=0.0; Sy=0.0; Sxy=0.0; Sz=0.0; Sp=0.0; Sm=0.0; ttt=0.0; sss=0.0; ppp=0.0; tmp=0.0; eee=0.0;

			// These accumulators may later be multiplied by quantities at the
			// opposite scale. Do not discard them using an absolute threshold.
			// For example, 1e-300 * 1e300 contributes at order one.
			for (Index i = 1; i <= _n_noise; i++) {
				ppp = 0.0;
				// Each of _elt._val[i] and (when in range) y._elt._val[i] was
				// re-read from the array up to four times per iteration below;
				// reading them once here gives the identical value every time
				// (nothing writes to either array within this loop).
				const double ai = _elt._val[i];
				const double bi = (i <= yn)? y._elt._val[i] : 0.0;

				if (i <= yn) {
					eee = _elt.twoProd(ai, bi, &ppp);
					ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));

					eee = _elt.twoSum(Sz, ppp, &tmp);
					ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
					Sz = tmp;

					eee = _elt.twoSum(Sxy, std::fabs(ppp), &tmp);
					ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
					Sxy = tmp;
				}

				eee = _elt.twoSum(Sx, std::fabs(ai), &tmp);
				ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
				Sx = tmp;

				if (i <= yn) {
					eee = _elt.twoSum(Sy, std::fabs(bi), &tmp);
					ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
					Sy = tmp;
				}
				// polarisation: u*v = ((u+v)^2-(u-v)^2)/4 needs the radii of
				// u+v and u-v, hence these two extra sums
				{
					eee = _elt.twoSum(Sp, std::fabs(ai+bi), &tmp);
					ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
					Sp = tmp;
					eee = _elt.twoSum(Sm, std::fabs(ai-bi), &tmp);
					ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
					Sm = tmp;
				}
			}

			const double xVal0 = _elt._val[0];
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
			for (Index i = 1; i <= yn; i++) {
				eee = _elt.twoProd(xVal0,y._elt._val[i], &ppp);
				ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
				xTmp[i] = ppp;

				if (std::fabs(xTmp[i]) < AF_EC) {
					sss = (1+2*AF_EM)*(sss+ std::fabs(xTmp[i]));
					xTmp[i] = 0.0;
				}

			}

			//RES =  RES + Xtmp = ( Y%(0) * X ) + ( X%T(0) * Y - X%T(0)*Y%(0) )
			for (Index i = 0; i <= yn; i++) {

				eee = _elt.twoSum(_elt._val[i],xTmp[i], &tmp);
				ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
				_elt._val[i] = tmp;

				if (std::fabs(_elt._val[i]) < AF_EC) {
					sss = (1+2*AF_EM)*(sss+ std::fabs(_elt._val[i]));
					_elt._val[i] = 0.0;
				}

			}

			// The quadratic remainder u*v, with u=sum(a_i.e_i) and
			// v=sum(b_i.e_i), admits two valid bands. Neither one dominates,
			// so the tightest is their intersection.
			//   B1 = Sz/2 +- (Sx.Sy - Sxy/2)            the AF2 band
			//   B2 = [-Sm^2/4 , Sp^2/4]                 by polarisation,
			//        u*v = ((u+v)^2-(u-v)^2)/4, each square lying in [0,S^2]
			// qbeta is the mid point of the intersection, qdelta its radius.
			// On a square, u=v gives Sm=0 and Sp=2.Sx, hence [0,Sx^2], the
			// Chebyshev band, which always beats B1.
			// The accumulated errors behave like two extra independent noise
			// symbols of coefficients _err, so they enter every radius below.
			// Every add/sub/mul below goes through twoSum/twoProd, exactly
			// like the loops above -- but their residuals are kept out of
			// the shared ttt (whose AF_EE penalty at the very end of the
			// function is calibrated for compounding over the _n_noise
			// loops above, not for this fixed handful of steps): they are
			// tracked in their own qerr instead, folded straight into
			// qdelta below, which already carries its own AF_EM margin at
			// its point of use.
			double qerr = 0.0;
			eee = _elt.twoSum(Sx, _elt._err, &tmp);
			qerr = (1+2*AF_EM)*(qerr+std::fabs(eee));
			const double sxe = tmp;

			eee = _elt.twoSum(Sy, y._elt._err, &tmp);
			qerr = (1+2*AF_EM)*(qerr+std::fabs(eee));
			const double sye = tmp;

			eee = _elt.twoSum(_elt._err, y._elt._err, &tmp);
			qerr = (1+2*AF_EM)*(qerr+std::fabs(eee));
			const double erry = tmp;   // _elt._err + y._elt._err, shared below

			eee = _elt.twoSum(Sp, erry, &tmp);
			qerr = (1+2*AF_EM)*(qerr+std::fabs(eee));
			const double spe = tmp;

			eee = _elt.twoSum(Sm, erry, &tmp);
			qerr = (1+2*AF_EM)*(qerr+std::fabs(eee));
			const double sme = tmp;

			double qlo, qhi;
			{
				eee = _elt.twoProd(sxe, sye, &ppp);
				qerr = (1+2*AF_EM)*(qerr+std::fabs(eee));
				const double sxy_hat = ppp;   // sxe*sye

				const double half_Sxy = 0.5*Sxy;   // exact: power-of-two scaling

				eee = _elt.twoSum(sxy_hat, -half_Sxy, &tmp);
				qerr = (1+2*AF_EM)*(qerr+std::fabs(eee));
				const double r1 = tmp;

				const double c1 = 0.5*Sz;   // exact

				eee = _elt.twoSum(c1, -r1, &tmp);
				qerr = (1+2*AF_EM)*(qerr+std::fabs(eee));
				qlo = tmp;

				eee = _elt.twoSum(c1, r1, &tmp);
				qerr = (1+2*AF_EM)*(qerr+std::fabs(eee));
				qhi = tmp;

				const double qlo1 = qlo, qhi1 = qhi;   // B1 band, kept for the empty-intersection fallback

				eee = _elt.twoProd(sme, sme, &ppp);
				qerr = (1+2*AF_EM)*(qerr+std::fabs(eee));
				const double lo2 = -0.25*ppp;   // exact: power-of-two scaling

				eee = _elt.twoProd(spe, spe, &ppp);
				qerr = (1+2*AF_EM)*(qerr+std::fabs(eee));
				const double hi2 = 0.25*ppp;   // exact: power-of-two scaling

				if (lo2 > qlo) qlo = lo2;
				if (hi2 < qhi) qhi = hi2;
				if (qhi < qlo) {   // B1 and B2 failed to overlap: fall back to B1 alone
					qlo = qlo1;
					qhi = qhi1;
				}
			}

			eee = _elt.twoSum(qlo, qhi, &tmp);
			qerr = (1+2*AF_EM)*(qerr+std::fabs(eee));
			const double qbeta = 0.5*tmp;   // exact: power-of-two scaling

			eee = _elt.twoSum(qbeta, -qlo, &tmp);
			qerr = (1+2*AF_EM)*(qerr+std::fabs(eee));
			const double d_lo = tmp;

			eee = _elt.twoSum(qhi, -qbeta, &tmp);
			qerr = (1+2*AF_EM)*(qerr+std::fabs(eee));
			const double d_hi = tmp;

			const double qdelta = std::max(d_lo, d_hi) + qerr;

			eee = _elt.twoSum(_elt._val[0],qbeta, &tmp);
			ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
			_elt._val[0] = tmp;

			if (std::fabs(_elt._val[0]) < AF_EC) {
				sss = (1+2*AF_EM)*(sss+ std::fabs(_elt._val[0]));
				_elt._val[0] = 0.0;
			}

			// qdelta already covers the cross terms, the accumulated errors
			// having been folded into the radii above.
			const double xerr = _elt._err;
			const double yerr = y._elt._err;

			_elt._err = (1+ 2*AF_EM) * (
					((1+ 2*AF_EM) *std::fabs(y._elt._val[0]) * xerr)  +
					((1+ 2*AF_EM) *std::fabs(xVal0) * yerr)  +
					((1+ 2*AF_EM) * qdelta)  +
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


// Kept for reference, not used: the AF2 product as it was before the
// polarisation bound was added to operator*=.
template<>
AffineMain<AF_fAF2>& AffineMain<AF_fAF2>::Amul_AF2(const AffineMain<AF_fAF2>& y) {

	if (is_active() && (y.is_active())) {
		if (y.is_degenerated()) {
			*this *= y._elt._val[0];
		}	else if (is_degenerated()) {
			const double tmp = _elt._val[0];
			*this = y;
			*this *= tmp;
		} else 	 {
			const Index yn = y.noise_count();   // y is const throughout: read it once, not once per iteration
			if (_n_noise < yn) {
				this->resize_noise(yn);
			}
			double Sx, Sy, Sxy, Sz, ttt, sss, ppp, tmp, eee;
			std::unique_ptr<double[]> xTmp;

			xTmp = std::make_unique<double[]>(_n_noise + 1);
			Sx=0.0; Sy=0.0; Sxy=0.0; Sz=0.0; ttt=0.0; sss=0.0; ppp=0.0; tmp=0.0; eee=0.0;

			// These accumulators may later be multiplied by quantities at the
			// opposite scale. Do not discard them using an absolute threshold.
			// For example, 1e-300 * 1e300 contributes at order one.
			for (Index i = 1; i <= _n_noise; i++) {
				ppp = 0.0;
				// _elt._val[i] and (when in range) y._elt._val[i] were each
				// re-read from the array up to twice per iteration below;
				// reading them once here gives the identical value every
				// time (nothing writes to either array within this loop).
				const double ai = _elt._val[i];
				const double bi = (i <= yn)? y._elt._val[i] : 0.0;

				if (i <= yn) {
					eee = _elt.twoProd(ai, bi, &ppp);
					ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));

					eee = _elt.twoSum(Sz, ppp, &tmp);
					ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
					Sz = tmp;

					eee = _elt.twoSum(Sxy, std::fabs(ppp), &tmp);
					ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
					Sxy = tmp;
				}

				eee = _elt.twoSum(Sx, std::fabs(ai), &tmp);
				ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
				Sx = tmp;

				if (i <= yn) {
					eee = _elt.twoSum(Sy, std::fabs(bi), &tmp);
					ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
					Sy = tmp;
				}
			}

			const double xVal0 = _elt._val[0];
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
			for (Index i = 1; i <= yn; i++) {
				eee = _elt.twoProd(xVal0,y._elt._val[i], &ppp);
				ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
				xTmp[i] = ppp;

				if (std::fabs(xTmp[i]) < AF_EC) {
					sss = (1+2*AF_EM)*(sss+ std::fabs(xTmp[i]));
					xTmp[i] = 0.0;
				}

			}

			//RES =  RES + Xtmp = ( Y%(0) * X ) + ( X%T(0) * Y - X%T(0)*Y%(0) )
			for (Index i = 0; i <= yn; i++) {

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
			const double half_Sz = ppp;

			eee = _elt.twoSum(_elt._val[0],half_Sz, &tmp);
			ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
			_elt._val[0] = tmp;

			if (std::fabs(_elt._val[0]) < AF_EC) {
				sss = (1+2*AF_EM)*(sss+ std::fabs(_elt._val[0]));
				_elt._val[0] = 0.0;
			}

			eee = _elt.twoSum(_elt._err,Sx, &tmp);
			ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
			const double xerr_Sx = tmp;   // _elt._err + Sx

			eee = _elt.twoSum(y._elt._err,Sy, &ppp);
			ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
			const double yerr_Sy = ppp;   // y._elt._err + Sy


			_elt._err = (1+ 2*AF_EM) * (
					((1+ 2*AF_EM) *std::fabs(y._elt._val[0]) * _elt._err)  +
					((1+ 2*AF_EM) *std::fabs(xVal0) * y._elt._err)  +
					((1+ 2*AF_EM) *(xerr_Sx * yerr_Sy)) +
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

		double Sx, ttt, sss, ppp, eee,tmp;
		Sx = 0; ttt = 0; sss = 0; ppp = 0; eee =0.0; tmp =0.0;

		// compute the error. The Chebyshev band below only needs sum|a_i|,
		// unlike the AF2 band of Asqr_AF2 which also needed sum(a_i^2).
		for (Index i = 1; i <= _n_noise; i++) {

			eee = _elt.twoSum(Sx,std::fabs(_elt._val[i]), &tmp);
			ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
			Sx = tmp;

			if (std::fabs(Sx) < AF_EC) {
				sss = (1+2*AF_EM)*(sss+ std::fabs(Sx));
				Sx = 0.0;
			}

		}
		// compute 2*_elt._val[0]*(*this)
		const double x0 = _elt._val[0];

		eee = _elt.twoProd(x0,x0, &ppp);
		ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
		_elt._val[0] = ppp;

		if (std::fabs(_elt._val[0]) < AF_EC) {
			sss = (1+2*AF_EM)*(sss+ std::fabs(_elt._val[0]));
			_elt._val[0] = 0.0;
		}

		// compute 2*_elt._val[0]*(*this)
		const double x0_2 = 2*x0;   // exact: power-of-two scaling, loop-invariant
		for (Index i = 1; i <= _n_noise; i++) {

			eee = _elt.twoProd(x0_2,_elt._val[i], &ppp);
			ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
			_elt._val[i] = ppp;

			if (std::fabs(_elt._val[i]) < AF_EC) {
				sss = (1+2*AF_EM)*(sss+ std::fabs(_elt._val[i]));
				_elt._val[i] = 0.0;
			}

		}

		// Chebyshev bound of the quadratic remainder: with u=sum(a_i.e_i),
		// u^2 lies in [0,Sx^2], which is a band centred on Sx^2/2 with radius
		// Sx^2/2. The AF2 band, centred on sum(a_i^2)/2 with radius
		// Sx^2-sum(a_i^2)/2, always contains it, so no intersection is needed.
		// This is what operator*= obtains by polarisation on x*x.
		// The accumulated error behaves like an extra independent noise symbol
		// of coefficient _err, so it belongs inside the radius. Like above,
		// the add and the square go through twoSum/twoProd; their residuals
		// are tracked in their own qerr rather than the shared ttt (whose
		// AF_EE penalty is calibrated for the _n_noise loop above, not this
		// fixed pair of steps) and folded straight into the radius term of
		// _elt._err below.
		double qerr = 0.0;
		eee = _elt.twoSum(Sx, _elt._err, &tmp);
		qerr = (1+2*AF_EM)*(qerr+std::fabs(eee));
		const double sxe = tmp;

		eee = _elt.twoProd(sxe, sxe, &ppp);
		qerr = (1+2*AF_EM)*(qerr+std::fabs(eee));
		const double qhalf = 0.5*ppp;   // exact: power-of-two scaling; mid and radius

		eee = _elt.twoSum(_elt._val[0],qhalf, &tmp);
		ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
		_elt._val[0] = tmp;

		if (std::fabs(_elt._val[0]) < AF_EC) {
			sss = (1+2*AF_EM)*(sss+ std::fabs(_elt._val[0]));
			_elt._val[0] = 0.0;
		}

		// qhalf already covers the cross terms, the accumulated error having
		// been folded into the radius above.
		const double xerr = _elt._err;

		_elt._err = (1+ 2*AF_EM) * (
				((1+ 2*AF_EM) *2*std::fabs(x0) * xerr)  +
				((1+ 2*AF_EM) * (qhalf + qerr))  +
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


// Kept for reference, not used: the AF2 square as it was before the
// Chebyshev bound was adopted in Asqr.
template<>
AffineMain<AF_fAF2>& AffineMain<AF_fAF2>::Asqr_AF2(const Interval& itv) {

	if (	(!is_active())||
			itv.is_empty()||
			itv.is_unbounded()||
			(itv.diam() < AF_EC)  ) {
		*this = pow(itv,2);

	} else  {

		double Sx, Sx2, ttt, sss, ppp, eee,tmp;
		Sx = 0; Sx2 = 0; ttt = 0; sss = 0; ppp = 0; eee =0.0; tmp =0.0;

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
		const double x0 = _elt._val[0];

		eee = _elt.twoProd(x0,x0, &ppp);
		ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
		_elt._val[0] = ppp;

		if (std::fabs(_elt._val[0]) < AF_EC) {
			sss = (1+2*AF_EM)*(sss+ std::fabs(_elt._val[0]));
			_elt._val[0] = 0.0;
		}

		// compute 2*_elt._val[0]*(*this)
		const double x0_2 = 2*x0;   // exact: power-of-two scaling, loop-invariant
		for (Index i = 1; i <= _n_noise; i++) {

			eee = _elt.twoProd(x0_2,_elt._val[i], &ppp);
			ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
			_elt._val[i] = ppp;

			if (std::fabs(_elt._val[i]) < AF_EC) {
				sss = (1+2*AF_EM)*(sss+ std::fabs(_elt._val[i]));
				_elt._val[i] = 0.0;
			}

		}

		eee = _elt.twoProd(0.5,Sx2, &ppp);
		ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
		const double half_Sx2 = ppp;

		eee = _elt.twoSum(_elt._val[0],half_Sx2, &tmp);
		ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
		_elt._val[0] = tmp;

		if (std::fabs(_elt._val[0]) < AF_EC) {
			sss = (1+2*AF_EM)*(sss+ std::fabs(_elt._val[0]));
			_elt._val[0] = 0.0;
		}

		eee = _elt.twoSum(_elt._err,Sx, &tmp);
		ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
		const double xerr_Sx = tmp;   // _elt._err + Sx

		_elt._err = (1+ 2*AF_EM) * (
				((1+ 2*AF_EM) *2*std::fabs(x0) * _elt._err)  +
				((1+ 2*AF_EM) *(xerr_Sx * xerr_Sx)) +
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
  			const double eee = _elt.twoSum(_elt._err,std::fabs(_elt._val[i]), &temp);
  			const double ttt = (1+2*AF_EM)*std::fabs(eee);
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





