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


const double MAX_DOUBLE = std::numeric_limits<double>::max();

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
template<>
AffineMain<AF_fAF2>::AffineMain() :
		 _actif (-2     ),
		 _n		(0		),
		 _elt	(nullptr	,oo)	{
 }


template<>
AffineMain<AF_fAF2>& AffineMain<AF_fAF2>::operator=(const Interval& x) {
	if (x.is_empty()) {
		_actif = -1;
		_elt._err = 0.0;
	} else if (x.ub()>= oo && x.lb()<= -oo ) {
		_actif = -2;
		_elt._err = 0.0;
	} else if (x.ub()>= oo ) {
		_actif = -3;
		if (x.lb()>= MAX_DOUBLE) this->_elt._err = MAX_DOUBLE;
		else this->_elt._err = x.lb();
	} else if (x.lb()<= -oo ) {
		_actif = -4;
		if (x.ub()<= -MAX_DOUBLE) this->_elt._err = -MAX_DOUBLE;
		else this->_elt._err = x.ub();
	} else  {
		if (_elt._val==nullptr) { _elt._val = new double[_n+1]; }
		_elt._val[0] = x.mid();
		for (int i=1; i<=_n;i++) {
			_elt._val[i] =0;
		}
		if ( x.is_degenerated()) {
			_actif=0;
			_elt._err	= 0;
		} else {
			_actif = 1;
			_elt._err	= x.rad();
		}
	}
	return *this;
}



template<>
AffineVarMain<AF_fAF2>& AffineVarMain<AF_fAF2>::operator=(const Interval& x) {
	assert(this->_var >= 0 && this->_n > this->_var);
	if (x.is_empty()) {
		this->_actif = -1;
		this->_elt._err = 0.0;
	} else if (x.ub()>= oo && x.lb()<= -oo ) {
		this->_actif = -2;
		this->_elt._err = 0.0;
	} else if (x.ub()>= oo ) {
		this->_actif = -3;
		if (x.lb()>= MAX_DOUBLE) this->_elt._err = MAX_DOUBLE;
		else this->_elt._err = x.lb();
	} else if (x.lb()<= -oo ) {
		this->_actif = -4;
		if (x.ub()<= -MAX_DOUBLE) this->_elt._err = -MAX_DOUBLE;
		else this->_elt._err = x.ub();
	} else  {
		if (this->_elt._val==nullptr) { this->_elt._val = new double[_n+1]; }
		this->_elt._val[0] = x.mid();
		for (int i=1; i<=_n;i++) {
			this->_elt._val[i] =0;
		}
		if ( x.is_degenerated()) {
			this->_actif = 0;
			this->_elt._err = 0.0;
		} else {
			this->_actif = 1;
			this->_elt._val[_var+1] = x.rad();
			this->_elt._err = 0.0;  // uncertainty fully captured in _val[_var+1]
		}
	} 
	return *this;
}





template<>
AffineMain<AF_fAF2>::AffineMain(int size, int var, const Interval& itv) :
			_actif	(0),
			_n 		(size),
			_elt	(nullptr,0.0)
{
	assert((size>=0) && (var>=0) && (var<=size));
	if (!(itv.is_unbounded()||itv.is_empty())) {
		_elt._val	=new double[size + 1];
		_elt._val[0] = itv.mid();
		for (int i = 1; i <= size; i++){
			_elt._val[i] = 0.0;
		}
		if (! itv.is_degenerated()) {
			_actif =1;
			_elt._val[var+1] = itv.rad();
		}
	} else if (itv.is_empty()) {
		_actif = -1;
		_elt._err = 0.0;
	} else if (itv.ub()>= oo && itv.lb()<= -oo ) {
		_actif = -2;
		_elt._err = 0.0;
	} else if (itv.ub()>= oo ) {
		_actif = -3;
		if (itv.lb()>= MAX_DOUBLE) this->_elt._err = MAX_DOUBLE;
		else this->_elt._err = itv.lb();
	} else if (itv.lb()<= -oo ) {
		_actif = -4;
		if (itv.ub()<= -MAX_DOUBLE) this->_elt._err = -MAX_DOUBLE;
		else this->_elt._err = itv.ub();
	}
}



template<>
AffineMain<AF_fAF2>::AffineMain(const AffineMain<AF_fAF2>& x) :
		_actif	(x._actif),
		_n		(x.size()),
		_elt	(nullptr	,x._elt._err ) {
	if (x.is_actif()) {
		_elt._val =new double[x.size() + 1];
		for (int i = 0; i <= x.size(); i++){
			_elt._val[i] = x._elt._val[i];
		}
	}
}



template<>
AffineMain<AF_fAF2>& AffineMain<AF_fAF2>::operator=(const AffineMain<AF_fAF2>& x) {
	if (this != &x) {
		_elt._err = x._elt._err;
		_actif = x._actif;
		if (x.is_actif()) {
			if (_n == x.size()) {
				if (_elt._val==nullptr) { _elt._val = new double[x.size()+1]; }
			} else {
				_n =x._n;
				if (_elt._val!=nullptr) { delete[] _elt._val; }
				_elt._val = new double[x.size()+1];
			}
			for (int i = 0; i <= x.size(); i++) {
				_elt._val[i] = x._elt._val[i];
			}
		} else {
			_n =x.size();
		}
	}
	return *this;
}


template<>
double AffineMain<AF_fAF2>::val(int i) const{
	assert(is_actif() &&(0<=i) && (i<size()));
	return _elt._val[i+1];
}

template<>
double AffineMain<AF_fAF2>::err() const{
	assert(is_actif() );
	return _elt._err;
}


template<>
const Interval AffineMain<AF_fAF2>::itv() const {

	switch(_actif) {
	case -1 : {
		return Interval::empty();
		break;
	}
	case -2 : {
		return Interval();
		break;
	}
	case -3 : {
		return Interval(_elt._err,oo);
		break;
	}
	case -4: {
		return Interval(-oo,_elt._err);
		break;
	}
	case 0: {
		return Interval(_elt._val[0]);
		break;
	}
	default: { // _actif==1
		Interval res(_elt._val[0]);
		Interval pmOne(-1.0, 1.0);
		for (int i = 1; i <= size(); i++){
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
	return (is_actif())? _elt._val[0] : itv().mid();
}




/**
 * Code for the particular case:
 * if the affine form is actif, _actif=1  and _n is the size of the affine form
 * if the set is degenerate, _actif = 0 and itv().diam()< AF_EC
 * if the set is empty, _actif = -1
 * if the set is ]-oo,+oo[, _actif = -2 and _err =]-oo,+oo[
 * if the set is [a, +oo[ , _actif = -3 and _err = [a, +oo[
 * if the set is ]-oo, a] , _actif = -4 and _err = ]-oo, a]
 *
 */

template<>
AffineMain<AF_fAF2>& AffineMain<AF_fAF2>::Aneg() {
	switch(_actif) {
	case -3 : {
		_elt._err = -_elt._err;
		_actif    = -4;
		break;
	}
	case -4 : {
		_elt._err = -_elt._err;
		_actif    = -3;
		break;
	}
	case 0 :{
		_elt._val[0] = (-_elt._val[0]);
		break;
	}
	case 1 : {
		for (int i = 0; i <= size(); i++) {
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
	if (_actif==1) {  // multiply by a scalar alpha
		if (alpha==0.0) {
			_actif = 0;
			for (int i=0; i<=_n;i++) {
				_elt._val[i]=0;
			}
			_elt._err = 0;
		} else if ( std::fabs(alpha) < oo) {
			double temp, ttt, sss, eee;
			ttt= 0.0;
			sss= 0.0;
			for (int i=0; i<=size();i++) {
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
			for (int i=0;i<=size();i++) {
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
		if (is_actif()) {
			*this = beta;
			return *this;
		}
		*this = this->itv() + beta;
		return *this;
	}

	if (_actif==1) {
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
		if (is_actif()) {
			if ((ddelta)<oo) {
				_actif=1;
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
void AffineMain<AF_fAF2>::resize(int n) {
	assert(n>=1);
	if (n>_n) {
		double * tmp= new double[n+1];
		int i=0;
		for (;i<=_n;i++) {
			tmp[i] = _elt._val[i];
		}
		for (;i<=n;i++) {
			tmp[i] = 0;
		}
		_n = n;
		delete[] _elt._val;
		_elt._val = tmp;
	} else if (n<_n) {
		double * tmp= new double[n+1];
		int i=0;
		for (;i<=n;i++) {
			tmp[i] = _elt._val[i];
		}
		_n = n;
		delete[] _elt._val;
		_elt._val = tmp;
	}
}





template<>
AffineMain<AF_fAF2>& AffineMain<AF_fAF2>::operator+=(const AffineMain<AF_fAF2>& y) {

	if (is_actif() && y.is_actif()) {
		if (y.is_degenerated()) {
			*this += y._elt._val[0];
		} else if (is_degenerated()) {
			double tmp = _elt._val[0];
			*this = y;
			*this += tmp;
		} else {
			if (_n < y.size()) {
				this->resize(y.size());
			}
			double temp, ttt, sss, eee;
			ttt=0.0;
			sss=0.0;
			for(int i=0;i<=y.size();i++) {
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
			for (int i=0;i<=_n;i++) {
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

	if (is_actif() && (y.is_actif())) {
		if (y.is_degenerated()) {
			*this *= y._elt._val[0];
		}	else if (is_degenerated()) {
			double tmp = _elt._val[0];
			*this = y;
			*this *= tmp;
		} else 	 {
			if (_n < y.size()) {
				this->resize(y.size());
			}
			double Sx, Sy, Sxy, Sz, ttt, sss, ppp, tmp, xVal0, eee;
			double * xTmp;

			xTmp = new double[_n + 1];
			Sx=0.0; Sy=0.0; Sxy=0.0; Sz=0.0; ttt=0.0; sss=0.0; ppp=0.0; tmp=0.0; xVal0=0.0; eee=0.0;

			for (int i = 1; i <= _n; i++) {
				ppp=0;
				if (i<=y.size()) {
					eee = _elt.twoProd(_elt._val[i],y._elt._val[i], &ppp);
					ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));

					eee = _elt.twoSum(Sz,ppp, &tmp);
					ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
					Sz = tmp;

					if (std::fabs(Sz) < AF_EC) {
						sss = (1+2*AF_EM)*(sss+ std::fabs(Sz));
						Sz = 0.0;
					}

					eee = _elt.twoSum(Sxy,std::fabs(ppp), &tmp);
					ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
					Sxy = tmp;

					if (std::fabs(Sxy) < AF_EC) {
						sss = (1+2*AF_EM)*(sss+ std::fabs(Sxy));
						Sxy = 0.0;
					}
				}
				eee = _elt.twoSum(Sx,std::fabs(_elt._val[i]), &tmp);
				ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
				Sx = tmp;

				if (std::fabs(Sx) < AF_EC) {
					sss = (1+2*AF_EM)*(sss+ std::fabs(Sx));
					Sx = 0.0;
				}

				if (i<=y.size()) {
					eee = _elt.twoSum(Sy,std::fabs(y._elt._val[i]), &tmp);
					ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
					Sy = tmp;

					if (std::fabs(Sy) < AF_EC) {
						sss = (1+2*AF_EM)*(sss+ std::fabs(Sy));
						Sy = 0.0;
					}
				}
			}

			xVal0 = _elt._val[0];
			// RES = X%T(0) * res
			for (int i = 0; i <= _n; i++) {
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
			for (int i = 1; i <= y.size(); i++) {
				eee = _elt.twoProd(xVal0,y._elt._val[i], &ppp);
				ttt = (1+2*AF_EM)*(ttt+std::fabs(eee));
				xTmp[i] = ppp;

				if (std::fabs(xTmp[i]) < AF_EC) {
					sss = (1+2*AF_EM)*(sss+ std::fabs(xTmp[i]));
					xTmp[i] = 0.0;
				}

			}

			//RES =  RES + Xtmp = ( Y%(0) * X ) + ( X%T(0) * Y - X%T(0)*Y%(0) )
			for (int i = 0; i <= y.size(); i++) {

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
			for (int i=0;i<=_n;i++) {
				b &= (std::fabs(_elt._val[i])<oo);
			}
			if (!b) *this = Interval();

			delete[] xTmp;
		}

	} else { // y or x is not a valid affine form. So we multiply y.itv() such as an interval
		*this = (itv() * y.itv());
	}

	//std::cout << "OUT *= "<<std::endl<< *this << std::endl;
	return *this;
}



template<>
AffineMain<AF_fAF2>& AffineMain<AF_fAF2>::operator*=(const Interval& y) {
	if (	(!is_actif())||
			y.is_empty()||
			y.is_unbounded() ) {
		*this = itv()*y;
	} else {
		if (y.is_degenerated()) {
			*this *= y.mid();
		} else {
			AffineMain<AF_fAF2> tmp;
			tmp= y.mid();	 // to check if it is the best way to do it
			tmp.inflate(y.rad());
			*this *= tmp;
		}
	}
	return *this;
}




template<>
AffineMain<AF_fAF2>& AffineMain<AF_fAF2>::Asqr(const Interval& itv) {

	if (	(!is_actif())||
			itv.is_empty()||
			itv.is_unbounded()||
			(itv.diam() < AF_EC)  ) {
		*this = pow(itv,2);

	} else  {

		double Sx, Sx2, ttt, sss, ppp, x0, eee,tmp;
		Sx = 0; Sx2 = 0; ttt = 0; sss = 0; ppp = 0; x0 = 0; eee =0.0; tmp =0.0;

		// compute the error
		for (int i = 1; i <= _n; i++) {

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
		for (int i = 1; i <= _n; i++) {

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
			for (int i=0;i<=_n;i++) {
				b &= (std::fabs(_elt._val[i])<oo);
			}
			if (!b) {
				*this = Interval();
			}
		}

	}

//	std::cout << "out sqr "<<std::endl;
	return *this;
}


template<>
void AffineMain<AF_fAF2>::compact(double tol){
	for (int i=1;i<=_n;i++) {
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





