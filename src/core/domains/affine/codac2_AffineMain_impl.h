/**
 *  \file codac2_AffineMain_impl.h
 *
 *  Inline implementations of the AffineMain class members declared in
 *  codac2_AffineMain.h.
 *
 * ----------------------------------------------------------------------------
 *  \date       2026
 *  \author     Jordan Ninin
 *  \copyright  Copyright 2026 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

// Inline functions

namespace codac2 {


template<class T>
inline AffineMain<T> AffineMain<T>::empty()
{
	AffineMain<T> result;
	result.set_empty();
	return result;
}


template<class T>
thread_local typename AffineMain<T>::Affine_Mode AffineMain<T>::mode = AffineMain<T>::AF_Lin_Chebyshev;// default mode is Chebyshev


template<class T>
inline double AffineMain<T>::lb() const { return this->itv().lb(); }

template<class T>
inline double AffineMain<T>::ub() const { return this->itv().ub(); }

template<class T>
inline double AffineMain<T>::rad() const { return this->itv().rad(); }

template<class T>
inline double AffineMain<T>::diam() const { return this->itv().diam(); }

template<class T>
inline double AffineMain<T>::mig() const { return this->itv().mig(); }

template<class T>
inline double AffineMain<T>::mag() const { return this->itv().mag(); }

template<class T>
inline double AffineMain<T>::smag() const { return this->itv().smag(); }

template<class T>
inline double AffineMain<T>::smig() const { return this->itv().smig(); }

template<class T>
inline double AffineMain<T>::volume() const { return this->itv().volume(); }

template<class T>
inline Index AffineMain<T>::size() const { return 1; }

template<class T>
inline AffineMain<T>& AffineMain<T>::init()
{
	*this = Interval(-oo,oo);
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::init(const Interval& x)
{
	*this = x;
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::init_from_list(const std::list<double>& l)
{
	if (l.size() == 1)
		*this = Interval(*l.begin());

	else if (l.size() == 2)
		*this = Interval(*l.begin(), *std::prev(l.end()));

	else
	{
		assert_release("'Affine' can only be defined by one or two 'double' values.");
	}

	return *this;
}

template<class T>
inline bool AffineMain<T>::is_subset(const Interval& x) const { return this->itv().is_subset(x); }

template<class T>
inline bool AffineMain<T>::is_strict_subset(const Interval& x) const { return this->itv().is_strict_subset(x); }

template<class T>
inline bool AffineMain<T>::is_interior_subset(const Interval& x) const { return this->itv().is_interior_subset(x); }

template<class T>
inline bool AffineMain<T>::is_strict_interior_subset(const Interval& x) const { return this->itv().is_strict_interior_subset(x); }

template<class T>
inline bool AffineMain<T>::is_superset(const Interval& x) const { return this->itv().is_superset(x); }

template<class T>
inline bool AffineMain<T>::is_strict_superset(const Interval& x) const { return this->itv().is_strict_superset(x); }


template<class T>
inline bool AffineMain<T>::is_subset(const AffineMain<T>& x) const { return this->itv().is_subset(x.itv()); }

template<class T>
inline bool AffineMain<T>::is_strict_subset(const AffineMain<T>& x) const { return this->itv().is_strict_subset(x.itv()); }

template<class T>
inline bool AffineMain<T>::is_interior_subset(const AffineMain<T>& x) const { return this->itv().is_interior_subset(x.itv()); }

template<class T>
inline bool AffineMain<T>::is_strict_interior_subset(const AffineMain<T>& x) const { return this->itv().is_strict_interior_subset(x.itv()); }

template<class T>
inline bool AffineMain<T>::is_superset(const AffineMain<T>& x) const { return this->itv().is_superset(x.itv()); }

template<class T>
inline bool AffineMain<T>::is_strict_superset(const AffineMain<T>& x) const { return this->itv().is_strict_superset(x.itv()); }


template<class T>
inline bool AffineMain<T>::contains(const double& d) const { return this->itv().contains(d); }

template<class T>
inline bool AffineMain<T>::interior_contains(const double& d) const { return this->itv().interior_contains(d); }

template<class T>
inline bool AffineMain<T>::intersects(const Interval &x) const { return this->itv().intersects(x); }

template<class T>
inline bool AffineMain<T>::overlaps(const Interval &x) const { return this->itv().overlaps(x); }

template<class T>
inline bool AffineMain<T>::is_disjoint(const Interval &x) const { return this->itv().is_disjoint(x); }

template<class T>
inline bool AffineMain<T>::intersects(const AffineMain<T> &x) const { return this->itv().intersects(x.itv()); }

template<class T>
inline bool AffineMain<T>::overlaps(const AffineMain<T> &x) const { return this->itv().overlaps(x.itv()); }

template<class T>
inline bool AffineMain<T>::is_disjoint(const AffineMain<T> &x) const { return this->itv().is_disjoint(x.itv()); }

/**
 * \brief Returns the additive inverse of an affine form.
 *
 * \return \f$-x\f$
 */
template<class T>
inline AffineMain<T> AffineMain<T>::operator-() const {
	AffineMain<T> res(*this);
	res.Aneg();
	return res;
}

template<class T>
inline void AffineMain<T>::change_mode(Affine_Mode tt) {
	switch (tt) {
	case Affine_Mode::AF_Lin_Chebyshev:
	case Affine_Mode::AF_Lin_MinRange:
		mode = tt;
		break;
	default:
		mode = Affine_Mode::AF_Lin_Chebyshev;
		break;
	}
}
template<class T>
inline AffineMain<T>::Affine_Mode AffineMain<T>::get_mode() {
	return mode;
}

template<class T>
inline void AffineMain<T>::compact(){	compact(AF_COMPAC_Tol); }


template<class T>
AffineMain<T>::AffineMain(double d) : AffineMain() {
	*this = Interval(d);
 }

template<class T>
AffineMain<T>::AffineMain(const Interval& x) : AffineMain() {
	*this = x;
 }

template<typename AF>
AffineMain<AF>::operator Interval() const {
    return this->itv();
}


template<class T>
inline AffineMain<T>& AffineMain<T>::operator=(double d) {
	*this = Interval(d);
	return *this;
}

template<class T>
inline bool AffineMain<T>::operator==(const AffineMain<T>& x) const{
	return (this->itv() == x.itv());
}

template<class T>
inline bool AffineMain<T>::operator==(const Interval& x) const{
	return (this->itv() == x);
}
template<class T>
inline bool AffineMain<T>::operator==(double x) const{
    return (this->itv() == x);
}

template<class T>
inline bool AffineMain<T>::operator!=(double x) const{
    return (this->itv() != x);
}

template<class T>
inline bool AffineMain<T>::operator!=(const AffineMain<T>& x) const{
	return (this->itv() != x.itv());
}

template<class T>
inline bool AffineMain<T>::operator!=(const Interval& x) const {
	return (this->itv() != x);
}

template<class T>
inline BoolInterval AffineMain<T>::operator<(const Interval& x) const
{
    return this->itv() < x;
}

template<class T>
inline BoolInterval AffineMain<T>::operator>(const Interval& x) const
{
    return this->itv() > x;
}

template<class T>
inline BoolInterval AffineMain<T>::operator<(const AffineMain& x) const
{
    return this->itv() < x.itv();
}

template<class T>
inline BoolInterval AffineMain<T>::operator>(const AffineMain& x) const
{
    return this->itv() > x.itv();
}

template<class T>
inline void AffineMain<T>::set_empty(){
	*this = Interval::empty();
}


template<class T>
inline Index AffineMain<T>::noise_count() const{
	return _n_noise;
}

template<class T>
inline bool AffineMain<T>::is_active() const{
	return (_status==AffineStatus::Active||_status==AffineStatus::Degenerate);
}

template<class T>
inline bool AffineMain<T>::is_empty() const{
	return (_status==AffineStatus::Empty);
}

template<class T>
inline bool AffineMain<T>::is_degenerated() const {
	//return (itv().diam() <	AF_EC);
	return (_status==AffineStatus::Degenerate);
}

template<class T>
inline bool AffineMain<T>::is_unbounded() const{
	return ((_status==AffineStatus::AllReals)||(_status==AffineStatus::UpperUnbounded)||(_status==AffineStatus::LowerUnbounded));
}

template<class T>
inline AffineMain<T>& AffineMain<T>::operator-=(double d){
	*this += (-d);
	return *this;
}

template<class T>
inline 	AffineMain<T>& AffineMain<T>::operator/=(double d) {
	return *this *= (1.0 / Interval(d)) ;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::operator+=(const Interval& x){
	if ( x.is_empty()) {
		*this = Interval::empty();
	} else if (x.is_unbounded() || (!is_active())) {
		*this = this->itv() + x;
	} else {
		*this += x.mid();
		this->inflate(x.rad());
	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::operator-=(const Interval& x){
	return *this += (-x);
}

template<class T>
inline AffineMain<T>& AffineMain<T>::operator/=(const Interval& x){
	return *this *= (1.0/x) ;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::operator-=(const AffineMain<T>& x){
	return *this += (AffineMain<T>(x).Aneg());
}

template<class T>
inline AffineMain<T>& AffineMain<T>::operator/=(const AffineMain<T>& x){
	return *this *= (AffineMain<T>(x).Ainv(x.itv()));
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Ainv(const Interval& itv){
	switch (mode) {
	case Affine_Mode::AF_Lin_Chebyshev:
		this->Ainv_CH(itv);
		break;
	case Affine_Mode::AF_Lin_MinRange:
		this->Ainv_MR(itv);
		break;
	default:
		*this = 1.0/itv;
		break;
	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Aexp(const Interval& itv){
	switch (mode) {
	case Affine_Mode::AF_Lin_Chebyshev:
		this->Aexp_CH(itv);
		break;
	case Affine_Mode::AF_Lin_MinRange:
		this->Aexp_MR(itv);
		break;
	default:
		*this = exp(itv);
		break;
	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Alog(const Interval& itv){
	switch (mode) {
	case Affine_Mode::AF_Lin_Chebyshev:
		this->Alog_CH(itv);
		break;
	case Affine_Mode::AF_Lin_MinRange:
		this->Alog_MR(itv);
		break;
	default:
		*this = log(itv);
		break;
	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Asqrt(const Interval& itv){
	switch (mode) {
	case Affine_Mode::AF_Lin_Chebyshev:
		this->Asqrt_CH(itv);
		break;
	case Affine_Mode::AF_Lin_MinRange:
		this->Asqrt_MR(itv);
		break;
	default:
		*this = sqrt(itv);
		break;
	}
	return *this;
}


// Chebyshev formula, for f monotonic and either convex or concave on x
// (so that g(t) = f(t) - alpha*t has a single extremum on x, at t=u):
// alpha = (f(ub(x)) - f(lb(x))) / diam(x)     (chord slope)
// u     = (f')^{-1}(alpha)                    (point where the tangent is parallel to the chord)
// d_a   = f(lb(x)) - alpha*lb(x)               (chord offset, endpoint a)
// d_b   = f(ub(x)) - alpha*ub(x)               (chord offset, endpoint b; d_a == d_b in exact arithmetic)
// d_tan = f(u) - alpha*u                       (tangent offset, interior point u)
//
// The chord lies on one side of the graph of f and the tangent at u on the
// other, so [d_a, d_b] (collapsed to a point) and d_tan bracket g over x:
// - f concave (f''<0, e.g. sqrt, log): chord is below the curve, tangent is
//   above it, so d_min = min(d_a, d_b) and d_max = d_tan.
// - f convex (f''>0, e.g. inv on one sign, exp, even powers): the roles
//   swap, chord is above the curve, tangent is below it, so
//   d_min = d_tan and d_max = max(d_a, d_b).
// In both cases:
// beta  = mid([d_min, d_max])
// zeta  = rad([d_min, d_max])
//
// Each routine below picks the ordering matching its own f: Asqrt_CH/Alog_CH
// (concave) place the chord offset first and d_tan last; Ainv_CH/Aexp_CH
// (convex) place d_tan first and the chord offset last.

template<class T>
AffineMain<T>& AffineMain<T>::Ainv_CH(const Interval& itv){
	Interval res_itv = 1.0/(itv);

	// Particular case
	if ((itv.is_unbounded()) || res_itv.is_empty() || res_itv.is_unbounded() || (!is_active()) || (itv.diam() < min_diam_for_slope(itv))) {
		*this = res_itv;
	}  else  {
		// General case
		Interval itv2 =abs(itv);
		if (itv2.diam() < min_diam_for_slope(itv2)) {
			*this = res_itv;
			return *this;
		}
		Interval dmm = (1.0/(itv2));
		const double alpha = -(dmm.diam()/itv2.diam());
		if (!std::isfinite(alpha) || alpha >= 0.0) {
			*this = res_itv;
			return *this;
		}

		// band = evaluate r(x)=f(x)-alpha*x at the bounds of itv2 and at the
		// point u=1/sqrt(-alpha) where f'(u)=alpha, with r(u)=2*sqrt(-alpha)
		Interval band = ((1.0/Interval((itv2).lb()))-alpha*Interval((itv2).lb()))
			 | ((1.0/Interval((itv2).ub()))-alpha*Interval((itv2).ub()));
		band = band | (2*sqrt(-Interval(alpha)));

		double beta = band.mid();
		const double t1 = (beta -band).ub();
		const double t2 = (band -beta).ub();
		const double ddelta = (t1>t2)? t1 : t2;

		if (itv.lb()<0.0) beta = -beta;

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);

	}


	return *this;
}


template<class T>
AffineMain<T>& AffineMain<T>::Asqrt_CH(const Interval& itv){
	Interval itv2 = itv & Interval(0,oo);
	Interval res_itv = sqrt(itv2);

	// Particular case
	if (res_itv.is_empty() || res_itv.is_unbounded() || (!is_active()) || (itv2.diam() < min_diam_for_slope(itv2))) {
		*this = res_itv;
	}  else  {
		// General case
		const double alpha = res_itv.diam()/itv2.diam();
		if (!std::isfinite(alpha) || alpha <= 0.0) {
			*this = res_itv;
			return *this;
		}

		// band = evaluate r(x)=f(x)-alpha*x at the bounds of itv2 and at the
		// point u=1/(4*alpha^2) where f'(u)=alpha, with r(u)=1/(4*alpha)
		Interval band = (sqrt(Interval(itv2.lb()))-alpha*Interval(itv2.lb()))
			 | (sqrt(Interval(itv2.ub()))-alpha*Interval(itv2.ub()));
		band = band | (1.0/(4*Interval(alpha)));

		const double beta = band.mid();
		const double t1 = (beta -band).ub();
		const double t2 = (band-beta).ub();
		const double ddelta = (t1>t2)? t1 : t2;

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);


	}
	return *this;
}

template<class T>
AffineMain<T>& AffineMain<T>::Aexp_CH(const Interval& itv){
	Interval res_itv = exp(itv);

	// Particular case
	if (res_itv.is_empty() || res_itv.is_unbounded() || (!is_active()) || (itv.diam() < min_diam_for_slope(itv))) {
		*this = res_itv;
	}  else  {
		// General case
		const double alpha = res_itv.diam()/itv.diam();
		if (!std::isfinite(alpha) || alpha <= 0.0) {
			*this = res_itv;
			return *this;
		}

		// band = evaluate r(x)=f(x)-alpha*x at the bounds of itv and at the
		// point u=log(alpha) where f'(u)=alpha, with r(u)=alpha*(1-log(alpha))
		Interval band = (exp(Interval(itv.lb()))-alpha*Interval(itv.lb()))
			 | (exp(Interval(itv.ub()))-alpha*Interval(itv.ub()));
		band = band | (alpha*(1-log(Interval(alpha))));

		const double beta = band.mid();
		const double t1 = (beta -band).ub();
		const double t2 = (band -beta).ub();
		const double ddelta = (t1>t2)? t1 : t2;

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);

	}
	return *this;
}

template<class T>
AffineMain<T>& AffineMain<T>::Alog_CH(const Interval& itv){
	Interval res_itv = log(itv);

	// Particular case
	if (res_itv.is_empty() || res_itv.is_unbounded() || (!is_active()) || (itv.diam() < min_diam_for_slope(itv))) {
		*this = res_itv;
	}  else  {
		// General case
		const double alpha = res_itv.diam()/itv.diam();
		if (!std::isfinite(alpha) || alpha <= 0.0) {
			*this = res_itv;
			return *this;
		}

		// band = evaluate r(x)=f(x)-alpha*x at the bounds of itv and at the
		// point u=1/alpha where f'(u)=alpha, with r(u)=-log(alpha)-1
		Interval band = (log(Interval(itv.lb()))-alpha*Interval(itv.lb()))
			 | (log(Interval(itv.ub()))-alpha*Interval(itv.ub()));
		band = band | (-log(Interval(alpha))-1);

		const double beta = band.mid();
		const double t1 = (beta -band).ub();
		const double t2 = (band -beta).ub();
		const double ddelta = (t1>t2)? t1 : t2;

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);


	}
	return *this;
}

template<class T>
AffineMain<T>& AffineMain<T>::Apow(const Interval& y, const Interval& itvx) {
  const Interval reference = pow(itvx, y);
  if (y.is_empty() || y.is_unbounded() || itvx.is_empty() ||
      itvx.is_unbounded()) {
    *this = reference;
    return *this;
  }

  // An interval singleton that is exactly representable as an int must use
  // the integer-power algorithm. This also preserves valid negative bases.
  if (y.is_degenerated() && std::isfinite(y.lb()) &&
      y.lb() >= static_cast<double>(std::numeric_limits<int>::min()) &&
      y.lb() <= static_cast<double>(std::numeric_limits<int>::max()) &&
      std::trunc(y.lb()) == y.lb()) {
    return this->Apow(static_cast<int>(y.lb()), itvx);
  }

  // exp(y*log(x)) is a real-valued identity only for a strictly positive
  // base. For zero or negative bases, use the Interval implementation,
  // whose domain conventions define the public result safely.
  if (itvx.lb() <= 0.0) {
    *this = reference;
    return *this;
  }

  this->Alog(itvx);
  *this *= y;
  this->Aexp(y*log(itvx));
  return *this;
}

template<class T>
AffineMain<T>& AffineMain<T>::Acos(const Interval& itv){
	Interval res_itv = cos(itv);

	// Particular case
	if (res_itv.is_empty() || res_itv.is_unbounded() || (!is_active()) || (itv.diam() < min_diam_for_slope(itv))) {
		*this = res_itv;
	}  else  {
		// General case
		if (itv.diam()>=Interval::two_pi().lb()) {
			*this = Interval(-1,1);
			return *this;
		}
		//  pour _itv = [a,b]
		// x0 = 1/sqrt(2)
		// x1= - x0
		// xb0 = 0.5*((b-a)*x0 +(a+b))
		// xb1 = 0.5*((b-a)*x1 +(a+b))
		// c0 = 0.5 (f(xb0)+f(xb1))
		// c1 = x0*f(xb0)+x1*f(xb1)
		// alpha = 2*c1/(b-a)
		//  old : beta = c0-c1*(a+b)/(b-a)
		//  old : ddelta = (b-a)^2 * f''(_itv)/16
		//  new : residual = evaluate the maximal error of f(x)-alpha*x at the bound and the points when f'(x)=alpha
		//  new : beta = evaluate the mid point of the residual
		//  new : ddelta = evaluate the radius of the residual

		const double x0 = 1.0/std::sqrt(2.);
		const double xb0 = (0.5)*(itv.diam()*x0 +itv.lb()+itv.ub());
		const double xb1 = (0.5)*(itv.diam()*(-x0) +itv.lb()+itv.ub());
		const double fxb0 = std::cos(xb0);
		const double fxb1 = std::cos(xb1);
		//c0 = (0.5)*(fxb0+fxb1);
		const double c1 = x0*fxb0-x0*fxb1;

		const double alpha = 2*c1/(itv.diam());

		// compute the first point such as f'(u) = alpha
		const Interval alpha_domain = Interval(alpha) & Interval(-1.0, 1.0);
		if (alpha_domain.is_empty()) {
			*this = res_itv;
			return *this;
		}
		const Interval u = asin(-alpha_domain);
		const Interval nb_period = (itv+Interval::half_pi()) / Interval::pi();

		// Large arguments cannot be converted safely to an integer period.
		const double period_limit =	static_cast<double>(std::numeric_limits<int>::max() - 2);
		if (!std::isfinite(nb_period.lb()) || !std::isfinite(nb_period.ub()) ||
			nb_period.lb() < -period_limit || nb_period.ub() > period_limit) {
			*this = res_itv;
			return *this;
		}

		// compute the residual r(x)=f(x)-alpha*x at the bounds of itv
		Interval residual = (cos(Interval(itv.lb()))-alpha*Interval(itv.lb()))
						  | (cos(Interval(itv.ub()))-alpha*Interval(itv.ub()));

		// evaluate the residual at the points such that f'(u) = alpha
		int p1 = static_cast<int>(nb_period.lb()) - 2;
		int p2 = static_cast<int>(nb_period.ub()) + 2;

		int i = p1;
		Interval TEMP1(0.0);
		while (i<=p2) { // looking for a point
			TEMP1 = (itv & (i%2==0? (u + i*Interval::pi()) : (i*Interval::pi() - u)));
			if (!(TEMP1.is_empty())) { // check if maximize the error
				residual = residual | (cos(TEMP1)-alpha*TEMP1);
			}
			i++;
		}

		const double beta = residual.mid();
		const double t1 = (beta-residual).ub();
		const double t2 = (residual-beta).ub();
		const double ddelta = (t1>t2) ? t1 : t2;

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);

	}
	return *this;
}

template<class T>
AffineMain<T>& AffineMain<T>::Asin(const Interval& itv){
	Interval res_itv = sin(itv);

	// Particular case
	if (res_itv.is_empty() || res_itv.is_unbounded() || (!is_active()) || (itv.diam() < min_diam_for_slope(itv))) {
		*this = res_itv;
	}  else  {
		// General case
		if (itv.diam()>=Interval::two_pi().lb()) {
			*this = Interval(-1,1);
			return *this;
		}
		//  pour _itv = [a,b]
		// x0 = 1/sqrt(2)
		// x1= - x0
		// xb0 = 0.5*((b-a)*x0 +(a+b))
		// xb1 = 0.5*((b-a)*x1 +(a+b))
		// c0 = 0.5 (f(xb0)+f(xb1))
		// c1 = x0*f(xb0)+x1*f(xb1)
		// alpha = 2*c1/(b-a)
		//  old : beta = c0-c1*(a+b)/(b-a)
		//  old : ddelta = (b-a)^2 * f''(_itv)/16
		//  new : residual = evaluate the maximal error of f(x)-alpha*x at the bound and the points when f'(x)=alpha
		//  new : beta = evaluate the mid point of the residual
		//  new : ddelta = evaluate the radius of the residual
		const double x0 = 1.0/std::sqrt(2.);
		const double xb0 = (0.5)*(itv.diam()*x0 +itv.lb()+itv.ub());
		const double xb1 = (0.5)*(itv.diam()*(-x0) +itv.lb()+itv.ub());
		const double fxb0 = std::sin(xb0);
		const double fxb1 = std::sin(xb1);
		//c0 = (0.5)*(fxb0+fxb1);
		const double c1 = x0*fxb0-x0*fxb1;

		const double alpha = 2*c1/(itv.diam());

		// compute the first point such as f'(u) = alpha
		const Interval alpha_domain =Interval(alpha) & Interval(-1.0, 1.0);
		if (alpha_domain.is_empty()) {
			*this = res_itv;
			return *this;
		}
		const Interval u = acos(alpha_domain);
		const Interval nb_period = itv / Interval::pi();

		const double period_limit =	static_cast<double>(std::numeric_limits<int>::max() - 2);
		if (!std::isfinite(nb_period.lb()) || !std::isfinite(nb_period.ub()) ||
			nb_period.lb() < -period_limit || nb_period.ub() > period_limit) {
			*this = res_itv;
			return *this;
		}

		// compute the residual r(x)=f(x)-alpha*x at the bounds of itv
		Interval residual = (sin(Interval(itv.lb()))-alpha*Interval(itv.lb()))
						  | (sin(Interval(itv.ub()))-alpha*Interval(itv.ub()));

		// evaluate the residual at the points such that f'(u) = alpha
		int p1 = static_cast<int>(nb_period.lb()) - 2;
		int p2 = static_cast<int>(nb_period.ub()) + 2;

		int i = p1;
		Interval TEMP1(0.0);
		while (i<=p2) { // looking for a point
			TEMP1 = (itv & (i%2==0? (u + i*Interval::pi()) : ((i+1)*Interval::pi() - u)));
			if (!(TEMP1.is_empty())) {
				residual = residual | (sin(TEMP1)-alpha*TEMP1);
			}
			i++;
		}

		const double beta = residual.mid();
		const double t1 = (beta-residual).ub();
		const double t2 = (residual-beta).ub();
		const double ddelta = (t1>t2) ? t1 : t2;

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);


	}
	return *this;
}

template<class T>
AffineMain<T>& AffineMain<T>::Atan(const Interval& itv){
	Interval res_itv = tan(itv);

	// Particular case
	if (res_itv.is_empty() || res_itv.is_unbounded() || (!is_active()) || (itv.diam() < min_diam_for_slope(itv))) {
		*this = res_itv;
	}  else  {
		// General case.
		// An interval wider than one period necessarily contains a pole, so
		// tan(itv) is unbounded and the case has normally been handled by the
		// guard above. Keep a defensive fallback rather than an assertion, so
		// that an interval library returning a bounded value here would only
		// lose precision, never soundness.
		if (itv.diam()>=Interval::two_pi().lb()) {
			*this = res_itv;
			return *this;
		}
		//  pour _itv = [a,b]
		// x0 = 1/sqrt(2)
		// x1= - x0
		// xb0 = 0.5*((b-a)*x0 +(a+b))
		// xb1 = 0.5*((b-a)*x1 +(a+b))
		// c0 = 0.5 (f(xb0)+f(xb1))
		// c1 = x0*f(xb0)+x1*f(xb1)
		// alpha = 2*c1/(b-a)
		//  old : beta = c0-c1*(a+b)/(b-a)
		//  old : ddelta = (b-a)^2 * f''(_itv)/16
		//  new : residual = evaluate the maximal error of f(x)-alpha*x at the bound and the points when f'(x)=alpha
		//  new : beta = evaluate the mid point of the residual
		//  new : ddelta = evaluate the radius of the residual
		const double x0 = 1.0/std::sqrt(2.);
		const double xb0 = (0.5)*(itv.diam()*  x0  +itv.lb()+itv.ub());
		const double xb1 = (0.5)*(itv.diam()*(-x0) +itv.lb()+itv.ub());
		const double fxb0 = std::tan(xb0);
		const double fxb1 = std::tan(xb1);
		//c0 = (0.5)*(fxb0+fxb1);
		const double c1 = x0*fxb0-x0*fxb1;

		const double alpha = 2*c1/(itv.diam());

		// compute the first point such as f'(u) = alpha
		if (!std::isfinite(alpha) || alpha < 1.0) {
			*this = res_itv;
			return *this;
		}
		const Interval cosine_domain =(Interval::one()/sqrt(Interval(alpha))) & Interval(0.0, 1.0);
		if (cosine_domain.is_empty()) {
			*this = res_itv;
			return *this;
		}
		const Interval u = acos(cosine_domain);
		const Interval nb_period = itv / Interval::pi();

		const double period_limit = static_cast<double>(std::numeric_limits<int>::max() - 2);
		if (!std::isfinite(nb_period.lb()) || !std::isfinite(nb_period.ub()) ||
			nb_period.lb() < -period_limit || nb_period.ub() > period_limit) {
			*this = res_itv;
			return *this;
		}

		// compute the residual r(x)=f(x)-alpha*x at the bounds of itv
		Interval residual = (tan(Interval(itv.lb()))-alpha*Interval(itv.lb()))
						  | (tan(Interval(itv.ub()))-alpha*Interval(itv.ub()));

		// evaluate the residual at the points such that f'(u) = alpha
		int p1 = static_cast<int>(nb_period.lb()) - 2;
		int p2 = static_cast<int>(nb_period.ub()) + 2;

		int i = p1;
		Interval TEMP1(0.0);
		while (i<=p2) { // looking for a point
			TEMP1 = (itv & ( i*Interval::pi() + u));
			if ((!(TEMP1.is_empty()))) {
				residual = residual | (tan(TEMP1)-alpha*TEMP1);
			}
			TEMP1 = (itv & ( i*Interval::pi() - u ));
			if ((!(TEMP1.is_empty()))) {
				residual = residual | (tan(TEMP1)-alpha*TEMP1);
			}
			i++;
		}

		const double beta = residual.mid();
		const double t1 = (beta-residual).ub();
		const double t2 = (residual-beta).ub();
		const double ddelta = (t1>t2) ? t1 : t2;

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);

	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Aacos(const Interval& itv){
	switch (mode) {
	case Affine_Mode::AF_Lin_Chebyshev:
		this->Aacos_CH(itv);
		break;
	case Affine_Mode::AF_Lin_MinRange:
		this->Aacos_MR(itv);
		break;
	default:
		*this = acos(itv);
		break;
	}
	return *this;
}

template<class T>
AffineMain<T>& AffineMain<T>::Aacos_CH(const Interval& itv) {
	const Interval domain = itv & Interval(-1.0, 1.0);
	const Interval res_itv = acos(domain);
	if (res_itv.is_empty() || res_itv.is_unbounded() || !is_active() ||
		domain.diam() < min_diam_for_slope(domain)) {
		*this = res_itv;
	} else {

		//  pour _itv = [a,b]
		// x0 = 1/sqrt(2)
		// x1= - x0
		// xb0 = 0.5*((b-a)*x0 +(a+b))
		// xb1 = 0.5*((b-a)*x1 +(a+b))
		// c0 = 0.5 (f(xb0)+f(xb1))
		// c1 = x0*f(xb0)+x1*f(xb1)
		// alpha = 2*c1/(b-a)
		//  old : beta = c0-c1*(a+b)/(b-a)
		//  old : ddelta = (b-a)^2 * f''(_itv)/16
		//  new : residual = evaluate the maximal error  of f(x)-alpha*x at the bound and the points when f'(x)=alpha
		//  new : beta = evaluate the mid point of the rseidual
		//  new : ddelta = evaluate the radius of the residual
		const double x0 = 1.0/std::sqrt(2.0);
		const double xb0 = 0.5*(domain.diam()*x0 + domain.lb() + domain.ub());
		const double xb1 = 0.5*(-domain.diam()*x0 + domain.lb() + domain.ub());
		const double fxb0 = std::acos(xb0);
		const double fxb1 = std::acos(xb1);
		//const double c0 = 0.5*(fxb0 + fxb1);
		const double c1 = x0*(fxb0 - fxb1);
		const double alpha = 2.0*c1/domain.diam();

		if (!std::isfinite(alpha) || alpha > -1.0) {
			*this = res_itv;
			return *this;
		}

		const Interval left  =  (acos(Interval(domain.lb())) - alpha*Interval(domain.lb()));
		const Interval right =  (acos(Interval(domain.ub())) - alpha*Interval(domain.ub())); 
		Interval residual = left | right;

		// acos'(u)=-1/sqrt(1-u^2)=alpha gives u=+-sqrt(1-1/alpha^2).
		const Interval critical = sqrt(1.0 - 1.0/sqr(Interval(alpha)));
		Interval TEMP1(0.0);
		for (const Interval& candidate : {critical, -critical}) {
			TEMP1 = candidate & domain;
			if (!TEMP1.is_empty()) {
				residual = residual | (acos(TEMP1) - alpha*TEMP1);
			}
		}

		const double beta = residual.mid();
		const double t1 = (beta-residual).ub();
		const double t2 = (residual-beta).ub();
		const double ddelta =  (t1>t2) ? t1 : t2;
		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Aasin(const Interval& itv){
	switch (mode) {
	case Affine_Mode::AF_Lin_Chebyshev:
		this->Aasin_CH(itv);
		break;
	case Affine_Mode::AF_Lin_MinRange:
		this->Aasin_MR(itv);
		break;
	default:
		*this = asin(itv);
		break;
	}
	return *this;
}

template<class T>
AffineMain<T>& AffineMain<T>::Aasin_CH(const Interval& itv) {
	const Interval domain = itv & Interval(-1.0, 1.0);
	const Interval res_itv = asin(domain);
	if (res_itv.is_empty() || res_itv.is_unbounded() || !is_active() ||
		domain.diam() < min_diam_for_slope(domain)) {
		*this = res_itv;
		return *this;
	}

	//  pour _itv = [a,b]
	// x0 = 1/sqrt(2)
	// x1= - x0
	// xb0 = 0.5*((b-a)*x0 +(a+b))
	// xb1 = 0.5*((b-a)*x1 +(a+b))
	// c0 = 0.5 (f(xb0)+f(xb1))
	// c1 = x0*f(xb0)+x1*f(xb1)
	// alpha = 2*c1/(b-a)
	//  old : beta = c0-c1*(a+b)/(b-a)
	//  old : ddelta = (b-a)^2 * f''(_itv)/16
	//  new : residual = evaluate the maximal error  of f(x)-alpha*x at the bound and the points when f'(x)=alpha
	//  new : beta = evaluate the mid point of the rseidual
	//  new : ddelta = evaluate the radius of the residual
	const double x0 = 1.0/std::sqrt(2.0);
	const double xb0 = 0.5*(domain.diam()*x0 + domain.lb() + domain.ub());
	const double xb1 = 0.5*(-domain.diam()*x0 + domain.lb() + domain.ub());
	const double fxb0 = std::asin(xb0);
	const double fxb1 = std::asin(xb1);
	//const double c0 = 0.5*(fxb0 + fxb1);
	const double c1 = x0*(fxb0 - fxb1);
	const double alpha = 2.0*c1/domain.diam();

	if (!std::isfinite(alpha) || alpha < 1.0) {
		*this = res_itv;
		return *this;
	}

	const Interval left  =  (asin(Interval(domain.lb())) - alpha*Interval(domain.lb()));
	const Interval right =  (asin(Interval(domain.ub())) - alpha*Interval(domain.ub()));
	Interval residual = left | right;

	// asin'(u)=1/sqrt(1-u^2)=alpha gives u=+-sqrt(1-1/alpha^2).
	const Interval critical = sqrt(1.0 - 1.0/sqr(Interval(alpha)));
	Interval TEMP1(0.0);
	for (const Interval& candidate : {critical, -critical}) {
		TEMP1  = candidate & domain;
		if (!TEMP1.is_empty()) {
			residual = residual | (asin(TEMP1) - alpha*TEMP1);
		}
	}

	const double beta = residual.mid();
	const double t1 = (beta-residual).ub();
	const double t2 = (residual-beta).ub();
	const double ddelta =  (t1>t2) ? t1 : t2;
	*this *= alpha;
	*this += beta;
	this->inflate(ddelta);
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Aatan(const Interval& itv){
	switch (mode) {
	case Affine_Mode::AF_Lin_Chebyshev:
		this->Aatan_CH(itv);
		break;
	case Affine_Mode::AF_Lin_MinRange:
		this->Aatan_MR(itv);
		break;
	default:
		*this = atan(itv);
		break;
	}
	return *this;
}

template<class T>
AffineMain<T>& AffineMain<T>::Aatan_CH(const Interval& itv) {
	const Interval res_itv = atan(itv);
	if (itv.is_unbounded() || res_itv.is_empty() || res_itv.is_unbounded() ||
		!is_active() || itv.diam() < min_diam_for_slope(itv)) {
		*this = res_itv;
		return *this;
	}

	//  pour _itv = [a,b]
	// x0 = 1/sqrt(2)
	// x1= - x0
	// xb0 = 0.5*((b-a)*x0 +(a+b))
	// xb1 = 0.5*((b-a)*x1 +(a+b))
	// c0 = 0.5 (f(xb0)+f(xb1))
	// c1 = x0*f(xb0)+x1*f(xb1)
	// alpha = 2*c1/(b-a)
	//  old : beta = c0-c1*(a+b)/(b-a)
	//  old : ddelta = (b-a)^2 * f''(_itv)/16
	//  new : residual = evaluate the maximal error  of f(x)-alpha*x at the bound and the points when f'(x)=alpha
	//  new : beta = evaluate the mid point of the rseidual
	//  new : ddelta = evaluate the radius of the residual
	const double x0 = 1.0/std::sqrt(2.0);
	const double xb0 = 0.5*(itv.diam()*x0 + itv.lb() + itv.ub());
	const double xb1 = 0.5*(-itv.diam()*x0 + itv.lb() + itv.ub());
	const double fxb0 = std::atan(xb0);
	const double fxb1 = std::atan(xb1);
	//const double c0 = 0.5*(fxb0 + fxb1);
	const double c1 = x0*(fxb0 - fxb1);
	const double alpha = 2.0*c1/itv.diam();

	if (!std::isfinite(alpha) || alpha <= 0.0 || alpha > 1.0) {
		*this = res_itv;
		return *this;
	}

	const Interval left  =  (atan(Interval(itv.lb())) - alpha*Interval(itv.lb()));
	const Interval right =  (atan(Interval(itv.ub())) - alpha*Interval(itv.ub()));
	Interval residual = left | right;

	// atan'(u)=1/(1+u^2)=alpha gives u=+-sqrt(1/alpha-1).
	const Interval critical = sqrt(1.0/Interval(alpha) - 1.0);
	Interval TEMP1(0.0);
	for (const Interval& candidate : {critical, -critical}) {
		TEMP1  = candidate & itv;
		if (!TEMP1.is_empty()) {
			residual = residual | (atan(TEMP1) - alpha*TEMP1);
		}
	}

	const double beta = residual.mid();
	const double t1 = (beta-residual).ub();
	const double t2 = (residual-beta).ub();
	const double ddelta =  (t1>t2) ? t1 : t2;
	*this *= alpha;
	*this += beta;
	this->inflate(ddelta);
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Acosh(const Interval& itv){
	switch (mode) {
	case Affine_Mode::AF_Lin_Chebyshev:
		this->Acosh_CH(itv);
		break;
	case Affine_Mode::AF_Lin_MinRange:
		this->Acosh_MR(itv);
		break;
	default:
		*this = cosh(itv);
		break;
	}
	return *this;
}

template<class T>
AffineMain<T>& AffineMain<T>::Acosh_CH(const Interval& itv){
	Interval res_itv = cosh(itv);

	// Particular case
	if (res_itv.is_empty() || res_itv.is_unbounded() || (!is_active()) || (itv.diam() < min_diam_for_slope(itv))) {
		*this = res_itv;
	}  else  {
		// General case
		const double alpha = ((cosh(Interval(itv.ub()))-cosh(Interval(itv.lb())))/itv.diam()).lb();

		// band = evaluate r(x)=f(x)-alpha*x at the bounds of itv and at the
		// point u=asinh(alpha) where f'(u)=alpha. As cosh(asinh(alpha)) =
		// sqrt(sqr(alpha)+1), r(u) = sqrt(sqr(alpha)+1)-alpha*asinh(alpha)
		Interval band = (cosh(Interval(itv.lb()))-alpha*Interval(itv.lb()))
					 | (cosh(Interval(itv.ub()))-alpha*Interval(itv.ub()));
		band = band | (sqrt(sqr(Interval(alpha))+1)-alpha*asinh(Interval(alpha)));

		const double beta = band.mid();
		const double t1 = (beta -band).ub();
		const double t2 = (band -beta).ub();
		const double ddelta = (t1>t2)? t1 : t2;

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);


	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Asinh(const Interval& itv){
	switch (mode) {
	case Affine_Mode::AF_Lin_Chebyshev:
		this->Asinh_CH(itv);
		break;
	case Affine_Mode::AF_Lin_MinRange:
		this->Asinh_MR(itv);
		break;
	default:
		*this = sinh(itv);
		break;
	}
	return *this;
}

template<class T>
AffineMain<T>& AffineMain<T>::Asinh_CH(const Interval& itv) {
	const Interval res_itv = sinh(itv);
	if (res_itv.is_empty() || res_itv.is_unbounded() || !is_active() ||
		itv.diam() < min_diam_for_slope(itv)) {
		*this = res_itv;
		return *this;
	}

	//  pour _itv = [a,b]
	// x0 = 1/sqrt(2)
	// x1= - x0
	// xb0 = 0.5*((b-a)*x0 +(a+b))
	// xb1 = 0.5*((b-a)*x1 +(a+b))
	// c0 = 0.5 (f(xb0)+f(xb1))
	// c1 = x0*f(xb0)+x1*f(xb1)
	// alpha = 2*c1/(b-a)
	//  old : beta = c0-c1*(a+b)/(b-a)
	//  old : ddelta = (b-a)^2 * f''(_itv)/16
	//  new : residual = evaluate the maximal error  of f(x)-alpha*x at the bound and the points when f'(x)=alpha
	//  new : beta = evaluate the mid point of the rseidual
	//  new : ddelta = evaluate the radius of the residual
	const double x0 = 1.0/std::sqrt(2.0);
	const double xb0 = 0.5*(itv.diam()*x0 + itv.lb() + itv.ub());
	const double xb1 = 0.5*(-itv.diam()*x0 + itv.lb() + itv.ub());
	const double fxb0 = std::sinh(xb0);
	const double fxb1 = std::sinh(xb1);
	//const double c0 = 0.5*(fxb0 + fxb1);
	const double c1 = x0*(fxb0 - fxb1);
	const double alpha = 2.0*c1/itv.diam();

	if (!std::isfinite(alpha) || alpha < 1.0) {
		*this = res_itv;
		return *this;
	}

	const Interval left  =  (sinh(Interval(itv.lb())) - alpha*Interval(itv.lb()));
	const Interval right =  (sinh(Interval(itv.ub())) - alpha*Interval(itv.ub()));
	Interval residual = left | right;

	// sinh'(u)=cosh(u)=alpha gives u=+-acosh(alpha).
	const Interval critical = acosh(Interval(alpha));
	Interval TEMP1(0.0);
	for (const Interval& candidate : {critical, -critical}) {
		TEMP1 = candidate & itv;
		if (!TEMP1.is_empty()) {
			residual = residual | (sinh(TEMP1) - alpha*TEMP1);
		}
	}

	const double beta = residual.mid();
	const double t1 = (beta-residual).ub();
	const double t2 = (residual-beta).ub();
	const double ddelta =  (t1>t2) ? t1 : t2;
	*this *= alpha;
	*this += beta;
	this->inflate(ddelta);
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Atanh(const Interval& itv){
	switch (mode) {
	case Affine_Mode::AF_Lin_Chebyshev:
		this->Atanh_CH(itv);
		break;
	case Affine_Mode::AF_Lin_MinRange:
		this->Atanh_MR(itv);
		break;
	default:
		*this = tanh(itv);
		break;
	}
	return *this;
}

template<class T>
AffineMain<T>& AffineMain<T>::Atanh_CH(const Interval& itv) {
	const Interval res_itv = tanh(itv);
	if (itv.is_unbounded() || res_itv.is_empty() || res_itv.is_unbounded() ||
		!is_active() || itv.diam() < min_diam_for_slope(itv)) {
		*this = res_itv;
		return *this;
	}

	//  pour _itv = [a,b]
	// x0 = 1/sqrt(2)
	// x1= - x0
	// xb0 = 0.5*((b-a)*x0 +(a+b))
	// xb1 = 0.5*((b-a)*x1 +(a+b))
	// c0 = 0.5 (f(xb0)+f(xb1))
	// c1 = x0*f(xb0)+x1*f(xb1)
	// alpha = 2*c1/(b-a)
	//  old : beta = c0-c1*(a+b)/(b-a)
	//  old : ddelta = (b-a)^2 * f''(_itv)/16
	//  new : residual = evaluate the maximal error  of f(x)-alpha*x at the bound and the points when f'(x)=alpha
	//  new : beta = evaluate the mid point of the rseidual
	//  new : ddelta = evaluate the radius of the residual
	const double x0 = 1.0/std::sqrt(2.0);
	const double xb0 = 0.5*(itv.diam()*x0 + itv.lb() + itv.ub());
	const double xb1 = 0.5*(-itv.diam()*x0 + itv.lb() + itv.ub());
	const double fxb0 = std::tanh(xb0);
	const double fxb1 = std::tanh(xb1);
	//const double c0 = 0.5*(fxb0 + fxb1);
	const double c1 = x0*(fxb0 - fxb1);
	const double alpha = 2.0*c1/itv.diam();

	if (!std::isfinite(alpha) || alpha <= 0.0 || alpha > 1.0) {
		*this = res_itv;
		return *this;
	}

	const Interval left  =  (tanh(Interval(itv.lb())) - alpha*Interval(itv.lb()));
	const Interval right =  (tanh(Interval(itv.ub())) - alpha*Interval(itv.ub()));
	Interval residual = left | right;

	// tanh'(u)=1/cosh(u)^2=alpha gives u=+-acosh(1/sqrt(alpha)).
	const Interval critical = acosh(1.0/sqrt(Interval(alpha)));
	Interval TEMP1(0.0);
	for (const Interval& candidate : {critical, -critical}) {
		TEMP1 = candidate & itv;
		if (!TEMP1.is_empty()) {
			residual = residual | (tanh(TEMP1) - alpha*TEMP1);
		}
	}

	const double beta = residual.mid();
	const double t1 = (beta-residual).ub();
	const double t2 = (residual-beta).ub();
	const double ddelta =  (t1>t2) ? t1 : t2;
	*this *= alpha;
	*this += beta;
	this->inflate(ddelta);
	return *this;
}


// acosh is concave on [1,+oo) (acosh''(x) = -x/(x^2-1)^{3/2} < 0), so it
// follows the same Chebyshev pattern as Asqrt_CH: band = [chord intercept
// (min), tangent intercept (max)]. Since acosh is strictly increasing,
// alpha = diam(res_itv)/diam(domain) is always > 0.
template<class T>
inline AffineMain<T>& AffineMain<T>::Aacosh(const Interval& itv){
	switch (mode) {
	case Affine_Mode::AF_Lin_Chebyshev:
		this->Aacosh_CH(itv);
		break;
	case Affine_Mode::AF_Lin_MinRange:
		this->Aacosh_MR(itv);
		break;
	default:
		*this = acosh(itv);
		break;
	}
	return *this;
}

template<class T>
AffineMain<T>& AffineMain<T>::Aacosh_CH(const Interval& itv ) {
	const Interval domain = itv & Interval(1.0, oo);
	const Interval res_itv = acosh(domain);

	// The model is built on the real domain intersection. Tests below check
	// that applying it to the original affine form still encloses all values
	// belonging to that valid domain when the input is only partly admissible.
	if (res_itv.is_empty() || res_itv.is_unbounded() || (!is_active()) || (domain.diam() < min_diam_for_slope(domain))) {
		*this = res_itv;
		return *this;
	}

	const double alpha = res_itv.diam()/domain.diam();
	if (!std::isfinite(alpha) || alpha <= 0.0) {
		*this = res_itv;
		return *this;
	}

	const Interval left = acosh(Interval(domain.lb())) - alpha*Interval(domain.lb());
	const Interval right = acosh(Interval(domain.ub())) - alpha*Interval(domain.ub());
	Interval band = left | right;

	// acosh'(u)=alpha gives u=sqrt(1+1/alpha^2).
	const Interval candidate = sqrt(1.0 + 1.0/sqr(Interval(alpha))) & domain;
	if (!candidate.is_empty()) {
		band = band | (acosh(candidate) - alpha*candidate);
	}

	const double beta = band.mid();
	const double t1 = (beta-band).ub();
	const double t2 = (band-beta).ub();
	const double ddelta =  (t1>t2) ? t1 : t2;
	*this *= alpha;
	*this += beta;
	this->inflate(ddelta);
	return *this;
}


template<class T>
inline AffineMain<T>& AffineMain<T>::Aasinh(const Interval& itv){
	switch (mode) {
	case Affine_Mode::AF_Lin_Chebyshev:
		this->Aasinh_CH(itv);
		break;
	case Affine_Mode::AF_Lin_MinRange:
		this->Aasinh_MR(itv);
		break;
	default:
		*this = asinh(itv);
		break;
	}
	return *this;
}

template<class T>
AffineMain<T>& AffineMain<T>::Aasinh_CH(const Interval& itv) {
	const Interval res_itv = asinh(itv);
	if (res_itv.is_empty() || res_itv.is_unbounded() || !is_active() ||
		itv.diam() < min_diam_for_slope(itv)) {
		*this = res_itv;
		return *this;
	}

	//  pour _itv = [a,b]
	// x0 = 1/sqrt(2)
	// x1= - x0
	// xb0 = 0.5*((b-a)*x0 +(a+b))
	// xb1 = 0.5*((b-a)*x1 +(a+b))
	// c0 = 0.5 (f(xb0)+f(xb1))
	// c1 = x0*f(xb0)+x1*f(xb1)
	// alpha = 2*c1/(b-a)
	//  old : beta = c0-c1*(a+b)/(b-a)
	//  old : ddelta = (b-a)^2 * f''(_itv)/16
	//  new : residual = evaluate the maximal error  of f(x)-alpha*x at the bound and the points when f'(x)=alpha
	//  new : beta = evaluate the mid point of the rseidual
	//  new : ddelta = evaluate the radius of the residual
	const double x0 = 1.0/std::sqrt(2.0);
	const double xb0 = 0.5*(itv.diam()*x0 + itv.lb() + itv.ub());
	const double xb1 = 0.5*(-itv.diam()*x0 + itv.lb() + itv.ub());
	const double fxb0 = std::asinh(xb0);
	const double fxb1 = std::asinh(xb1);
	//const double c0 = 0.5*(fxb0 + fxb1);
	const double c1 = x0*(fxb0 - fxb1);
	const double alpha = 2.0*c1/itv.diam();

	if (!std::isfinite(alpha) || alpha <= 0.0 || alpha > 1.0) {
		*this = res_itv;
		return *this;
	}

	const Interval left  =  (asinh(Interval(itv.lb())) - alpha*Interval(itv.lb()));
	const Interval right =  (asinh(Interval(itv.ub())) - alpha*Interval(itv.ub()));
	Interval residual = left | right;

	// asinh'(u)=alpha gives u=+-sqrt(1/alpha^2-1).
	const Interval critical = sqrt(1.0/sqr(Interval(alpha)) - 1.0);
	Interval TEMP1(0.0);
	for (const Interval& candidate : {critical, -critical}) {
		TEMP1 = candidate & itv;
		if (!TEMP1.is_empty()) {
			residual = residual | (asinh(TEMP1) - alpha*TEMP1);
		}
	}

	const double beta = residual.mid();
	const double t1 = (beta-residual).ub();
	const double t2 = (residual-beta).ub();
	const double ddelta =  (t1>t2) ? t1 : t2;
	*this *= alpha;
	*this += beta;
	this->inflate(ddelta);
	return *this;
}


// atanh also has an inflection point at x=0 (same S-shape as tanh/asinh),
// but unlike them its derivative 1/(1-x^2) is minimal (=1) at x=0 and grows
// without bound towards the domain edges +-1, so alpha is always >= 1
// (mean value theorem), guaranteeing 1-1/alpha >= 0 below.
template<class T>
inline AffineMain<T>& AffineMain<T>::Aatanh(const Interval& itv){
	switch (mode) {
	case Affine_Mode::AF_Lin_Chebyshev:
		this->Aatanh_CH(itv);
		break;
	case Affine_Mode::AF_Lin_MinRange:
		this->Aatanh_MR(itv);
		break;
	default:
		*this = atanh(itv);
		break;
	}
	return *this;
}

template<class T>
inline AffineMain<T>&
AffineMain<T>::Aatanh_CH(const Interval& itv) {
	const Interval domain = itv & Interval(-1.0, 1.0);
	const Interval res_itv = atanh(itv);
	if (domain != itv || res_itv.is_empty() || res_itv.is_unbounded() ||
		!is_active() || itv.diam() < min_diam_for_slope(itv)) {
		*this = res_itv;
		return *this;
	}

	const double x0 = 1.0/std::sqrt(2.0);
	const double xb0 = 0.5*(itv.diam()*x0 + itv.lb() + itv.ub());
	const double xb1 = 0.5*(-itv.diam()*x0 + itv.lb() + itv.ub());
	const double fxb0 = std::atanh(xb0);
	const double fxb1 = std::atanh(xb1);
	//const double c0 = 0.5*(fxb0 + fxb1);
	const double c1 = x0*(fxb0 - fxb1);
	const double alpha = 2.0*c1/itv.diam();

	if (!std::isfinite(alpha) || alpha < 1.0) {
		*this = res_itv;
		return *this;
	}

	const Interval left  =  (atanh(Interval(itv.lb())) - alpha*Interval(itv.lb()));
	const Interval right =  (atanh(Interval(itv.ub())) - alpha*Interval(itv.ub()));
	Interval residual = left | right;

	// atanh'(u)=alpha gives u=+-sqrt(1-1/alpha).
	const Interval critical = sqrt(1.0 - 1.0/Interval(alpha));
	Interval TEMP1(0.0);
	for (const Interval& candidate : {critical, -critical}) {
		TEMP1 = candidate & itv;
		if (!TEMP1.is_empty()) {
			residual = residual | (atanh(TEMP1) - alpha*TEMP1);
		}
	}

	const double beta = residual.mid();
	const double t1 = (beta-residual).ub();
	const double t2 = (residual-beta).ub();
	const double ddelta =  (t1>t2) ? t1 : t2;
	*this *= alpha;
	*this += beta;
	this->inflate(ddelta);
	return *this;
}


template<class T>
inline AffineMain<T>& AffineMain<T>::Aatan2(const AffineMain<T>& x, const Interval& itvY, const Interval& itvX){
    if(itvY.is_empty() || itvX.is_empty()){
		*this = Interval::empty();
		return *this;
	}
    // We handle the special case x=[0,0] separately
    else if(itvX == Interval::zero())  {
      if(itvY.lb() >= 0) {
        if(itvY.ub() == 0)
          *this = Interval::empty(); // atan2(0,0) is undefined
        else 
          *this = Interval::half_pi();
      }
      else if(itvY.ub() <= 0) 
        *this = (-Interval::half_pi());
      else
        *this = Interval(-1,1)*Interval::half_pi();
	  return *this;
    }
    else if(itvX.lb() >= 0) {
		// Zone sûre (x>0 sur toute la boîte): atan2(y,x) == atan(y/x)
    	// exactement, donc on réutilise la division affine et l'atan affine
    	// déjà validés, sans nouvelle linéarisation.
		*this /= x;                 // *this contient AF[y] -> devient AF[y]/AF[x]
    	this->Aatan(this->itv());
      	return *this; // now, x.ub()>0 -> atan does not give an empty set
	}
    else if(itvX.ub() <= 0) {
      if(itvY.lb() >= 0) {
		*this /= x;                 // *this contient AF[y] -> devient AF[y]/AF[x]
    	this->Aatan(this->itv());
		*this += Interval::pi(); // x.ub()<0
      	return *this;
		}
      else if(itvY.ub() < 0) {
		*this /= x;         
    	this->Aatan(this->itv());
		*this -= Interval::pi(); 
      	return *this;
		}
      else {
        *this = Interval(-1,1)*Interval::pi();
		return *this;
	  }
    }
    else  {     
	// Cas singulier ou coupure de branche potentielle: on ne tente pas de
    // préserver la corrélation affine, on redescend sur l'Interval.
		Interval res_itv = atan2(itvY, itvX);
        *this = res_itv;
        return *this;
    }
  }


template<class T>
inline AffineMain<T>& AffineMain<T>::Aabs(const Interval& itv){
	Interval res_itv = abs(itv);

	// Particular case
	if (0<=itv.lb()) return *this;
	if (itv.ub()<=0) {
		this->Aneg();
		return *this;
	}
	if (res_itv.is_empty() || res_itv.is_unbounded() || (!is_active()) || (itv.diam() < min_diam_for_slope(itv))) {
		*this = res_itv;
	}  else  {
		// General case

		const double alpha = ((abs(Interval(itv.ub()))-abs(Interval(itv.lb())))/itv.diam()).ub();

		// band = evaluate r(x)=f(x)-alpha*x at the bounds of itv and at the
		// only non-differentiable point u=0, which lies in itv here and
		// gives the minimum r(0)=0
		Interval band = (abs(Interval(itv.lb()))-alpha*Interval(itv.lb()))
			 | (abs(Interval(itv.ub()))-alpha*Interval(itv.ub()));
		band = band | Interval(0.0);

		const double beta = band.mid();
		const double t1 = (beta-band).ub();
		const double t2 = (band-beta).ub();
		const double ddelta =  (t1>t2) ? t1 : t2;
		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);
	}
	return *this;
}


template<class T>
inline AffineMain<T>& AffineMain<T>::Apow(int n, const Interval& itv) {
	//	std::cout << "in power "<<std::endl;
	// Particular case
	if (itv.is_empty()||itv.is_unbounded()) {
		*this = pow(itv,n);
	} else if (!is_active()) {
		*this = pow(itv,n);
	} else if (itv.diam() < min_diam_for_slope(itv)) {
		*this = pow(itv,n);
	} else {
		// General Case
		if (n == 0) {
			*this = Interval::one();
		} else if (n == 1)
			return *this;
		else if (n == 2)
			this->Asqr(itv);
		else if (n<0) {
			if (n == std::numeric_limits<int>::min()) {
				*this = pow(itv,n);
			} else {
				this->Apow(-n,itv);
				this->Ainv(pow(itv,-n));
			}

		} else {
			// General case, n>=3: the linearization depends on the mode
			switch (mode) {
			case Affine_Mode::AF_Lin_Chebyshev:
				this->Apow_CH(n,itv);
				break;
			case Affine_Mode::AF_Lin_MinRange:
				this->Apow_MR(n,itv);
				break;
			default:
				*this = pow(itv,n);
				break;
			}
		}

	}
	//	std::cout << "out power "<<std::endl;
	return *this;
}


template<class T>
inline AffineMain<T>& AffineMain<T>::Apow_CH(int n, const Interval& itv) {
	if (n % 2 == 0) {
		// alpha = (f(sup(x)) - f(inf(x)))/diam(x)
		// f is convex for n even, so the residual r(x)=f(x)-alpha*x is
		// maximal at a bound of x and minimal at the point u where
		// f'(u)=alpha, i.e. u = (alpha/n)^(1/(n-1))
		// band   = evaluate r at the bounds of x and at u
		// beta   = evaluate the mid point of the band
		// ddelta = evaluate the radius of the band
		const double f_lb = std::pow(itv.lb(), static_cast<unsigned int>(n));
		const double f_ub = std::pow(itv.ub(), static_cast<unsigned int>(n));
		const double alpha = (f_ub - f_lb)/itv.diam();
		if (!std::isfinite(f_lb) || !std::isfinite(f_ub) ||	!std::isfinite(alpha)) {
			*this = pow(itv, n);
			return *this;
		}

		Interval band =	(pow(Interval(itv.lb()), n) - alpha*Interval(itv.lb())) |
						(pow(Interval(itv.ub()), n) - alpha*Interval(itv.ub()));

		// u = (alpha/n)^(1/(n-1)) always belongs to itv (mean value
		// theorem), and r(u) = u^n - alpha*u = (1-n)*(alpha/n)*u
		const Interval critical = Interval(alpha) / n;
		band = band | ((1 - n) * critical * (root(critical, n - 1)));

		const double beta = band.mid();
		const double t1 = (beta-band).ub();
		const double t2 = (band-beta).ub();
		const double ddelta =  (t1>t2) ? t1 : t2;
		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);

	} else {
		// for _itv = [a,b]
		// x0 = 1/sqrt(2)
		// x1= - x0
		// xb0 = 0.5*((b-a)*x0 +(a+b))
		// xb1 = 0.5*((b-a)*x1 +(a+b))
		// c0 = 0.5 (f(xb0)+f(xb1))
		// c1 = x0*f(xb0)+x1*f(xb1)
		// alpha = 2*c1/(b-a)
		//  old : beta = c0-c1*(a+b)/(b-a)
		//  old : ddelta = (b-a)^2 * f''(_itv)/16
		//  new : residual = evaluate the maximal error of f(x)-alpha*x at the bound and the points when f'(x)=alpha
		//  new : beta = evaluate the mid point of the residual
		//  new : ddelta = evaluate the radius of the residual

		const double x0  = 1.0/std::sqrt(2.);
		const double xb0 = (0.5) * (itv.diam() * ( x0) + itv.lb() + itv.ub());
		const double xb1 = (0.5) * (itv.diam() * (-x0) + itv.lb() + itv.ub());
		const double fxb0 = std::pow(xb0, static_cast<unsigned int>(n));
		const double fxb1 = std::pow(xb1, static_cast<unsigned int>(n));
		//const double c0 = (0.5) * (fxb0 + fxb1);
		const double c1 = x0 * fxb0 - x0 * fxb1;

		const double alpha = 2 * c1 / (itv.diam());
		if (!std::isfinite(fxb0) || !std::isfinite(fxb1) ||	!std::isfinite(alpha)) {
			*this = pow(itv, n);
			return *this;
		}

		// compute the residual r(x)=f(x)-alpha*x at the bounds of itv
		Interval residual =	(pow(Interval(itv.lb()), n) - alpha*Interval(itv.lb())) |
							(pow(Interval(itv.ub()), n) - alpha*Interval(itv.ub()));

		// f'(u)=n*u^(n-1)=alpha gives u=+-(alpha/n)^(1/(n-1)).
		const Interval critical = pow(Interval(alpha) / n, 1.0 / Interval(n - 1));
		Interval TEMP1(0.0);
		for (const Interval& candidate : {critical, -critical}) {
			TEMP1 = candidate & itv;
			if (!TEMP1.is_empty()) {
				residual = residual | (pow(TEMP1, n) - alpha*TEMP1);
			}
		}

		const double beta = residual.mid();
		const double t1 = (beta-residual).ub();
		const double t2 = (residual-beta).ub();
		const double ddelta =  (t1>t2) ? t1 : t2;
		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);
	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Apow_MR(int n, const Interval& itv) {

	Interval res_itv = pow(itv,n);
	if (itv.is_unbounded() || res_itv.is_empty() || res_itv.is_unbounded() || (!is_active()) || (itv.diam() < min_diam_for_slope(itv))) {
		*this = res_itv;
	}  else  {  // _actif && b
		double alpha = 0.0;
		Interval band(0.0);
		// f'(x)=n*x^(n-1)
		if (n % 2 == 1) {
			// n-1 is even, so f is increasing and its slope is smallest at the
			// point of itv closest to zero, rounded downward
			alpha = (n*pow(Interval(itv.mig()),n-1)).lb();
		} else if (0.0 <= itv.lb()) {
			// f is increasing on itv: the smallest slope is at the lower bound
			alpha = (n*pow(Interval(itv.lb()),n-1)).lb();
		} else if (itv.ub() <= 0.0) {
			// f is decreasing on itv: the slope of smallest magnitude is at the
			// upper bound and is negative, rounded upward
			alpha = (n*pow(Interval(itv.ub()),n-1)).ub();
		} else {
			// the minimum of an even power is reached inside itv, so the
			// smallest slope is zero
			alpha = 0.0;
		}
		if (!std::isfinite(alpha) || alpha == 0.0) {
			alpha = 0.0;
			band = res_itv;
		} else {
			// alpha is the slope of smallest magnitude on itv, so r(x)=f(x)-alpha*x
			// is monotonic there: the band is the hull of r at both bounds
			band = (pow(Interval(itv.lb()),n) -alpha*Interval(itv.lb()))
				 | (pow(Interval(itv.ub()),n) -alpha*Interval(itv.ub()));
		}

		const double beta = band.mid();
		const double t1 = (beta -band).ub();
		const double t2 = (band -beta).ub();
		const double ddelta = (t1>t2)? t1 : t2;

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);

	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Aroot_MR(int n, const Interval& itv) {

	Interval res_itv = root(itv,n);
	if (itv.is_unbounded() || res_itv.is_empty() || res_itv.is_unbounded() || (!is_active()) || (itv.diam() < min_diam_for_slope(itv))) {
		*this = res_itv;
	}  else  {  // _actif && b
		double alpha = 0.0;
		Interval band(0.0);
		// root'(x)=|x|^(1/n-1)/n decreases with |x| for an odd index, so the
		// smallest slope is reached at the bound of itv farthest from zero,
		// rounded downward. It is written root(m,n)/(n*m) to avoid a negative
		// exponent.
		const Interval m(itv.mag());
		alpha = (root(m,n)/(n*m)).lb();  // compute the derivative
		if (!std::isfinite(alpha) || alpha <= 0.0) {
			alpha = 0.0;
			band = res_itv;
		} else {
			// alpha is the smallest slope on itv, so r(x)=f(x)-alpha*x is
			// monotonic there: the band is the hull of r at both bounds
			band = (root(Interval(itv.lb()),n) -alpha*Interval(itv.lb()))
				 | (root(Interval(itv.ub()),n) -alpha*Interval(itv.ub()));
		}

		const double beta = band.mid();
		const double t1 = (beta -band).ub();
		const double t2 = (band -beta).ub();
		const double ddelta = (t1>t2)? t1 : t2;

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);

	}
	return *this;
}

template<class T>
inline AffineMain<T>& AffineMain<T>::Apow(double d, const Interval& itv) {
	if (!std::isfinite(d)) {
		*this = pow(itv, d);
		return *this;
	}
	if (d >= static_cast<double>(std::numeric_limits<int>::min()) &&
	    d <= static_cast<double>(std::numeric_limits<int>::max()) &&
	    std::trunc(d) == d) {
	  this->Apow(static_cast<int>(d), itv);
	  return *this;
	} else if (d<0) {
		this->Apow(Interval(-d), itv);
		this->Ainv(pow(itv,-d));
		return *this;
	} else {
		this->Apow(Interval(d),itv);
		return *this;
	}
}


template<class T>
inline AffineMain<T>& AffineMain<T>::Aroot(int n, const Interval& itv) {

	if (is_empty()) return *this;
	else if (n==0)  return *this = Interval::empty();
	else if (n==1)  return *this;
	else if (is_degenerated()) {
		return *this = root(itv,n);
	}
	else if (n<0) {
		if (n == std::numeric_limits<int>::min()) {
			// |n| is not representable as an int, and root(itv,n) computes -n,
			// which overflows. As n is even here, the root is only defined on
			// the non-negative part of itv, where the identity
			// root(x,n) = exp(log(x)/n) is rigorous.
			return *this = exp(log(itv & Interval(0.0,oo))
					/ Interval(static_cast<double>(n)));
		}
		this->Aroot(-n,itv);
		this->Ainv(root(itv,-n));
		return *this;
	}
	else if (n % 2 == 0) {
		this->Apow(Interval::one()/n,itv);
		return *this; // the negative part of x should be removed
	}
	else if (mode == Affine_Mode::AF_Lin_MinRange) {
		// an odd root is increasing on the whole of R, so a single
		// MinRange linearization covers the three cases below
		return this->Aroot_MR(n,itv);
	}
	else if (0 <= itv.lb()) {
		this->Apow(Interval::one()/n,itv);
		return  *this;
	}
	else if (itv.ub() <= 0) {
		this->Aneg();
		this->Apow(Interval::one()/n,-itv);
		this->Aneg();
		return  *this;
	}
	else {
		// Odd root across zero. The derivative is singular at zero, but the
		// extrema of r(x)=root(x,n)-alpha*x are the endpoints, zero,
		// and the two points where root'(x)=alpha.
		const Interval root_itv = root(itv, n);
		const double alpha = root_itv.diam()/itv.diam();
		if (!std::isfinite(alpha) || alpha <= 0.0) {
			return *this = root_itv;
		}

		Interval residual = (root(Interval(itv.lb()), n) - alpha*Interval(itv.lb())) |
							(root(Interval(itv.ub()), n) - alpha*Interval(itv.ub()));
		residual = residual | Interval(0.0); // residual at x=0

		// root'(u)=alpha gives u=+-(1/(n*alpha))^(n/(n-1))
		const double exponent = static_cast<double>(n)/static_cast<double>(n-1);
		const Interval critical = pow(Interval::one()/(static_cast<double>(n)*alpha), exponent);

		Interval TEMP1(0.0);
		for (const Interval& candidate : {critical, -critical}) {
			TEMP1 = candidate & itv;
			if (!TEMP1.is_empty()) {
				residual = residual | (root(TEMP1, n) - alpha*TEMP1);
			}
		}

		const double beta = residual.mid();
		const double t1 = (beta-residual).ub();
		const double t2 = (residual-beta).ub();
		const double ddelta =  (t1>t2) ? t1 : t2;

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		return *this;
	}

}

//////////////////////////////////////////////////////////////////////////////////
// MinRange formula


// For x = [a,b]:
// if f is increasing:
//   alpha = f'(a)
//   band  = [f(a)-alpha*a, f(b)-alpha*b]
// if f is decreasing:
//   alpha = f'(b)
//   band  = [f(a)-alpha*b, f(b)-alpha*a]
// then:
//   beta   = mid(band)
//   ddelta = rad(band)

template<class T>
inline AffineMain<T>& AffineMain<T>::Asqrt_MR(const Interval& itv) {

	Interval res_itv = sqrt(itv);
	if (itv.is_unbounded() || res_itv.is_empty() || res_itv.is_unbounded() || (!is_active()) || (itv.diam() < min_diam_for_slope(itv))) {
		*this = res_itv;
	}  else  {  // _actif && b
		const Interval itv2 = itv & Interval(0.0, oo);
		const Interval dmm = sqrt(itv2);
		Interval band(0.0);
		double alpha;
		if (itv2.diam() < min_diam_for_slope(itv2)) {
			alpha = 0.0;
			band =dmm;
		} else {
			alpha = (1.0/(2.0*dmm)).lb();  // compute the derivative
			if (alpha<=0) {
				alpha = 0.0;
				band = dmm;
			} else {
				// alpha is the smallest slope on itv2, so r(x)=f(x)-alpha*x is
				// monotonic there: the band is the hull of r at both bounds
				band = (sqrt(Interval(itv2.lb())) -alpha*Interval(itv2.lb()))
					 | (sqrt(Interval(itv2.ub())) -alpha*Interval(itv2.ub()));
			}
		}

		const double beta = band.mid();
		const double t1 = (beta -band).ub();
		const double t2 = (band -beta).ub();
		const double ddelta = (t1>t2)? t1 : t2;

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);

	}
	return *this;
}

template<class T>
AffineMain<T>& AffineMain<T>::Aexp_MR(const Interval& itv) {

	Interval res_itv = exp(itv);
	if (itv.is_unbounded() || res_itv.is_empty() || res_itv.is_unbounded() || (!is_active()) || (itv.diam() < min_diam_for_slope(itv))) {
		*this = res_itv;
	}  else  {  // _actif && b
		double alpha = 0.0;
		Interval band(0.0);
		if (itv.diam() < min_diam_for_slope(itv)) {
			alpha = 0.0;
			band = res_itv;
		} else {
			alpha = res_itv.lb();  // compute the derivative
			if (alpha<=0) {
				alpha = 0.0;
				band = res_itv;
			} else {
				// alpha is the smallest slope on itv, so r(x)=f(x)-alpha*x is
				// monotonic there: the band is the hull of r at both bounds
				band = (exp(Interval(itv.lb())) -alpha*Interval(itv.lb()))
					 | (exp(Interval(itv.ub())) -alpha*Interval(itv.ub()));
			}
		}

		const double beta = band.mid();
		const double t1 = (beta -band).ub();
		const double t2 = (band -beta).ub();
		const double ddelta = (t1>t2)? t1 : t2;

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);

	}
	return *this;
}

template<class T>
AffineMain<T>& AffineMain<T>::Alog_MR(const Interval& itv) {

	Interval res_itv = log(itv);
	if (itv.is_unbounded() || res_itv.is_empty() || res_itv.is_unbounded() || (!is_active()) || (itv.diam() < min_diam_for_slope(itv))) {
		*this = res_itv;
	}  else  {  // _actif && b

		double alpha = 0.0;
		Interval band(0.0);
		if (itv.diam() < min_diam_for_slope(itv)) {
			alpha = 0.0;
			band =res_itv;
		} else {
			alpha = (1.0/itv).lb();  // compute the derivative
			if (alpha<=0) {
				alpha = 0.0;
				band = res_itv;
			} else {
				// alpha is the smallest slope on itv, so r(x)=f(x)-alpha*x is
				// monotonic there: the band is the hull of r at both bounds
				band = (log(Interval(itv.lb())) -alpha*Interval(itv.lb()))
					 | (log(Interval(itv.ub())) -alpha*Interval(itv.ub()));
			}
		}
		const double beta = band.mid();
		const double t1 = (beta -band).ub();
		const double t2 = (band -beta).ub();
		const double ddelta = (t1>t2)? t1 : t2;

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);
		//saxpy(alpha, AffineMain<T>(), beta, ddelta, true,false,true,true);

		//}

	}
	return *this;
}

template<class T>
AffineMain<T>& AffineMain<T>::Acosh_MR(const Interval& itv) {

	Interval res_itv = cosh(itv);
	if (itv.is_unbounded() || res_itv.is_empty() || res_itv.is_unbounded() || (!is_active()) || (itv.diam() < min_diam_for_slope(itv))) {
		*this = res_itv;
	}  else  {  // _actif && b
		double alpha = 0.0;
		Interval band(0.0);
		if (itv.diam() < min_diam_for_slope(itv)) {
			alpha = 0.0;
		} else if (0.0 <= itv.lb()) {
			// cosh is increasing on itv: the smallest slope is at the lower
			// bound, rounded downward so that sinh(x)>=alpha holds on all itv
			alpha = sinh(Interval(itv.lb())).lb();  // compute the derivative
		} else if (itv.ub() <= 0.0) {
			// cosh is decreasing on itv: the slope of smallest magnitude is at
			// the upper bound and is negative, rounded upward so that
			// sinh(x)<=alpha holds on all itv
			alpha = sinh(Interval(itv.ub())).ub();  // compute the derivative
		} else {
			// the minimum of cosh is reached inside itv, so the smallest slope
			// is zero and the affine form degenerates into the interval image
			alpha = 0.0;
		}
		if (alpha == 0.0) {
			band = res_itv;
		} else {
			// alpha is the smallest slope on itv, so r(x)=f(x)-alpha*x is
			// monotonic there: the band is the hull of r at both bounds
			band = (cosh(Interval(itv.lb())) -alpha*Interval(itv.lb()))
				 | (cosh(Interval(itv.ub())) -alpha*Interval(itv.ub()));
		}

		const double beta = band.mid();
		const double t1 = (beta -band).ub();
		const double t2 = (band -beta).ub();
		const double ddelta = (t1>t2)? t1 : t2;

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);

	}
	return *this;
}

template<class T>
AffineMain<T>& AffineMain<T>::Atanh_MR(const Interval& itv) {

	Interval res_itv = tanh(itv);
	if (itv.is_unbounded() || res_itv.is_empty() || res_itv.is_unbounded() || (!is_active()) || (itv.diam() < min_diam_for_slope(itv))) {
		*this = res_itv;
	}  else  {  // _actif && b
		double alpha = 0.0;
		Interval band(0.0);
		if (itv.diam() < min_diam_for_slope(itv)) {
			alpha = 0.0;
			band = res_itv;
		} else {
			// tanh'(x)=1-tanh(x)^2 decreases with |x|, so the smallest slope is
			// reached at the bound of itv farthest from zero, rounded downward
			// so that tanh'(x)>=alpha holds on all itv
			alpha = (1.0-sqr(tanh(Interval(itv.mag())))).lb();  // compute the derivative
			if (!std::isfinite(alpha) || alpha <= 0.0) {
				// a slope rounded down to zero or below leaves a constant band
				alpha = 0.0;
				band = res_itv;
			} else {
				// alpha is the smallest slope on itv, so r(x)=f(x)-alpha*x is
				// monotonic there: the band is the hull of r at both bounds
				band = (tanh(Interval(itv.lb())) -alpha*Interval(itv.lb()))
					 | (tanh(Interval(itv.ub())) -alpha*Interval(itv.ub()));
			}
		}

		const double beta = band.mid();
		const double t1 = (beta -band).ub();
		const double t2 = (band -beta).ub();
		const double ddelta = (t1>t2)? t1 : t2;

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);

	}
	return *this;
}

template<class T>
AffineMain<T>& AffineMain<T>::Aatan_MR(const Interval& itv) {

	Interval res_itv = atan(itv);
	if (itv.is_unbounded() || res_itv.is_empty() || res_itv.is_unbounded() || (!is_active()) || (itv.diam() < min_diam_for_slope(itv))) {
		*this = res_itv;
	}  else  {  // _actif && b
		double alpha = 0.0;
		Interval band(0.0);
		if (itv.diam() < min_diam_for_slope(itv)) {
			alpha = 0.0;
			band = res_itv;
		} else {
			// atan'(x)=1/(1+x^2) decreases with |x|, so the smallest slope is
			// reached at the bound of itv farthest from zero, rounded downward
			// so that atan'(x)>=alpha holds on all itv
			alpha = (1.0/(1.0+sqr(Interval(itv.mag())))).lb();  // compute the derivative
			if (!std::isfinite(alpha) || alpha <= 0.0) {
				// a slope rounded down to zero or below leaves a constant band
				alpha = 0.0;
				band = res_itv;
			} else {
				// alpha is the smallest slope on itv, so r(x)=f(x)-alpha*x is
				// monotonic there: the band is the hull of r at both bounds
				band = (atan(Interval(itv.lb())) -alpha*Interval(itv.lb()))
					 | (atan(Interval(itv.ub())) -alpha*Interval(itv.ub()));
			}
		}

		const double beta = band.mid();
		const double t1 = (beta -band).ub();
		const double t2 = (band -beta).ub();
		const double ddelta = (t1>t2)? t1 : t2;

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);

	}
	return *this;
}

template<class T>
AffineMain<T>& AffineMain<T>::Aasin_MR(const Interval& itv) {

	const Interval domain = itv & Interval(-1.0,1.0);
	Interval res_itv = asin(domain);
	if (res_itv.is_empty() || res_itv.is_unbounded() || (!is_active()) || (domain.diam() < min_diam_for_slope(domain))) {
		*this = res_itv;
	}  else  {  // _actif && b
		double alpha = 0.0;
		Interval band(0.0);
		if (domain.diam() < min_diam_for_slope(domain)) {
			alpha = 0.0;
			band = res_itv;
		} else {
			// asin'(x)=1/sqrt(1-x^2) grows with |x|, so the smallest slope is
			// reached at the point of itv closest to zero, rounded downward
			// so that asin'(x)>=alpha holds on all itv
			alpha = (1.0/sqrt(1.0-sqr(Interval(domain.mig())))).lb();  // compute the derivative
			if (!std::isfinite(alpha) || alpha == 0.0) {
				alpha = 0.0;
				band = res_itv;
			} else {
				// alpha is the slope of smallest magnitude on itv, so
				// r(x)=f(x)-alpha*x is monotonic there: the band is the hull
				// of r at both bounds
				band = (asin(Interval(domain.lb())) -alpha*Interval(domain.lb()))
					 | (asin(Interval(domain.ub())) -alpha*Interval(domain.ub()));
			}
		}

		const double beta = band.mid();
		const double t1 = (beta -band).ub();
		const double t2 = (band -beta).ub();
		const double ddelta = (t1>t2)? t1 : t2;

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);

	}
	return *this;
}

template<class T>
AffineMain<T>& AffineMain<T>::Aacos_MR(const Interval& itv) {

	const Interval domain = itv & Interval(-1.0,1.0);
	Interval res_itv = acos(domain);
	if (res_itv.is_empty() || res_itv.is_unbounded() || (!is_active()) || (domain.diam() < min_diam_for_slope(domain))) {
		*this = res_itv;
	}  else  {  // _actif && b
		double alpha = 0.0;
		Interval band(0.0);
		if (domain.diam() < min_diam_for_slope(domain)) {
			alpha = 0.0;
			band = res_itv;
		} else {
			// acos'(x)=-1/sqrt(1-x^2) is negative and its magnitude grows with
			// |x|, so the slope of smallest magnitude is reached at the point of
			// itv closest to zero, rounded upward so that acos'(x)<=alpha on all itv
			alpha = (-1.0/sqrt(1.0-sqr(Interval(domain.mig())))).ub();  // compute the derivative
			if (!std::isfinite(alpha) || alpha == 0.0) {
				alpha = 0.0;
				band = res_itv;
			} else {
				// alpha is the slope of smallest magnitude on itv, so
				// r(x)=f(x)-alpha*x is monotonic there: the band is the hull
				// of r at both bounds
				band = (acos(Interval(domain.lb())) -alpha*Interval(domain.lb()))
					 | (acos(Interval(domain.ub())) -alpha*Interval(domain.ub()));
			}
		}

		const double beta = band.mid();
		const double t1 = (beta -band).ub();
		const double t2 = (band -beta).ub();
		const double ddelta = (t1>t2)? t1 : t2;

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);

	}
	return *this;
}

template<class T>
AffineMain<T>& AffineMain<T>::Asinh_MR(const Interval& itv) {

	Interval res_itv = sinh(itv);
	if (itv.is_unbounded() || res_itv.is_empty() || res_itv.is_unbounded() || (!is_active()) || (itv.diam() < min_diam_for_slope(itv))) {
		*this = res_itv;
	}  else  {  // _actif && b
		double alpha = 0.0;
		Interval band(0.0);
		// sinh'(x)=cosh(x) grows with |x|, so the smallest slope is reached at
		// the point of itv closest to zero, rounded downward so that
		// sinh'(x)>=alpha holds on all itv
		alpha = cosh(Interval(itv.mig())).lb();  // compute the derivative
		if (!std::isfinite(alpha) || alpha <= 0.0) {
			// a slope rounded down to zero or below leaves a constant band
			alpha = 0.0;
			band = res_itv;
		} else {
			// alpha is the smallest slope on the domain, so r(x)=f(x)-alpha*x
			// is monotonic there: the band is the hull of r at both bounds
			band = (sinh(Interval(itv.lb())) -alpha*Interval(itv.lb()))
				 | (sinh(Interval(itv.ub())) -alpha*Interval(itv.ub()));
		}

		const double beta = band.mid();
		const double t1 = (beta -band).ub();
		const double t2 = (band -beta).ub();
		const double ddelta = (t1>t2)? t1 : t2;

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);

	}
	return *this;
}

template<class T>
AffineMain<T>& AffineMain<T>::Aasinh_MR(const Interval& itv) {

	Interval res_itv = asinh(itv);
	if (itv.is_unbounded() || res_itv.is_empty() || res_itv.is_unbounded() || (!is_active()) || (itv.diam() < min_diam_for_slope(itv))) {
		*this = res_itv;
	}  else  {  // _actif && b
		double alpha = 0.0;
		Interval band(0.0);
		// asinh'(x)=1/sqrt(1+x^2) decreases with |x|, so the smallest slope is
		// reached at the bound of itv farthest from zero, rounded downward
		alpha = (1.0/sqrt(1.0+sqr(Interval(itv.mag())))).lb();  // compute the derivative
		if (!std::isfinite(alpha) || alpha <= 0.0) {
			// a slope rounded down to zero or below leaves a constant band
			alpha = 0.0;
			band = res_itv;
		} else {
			// alpha is the smallest slope on the domain, so r(x)=f(x)-alpha*x
			// is monotonic there: the band is the hull of r at both bounds
			band = (asinh(Interval(itv.lb())) -alpha*Interval(itv.lb()))
				 | (asinh(Interval(itv.ub())) -alpha*Interval(itv.ub()));
		}

		const double beta = band.mid();
		const double t1 = (beta -band).ub();
		const double t2 = (band -beta).ub();
		const double ddelta = (t1>t2)? t1 : t2;

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);

	}
	return *this;
}

template<class T>
AffineMain<T>& AffineMain<T>::Aacosh_MR(const Interval& itv) {

	const Interval domain = itv & Interval(1.0, oo);
	Interval res_itv = acosh(domain);
	if (res_itv.is_empty() || res_itv.is_unbounded() || (!is_active()) || (domain.diam() < min_diam_for_slope(domain))) {
		*this = res_itv;
	}  else  {  // _actif && b
		double alpha = 0.0;
		Interval band(0.0);
		// acosh'(x)=1/sqrt(x^2-1) decreases with x, so the smallest slope is
		// reached at the upper bound of the domain, rounded downward
		alpha = (1.0/sqrt(sqr(Interval(domain.ub()))-1.0)).lb();  // compute the derivative
		if (!std::isfinite(alpha) || alpha <= 0.0) {
			// a slope rounded down to zero or below leaves a constant band
			alpha = 0.0;
			band = res_itv;
		} else {
			// alpha is the smallest slope on the domain, so r(x)=f(x)-alpha*x
			// is monotonic there: the band is the hull of r at both bounds
			band = (acosh(Interval(domain.lb())) -alpha*Interval(domain.lb()))
				 | (acosh(Interval(domain.ub())) -alpha*Interval(domain.ub()));
		}

		const double beta = band.mid();
		const double t1 = (beta -band).ub();
		const double t2 = (band -beta).ub();
		const double ddelta = (t1>t2)? t1 : t2;

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);

	}
	return *this;
}

template<class T>
AffineMain<T>& AffineMain<T>::Aatanh_MR(const Interval& itv) {

	const Interval domain = itv & Interval(-1.0, 1.0);
	Interval res_itv = atanh(domain);
	if (res_itv.is_empty() || res_itv.is_unbounded() || (!is_active()) || (domain.diam() < min_diam_for_slope(domain))) {
		*this = res_itv;
	}  else  {  // _actif && b
		double alpha = 0.0;
		Interval band(0.0);
		// atanh'(x)=1/(1-x^2) grows with |x|, so the smallest slope is reached
		// at the point of the domain closest to zero, rounded downward
		alpha = (1.0/(1.0-sqr(Interval(domain.mig())))).lb();  // compute the derivative
		if (!std::isfinite(alpha) || alpha <= 0.0) {
			// a slope rounded down to zero or below leaves a constant band
			alpha = 0.0;
			band = res_itv;
		} else {
			// alpha is the smallest slope on the domain, so r(x)=f(x)-alpha*x
			// is monotonic there: the band is the hull of r at both bounds
			band = (atanh(Interval(domain.lb())) -alpha*Interval(domain.lb()))
				 | (atanh(Interval(domain.ub())) -alpha*Interval(domain.ub()));
		}

		const double beta = band.mid();
		const double t1 = (beta -band).ub();
		const double t2 = (band -beta).ub();
		const double ddelta = (t1>t2)? t1 : t2;

		*this *= alpha;
		*this += beta;
		this->inflate(ddelta);

	}
	return *this;
}

template<class T>
AffineMain<T>& AffineMain<T>::Ainv_MR(const Interval& itv) {
	const Interval res_itv = 1.0/itv;
	if (itv.is_unbounded() || res_itv.is_empty() || res_itv.is_unbounded() || !is_active() || itv.diam() < min_diam_for_slope(itv)) {
		*this = res_itv;
		return *this;
	}

	// On either side of zero, 1/x is decreasing. The endpoint farthest from
	// zero has the derivative of smallest magnitude and gives the MinRange
	// slope. Keep the negative slope instead of collapsing it to zero.
	const double endpoint = (itv.lb() > 0.0) ? itv.ub() : itv.lb();
	const double alpha = (-1.0/sqr(Interval(endpoint))).mid();
	if (!std::isfinite(alpha) || alpha >= 0.0) {
		*this = res_itv;
		return *this;
	}

	Interval band = (1.0/Interval(itv.lb()) - alpha*Interval(itv.lb())) |
					(1.0/Interval(itv.ub()) - alpha*Interval(itv.ub()));

	// If rounding moved the stationary point slightly inside the interval,
	// include it explicitly in the residual range.
	const Interval stationary_abs = sqrt(-1.0/Interval(alpha));
	Interval TEMP1(0.0);
	for (const Interval& candidate : {stationary_abs, -stationary_abs}) {
		TEMP1 = candidate & itv;
		if (!TEMP1.is_empty()) {
			band = band | (1.0/TEMP1 - alpha*TEMP1);
		}
	}

	const double beta = band.mid();
	const double t1 = (beta -band).ub();
	const double t2 = (band -beta).ub();
	const double ddelta = (t1>t2)? t1 : t2;

	*this *= alpha;
	*this += beta;
	this->inflate(ddelta);
	return *this;
}


} // end namespace codac2
