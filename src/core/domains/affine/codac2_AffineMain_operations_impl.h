/**
 *  \file codac2_AffineMain_operations_impl.h
 *
 *  Inline implementations of the arithmetic operators and mathematical
 *  functions on AffineMain declared in codac2_AffineMain.h.
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
inline Interval operator&(const AffineMain<T>& x1, const AffineMain<T>& x2) {
	return ((x1.itv())&(x2.itv()));
}

template<class T>
inline Interval operator&(const Interval& x1, const AffineMain<T>& x2) {
	return ((x1)&(x2.itv()));
}

template<class T>
inline Interval operator&(const AffineMain<T>& x1, const Interval& x2) {
	return ((x1.itv())&(x2));
}

template<class T>
inline Interval operator|(const AffineMain<T>& x1, const AffineMain<T>& x2) {
	return ((x1.itv())|(x2.itv()));
}

template<class T>
inline Interval operator|(const Interval& x1, const AffineMain<T>& x2) {
	return ((x1)|(x2.itv()));
}

template<class T>
inline Interval operator|(const AffineMain<T>& x1, const Interval& x2) {
	return ((x1.itv())|(x2));
}

template<class T>
inline BoolInterval operator<(const Interval& x, const AffineMain<T>& y) {
    return x < y.itv();
}

template<class T>
inline BoolInterval operator>(const Interval& x, const AffineMain<T>& y){
    return x > y.itv();
}


template<class T>
inline AffineMain<T> operator+(const AffineMain<T>& x1, const AffineMain<T>& x2){
	return AffineMain<T>(x1) += x2;
}

template<class T>
inline AffineMain<T> operator+(const AffineMain<T>& x, double d){
	return AffineMain<T>(x) += d;
}

template<class T>
inline AffineMain<T> operator+(double d, const AffineMain<T>& x){
	return AffineMain<T>(x) += d;
}

template<class T>
inline AffineMain<T> operator+(const AffineMain<T>& x1, const Interval& x2){
	return AffineMain<T>(x1) += x2;
}

template<class T>
inline AffineMain<T> operator+(const Interval& x1, const AffineMain<T>& x2){
	return AffineMain<T>(x2) += x1;
}

template<class T>
inline AffineMain<T> operator-(const AffineMain<T>& x1, const AffineMain<T>& x2){
	return AffineMain<T>(x1) += (-x2);
}

template<class T>
inline AffineMain<T> operator-(const AffineMain<T>& x, double d){
	return AffineMain<T>(x) -= d;
}

template<class T>
inline AffineMain<T> operator-(double d, const AffineMain<T>& x){
	AffineMain<T> res = (-x);
	return res += d;
}

template<class T>
inline AffineMain<T> operator-(const AffineMain<T>& x1, const Interval& x2) {
	return AffineMain<T>(x1) -= x2;
}

template<class T>
inline AffineMain<T> operator-(const Interval& x1, const AffineMain<T>& x2) {
	AffineMain<T> res = (- x2);
	return res += x1;
}

template<class T>
inline AffineMain<T> operator*(const AffineMain<T>& x1, const AffineMain<T>& x2) {
	return AffineMain<T>(x1) *= x2;
}

template<class T>
inline AffineMain<T> operator*(const AffineMain<T>& x, double d){
	return AffineMain<T>(x) *= d;
}

template<class T>
inline AffineMain<T> operator*(double d, const AffineMain<T>& x){
	return AffineMain<T>(x) *= d;
}

template<class T>
inline AffineMain<T> operator*(const AffineMain<T>& x1, const Interval& x2){
	return AffineMain<T>(x1) *= x2;
}

template<class T>
inline AffineMain<T> operator*(const Interval& x1, const AffineMain<T>& x2){
	return AffineMain<T>(x2) *= x1;
}

template<class T>
inline AffineMain<T> operator/(const AffineMain<T>& x1, const AffineMain<T>& x2){
	return AffineMain<T>(x1) /= x2;
}

template<class T>
inline AffineMain<T> operator/(const AffineMain<T>& x, double d){
	return AffineMain<T>(x) /= d;
}

template<class T>
inline AffineMain<T> operator/(double d, const AffineMain<T>& x){
	AffineMain<T> out;
	out = d;
	return out *= (AffineMain<T>(x).Ainv(x.itv()));
}

template<class T>
inline AffineMain<T> operator/(const AffineMain<T>& x1, const Interval& x2){
	return AffineMain<T>(x1) /= x2;
}

template<class T>
inline AffineMain<T> operator/(const Interval& x1, const AffineMain<T>& x2){
	return AffineMain<T>(x1) *= (AffineMain<T>(x2).Ainv(x2.itv()));
}

template<class T>
inline AffineMain<T> inv(const AffineMain<T>& x){
	AffineMain<T> out(x);
	out.Ainv(x.itv());
	return out;
}

template<class T>
inline AffineMain<T> sqr(const AffineMain<T>& x){
	AffineMain<T> out(x);
	out.Asqr(x.itv());
	return out;
}

template<class T>
inline AffineMain<T> sqrt(const AffineMain<T>& x){
	AffineMain<T> out(x);
	out.Asqrt(x.itv());
	return out;
}

template<class T>
inline AffineMain<T> exp(const AffineMain<T>& x){
	AffineMain<T> out(x);
	out.Aexp(x.itv());
	return out;
}

template<class T>
inline AffineMain<T> log(const AffineMain<T>& x){
	AffineMain<T> out(x);
	out.Alog(x.itv());
	return out;
}

template<class T>
inline AffineMain<T> pow(const AffineMain<T>& x, int n) {
	AffineMain<T> out(x);
	out.Apow(n,x.itv());
	return out;
}
template<class T>
inline AffineMain<T> pow(const AffineMain<T>& x, double d){
	AffineMain<T> out(x);
	out.Apow(d,x.itv());
	return out;
}
template<class T>
inline AffineMain<T> pow(const AffineMain<T> &x, const Interval &y){
	// return exp(y * log(x));
	AffineMain<T> out(x);
	out.Apow(y,x.itv());
	return out;
}
template<class T>
inline AffineMain<T> pow(const AffineMain<T>& x, const AffineMain<T>& y)
{
  // Deliberate design choice: the affine dependency structure of the
  // exponent y is not preserved. The exponent is evaluated through y.itv()
  // as an independent interval. The affine structure of the base x is kept
  // whenever Apow(const Interval&, const Interval&) can safely preserve it.
  AffineMain<T> out(x);
  out.Apow(y.itv(), x.itv());
  return out;
}
template<class T>
inline AffineMain<T> pow(const Interval& x, const AffineMain<T>& y)
{
  // Deliberate design choice: the affine dependency structure of the
  // exponent y is not preserved. The exponent is evaluated through y.itv()
  // as an independent interval. The affine structure of the base x is kept
  // whenever Apow(const Interval&, const Interval&) can safely preserve it.
  AffineMain<T> out(x);
  out.Apow(y.itv(), x);
  return out;
}

template<class T>
inline AffineMain<T> pow(double x, const AffineMain<T>& y)
{
  // Deliberate design choice: the affine dependency structure of the
  // exponent y is not preserved. The exponent is evaluated through y.itv()
  // as an independent interval. The affine structure of the base x is kept
  // whenever Apow(const Interval&, const Interval&) can safely preserve it.
  AffineMain<T> out(x);
  out.Apow(y.itv(), Interval(x));
  return out;
}

template<class T>
inline AffineMain<T> root(const AffineMain<T>& x, int n) {
	AffineMain<T> out(x);
	out.Aroot(n,x.itv());
	return out;
}

template<class T>
inline AffineMain<T> cos(const AffineMain<T>& x){
	AffineMain<T> out(x);
	out.Acos(x.itv());
	return out;
}

template<class T>
inline AffineMain<T> sin(const AffineMain<T>& x){
	AffineMain<T> out(x);
	out.Asin(x.itv());
	return out;
}

template<class T>
inline AffineMain<T> tan(const AffineMain<T>& x){
	AffineMain<T> out(x);
	out.Atan(x.itv());
	return out;
}

template<class T>
inline AffineMain<T> acos(const AffineMain<T>& x){
	AffineMain<T> out(x);
	out.Aacos(x.itv());
	return out;
}

template<class T>
inline AffineMain<T> asin(const AffineMain<T>& x){
	AffineMain<T> out(x);
	out.Aasin(x.itv());
	return out;
}

template<class T>
inline AffineMain<T> atan(const AffineMain<T>& x){
	AffineMain<T> out(x);
	out.Aatan(x.itv());
	return out;
}

template<class T>
inline AffineMain<T> cosh(const AffineMain<T>& x){
	AffineMain<T> out(x);
	out.Acosh(x.itv());
	return out;
}

template<class T>
inline AffineMain<T> sinh(const AffineMain<T>& x){
	AffineMain<T> out(x);
	out.Asinh(x.itv());
	return out;
}

template<class T>
inline AffineMain<T> tanh(const AffineMain<T>& x){
	AffineMain<T> out(x);
	out.Atanh(x.itv());
	return out;
}

template<class T>
inline AffineMain<T> acosh(const AffineMain<T>& x){
	AffineMain<T> out(x);
	out.Aacosh(x.itv());
	return out;
}

template<class T>
inline AffineMain<T> asinh(const AffineMain<T>& x){
	AffineMain<T> out(x);
	out.Aasinh(x.itv());
	return out;
}

template<class T>
inline AffineMain<T> atanh(const AffineMain<T>& x){
	AffineMain<T> out(x);
	out.Aatanh(x.itv());
	return out;
}

template<class T>
inline AffineMain<T> atan2(const AffineMain<T>& y, const AffineMain<T>& x)
{
    AffineMain<T> res(y);
    res.Aatan2(x, y.itv(), x.itv());
    return res;
}


template<class T>
inline AffineMain<T> abs(const AffineMain<T> &x){
	AffineMain<T> out(x);
	out.Aabs(x.itv());
	return out;
}


template<class T>
inline Interval max(const AffineMain<T>& x, const AffineMain<T>& y) {
	return max(x.itv(), y.itv());
}
template<class T>
inline Interval max(const Interval& x, const AffineMain<T>& y) {
	return max(x, y.itv());
}
template<class T>
inline Interval max(const AffineMain<T>& x, const Interval& y) {
	return max(x.itv(), y);
}


template<class T>
inline Interval min(const AffineMain<T>& x, const AffineMain<T>& y) {
	return min(x.itv(), y.itv());
}
template<class T>
inline Interval min(const Interval& x, const AffineMain<T>& y) {
	return min(x, y.itv());
}
template<class T>
inline Interval min(const AffineMain<T>& x, const Interval& y) {
	return min(x.itv(), y);
}


template<class T>
inline Interval integer(const AffineMain<T>& x)
{
  return integer(x.itv());
}

template<class T>
inline Interval floor(const AffineMain<T>& x)
{
  return floor(x.itv());
}

template<class T>
inline Interval ceil(const AffineMain<T>& x)
{
  return ceil(x.itv());
}

template<class T>
inline AffineMain<T> sign(const AffineMain<T>& x) {
	AffineMain<T> out;
	out =sign(x.itv());
	return out;
}


template<class T>
inline AffineMain<T> chi(const AffineMain<T>&  a,const AffineMain<T>&  b,const AffineMain<T>&  c) {
	return chi(a.itv(),b,c);
}
template<class T>
inline AffineMain<T> chi(const Interval&  a,const Interval&  b,const AffineMain<T>&  c) {
	AffineMain<T> out(c);
	out = b;
	return chi(a,out,c);
}
template<class T>
inline AffineMain<T> chi(const Interval&  a,const AffineMain<T>&  b,const Interval&  c) {
	AffineMain<T> out(b);
	out = c;
	return chi(a,b,out);
}
template<class T>
inline AffineMain<T> chi(const AffineMain<T>&  a,const Interval&  b,const AffineMain<T>&  c) {
	AffineMain<T> out(c);
	out = b;
	return chi(a.itv(),out,c);
}
template<class T>
inline AffineMain<T> chi(const AffineMain<T>&  a,const AffineMain<T>&  b,const Interval&  c){
	AffineMain<T> out(b);
	out = c;
	return chi(a.itv(),b,out);
}
template<class T>
inline AffineMain<T> chi(const Interval& a,
                         const AffineMain<T>& b,
                         const AffineMain<T>& c)
{
  if (a.is_empty()) {
    AffineMain<T> out;
    out = Interval::empty();
    return out;
  }
  if (a.ub() <= 0.0) {
    return AffineMain<T>(b);
  }
  if (a.lb() > 0.0) {
    return AffineMain<T>(c);
  }

  // Both branches are possible. Returning the interval hull deliberately
  // drops either branch's dependency rather than selecting an arbitrary one.
		Interval b_itv = b.itv();
		Interval c_itv = c.itv();
		if (b_itv.is_subset(c_itv)) {
			return AffineMain<T>(c);
		} else if (c_itv.is_subset(b_itv)) {
			return AffineMain<T>(b);
		} else {
			AffineMain<T> out;
			out = b_itv|c_itv;
			return  out;
		}
}

template<class T>
inline std::ostream& operator<<(std::ostream& os, const AffineMain<T>& x) {
	{
		os << x.itv() << " : ";
		if (x.is_active()) {
			os << x.mid();
			for (int i = 0; i < x.noise_count(); i++) {
				os << " + " << x.noise(i) << " eps_" << i;
			}
			os << " + " << x.err() << " [-1,1] ";
		} else {
			os << "Affine Form is not enabled. ";
		}
		return os;

	}
}


} // end namespace codac2
