/** 
 *  \file codac2_SlicedTube_operations.h
 * ----------------------------------------------------------------------------
 *  \date       2026
 *  \author     Simon Rohou
 *  \copyright  Copyright 2026 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include "codac2_SlicedTube.h"
#include "codac2_math.h"

namespace codac2
{
  #define macro_unary_tube(f) \
  { \
    auto y = x1; \
    for(auto it = y.tdomain()->begin() ; it != y.tdomain()->end() ; it++) \
    { \
      auto sy = y.slice(it); \
      sy->codomain() = f(sy->codomain()); \
    } \
    return y; \
  }; \

  #define macro_binary_tube_tube(f) \
  { \
    assert_release(x1.tdomain() == x2.tdomain()); \
    auto y = x1; \
    for(auto it = y.tdomain()->begin() ; it != y.tdomain()->end() ; it++) \
    { \
      auto sy = y.slice(it); \
      sy->codomain() = f(sy->codomain(),x2.slice(it)->codomain()); \
    } \
    return y; \
  } \

  #define macro_binary_real_tube(f) \
  { \
    auto y = x2; \
    for(auto it = y.tdomain()->begin() ; it != y.tdomain()->end() ; it++) \
    { \
      auto sy = y.slice(it); \
      sy->codomain() = f(x1,sy->codomain()); \
    } \
    return y; \
  } \

  #define macro_binary_tube_real(f) \
  { \
    auto y = x1; \
    for(auto it = y.tdomain()->begin() ; it != y.tdomain()->end() ; it++) \
    { \
      auto sy = y.slice(it); \
      sy->codomain() = f(sy->codomain(),x2); \
    } \
    return y; \
  } \

  #define macro_member_binary_tube_tube(f) \
  { \
    assert_release(x1.tdomain() == x2.tdomain()); \
    for(auto it = x1.tdomain()->begin() ; it != x1.tdomain()->end() ; it++) \
    { \
      auto sx1 = x1.slice(it); \
      sx1->codomain() = f(sx1->codomain(),x2.slice(it)->codomain()); \
    } \
    return x1; \
  } \

  #define macro_member_binary_tube_real(f) \
  { \
    for(auto it = x1.tdomain()->begin() ; it != x1.tdomain()->end() ; it++) \
    { \
      auto sx1 = x1.slice(it); \
      sx1->codomain() = f(sx1->codomain(),x2); \
    } \
    return x1; \
  } \

  template<typename T,typename X1,typename X2>
  inline T operator_tube_add(const X1& x1, const X2& x2) { return x1 + x2; }

  template<typename T,typename X1,typename X2>
  inline T operator_tube_sub(const X1& x1, const X2& x2) { return x1 - x2; }

  template<typename T,typename X1,typename X2>
  inline T operator_tube_mul(const X1& x1, const X2& x2) { return x1 * x2; }

  template<typename T,typename X2>
    requires (!std::is_same_v<X2,Interval>)
  inline T operator_tube_mul_scal(const Interval& x1, const X2& x2) { return x1 * x2; }

  template<typename T,typename X1>
  inline T operator_tube_mul_scal(const X1& x1, const Interval& x2) { return x1 * x2; }

  inline IntervalMatrix operator_tube_mul_vec(const IntervalMatrix& x1, const IntervalVector& x2) { return x1 * x2; }

  template<typename T,typename X1,typename X2>
  inline T operator_tube_div(const X1& x1, const X2& x2) { return x1 / x2; }

  template<typename T,typename X1>
  inline T operator_tube_div_scal(const X1& x1, const Interval& x2) { return x1 / x2; }

  /** \brief \f$x1(\cdot)\f$
    * \param x1
    * \return tube output
    */
  template<typename T>
  inline const SlicedTube<T>& operator+(const SlicedTube<T>& x1) {
    return x1;
  }

  /** \brief \f$x_1(\cdot)+x_2(\cdot)\f$
    * \param x1
    * \param x2
    * \return tube output
    */
  template<typename T>
  inline SlicedTube<T> operator+(const SlicedTube<T>& x1, const SlicedTube<T>& x2)
    macro_binary_tube_tube(operator_tube_add<T>);

  /** \brief \f$x_1(\cdot)+x_2\f$
    * \param x1
    * \param x2
    * \return tube output
    */
  template<typename T,typename Q>
  inline SlicedTube<T> operator+(const SlicedTube<T>& x1, const Q& x2)
    macro_binary_tube_real(operator_tube_add<T>);

  /** \brief \f$x+x_2(\cdot)\f$
    * \param x1
    * \param x2
    * \return tube output
    */
  template<typename T,typename Q>
  inline SlicedTube<T> operator+(const Q& x1, const SlicedTube<T>& x2)
    macro_binary_real_tube(operator_tube_add<T>);

  /**
   * \brief Operates +=
   * \param x1
   * \param x2
   * \return updated output
   */
  template<typename T,typename Q>
  inline SlicedTube<T>& operator+=(SlicedTube<T>& x1, const Q& x2)
    macro_member_binary_tube_real(operator_tube_add<T>);

  /**
   * \brief Operates +=
   * \param x1
   * \param x2
   * \return updated output
   */
  template<typename T>
  inline SlicedTube<T>& operator+=(SlicedTube<T>& x1, const SlicedTube<T>& x2)
    macro_member_binary_tube_tube(operator_tube_add<T>);

  /** \brief \f$-x_1(\cdot)\f$
    * \param x1
    * \return tube output
    */
  template<typename T>
  inline SlicedTube<T> operator-(const SlicedTube<T>& x1) {
    return -1.*x1;
  }

  /** \brief \f$x_1(\cdot)-x_2(\cdot)\f$
    * \param x1
    * \param x2
    * \return tube output
    */
  template<typename T>
  inline SlicedTube<T> operator-(const SlicedTube<T>& x1, const SlicedTube<T>& x2)
    macro_binary_tube_tube(operator_tube_sub<T>);

  /** \brief \f$x_1(\cdot)-x_2\f$
    * \param x1
    * \param x2
    * \return tube output
    */
  template<typename T,typename Q>
  inline SlicedTube<T> operator-(const SlicedTube<T>& x1, const Q& x2)
    macro_binary_tube_real(operator_tube_sub<T>);

  /** \brief \f$x-x_2(\cdot)\f$
    * \param x1
    * \param x2
    * \return tube output
    */
  template<typename T,typename Q>
  inline SlicedTube<T> operator-(const Q& x1, const SlicedTube<T>& x2)
    macro_binary_real_tube(operator_tube_sub<T>);

  /**
   * \brief Operates -=
   * \param x1
   * \param x2
   * \return updated output
   */
  template<typename T,typename Q>
  inline SlicedTube<T>& operator-=(SlicedTube<T>& x1, const Q& x2)
    macro_member_binary_tube_real(operator_tube_sub<T>);

  /**
   * \brief Operates -=
   * \param x1
   * \param x2
   * \return updated output
   */
  template<typename T>
  inline SlicedTube<T>& operator-=(SlicedTube<T>& x1, const SlicedTube<T>& x2)
    macro_member_binary_tube_tube(operator_tube_sub<T>);

  /** \brief \f$x_1\cdot x_2(\cdot)\f$
    * \param x1
    * \param x2
    * \return tube output
    */
  template<typename T>
    requires (!std::is_same_v<T,double>)
  inline SlicedTube<T> operator*(const Interval& x1, const SlicedTube<T>& x2)
    macro_binary_real_tube(operator_tube_mul_scal<T>);

  /** \brief \f$x_1(\cdot)\cdot x_2\f$
    * \param x1
    * \param x2
    * \return tube output
    */
  template<typename T>
    requires (!std::is_same_v<T,double>)
  inline SlicedTube<T> operator*(const SlicedTube<T>& x1, const Interval& x2)
    macro_binary_tube_real(operator_tube_mul_scal<T>);

  /** \brief \f$x_1(\cdot)\cdot x_2(\cdot)\f$
    * \param x1
    * \param x2
    * \return tube output
    */
  template<typename T>
  inline SlicedTube<T> operator*(const SlicedTube<T>& x1, const SlicedTube<T>& x2)
    macro_binary_tube_tube(operator_tube_mul<T>);

  /** \brief \f$x_1(\cdot)\cdot x_2\f$
    * \param x1
    * \param x2
    * \return tube output
    */
  template<typename T,typename Q>
  inline SlicedTube<T> operator*(const SlicedTube<T>& x1, const Q& x2)
    macro_binary_tube_real(operator_tube_mul<T>);

  /** \brief \f$x\cdot x_2(\cdot)\f$
    * \param x1
    * \param x2
    * \return tube output
    */
  template<typename T,typename Q>
  inline SlicedTube<T> operator*(const Q& x1, const SlicedTube<T>& x2)
    macro_binary_real_tube(operator_tube_mul<T>);

  /** \brief \f$x_1(\cdot)\cdot x_2\f$
    * \param x1
    * \param x2
    * \return tube output
    */
  inline SlicedTube<IntervalVector> operator*(const SlicedTube<IntervalMatrix>& x1, const SlicedTube<IntervalVector>& x2)
    macro_binary_tube_tube(operator_tube_mul_vec);

  /**
   * \brief Operates *=
   * \param x1
   * \param x2
   * \return updated output
   */
  template<typename T,typename Q>
  inline SlicedTube<T>& operator*=(SlicedTube<T>& x1, const Q& x2)
    macro_member_binary_tube_real(operator_tube_mul<T>);

  /**
   * \brief Operates *=
   * \param x1
   * \param x2
   * \return updated output
   */
  template<typename T>
  inline SlicedTube<T>& operator*=(SlicedTube<T>& x1, const SlicedTube<T>& x2)
    macro_member_binary_tube_tube(operator_tube_mul<T>);

  /** \brief \f$x_2(\cdot)/x_1\f$
    * \param x1
    * \param x2
    * \return tube output
    */
  template<typename T>
    requires (!std::is_same_v<T,double>)
  inline SlicedTube<T> operator/(const SlicedTube<T>& x1, const Interval& x2)
    macro_binary_tube_real(operator_tube_div_scal<T>);

  /** \brief \f$x_1(\cdot)/x_2(\cdot)\f$
    * \param x1
    * \param x2
    * \return tube output
    */
  template<typename T>
  inline SlicedTube<T> operator/(const SlicedTube<T>& x1, const SlicedTube<T>& x2)
    macro_binary_tube_tube(operator_tube_div<T>);

  /** \brief \f$x_1(\cdot)/x_2\f$
    * \param x1
    * \param x2
    * \return tube output
    */
  template<typename T,typename Q>
  inline SlicedTube<T> operator/(const SlicedTube<T>& x1, const Q& x2)
    macro_binary_tube_real(operator_tube_div<T>);

  /** \brief \f$x/x_2(\cdot)\f$
    * \param x1
    * \param x2
    * \return tube output
    */
  template<typename T,typename Q>
  inline SlicedTube<T> operator/(const Q& x1, const SlicedTube<T>& x2)
    macro_binary_real_tube(operator_tube_div<T>);

  /**
   * \brief Operates /=
   * \param x1
   * \param x2
   * \return updated output
   */
  template<typename T,typename Q>
  inline SlicedTube<T>& operator/=(SlicedTube<T>& x1, const Q& x2)
    macro_member_binary_tube_real(operator_tube_div<T>);

  /**
   * \brief Operates /=
   * \param x1
   * \param x2
   * \return updated output
   */
  template<typename T>
  inline SlicedTube<T>& operator/=(SlicedTube<T>& x1, const SlicedTube<T>& x2)
    macro_member_binary_tube_tube(operator_tube_div<T>);

  /** \brief \f$x^2(\cdot)\f$
    * \param x1
    * \return tube output
    */
  SlicedTube<Interval> sqr(const SlicedTube<Interval>& x1);

  /** \brief \f$\sqrt{x_1(\cdot)}\f$
    * \param x1
    * \return tube output
    */
  SlicedTube<Interval> sqrt(const SlicedTube<Interval>& x1);

  /** \brief \f$x^x_2(\cdot)\f$
    * \param x1
    * \param x2
    * \return tube output
    */
  SlicedTube<Interval> pow(const SlicedTube<Interval>& x1, const Interval& x2);

  /** \brief \f$\exp(x_1(\cdot))\f$
    * \param x1
    * \return tube output
    */
  SlicedTube<Interval> exp(const SlicedTube<Interval>& x1);

  /** \brief \f$\log(x_1(\cdot))\f$
    * \param x1
    * \return tube output
    */
  SlicedTube<Interval> log(const SlicedTube<Interval>& x1);

  /** \brief \f$\cos(x_1(\cdot))\f$
    * \param x1
    * \return tube output
    */
  SlicedTube<Interval> cos(const SlicedTube<Interval>& x1);

  /** \brief \f$\sin(x_1(\cdot))\f$
    * \param x1
    * \return tube output
    */
  SlicedTube<Interval> sin(const SlicedTube<Interval>& x1);

  /** \brief \f$\tan(x_1(\cdot))\f$
    * \param x1
    * \return tube output
    */
  SlicedTube<Interval> tan(const SlicedTube<Interval>& x1);

  /** \brief \f$\arccos(x_1(\cdot))\f$
    * \param x1
    * \return tube output
    */
  SlicedTube<Interval> acos(const SlicedTube<Interval>& x1);

  /** \brief \f$\arcsin(x_1(\cdot))\f$
    * \param x1
    * \return tube output
    */
  SlicedTube<Interval> asin(const SlicedTube<Interval>& x1);

  /** \brief \f$\arctan(x_1(\cdot))\f$
    * \param x1
    * \return tube output
    */
  SlicedTube<Interval> atan(const SlicedTube<Interval>& x1);

  /** \brief \f$\mathrm{arctan2}(x_1(\cdot),x_2(\cdot))\f$
    * \param x1
    * \param x2
    * \return tube output
    */
  SlicedTube<Interval> atan2(const SlicedTube<Interval>& x1, const SlicedTube<Interval>& x2);

  /** \brief \f$\mathrm{arctan2}(x_1(\cdot),x_2)\f$
    * \param x1
    * \param x2
    * \return tube output
    */
  SlicedTube<Interval> atan2(const SlicedTube<Interval>& x1, const Interval& x2);

  /** \brief \f$\mathrm{arctan2}(x_1, x_2(\cdot))\f$
    * \param x1
    * \param x2
    * \return tube output
    */
  SlicedTube<Interval> atan2(const Interval& x1, const SlicedTube<Interval>& x2);

  /** \brief \f$\cosh(x_1(\cdot))\f$
    * \param x1
    * \return tube output
    */
  SlicedTube<Interval> cosh(const SlicedTube<Interval>& x1);

  /** \brief \f$\sinh(x_1(\cdot))\f$
    * \param x1
    * \return tube output
    */
  SlicedTube<Interval> sinh(const SlicedTube<Interval>& x1);

  /** \brief \f$\tanh(x_1(\cdot))\f$
    * \param x1
    * \return tube output
    */
  SlicedTube<Interval> tanh(const SlicedTube<Interval>& x1);

  /** \brief \f$\mathrm{arccosh}(x_1(\cdot))\f$
    * \param x1
    * \return tube output
    */
  SlicedTube<Interval> acosh(const SlicedTube<Interval>& x1);

  /** \brief \f$\mathrm{arcsinh}(x_1(\cdot))\f$
    * \param x1
    * \return tube output
    */
  SlicedTube<Interval> asinh(const SlicedTube<Interval>& x1);

  /** \brief \f$\mathrm{arctanh}(x_1(\cdot))\f$
    * \param x1
    * \return tube output
    */
  SlicedTube<Interval> atanh(const SlicedTube<Interval>& x1);

  /** \brief \f$\mid x_1(\cdot)\mid\f$
    * \param x1
    * \return tube output
    */
  SlicedTube<Interval> abs(const SlicedTube<Interval>& x1);

  /** \brief \f$\min(x_1(\cdot),x_2(\cdot))\f$
    * \param x1
    * \param x2
    * \return tube output
    */
  SlicedTube<Interval> min(const SlicedTube<Interval>& x1, const SlicedTube<Interval>& x2);

  /** \brief \f$\min(x_1(\cdot),x_2)\f$
    * \param x1
    * \param x2
    * \return tube output
    */
  SlicedTube<Interval> min(const SlicedTube<Interval>& x1, const Interval& x2);

  /** \brief \f$\min(x_1, x_2(\cdot))\f$
    * \param x1
    * \param x2
    * \return tube output
    */
  SlicedTube<Interval> min(const Interval& x1, const SlicedTube<Interval>& x2);

  /** \brief \f$\max(x_1(\cdot),x_2(\cdot))\f$
    * \param x1
    * \param x2
    * \return tube output
    */
  SlicedTube<Interval> max(const SlicedTube<Interval>& x1, const SlicedTube<Interval>& x2);

  /** \brief \f$\max(x_1(\cdot),x_2)\f$
    * \param x1
    * \param x2
    * \return tube output
    */
  SlicedTube<Interval> max(const SlicedTube<Interval>& x1, const Interval& x2);

  /** \brief \f$\max(x_1, x_2(\cdot))\f$
    * \param x1
    * \param x2
    * \return tube output
    */
  SlicedTube<Interval> max(const Interval& x1, const SlicedTube<Interval>& x2);

  SlicedTube<Interval> sign(const SlicedTube<Interval>& x1);

  SlicedTube<Interval> integer(const SlicedTube<Interval>& x1);

  SlicedTube<Interval> floor(const SlicedTube<Interval>& x1);

  SlicedTube<Interval> ceil(const SlicedTube<Interval>& x1);
}