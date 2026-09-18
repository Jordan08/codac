/** 
 *  \file codac2_math.h
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include <limits>
#include <numbers>
#include <cmath>

namespace codac2
{
  const double oo = std::numeric_limits<double>::infinity();

  constexpr double PI = std::numbers::pi; // Need C++20

  template<typename T>
  requires std::is_arithmetic_v<T>
  inline constexpr int sign(T x)
  {
    return (x > T(0)) - (x < T(0));
  }

  inline int integer(double x)
  {
    return (int)x;
  }
  
  inline double pow(double v, double p)
  {
    return v < 0 ? -std::pow(-v,p) : std::pow(v,p);
  }
  
  inline double root(double v, double p)
  {
    return pow(v, 1./p);
  }

  template<typename T>
  inline bool is_nan(const T& x)
  {
    if constexpr(std::is_arithmetic_v<T>)
      return std::isnan(x);
    else
      return x.is_nan();
  }
}