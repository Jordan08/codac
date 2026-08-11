/** 
 *  \file codac2_MatrixBase_addons_AffineMainVector.h
 * 
 *  This file is included in the declaration of Eigen::MatrixBase,
 *  thanks to the preprocessor token EIGEN_MATRIXBASE_PLUGIN.
 *  See: https://eigen.tuxfamily.org/dox/TopicCustomizing_Plugins.html
 *  and the file codac2_matrices.h
 * 
 * ----------------------------------------------------------------------------
 *  \date       2026
 *  \author     Jordan Ninin
 *  \copyright  Copyright 2026 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */


template<typename U = Scalar, typename OtherDerived>
  requires IsIntervalDomain<U>
	&& (IsIntervalDomain<OtherDerived> ||std::is_arithmetic_v<OtherDerived>)
	&& requires(U& x, const OtherDerived& d) { x *= d; }
inline Derived& operator*=(const OtherDerived& d)
{
  for(Index i = 0; i < this->size(); ++i)
    this->derived().coeffRef(i) *= d;

  return this->derived();
}

template<typename U = Scalar, typename OtherDerived>
  requires IsIntervalDomain<U>
	&& (IsIntervalDomain<OtherDerived> ||std::is_arithmetic_v<OtherDerived>)
	&& requires(U& x, const OtherDerived& d) { x /= d; }
inline Derived& operator/=(OtherDerived d)
{
  for(Index i = 0; i < this->size(); ++i)
    this->derived().coeffRef(i) /= d;

  return this->derived();
}

template<typename U = Scalar, typename OtherDerived>
  requires  IsIntervalDomain<U>
	&& (IsIntervalDomain<OtherDerived> ||std::is_arithmetic_v<OtherDerived>)
	&& requires(U& x, const OtherDerived& d) { x += d; }
inline Derived& operator+=(const OtherDerived& d)
{
  for(Index i = 0; i < this->size(); ++i)
    this->derived().coeffRef(i) += d;

  return this->derived();
}

template<typename U = Scalar, typename OtherDerived>
  requires  IsIntervalDomain<U>
	&& (IsIntervalDomain<OtherDerived> ||std::is_arithmetic_v<OtherDerived>)
	&& requires(U& x, const OtherDerived& d) { x -= d; }
inline Derived& operator-=(const OtherDerived& d)
{
  for(Index i = 0; i < this->size(); ++i)
    this->derived().coeffRef(i) -= d;

  return this->derived();
}


/**
 * \brief Returns the interval hull of a vector whose scalar type provides
 *        an ``itv()`` member.
 *
 * This method is intentionally implemented on Eigen::MatrixBase so that
 * domain vector aliases such as ``AffineMainVector`` do not need to derive
 * from Eigen::Matrix just to expose this operation.
 */
template<typename U=Scalar,int R=RowsAtCompileTime,int C=ColsAtCompileTime>
  requires IsAffineDomain<U> && requires(const U& x) { x.itv(); }
inline auto itv() const
{
  Matrix<codac2::Interval,R,C> m(this->rows(),this->cols());

  for(Index i = 0 ; i < this->rows() ; i++)
    for(Index j = 0 ; j < this->cols() ; j++)
      m(i,j) = (*this)(i,j).itv();

  return m;
}



