/**
 *  \file codac2_AffineVector.h
 *
 *  Dynamic-size vector of affine forms.
 *
 * ----------------------------------------------------------------------------
 *  \date       2026
 *  \author     Jordan Ninin
 *  \copyright  Copyright 2026 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include <ostream>
#include "codac2_matrices.h"
#include "codac2_Vector.h"
#include "codac2_IntervalVector.h"
#include "codac2_AffineMain.h"

namespace codac2 {

template<class T> class AffineVarMainVector;

/**
 *
 * \brief Dynamic-size column vector of affine forms.
 *
 * This is deliberately an alias of Eigen::Matrix rather than a derived
 * wrapper class. Consequently all Eigen constructors, assignments and
 * expression-template operations remain available without having to
 * duplicate Eigen's API in a derived type.
 */
template<class T>
using AffineMainVector = Eigen::Matrix<AffineMain<T>, -1, 1>;

/**
 * \brief Streams out x
 *
 * \param os the stream to be updated
 * \param x the affine vector to stream out
 * \return a reference to the updated stream
 */
template<class T>
inline std::ostream& operator<<(std::ostream& os, const AffineMainVector<T>& x)
{
  if(x.is_empty())
    return os << "[ empty " << x.size() << "d box ]";

  return os << x.format(codac_vector_fmt());
}

} // namespace codac2
