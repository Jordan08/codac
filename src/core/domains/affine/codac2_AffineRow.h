/**
 *  \file codac2_AffineRow.h
 *
 *  Dynamic-size row vector of affine forms.
 *
 * ----------------------------------------------------------------------------
 *  \date       2026
 *  \author     Jordan Ninin
 *  \copyright  Copyright 2026 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include "codac2_matrices.h"
#include "codac2_Vector.h"
#include "codac2_IntervalVector.h"
#include "codac2_AffineMain.h"

namespace codac2
{
  /**
   * \brief Alias for a dynamic-size row vector of affine forms.
   *
   * Represents a row vector with a dynamic number of columns,
   * where each element is an \c AffineMain<T> object.
   *
   * This type alias is based on Eigen's matrix template and corresponds to ``Eigen::Matrix<AffineMain<T>,1,-1>``.
   */
  template<class T>
  using AffineMainRow = Eigen::Matrix<AffineMain<T>,1,-1>;

  /**
   * \brief Streams out x
   *
   * \param os the stream to be updated
   * \param x the affine row vector to stream out
   * \return a reference to the updated stream
   */
  template<class T>
  inline std::ostream& operator<<(std::ostream& os, const AffineMainRow<T>& x)
  {
    if(x.is_empty())
      return os << "[ empty row ]";

    else
    {
      os << x.format(codac_row_fmt());
      return os;
    }
  }
}
