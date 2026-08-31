/**
 *  \file codac2_AffineMatrix.h
 *
 *  Dynamic-size matrix of affine forms.
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
#include "codac2_AffineMain.h"

namespace codac2
{
  /**
   * \brief Alias for a dynamic-size matrix of affine forms.
   *
   * Represents a matrix with a dynamic number of rows and columns,
   * where each element is an \c AffineMain<T> object.
   *
   * This type alias is based on Eigen's matrix template and corresponds to ``Eigen::Matrix<AffineMain<T>,-1,-1>``.
   */
  template<class T>
  using AffineMainMatrix = Eigen::Matrix<AffineMain<T>,-1,-1>;

  /**
   * \brief Streams out x
   *
   * \param os the stream to be updated
   * \param x the affine matrix to stream out
   * \return a reference to the updated stream
   */
  template<class T>
  inline std::ostream& operator<<(std::ostream& os, const AffineMainMatrix<T>& x)
  {
    if(x.is_empty())
      return os << "[ empty " << x.rows() << "x" << x.cols() << " mat ]";

    else
    {
      os << x.format(codac_matrix_fmt());
      return os;
    }
  }
}
