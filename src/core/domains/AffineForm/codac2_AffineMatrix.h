/**
 *  \file codac2_AffineMatrix.h
 *
 *  Dynamic-size matrix wrapper for affine forms.
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


namespace codac2
{
  /**
   * \brief Alias for a dynamic-size matrix of affine forms.
   *
   * This alias corresponds to ``Eigen::Matrix<AffineMain<T>,-1,-1>``.
   */
  template<class T>
  using AffineMainMatrix = Eigen::Matrix< AffineMain<T>,-1,-1>;


  /**
   * \brief Streams out an affine matrix.
   *
   * Prints ``[ empty RxC mat ]`` for empty matrices, otherwise uses
   * ``codac_matrix_fmt()`` for Eigen-style formatting.
   *
   * \param os output stream
   * \param x affine matrix to print
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
