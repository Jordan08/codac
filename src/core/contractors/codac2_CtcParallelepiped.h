/** 
 *  \file codac2_CtcParallelepiped.h
 * ----------------------------------------------------------------------------
 *  \date       2026
 *  \author     Simon Rohou, Benoit Desrochers
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include "codac2_Parallelepiped.h"
#include "codac2_Ctc.h"

namespace codac2
{
  /**
   * \class CtcParallelepiped
   * \brief Implements the contractor associated with a parallelepiped constraint.
   * 
   * This contractor requires a parallelepiped with a square shape matrix \f$A\f$.
   */
  class CtcParallelepiped : public Ctc<CtcParallelepiped,IntervalVector>
  {
    public:

      /**
       * \brief Creates the contractor associated with a parallelepiped constraint.
       */
      CtcParallelepiped(const Parallelepiped& p)
        : Ctc<CtcParallelepiped,IntervalVector>(p.c.size()), _p(p)
      { 
        assert_release(p.A.rows()==p.A.cols() && "Parallelepiped's matrix A must be square");
      }

      /**
       * \brief Applies \f$\mathcal{C}_{\textrm{parallelepiped}}\big([\mathbf{x}]\big)\f$.
       */
      void contract(IntervalVector& x) const;

    protected:

      const Parallelepiped _p;
  };
}