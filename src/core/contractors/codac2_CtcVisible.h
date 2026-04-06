/** * \file codac2_CtcVisible.h
 * ----------------------------------------------------------------------------
 * \date       2026
 * \author     Quentin Brateau
 * \copyright  Copyright 2026 Codac Team
 * \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include "codac2_Ctc.h"
#include "codac2_Segment.h"

namespace codac2
{
  class CtcVisible : public Ctc<CtcVisible, IntervalVector>
  {
    public:
      /**
       * \brief Constructor for visibility from point 'a' relative to segment 's'.
       */
      CtcVisible(const IntervalVector& a, const Segment& s);

      void contract(IntervalVector& x) const;

    private:
      const IntervalVector _a;
      const Segment _s;
      
      // Pre-calculated constants for the obstacle
      const IntervalVector _e1, _e2;
      const IntervalVector _v_e2e1; // e2 - e1
      const IntervalVector _v_ae1;  // a - e1
      const IntervalVector _v_ae2;  // a - e2
      const IntervalVector _s_box;  // Bounding box of the segment
      
      double _k; // Orientation sign (ksi)

      // Internal helpers for the 4 conditions
      void contract_det(IntervalVector& x, const IntervalVector& p, const IntervalVector& v, double sign) const;
      void contract_aabb(IntervalVector& x) const;
  };

  class CtcNoVisible : public Ctc<CtcNoVisible, IntervalVector>
  {
    public:
      CtcNoVisible(const IntervalVector& a, const Segment& s);

      void contract(IntervalVector& x) const;

    private:
      const IntervalVector _a;
      const Segment _s;
      const IntervalVector _v_e2e1, _v_ae1, _v_ae2;
      double _k;

      void contract_det(IntervalVector& x, const IntervalVector& p, const IntervalVector& v, double sign) const;
      void contract_aabb(IntervalVector& x) const;
  };
}