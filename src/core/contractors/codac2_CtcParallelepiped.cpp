/** 
 *  CtcParallelepiped.cpp
 * ----------------------------------------------------------------------------
 *  \date       2026
 *  \author     Maël Godard
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include "codac2_CtcParallelepiped.h"
#include "codac2_inversion.h"

// TO DELETE

#include <iostream>

using namespace std;

namespace codac2
{
  void CtcParallelepiped::contract(IntervalVector& x) const
  {
    assert_release(x.size()==_p.c.size());
    IntervalVector x_p = inverse_enclosure(_p.A)*(x - _p.c);  // projection in the parallelepiped's coordinate system
    x &= _p.c + _p.A*(x_p & IntervalVector::constant(x.size(),{-1,1}));
  }
}