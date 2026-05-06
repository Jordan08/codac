#!/usr/bin/env python

#  Codac tests
# ----------------------------------------------------------------------------
#  \date       2026
#  \author     Maël Godard
#  \copyright  Copyright 2024 Codac Team
#  \license    GNU Lesser General Public License (LGPL)

import unittest
from codac import *

class TestCtcParallelepiped(unittest.TestCase):
  
  def test_CtcParallelepiped(self):

    p1 = Parallelepiped([1.5,2.8],Matrix([[0.5,0.4],[0,0.2]]))
    ctc_par1 = CtcParallelepiped(p1)

    x0 = IntervalVector([[0,5],[0,5]])
    x1 = IntervalVector([[0,0.5],[0,0.5]])

    ctc_par1.contract(x0)
    ctc_par1.contract(x1)

    self.assertTrue(Approx(x0) == IntervalVector([[0.6,2.4],[2.6,3]]))
    self.assertTrue(x1.is_empty())

    p2 = Parallelepiped([-5,-4,3,8],Matrix([[0.5,0.4,0.1,0],[0,0.2,0.3,0],[0,0,0.1,0.4],[0,0,0,0.2]]))
    ctc_par2 = CtcParallelepiped(p2)

    x2 = IntervalVector([[-10,10],[-10,10],[-10,10],[-10,10]])
    x3 = IntervalVector([[-0.5,0.5],[-0.5,0.5],[-0.5,0.5],[-0.5,0.5]])

    ctc_par2.contract(x2)
    ctc_par2.contract(x3)

    self.assertTrue(Approx(x2) == IntervalVector([[-6,-4],[-4.5,-3.5],[2.5,3.5],[7.8,8.2]]))
    self.assertTrue(x3.is_empty())


if __name__ ==  '__main__':
  unittest.main()