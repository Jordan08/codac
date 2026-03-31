#!/usr/bin/env python

#  Codac tests - Visibility
# ----------------------------------------------------------------------------
#  \date       2026
#  \author     Quentin Brateau
#  \copyright  Copyright 2026 Codac Team
#  \license    GNU Lesser General Public License (LGPL)

import unittest
from codac import *

class TestVisibility(unittest.TestCase):
  
  def test_CtcVisibility(self):
    # Observer at origin, obstacle vertical at x=2 from y=-1 to y=1
    a = [0.0, 0.0]
    s = Segment([2.0, -1.0], [2.0, 1.0])
    
    ctc_vis = CtcVisible(a, s)
    ctc_nvis = CtcNoVisible(a, s)

    # 1. Box fully in visible zone (in front of obstacle)
    x = IntervalVector([[0.5, 1.5], [-0.5, 0.5]])
    x_orig = IntervalVector(x)
    ctc_vis.contract(x)
    self.assertTrue(x == x_orig)

    x_test_no = IntervalVector(x_orig)
    ctc_nvis.contract(x_test_no)
    self.assertTrue(x_test_no.is_empty())

    # 2. Box fully in hidden zone (shadow)
    x = IntervalVector([[3.0, 4.0], [-0.2, 0.2]])
    x_orig = IntervalVector(x)
    ctc_vis.contract(x)
    self.assertTrue(x.is_empty())

    x_test_no = IntervalVector(x_orig)
    ctc_nvis.contract(x_test_no)
    self.assertTrue(x_test_no == x_orig)

    # 3. Box behind the observer
    x = IntervalVector([[-2.0, -1.0], [-1.0, 1.0]])
    x_orig = IntervalVector(x)
    ctc_vis.contract(x)
    self.assertTrue(x == x_orig)

    x_test_no = IntervalVector(x_orig)
    ctc_nvis.contract(x_test_no)
    self.assertTrue(x_test_no.is_empty())

    # 4. Box on the side (outside the angular cone)
    x = IntervalVector([[1.0, 4.0], [2.0, 3.0]])
    x_orig = IntervalVector(x)
    ctc_vis.contract(x)
    self.assertTrue(x == x_orig)

    x_test_no = IntervalVector(x_orig)
    ctc_nvis.contract(x_test_no)
    self.assertTrue(x_test_no.is_empty())

    # 5. Straddling the shadow edge (angular boundary)
    # Boundary at x=4 is y=2. Visible: [2, 2.5], Hidden: [1.5, 2]
    x_vis = IntervalVector([[4.0, 4.0], [1.5, 2.5]])
    ctc_vis.contract(x_vis)
    self.assertGreaterEqual(x_vis[1].lb(), 1.99)

    x_hid = IntervalVector([[4.0, 4.0], [1.5, 2.5]])
    ctc_nvis.contract(x_hid)
    self.assertLessEqual(x_hid[1].ub(), 2.01)

    # 6. AABB boundary test (Sight-line doesn't reach obstacle)
    s2 = Segment([2.0, 0.0], [3.0, 0.0])
    ctc_vis2 = CtcVisible([0.0, 0.0], s2)
    x = IntervalVector([[1.0, 1.5], [-0.5, 0.5]])
    x_orig = IntervalVector(x)
    ctc_vis2.contract(x)
    self.assertTrue(x == x_orig)

if __name__ == '__main__':
  unittest.main()