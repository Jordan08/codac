/** * Codac tests - Visibility Contractors
 * ----------------------------------------------------------------------------
 * \date       2026
 * \author     Quentin Brateau
 * \copyright  Copyright 2026 Codac Team
 * \license    GNU Lesser General Public License (LGPL)
 */

#include <catch2/catch_test_macros.hpp>
#include <codac2_CtcVisible.h>
#include <codac2_Vector.h>

using namespace codac2;

TEST_CASE("CtcVisible & CtcNoVisible - Point/Segment visibility")
{
  // Observer at origin, obstacle horizontal from (2, -1) to (2, 1)
  Vector a({0.0, 0.0});
  Segment s({2.0, -1.0}, {2.0, 1.0});
  
  CtcVisible ctc_vis(a, s);
  CtcNoVisible ctc_no_vis(a, s);

  SECTION("Box fully in visible zone (in front of obstacle)")
  {
    IntervalVector x({{0.5, 1.5}, {-0.5, 0.5}});
    IntervalVector x_orig = x;

    ctc_vis.contract(x);
    CHECK(x == x_orig); // Should not be contracted

    ctc_no_vis.contract(x_orig);
    CHECK(x_orig.is_empty()); // Cannot be hidden if in front
  }

  SECTION("Box fully in hidden zone (shadow)")
  {
    IntervalVector x({{3.0, 4.0}, {-0.2, 0.2}});
    IntervalVector x_orig = x;

    ctc_vis.contract(x);
    CHECK(x.is_empty()); // Fully in shadow -> not visible

    ctc_no_vis.contract(x_orig);
    CHECK(x_orig == IntervalVector({{3.0, 4.0}, {-0.2, 0.2}})); // Fully hidden
  }

  SECTION("Box behind the observer")
  {
    IntervalVector x({{-2.0, -1.0}, {-1.0, 1.0}});
    IntervalVector x_orig = x;

    ctc_vis.contract(x);
    CHECK(x == x_orig); // Visible (obstacle is far away in the other direction)

    ctc_no_vis.contract(x_orig);
    CHECK(x_orig.is_empty());
  }

  SECTION("Box on the side (outside the angular cone)")
  {
    IntervalVector x({{1.0, 4.0}, {2.0, 3.0}}); // High above the obstacle
    IntervalVector x_orig = x;

    ctc_vis.contract(x);
    CHECK(x == x_orig); 

    ctc_no_vis.contract(x_orig);
    CHECK(x_orig.is_empty());
  }

  SECTION("Straddling the shadow edge (angular boundary)")
  {
    // Obstacle is y in [-1, 1] at x=2. Upper shadow boundary is line y = 0.5*x (roughly)
    // We place a box at x=4, y in [1.5, 2.5]. 
    // Boundary at x=4 is y=2. So [1.5, 2] is hidden, [2, 2.5] is visible.
    
    IntervalVector x_vis({{4.0, 4.0}, {1.5, 2.5}});
    ctc_vis.contract(x_vis);
    // Note: Depends on precision/ksi, but x_vis[1] should be pruned to [2, 2.5]
    CHECK(x_vis[1].lb() >= 1.99); 

    IntervalVector x_hid({{4.0, 4.0}, {1.5, 2.5}});
    ctc_no_vis.contract(x_hid);
    CHECK(x_hid[1].ub() <= 2.01);
  }

  SECTION("AABB boundary test (Cinterseg logic)")
  {
    // Box is angularly behind the segment, but outside its X-AABB
    // Observer (0,0), Segment x in [2,3], y=0.
    // Box at x=1.5 (between observer and obstacle)
    Segment s2({2.0, 0.0}, {3.0, 0.0});
    Vector a2({0.0, 0.0});
    CtcVisible c_vis2(a2, s2);
    
    IntervalVector x({{1.0, 1.5}, {-0.5, 0.5}});
    IntervalVector x_orig = x;
    c_vis2.contract(x);
    CHECK(x == x_orig); // Visible because it's in front of the obstacle's X-range
  }
}