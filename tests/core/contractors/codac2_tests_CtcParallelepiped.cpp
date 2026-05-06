/** 
 *  Codac tests
 * ----------------------------------------------------------------------------
 *  \date       2026
 *  \author     Maël Godard
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <catch2/catch_test_macros.hpp>
#include <codac2_CtcParallelepiped.h>
#include <codac2_Approx.h>

using namespace std;
using namespace codac2;


TEST_CASE("CtcParallelepiped")
{
  CHECK(true);
  Parallelepiped p1 ({1.5,2.8},Matrix({{0.5,0.4},{0,0.2}}));
  CtcParallelepiped ctc_par1 (p1);

  IntervalVector x0 ({{0,5},{0,5}});
  IntervalVector x1 ({{0,0.5},{0,0.5}});

  ctc_par1.contract(x0);
  ctc_par1.contract(x1);

  CHECK(Approx(x0) == IntervalVector({{0.6,2.4},{2.6,3}}));
  CHECK(x1.is_empty());

  Parallelepiped p2 ({-5,-4,3,8},Matrix({{0.5,0.4,0.1,0},{0,0.2,0.3,0},{0,0,0.1,0.4},{0,0,0,0.2}}));
  CtcParallelepiped ctc_par2 (p2);

  IntervalVector x2 ({{-10,10},{-10,10},{-10,10},{-10,10}});
  IntervalVector x3 ({{-0.5,0.5},{-0.5,0.5},{-0.5,0.5},{-0.5,0.5}});

  ctc_par2.contract(x2);
  ctc_par2.contract(x3);

  CHECK(Approx(x2) == IntervalVector({{-6,-4},{-4.5,-3.5},{2.5,3.5},{7.8,8.2}}));
  CHECK(x3.is_empty());
}