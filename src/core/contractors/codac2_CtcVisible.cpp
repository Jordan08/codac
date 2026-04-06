/** * codac2_CtcVisible.cpp
 * ----------------------------------------------------------------------------
 * \date       2026
 * \author     Quentin Brateau
 * \copyright  Copyright 2026 Codac Team
 * \license    GNU Lesser General Public License (LGPL)
 */


 /** * codac2_CtcVisible.cpp
 */

#include "codac2_CtcVisible.h"
#include "codac2_det.h"
#include "codac2_arith_sub.h"
#include "codac2_max.h"
#include "codac2_min.h"
#include "codac2_IntervalVector.h"


using namespace std;
using namespace codac2;

CtcVisible::CtcVisible(const IntervalVector& a, const Segment& s)
  : Ctc<CtcVisible, IntervalVector>(2), 
    _a(a), _s(s),
    _v_e2e1(s[1] - s[0]), 
    _v_ae1(a - s[0]), 
    _v_ae2(a - s[1])
{
    // Compute orientation (ksi)
    double det_val = (_a[0].mid() - _s[0][0].mid()) * (_v_e2e1[1].mid()) - 
                     (_a[1].mid() - _s[0][1].mid()) * (_v_e2e1[0].mid());
    _k = (det_val > 0) ? 1.0 : -1.0;
}

void CtcVisible::contract(IntervalVector& x) const
{
  IntervalVector x1(x), x2(x), x3(x), x4(x);

  contract_det(x1, _s[0], _v_e2e1, _k);
  contract_det(x2, _s[0], _v_ae1, _k);
  contract_det(x3, _s[1], _v_ae2, -_k);
  contract_aabb(x4);

  x &= (x1 | x2 | x3 | x4);
}

void CtcVisible::contract_det(IntervalVector& x, const IntervalVector& p, const IntervalVector& v, double sign) const
{
  IntervalVector v_xp = x - p;
  IntervalVector v_fixed = v; 
  
  Interval target = (sign > 0) ? Interval(0, oo) : Interval(-oo, 0);
  
  DetOp::bwd(target, v_xp, v_fixed);
  x &= v_xp + p;
}

void CtcVisible::contract_aabb(IntervalVector& x) const
{
  auto contract_1dim = [](double a, Interval& x_val, double c, double d) {
    double min_cd = std::min(c, d);
    double max_cd = std::max(c, d);

    // Forward contractions
    Interval i1 = MinOp::fwd(Interval(a), x_val);
    Interval i2 = MaxOp::fwd(i1, Interval(min_cd));
    Interval i3 = MaxOp::fwd(Interval(a), x_val);
    Interval i4 = MinOp::fwd(i3, Interval(max_cd));
    Interval i5 = i2 - i4;

    // Top of the DAG
    if ((i5 &= Interval(0, oo)).is_empty()) {
        x_val.set_empty();
        return;
    }

    // Backward contractions
    i2 &= i5 + i4;
    i4 &= i2 - i5;

    Interval tmp_max_cd = Interval(max_cd);
    Interval tmp_min_cd = Interval(min_cd);
    Interval tmp_a = Interval(a);

    MinOp::bwd(i4, i3, tmp_max_cd);

    tmp_a = Interval(a); // reset
    MaxOp::bwd(i3, tmp_a, x_val);

    MaxOp::bwd(i2, i1, tmp_min_cd);

    tmp_a = Interval(a); // reset
    MinOp::bwd(i1, tmp_a, x_val);
  };

  // Apply the 1D contraction to each dimension
  // Note: _a and _s are assumed to be Point-like (degenerate intervals) 
  // so we use .mid() to get the double values c, d, and a.
  contract_1dim(_a[0].mid(), x[0], _s[0][0].mid(), _s[1][0].mid());
  contract_1dim(_a[1].mid(), x[1], _s[0][1].mid(), _s[1][1].mid());
}

CtcNoVisible::CtcNoVisible(const IntervalVector& a, const Segment& s)
  : Ctc<CtcNoVisible, IntervalVector>(2), 
    _a(a), _s(s),
    _v_e2e1(s[1] - s[0]), 
    _v_ae1(a - s[0]), 
    _v_ae2(a - s[1])
{
    double det_val = (_a[0].mid() - _s[0][0].mid()) * (_v_e2e1[1].mid()) - 
                     (_a[1].mid() - _s[0][1].mid()) * (_v_e2e1[0].mid());
    _k = (det_val > 0) ? 1.0 : -1.0;
}

void CtcNoVisible::contract(IntervalVector& x) const
{
  IntervalVector xi(x);
  contract_det(xi, _s[0], _v_e2e1, -_k);
  contract_det(xi, _s[0], _v_ae1, -_k);
  contract_det(xi, _s[1], _v_ae2, _k);
  contract_aabb(xi);

  x &= xi;
}

void CtcNoVisible::contract_det(IntervalVector& x, const IntervalVector& p, const IntervalVector& v, double sign) const
{
  IntervalVector v_xp = x - p;
  IntervalVector v_fixed = v; 
  Interval target = (sign > 0) ? Interval(0, oo) : Interval(-oo, 0);
  
  DetOp::bwd(target, v_xp, v_fixed);
  x &= v_xp + p;
}

void CtcNoVisible::contract_aabb(IntervalVector& x) const
{
  auto contract_1dim = [](double a, Interval& x_val, double c, double d) {
    double min_cd = std::min(c, d);
    double max_cd = std::max(c, d);

    Interval i1 = MinOp::fwd(Interval(a), x_val);
    Interval i2 = MaxOp::fwd(i1, Interval(min_cd));
    Interval i3 = MaxOp::fwd(Interval(a), x_val);
    Interval i4 = MinOp::fwd(i3, Interval(max_cd));
    Interval i5 = i2 - i4;

    if ((i5 &= Interval(-oo, 0)).is_empty()) {
        x_val.set_empty();
        return;
    }

    i2 &= i5 + i4;
    i4 &= i2 - i5;

    Interval tmp_max_cd = Interval(max_cd);
    Interval tmp_min_cd = Interval(min_cd);
    Interval tmp_a = Interval(a);

    MinOp::bwd(i4, i3, tmp_max_cd);
    tmp_a = Interval(a); 
    MaxOp::bwd(i3, tmp_a, x_val);
    MaxOp::bwd(i2, i1, tmp_min_cd);
    tmp_a = Interval(a);
    MinOp::bwd(i1, tmp_a, x_val);
  };

  contract_1dim(_a[0].mid(), x[0], _s[0][0].mid(), _s[1][0].mid());
  contract_1dim(_a[1].mid(), x[1], _s[0][1].mid(), _s[1][1].mid());
}
