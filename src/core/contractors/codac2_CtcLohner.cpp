/** 
 *  codac2_CtcLohner.cpp
 * ----------------------------------------------------------------------------
 *  \date       2026
 *  \author     Auguste Bourgois, Maël Godard
 *  \copyright  Copyright 2025 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include "codac2_CtcLohner.h"
#include "codac2_Slice.h"

// TO DELETE
#include <iostream>

using namespace std;
using namespace codac2;

LohnerAlgorithm::LohnerAlgorithm(const AnalyticFunction<VectorType> *f,
                                 double h,
                                 bool forward,
                                 const IntervalVector &u0,
                                 int contractions,
                                 double eps)
    : dim(f->input_size()),
      h(h),
      direction((forward) ? FWD : BWD),
      eps(eps),
      contractions(contractions),
      u(u0),
      z(u0 - u0.mid()),
      r(z),
      B(Matrix::eye(dim,dim)),
      Binv(Matrix::eye(dim,dim)),
      u_hat(u0.mid()),
      f(f) {}

const IntervalVector &LohnerAlgorithm::integrate(unsigned int steps, double H) {
  if (H > 0) h = H;
  for (unsigned int i = 0; i < steps; ++i) {
    IntervalVector z1 = z, r1 = r, u1 = u;
    Matrix B1 = B, B1inv = Binv;
    Vector u_hat1 = u_hat;
    IntervalVector u_t = globalEnclosure(u, FWD);
    for (int j = 0; j < contractions; ++j) {
      z1 = 0.5 * h * h * f->diff(u_t) * f->eval(u_t);
      Vector m1 = z1.mid();
      IntervalMatrix A = Matrix::eye(dim,dim) + h * direction * f->diff(u);
      Eigen::HouseholderQR<Eigen::MatrixXd> qr((A * B.template cast<Interval>()).mid());
      B1 = qr.householderQ();
      B1inv = B1.inverse();
      r1 = (B1inv.template cast<Interval>() * A * B.template cast<Interval>()) * r + B1inv.template cast<Interval>() * (z1 - m1);
      IntervalVector inter = f->eval(u_hat);
      u_hat1 = u_hat + h * direction * inter.mid() + m1;
      u1 = u_hat1 + B1 * r1;
      if (j < contractions - 1) {
        u_t = u_t & globalEnclosure(u1, BWD);
      }
    }
    z = z1, r = r1, u = u1, B = B1, Binv = B1inv, u_hat = u_hat1, u_tilde = u_t;
  }
  return u;
}

IntervalVector LohnerAlgorithm::globalEnclosure(const IntervalVector &initialGuess, double dir) {
  IntervalVector u_0 = initialGuess;
  for (unsigned int i = 0; i < 30; ++i) {
    IntervalVector u_1 = initialGuess + dir * direction * Interval(0, h) * f->eval(u_0);
    if (u_0.is_superset(u_1)) {
      return u_0;
    } else {
      u_0 = (1 + eps) * u_1 - eps * u_1;
    }
  }
  throw GlobalEnclosureError();
}

void LohnerAlgorithm::contractStep(const IntervalVector &x) {
  u = x & u;
  u_hat = u.mid();
  r = r & (Binv * (u - u_hat));
}

const IntervalVector &LohnerAlgorithm::getLocalEnclosure() const {
  return u;
}

const IntervalVector &LohnerAlgorithm::getGlobalEnclosure() const {
  return u_tilde;
}

CtcLohner::CtcLohner(const AnalyticFunction<VectorType>& f, int contractions, double eps)
  :   m_f(f),
      contractions(contractions),
      dim(f.input_size()),
      eps(eps) {}

void CtcLohner::contract(codac2::SlicedTube<IntervalVector> &tube, TimePropag t_propa) const
{
  assert((!tube.is_empty()) && (tube.size() == dim));
  IntervalVector input_gate = IntervalVector::constant(dim, Interval(0));
  IntervalVector output_gate = IntervalVector::constant(dim, Interval(0));
  IntervalVector slice = IntervalVector::constant(dim, Interval(0));
  double h;
  if ((t_propa & TimePropag::FWD) != static_cast<TimePropag>(0)) 
  {
    input_gate = tube.first_slice()->input_gate();
    LohnerAlgorithm lo(&m_f, 0.1, true, input_gate, contractions, eps);
    // Forward loop
    for (auto it = tube[0].tdomain()->begin(); it != tube[0].tdomain()->end(); ++it) 
    {
      if (!tube.slice(it)->is_gate())
      {
        h = (*it).diam();
        output_gate = tube.slice(it)->output_gate();
        slice = tube.slice(it)->codomain();
        lo.integrate(1, h);
        lo.contractStep(output_gate);
        tube.slice(it)->set(slice & lo.getGlobalEnclosure());
      }
    }
    tube.set(output_gate & lo.getLocalEnclosure(), tube.tdomain()->t0_tf().ub());
  }
  if ((t_propa & TimePropag::BWD) != static_cast<TimePropag>(0)) 
  {
    input_gate = tube.last_slice()->output_gate();
    LohnerAlgorithm lo2(&m_f, 0.1, false, input_gate, contractions, eps);
    // Backward loop
    for (auto it = tube[0].tdomain()->rbegin(); it != tube[0].tdomain()->rend(); ++it) 
    {
      if (!tube.slice(it)->is_gate())
      {
        h = (*it).diam();
        output_gate = tube.slice(it)->input_gate();
        slice = tube.slice(it)->codomain();
        lo2.integrate(1, h);
        lo2.contractStep(output_gate);
        tube.slice(it)->set(slice & lo2.getGlobalEnclosure());
      }
    }
    tube.set(output_gate & lo2.getLocalEnclosure(), tube.tdomain()->t0_tf().lb());
  }
}