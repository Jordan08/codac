/** 
 *  \file codac2_peibos_capd.h
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Maël Godard
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include <vector>
#include "codac2_peibos.h"
#include "codac2_capd.h"
#include "codac2_OctaSym.h"
#include <capd/poincare/lib.h>

namespace codac2
{
  struct PEIBOS_CAPD_Key
  {
    IntervalVector box;
    AnalyticFunction<VectorType> psi_0;
    OctaSym sigma;
    Vector offset;
  };

  using T = std::pair<PEIBOS_CAPD_Key,std::pair<IntervalVector,IntervalMatrix>>;
  std::map<double, std::vector<T>> PEIBOS(const capd::IMap& i_map, double tf, double dt, const AnalyticFunction<VectorType>& psi_0, const std::vector<OctaSym>& Sigma, double epsilon, bool verbose = false);
  using T = std::pair<PEIBOS_CAPD_Key,std::pair<IntervalVector,IntervalMatrix>>;
  std::map<double, std::vector<T>> PEIBOS(const capd::IMap& i_map, double tf, double dt, const AnalyticFunction<VectorType>& psi_0, const std::vector<OctaSym>& Sigma, double epsilon, const Vector& offset, bool verbose = false);

  using T = std::pair<PEIBOS_CAPD_Key,std::pair<IntervalVector,IntervalMatrix>>;
  std::map<double, std::vector<Parallelepiped>> reach_set(const std::map<double, std::vector<T>>& peibos_output);
}