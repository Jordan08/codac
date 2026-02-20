/** 
 *  codac2_CtcLohner.cpp
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Simon Rohou
 *  \copyright  Copyright 2025 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include "codac2_CtcLohner.h"
#include "codac2_Slice.h"
#include "codac2_Interval.h"

using namespace std;
using namespace codac2;

CtcLohner::CtcLohner(const Function &f, int contractions, double eps)
    : Ctc(),
      m_f(f),
      contractions(contractions),
      dim(f.nb_var()),
      eps(eps) {}
