/** 
 *  \file codac2_CtcLohner.h
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Simon Rohou
 *  \copyright  Copyright 2025 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include "codac2_Ctc.h"
#include "codac2_AnalyticFunction.h"
#include "codac2_TimePropag.h"
#include "codac2_TDomain.h"

namespace codac2
{
  template<class T>
  class Slice;

  template<typename T>
  class SlicedTube;
  
  class CtcLohner
  {
    public:

      CtcLohner(const AnalyticFunction& f, int contractions = 5, double eps = 0.1);

    protected:

      AnalyticFunction m_f; //!< forward function
      int contractions; //!< number of contractions of the global enclosure by the estimated local enclosure
      int dim; //!< dimension of the state vector
      double eps; //!< inflation parameter for the global enclosure

      static const std::string m_ctc_name; //!< class name (mainly used for CN Exceptions)
      static std::vector<std::string> m_str_expected_doms;
  };
}