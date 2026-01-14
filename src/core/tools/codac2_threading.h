/** 
 *  \file codac2_threading.h
 * ----------------------------------------------------------------------------
 *  \date       2026
 *  \author     Maël Godard
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include <thread>

namespace codac2_threading
{
  inline int threads_used = 1;
}

namespace codac2
{
  int max_threads();
  void set_threads_used(int n);
  int get_threads_used();
  
}