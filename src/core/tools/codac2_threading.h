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

namespace codac2
{
  inline int threads_used = 1;

  // FOR SIMON : when is this resolved ? 
  // i.e. if I compile the library on a machine (say a github one), and used it on my personnal machine
  // will I get the hardware_concurrency of my machine or the github one ?
  // And what if I just write const int max_threads = std::thread::hardware_concurrency(); ?
  inline int max_threads()
  {
    return std::thread::hardware_concurrency(); 
  }
  
  inline void set_threads_used(int n)
  {
    threads_used = (n > 1 ? std::min(n,max_threads()) : 1);
  }

}