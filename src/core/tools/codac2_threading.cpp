/** 
 *  codac2_threading.cpp
 * ----------------------------------------------------------------------------
 *  \date       2026
 *  \author     Maël Godard
 *  \copyright  Copyright 2025 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include "codac2_threading.h"

using namespace std;
using namespace codac2;

namespace
{
  int threads_used = 1;
}

namespace codac2
{
  // FOR SIMON : when is this resolved ? 
  // i.e. if I compile the library on a machine (say a github one), and used it on my personnal machine
  // will I get the hardware_concurrency of my machine or the github one ?
  // And what if I just write const int max_threads = std::thread::hardware_concurrency(); ?
  int max_threads()
  {
    return std::thread::hardware_concurrency(); 
  }
  
  void set_threads_used(int n)
  {
    threads_used = (n > 1 ? std::min(n,max_threads()) : 1);
  }

  int get_threads_used()
  {
    return threads_used;
  }
}