/* src/autoconfig.hpp.  Generated from autoconfig.hpp.in by configure.  */
//
//  autoconfig.hpp
//  pense
//
//  Created by David Kepplinger on 2019-04-03.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifndef AUTOCONFIG_HPP_
#define AUTOCONFIG_HPP_

/* OpenMP is optional. Enable it only when the compiler can see omp.h. */
#if defined(__has_include)
#  if __has_include(<omp.h>)
#    define PENSE_ENABLE_OPENMP 1
#  else
#    define PENSE_DISABLE_OPENMP 1
#  endif
#else
#  define PENSE_ENABLE_OPENMP 1
#endif
/* #undef PENSE_OPENMP_ADD_CONST_SHARED */
#define NSOPTIM_METRICS_DISABLED 1
/* #undef NSOPTIM_METRICS_ENABLED */
/* #undef NSOPTIM_METRICS_DETAILED */

#endif  // AUTOCONFIG_HPP_
