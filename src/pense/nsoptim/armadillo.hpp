//
//  armadillo.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2018-11-30.
//  Copyright © 2018 David Kepplinger. All rights reserved.
//

#ifndef NSOPTIM_ARMADILLO_HPP_
#define NSOPTIM_ARMADILLO_HPP_

#ifndef ARMA_USE_CXX11
#  define ARMA_USE_CXX11 1
#endif

#ifndef ARMA_DONT_USE_OPENMP
#  define ARMA_DONT_USE_OPENMP 1
#endif

#ifdef __clang__
#  pragma clang diagnostic push
#  pragma clang diagnostic ignored "-Weverything"
#endif

# include <RcppArmadillo.h>


#ifdef __clang__
#  pragma clang diagnostic pop
#endif

#endif  // NSOPTIM_ARMADILLO_HPP_
