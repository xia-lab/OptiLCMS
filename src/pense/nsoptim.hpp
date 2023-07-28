//
//  nsoptim.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2018-11-30.
//  Copyright © 2018 David Kepplinger. All rights reserved.
//

#ifndef NSOPTIM_HPP_
#define NSOPTIM_HPP_

#if defined(Rcpp_hpp)
    #error "The file 'Rcpp.h' should not be included. Please correct to include only 'nsoptim.hpp'."
#endif
#if defined(RcppArmadillo__RcppArmadillo__h)
    #error "The file 'RcppArmadillo.h' should not be included. Please correct to include only 'nsoptim.hpp'."
#endif

#include "autoconfig.hpp"
#include "nsoptim_forward.hpp"

#include "nsoptim/armadillo.hpp"
#include "nsoptim/rcpp_integration.hpp"
#include "nsoptim/objective.hpp"
#include "nsoptim/container.hpp"
#include "nsoptim/optimizer.hpp"

#endif  // NSOPTIM_HPP_
