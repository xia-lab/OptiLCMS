//
//  nsoptim_forward.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2018-11-30.
//  Copyright © 2018 David Kepplinger. All rights reserved.
//

#ifndef NSOPTIM_FORWARD_HPP_
#define NSOPTIM_FORWARD_HPP_

#include "autoconfig.hpp"
#include "nsoptim/config.hpp"
#include "nsoptim/armadillo_forward.hpp"
#include "nsoptim/container/forward.hpp"
#include "nsoptim/objective/forward.hpp"

namespace Rcpp {
// Specialize Rcpp::wrap for RegressionCoefficients from <nsoptim/container/regression_coefficients.h>
template<typename T>
SEXP wrap(const nsoptim::RegressionCoefficients<T>&);

//! Specialize Rcpp::wrap for armadillo sparse vectors
template<typename T>
SEXP wrap(const arma::SpCol<T>&);

namespace traits {
//! Specialize Rcpp::as for armadillo sparse vectors
template<typename T>
class Exporter<arma::SpCol<T>>;

// Specialize Rcpp::as for RegressionCoefficients from <nsoptim/container/regression_coefficients.h>
template<typename T>
class Exporter<nsoptim::RegressionCoefficients<T>>;
}  // namespace traits
}  // namespace Rcpp

#endif  // NSOPTIM_FORWARD_HPP_
