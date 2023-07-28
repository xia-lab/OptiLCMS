#ifndef LINEAR_REGRESSION_HPP_
#define LINEAR_REGRESSION_HPP_

#include <RcppArmadillo.h>
#include <type_traits>
#include "pense/nsoptim_forward.hpp"
#include "pense/r_en_regression.hpp"

using namespace Rcpp;
using namespace arma;
using namespace pense::r_interface;

NumericVector PerformLinearRegress(NumericMatrix X, NumericVector Y, NumericVector penalty_loadings);


#endif  // LINEAR_REGRESSION_HPP_
