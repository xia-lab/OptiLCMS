//
//  r_en_regression.hpp
//  pense
//
//  Created by David Kepplinger on 2019-04-03.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifndef R_EN_REGRESSION_HPP_
#define R_EN_REGRESSION_HPP_

#include "nsoptim_forward.hpp"

namespace pense {
namespace r_interface {
//! Compute the EN Regularization Path.
//!
//! @param x numeric predictor matrix with `n` rows and `p` columns.
//! @param y numeric response vector with `n` elements.
//! @param penalties a list of EN penalties with decreasing values of the lambda hyper-parameter.
//! @param include_intercept include an intercept in the loss function?
//! @param optional_args a list containing the following named items:
//!                       `en_options` ... control options for the EN algorithm
//!                       `obs_weights` ... optional vector of length `n` with non-negative observation weights.
//!                       `pen_loadings` ... optional vector of length `p` with non-negative penalty loadings
SEXP LsEnRegression(SEXP x, SEXP y, SEXP penalties, SEXP include_intercept, SEXP optional_args) noexcept;
}  // namespace r_interface
}  // namespace pense

#endif  // R_EN_REGRESSION_HPP_
