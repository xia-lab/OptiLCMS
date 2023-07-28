//
//  r_pensem_regression.hpp
//  pense
//
//  Created by David Kepplinger on 2020-06-08
//  Copyright © 2020 David Kepplinger. All rights reserved.
//

#ifndef R_MESTEN_REGRESSION_HPP_
#define R_MESTEN_REGRESSION_HPP_

#include "nsoptim_forward.hpp"

namespace pense {
namespace r_interface {
//! Compute the (Adaptive) M-EN Regularization Path.
//!
//! @param x numeric predictor matrix with `n` rows and `p` columns.
//! @param y numeric response vector with `n` elements.
//! @param scale scalar, numeric, auxiliary scale.
//! @param penalties a list of EN penalties with decreasing values of the lambda hyper-parameter.
//! @param mest_opts a list of options for the M-estimation algorithm.
//! @param optional_args a list containing the following named items:
//!                      `shared_starts` ... optional list of coefficients to start at every penalty.
//!                      `individual_starts` ... optional list the same length as `penalties` with a list coefficients
//!                                               to start at the corresponding penalty.
//!                      `pen_loadings` ... optional vector of length `p` with non-negative penalty loadings.
SEXP MestEnRegression(SEXP x, SEXP y, SEXP scale, SEXP penalties, SEXP mest_opts, SEXP optional_args) noexcept;

//! Get the smallest lambda such that the (Adaptive) M-EN-estimate gives the empty model.
//!
//! @param x numeric predictor matrix with `n` rows and `p` columns.
//! @param y numeric response vector with `n` elements.
//! @param scale scalar, numeric, auxiliary scale.
//! @param mest_opts a list of options for the M-estimation algorithm.
//! @param optional_args a list containing the following named items:
//!                       `pen_loadings` ... optional vector of length `p` with non-negative penalty loadings.
SEXP MestEnMaxLambda(SEXP x, SEXP y, SEXP scale, SEXP mest_opts, SEXP optional_args) noexcept;

}  // namespace r_interface
}  // namespace pense

#endif  // R_MESTEN_REGRESSION_HPP_
