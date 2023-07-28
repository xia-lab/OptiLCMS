//
//  r_pense_regression.hpp
//  pense
//
//  Created by David Kepplinger on 2019-04-03.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifndef R_PENSE_REGRESSION_HPP_
#define R_PENSE_REGRESSION_HPP_

#include "nsoptim_forward.hpp"

namespace pense {
namespace r_interface {
//! Compute the (Adaptive) PENSE Regularization Path.
//!
//! @param x numeric predictor matrix with `n` rows and `p` columns.
//! @param y numeric response vector with `n` elements.
//! @param penalties a list of EN penalties with decreasing values of the lambda hyper-parameter.
//! @param enpy_inds a vector of 1-based indices for the `penalties` list, at which initial ENPY estimates should be
//!                  computed.
//! @param pense_opts a list of options for the PENSE algorithm.
//! @param enpy_opts a list of options for the ENPY algorithm.
//! @param optional_args a list containing the following named items:
//!                       `shared_starts` ... optional list of coefficients to start at every penalty.
//!                       `individual_starts` ... optional list the same length as `penalties` with a list coefficients
//!                                               to start at the corresponding penalty.
//!                       `pen_loadings` ... optional vector of length `p` with non-negative penalty loadings.
SEXP PenseEnRegression(SEXP x, SEXP y, SEXP penalties, SEXP enpy_inds, SEXP pense_opts, SEXP enpy_opts,
                       SEXP optional_args) noexcept;

//! Get the smallest lambda such that the (Adaptive) PENSE estimate gives the empty model.
//!
//! @param x numeric predictor matrix with `n` rows and `p` columns.
//! @param y numeric response vector with `n` elements.
//! @param pense_opts a list of options for the PENSE algorithm.
//! @param optional_args a list containing the following named items:
//!                       `pen_loadings` ... optional vector of length `p` with non-negative penalty loadings.
SEXP PenseMaxLambda(SEXP x, SEXP y, SEXP pense_opts, SEXP optional_args) noexcept;

}  // namespace r_interface
}  // namespace pense

#endif  // R_PENSE_REGRESSION_HPP_
