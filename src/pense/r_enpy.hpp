//
//  r_enpy.hpp
//  pense
//
//  Created by David Kepplinger on 2019-04-03.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifndef R_ENPY_HPP_
#define R_ENPY_HPP_

#include "nsoptim_forward.hpp"

namespace pense {
namespace r_interface {
//! Compute the (Adaptive) penalized PY Initial Estimators.
//!
//! @param x numeric predictor matrix with `n` rows and `p` columns.
//! @param y numeric response vector with `n` elements.
//! @param penalties a list of EN penalties with decreasing values of the lambda hyper-parameter.
//! @param sloss_params parameters for the M-scale.
//! @param enpy_opts a list of options for the ENPY algorithm.
//! @param optional_args a list containing the following named items:
//!                       `pen_loadings` ... optional vector of length `p` with non-negative penalty loadings.
SEXP PenPyInitialEstimator(SEXP x, SEXP y, SEXP penalties, SEXP sloss_params, SEXP enpy_opts,
                           SEXP optional_args) noexcept;

//! Compute the Principal Sensitivity Components.
//!
//! @param x numeric predictor matrix with `n` rows and `p` columns.
//! @param y numeric response vector with `n` elements.
//! @param penalties a list of EN penalties with decreasing values of the lambda hyper-parameter.
//! @param en_options options for the EN algorithm.
//! @param optional_args a list containing the following named items:
//!                       `intercept` ... boolean determining if an intercept should be included.
//!                       `num_threads` ... number of threads.
//!                       `pen_loadings` ... optional vector of length `p` with non-negative penalty loadings.
SEXP PrincipalSensitivityComponents(SEXP x, SEXP y, SEXP penalties, SEXP en_options, SEXP optional_args) noexcept;

}  // namespace r_interface
}  // namespace pense

#endif  // R_ENPY_HPP_
