//
//  robust_scale_location.cc
//  pense
//
//  Created by David Kepplinger on 2019-01-30.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#include <string>
#include "nsoptim.hpp"
#include "rcpp_utils.hpp"
#include "robust_scale_location.hpp"
#include "constants.hpp"

using arma::vec;
using arma::median;
using arma::abs;
using arma::uword;

namespace {
constexpr double kTauSizeC2Squared = 9.;
constexpr double kTauSizeConsistencyConstant = 1. / 0.961;

constexpr double kMadScaleConsistencyConstant = 1.4826;
}  // namespace

namespace pense {
double TauSize(const vec& values) noexcept {
  const vec abs_values(abs(values));
  const double sigma_0 = median(abs_values);

  if (sigma_0 < kNumericZero) {
    return 0.;
  }

  const double tau_size = arma::mean(arma::clamp(arma::square(abs_values / sigma_0),
                                                 0, kTauSizeC2Squared));
  return sigma_0 * kTauSizeConsistencyConstant * sqrt(tau_size);
}

namespace robust_scale_location {
double InitialScaleEstimate(const vec& values, const double delta, const double eps) {
  // Try the MAD of the uncentered values.
  const double mad = kMadScaleConsistencyConstant * median(abs(values));
  if (mad > eps) {
    return mad;
  } else if (static_cast<uword>((1 - delta) * values.n_elem) > values.n_elem / 2) {
    // If the MAD is also (almost) 0, but the M-scale takes into account more observations than the MAD,
    // compute the variance of the additional elements (i.e., the variance without considering the smallest
    // 50% of the observations)
    const uword lower_index = values.n_elem / 2;
    const uword upper_index = static_cast<uword>((1 - delta) * values.n_elem);
    const vec ordered_values = arma::sort(abs(values));
    const double scale = arma::var(ordered_values.rows(lower_index, upper_index));
    if (scale > eps) {
      return scale;
    }
  }
  return 0.;
}

}  // namespace robust_scale_location

}  // namespace pense
