//
//  robust_scale_location.hpp
//  pense
//
//  Created by David Kepplinger on 2019-01-30.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifndef ROBUST_SCALE_LOCATION_HPP_
#define ROBUST_SCALE_LOCATION_HPP_

#include <exception>
#include <string>

#include "nsoptim.hpp"
#include "rho.hpp"
#include "constants.hpp"
#include "rcpp_utils.hpp"

namespace pense {
namespace robust_scale_location {
//! Default breakdown point for the M-scale equation.
constexpr double kDefaultMscaleDelta = 0.5;
//! Default number of iterations for the M-scale algorithm.
constexpr int kDefaultMscaleMaxIt = 100;

template <typename T>
struct DefaultMscaleConstant {
  constexpr static double value = 0;
};

template <>
struct DefaultMscaleConstant<RhoBisquare> {
  constexpr static double value = kDefaultBisquareMscaleCc;
};

//! Compute an initial estimate of the scale of `values`.
//!
//! @param values a vector of values.
//! @param delta the right-hand-side in the M-scale equation.
//! @param eps numerical tolerance.
//! @return an initial, inaccurate, estimate of the scale of `values`.
double InitialScaleEstimate(const arma::vec& values, const double delta, const double eps);

}  // namespace robust_scale_location

//! The result of simultaneous estimation of the M-location and M-scale.
struct LocationScaleEstimate {
  double location;
  double scale;
};

//! Exception class for when all the weights become zero.
class ZeroWeightsException : public std::runtime_error {
 public:
  ZeroWeightsException() : std::runtime_error("all weights are zero") {}
};

//! Compute the tau-scale of the uncentered values.
//! This function fixes the value of C_2 to 3!
//!
//! @param values values
double TauSize(const arma::vec& values) noexcept;

//! Functor to compute the Mscale using the specified `rho` function.
//! The M-scale `s` is defined by the M-estimation equation:
//!  (1/n) sum_{i = 1}^n rho(value[i] / s) == delta
template <class RhoFunction>
class Mscale {
 public:
  //! Construct the M-scale function.
  //!
  //! @param user_options an R list of user options.
  explicit Mscale(const Rcpp::List& user_options) noexcept
    : rho_(GetFallback(user_options, "cc",
        robust_scale_location::DefaultMscaleConstant<RhoFunction>::value)),
      delta_(GetFallback(user_options, "delta",
        robust_scale_location::kDefaultMscaleDelta)),
      max_it_(GetFallback(user_options, "max_it",
        robust_scale_location::kDefaultMscaleMaxIt)),
      eps_(GetFallback(user_options, "eps", kDefaultConvergenceTolerance)),
      scale_(-1) {}

  //! Construct the M-scale function.
  //!
  //! @param rho rho function to use.
  //! @param delta right-hand side of the M-estimation equation.
  //! @param max_it maximum number of iterations.
  //! @param eps numerical tolerance for convergence.
  Mscale(const RhoFunction& rho, const double delta, const int max_it,
         const double eps) noexcept
      : rho_(rho), delta_(delta), max_it_(max_it), eps_(eps), scale_(-1) {}

  Mscale(const Mscale&) = default;
  Mscale& operator=(const Mscale&) = default;
  Mscale(Mscale&&) = default;
  Mscale& operator=(Mscale&&) = default;

  //! Reset the M-scale estimate.
  void Reset() noexcept {
    scale_ = -1;
  }

  //! Compute the M-scale of the given values. The initial guess is either the given scale (if positive),
  //! the previous scale estimate (if availalbe), or the median of the absolute values (MAD).
  //!
  //! @param values a vector of values.
  //! @return the M-scale of the given values.
  double operator()(const arma::vec& values) const {
    return ComputeMscale(values, InitialEstimate(values));
  }

  //! Compute the M-scale of the given values. The initial guess is either the given scale (if positive),
  //! the previous scale estimate (if availalbe), or the median of the absolute values (MAD).
  //!
  //! @param values a vector of values.
  //! @return the M-scale of the given values.
  double operator()(const arma::vec& values, double scale) const {
    return ComputeMscale(values, scale);
  }

  //! Compute the M-scale of the given values. The previous scale estimate (if present) is used as an initial guess.
  //! If no previous estimate is present (or the M-scale object was reset), the median of the absolute values (MAD)
  //! is used as an initial guess.
  //!
  //! @param values a vector of values.
  //! @return the M-scale of the given values.
  double operator()(const arma::vec& values) {
    scale_ = ComputeMscale(values, InitialEstimate(values));
    return scale_;
  }

  //! Set the initial estimate of the scale.
  //!
  //! @param scale initial scale estimate (> 0).
  void SetInitial(const double scale) noexcept {
    scale_ = scale;
  }

  //! Get the number of iterations required for the last scale estimate.
  //!
  //! @return number of iterations, or -1 if none have been performed yet.
  int LastIterations() const noexcept {
    return it_;
  }

  //! Compute the 1st derivative of the M-scale function with respect to each element.
  //!
  //! @param values vector of values
  //! @return a vector of derivatives, one for each element in `values`. If the scale is 0 or the M-scale equation
  //!   is violated, an empty vector is returned.
  arma::vec Derivative(const arma::vec& values) const {
    const double scale = ComputeMscale(values, InitialEstimate(values));
    if (scale < eps_) {
      return arma::vec();
    }

    const auto deriv_rho = rho_.Derivative(values, scale);
    const auto denom = sum(deriv_rho % values) / scale;
    if (denom < eps_) {
      return arma::vec(values.n_elem, arma::fill::value(R_PosInf));
    } else {
      return deriv_rho / denom;
    }
  }

  //! Compute the gradient and Hessian of the M-scale
  //! function evaluated at the given vector.
  //!
  //! @param values vector of values
  //! @return a matrix of dimension n x n + 1, where the 1st column is the
  //!    gradient and the other columns are the Hessian matrix.
  //!    If the scale is 0 or the M-scale equation
  //!    is violated, an empty matrix is returned.
  arma::mat GradientHessian(const arma::vec& values) {
    const double scale = this->operator()(values);
    if (scale < eps_) {
      return arma::mat(1, 1, arma::fill::value(scale));
    }
    const auto violation = rho_.SumStd(values, scale) / values.n_elem - delta_;

    arma::mat grad_hess(values.n_elem, values.n_elem + 2, arma::fill::zeros);

    // Compute the gradient and its maximum
    grad_hess.col(0) = rho_.Derivative(values, scale);
    const auto denom = arma::sum(grad_hess.col(0) % values);

    grad_hess.at(1, 2) = denom;
    grad_hess.at(2, 2) = scale;
    grad_hess.at(3, 2) = violation;

    // Compute the Hessian and its maximum
    const auto rho_2nd = rho_.SecondDerivative(values, scale);
    const auto sum_2nd = arma::sum(rho_2nd % values % values) / denom;
    grad_hess.col(1) = rho_2nd;
    double diag_offset;
    for (int i = 0; i < values.n_elem; ++i) {
      diag_offset = denom * rho_2nd[i];
      for (int k = i; k < values.n_elem; ++k) {
        grad_hess(i, k + 2) = HessianElementUnscaled(
          i, k, grad_hess.unsafe_col(0), rho_2nd, values, sum_2nd, diag_offset);

        grad_hess(i, k + 2) *= scale / (denom * denom);
        diag_offset = 0;
      }
    }

    // Final pass to get gradient right
    grad_hess.col(0) *= scale / denom;

    return grad_hess;
  }

  //! Compute the maximum of the 1st and 2nd derivatives of the M-scale
  //! function evaluated at all elements in the given vector.
  //!
  //! @param values vector of values
  //! @return a vector with 3 elements: the M-scale,
  //!    the maximum element of the gradient
  //!    and the maximum element in the Hessian.
  //!    If the scale is 0 or the M-scale equation
  //!    is violated, an empty vector is returned.
  arma::vec::fixed<3> MaxGradientHessian(const arma::vec& values) {
    arma::vec::fixed<3> maxima(arma::fill::zeros);
    maxima[0] = this->operator()(values, InitialEstimate(values));
    if (maxima[0] < eps_) {
      return maxima;
    }
    const auto violation = rho_.SumStd(values, maxima[0]) -
      values.n_elem * delta_;

    if (violation * violation > values.n_elem * values.n_elem * eps_ * eps_) {
      return maxima;
    }

    // Compute the gradient and its maximum
    const auto rho_1st = rho_.Derivative(values, maxima[0]);
    const auto denom = sum(rho_1st % values);
    maxima[1] = (denom < eps_) ? R_PosInf :
      (arma::max(rho_1st) * maxima[0] / denom);

    // Compute the Hessian and its maximum
    const auto rho_2nd = rho_.SecondDerivative(values, maxima[0]);
    const auto sum_2nd = sum(rho_2nd % values % values) / denom;
    double diag_offset;

    for (int i = 0; i < values.n_elem; ++i) {
      diag_offset = denom * rho_2nd[i];
      for (int k = i; k < values.n_elem; ++k) {
        const auto tmp = std::abs(HessianElementUnscaled(
          i, k, rho_1st, rho_2nd, values, sum_2nd, diag_offset));

        diag_offset = 0;
        if (tmp > maxima[2]) {
          maxima[2] = tmp;
        }
      }
    }
    maxima[2] *= maxima[0] / (denom * denom);

    return maxima;
  }

  //! Get the rho function object.
  const RhoFunction& rho() const noexcept {
    return rho_;
  }

  //! Get the maximum number of iterations.
  int max_it() const noexcept {
    return max_it_;
  }

  //! Get the numerical tolerance.
  double eps() const noexcept {
    return eps_;
  }

  //! Get delta, i.e., the right-hand-side of the M-scale equation
  double delta() const noexcept {
    return delta_;
  }

 private:
  double ComputeMscale(const arma::vec& values, double scale) const {
    const double rho_denom = 1. / (delta_ * values.n_elem);
    if (scale < kNumericZero) {
      return 0;
    }

    int iter = 0;
    double err = eps_;
    // Start iterations
    do {
      const double rho_sum = rho_.SumStd(values, scale);
      const double new_scale = scale * std::sqrt(rho_sum * rho_denom);
      err = std::abs(new_scale - scale);
      scale = new_scale;
    } while (++iter < max_it_ && err > eps_ * scale);

    return scale;
  }

  double ComputeMscale(const arma::vec& values, double scale) {
    const double rho_denom = 1. / (delta_ * values.n_elem);
    if (scale < kNumericZero) {
      return 0;
    }

    it_ = 0;
    double err = eps_;
    // Start iterations
    do {
      const double rho_sum = rho_.SumStd(values, scale);
      const double new_scale = scale * std::sqrt(rho_sum * rho_denom);
      err = std::abs(new_scale - scale);
      scale = new_scale;
    } while (++it_ < max_it_ && err > eps_ * scale);

    return scale;
  }

  double InitialEstimate(const arma::vec& values) const {
    // If the internal scale is already set, use it as initial estimate.
    if (scale_ > eps_) {
      return scale_;
    }
    return robust_scale_location::InitialScaleEstimate(values, delta_, eps_);
  }

  double HessianElementUnscaled(const int i, const int k,
                                const arma::vec& rho_1st,
                                const arma::vec& rho_2nd,
                                const arma::vec& values,
                                const double sum_2nd,
                                const double diag_offset) const {
    return diag_offset +
      rho_1st[i] * rho_1st[k] * sum_2nd -
      rho_1st[i] * rho_2nd[k] * values[k] -
      rho_1st[k] * rho_2nd[i] * values[i];
  }

  RhoFunction rho_;
  double delta_;
  int max_it_;
  int it_ = -1;
  double eps_;
  double scale_;
};

//! Computation of the M-location of the given vector.
//!
//! @param values values to compute the location and scale from.
//! @param rho rho-function for the M-location.
//! @param scale the scale of the values.
//! @param convergence_tol numeric convergence tolerance.
//! @param max_it maximum number of iterations.
//! @return location of the given values.
template <class RhoFunction>
double MLocation(const arma::vec& values, const RhoFunction& rho, const double scale,
                 const double convergence_tol, const int max_it) {
  const double scaled_conv_tol = convergence_tol * scale;
  int it = 0;
  double location = arma::median(values);

  arma::vec residuals(values.n_elem);
  arma::vec w_loc(values.n_elem);
  while (it++ < max_it) {
    residuals = values - location;
    rho.Weight(residuals, scale, &w_loc);
    const double prev_location = location;
    const double w_loc_sum = arma::accu(w_loc);

    if (w_loc_sum < convergence_tol) {
      throw ZeroWeightsException();
    }

    location = arma::accu(w_loc % values) / w_loc_sum;

    if (std::abs(prev_location - location) < scaled_conv_tol) {
      break;
    }
  }

  return location;
}

//! Simultaneous computation of the M-location and M-scale of the given vector.
//!
//! @param values values to compute the location and scale from.
//! @param mscale M-scale definition.
//! @param location_rho rho-function for the M-location.
//! @return location (first) and scale (second) of the given values.
template <class ScaleRhoFunction, class LocationRhoFunction>
LocationScaleEstimate MLocationScale(const arma::vec& values, const Mscale<ScaleRhoFunction>& mscale,
                                     const LocationRhoFunction& location_rho) {
  int it = 0;
  LocationScaleEstimate est {arma::median(values)};
  est.scale = robust_scale_location::InitialScaleEstimate(values - est.location, mscale.delta(), mscale.eps());

  if (est.scale < mscale.eps()) {
    est.scale = 0;
    return est;
  }

  const double convergence_tol = est.scale * mscale.eps();
  const double recip_sqrt_delta = 1. / sqrt(mscale.delta());

  arma::vec residuals(values.n_elem);
  arma::vec w_loc(values.n_elem);
  while (it++ < mscale.max_it()) {
    residuals = values - est.location;
    location_rho.Weight(residuals, est.scale, &w_loc);
    const double w_scale_mean = mscale.rho().SumStd(residuals, est.scale) / residuals.n_elem;
    const double w_loc_sum = arma::accu(w_loc);

    if (w_loc_sum < convergence_tol) {
      throw ZeroWeightsException();
    }

    const LocationScaleEstimate prev_est = est;
    est.location = arma::accu(w_loc % values) / w_loc_sum;
    est.scale = prev_est.scale * std::sqrt(w_scale_mean) * recip_sqrt_delta;

    if (std::abs(prev_est.location - est.location) < convergence_tol &&
        std::abs(prev_est.scale - est.scale) < convergence_tol) {
      break;
    }
  }

  return est;
}
}  // namespace pense

#endif  // ROBUST_SCALE_LOCATION_HPP_
