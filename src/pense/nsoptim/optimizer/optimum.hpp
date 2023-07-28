//
//  optimum.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2018-11-30.
//  Copyright © 2018 David Kepplinger. All rights reserved.
//

#ifndef NSOPTIM_OPTIMIZER_OPTIMUM_HPP_
#define NSOPTIM_OPTIMIZER_OPTIMUM_HPP_

#include <string>
#include <limits>
#include <memory>
#include <type_traits>

#include "../container/metrics.hpp"

namespace nsoptim {

enum class OptimumStatus { kOk, kWarning, kError };

namespace optimum_internal {
using MetricsPtr = std::unique_ptr<Metrics>;

//! Wrapper around the information at an optimum point.
template <typename T, typename U, typename V>
struct Optimum {
  using LossFunction = T;
  using PenaltyFunction = U;
  using Coefficients = V;
  using ResidualType = typename LossFunction::ResidualType;

  Optimum(const LossFunction& _loss, const PenaltyFunction& _penalty) noexcept : loss(_loss), penalty(_penalty) {}

  Optimum(const LossFunction& _loss, const PenaltyFunction& _penalty, const Coefficients& _coefs,
          const arma::vec& _residuals, const double _objf_value, MetricsPtr _metrics,
          const OptimumStatus _status, const std::string& _message) noexcept
    : loss(_loss), penalty(_penalty), coefs(_coefs), residuals(_residuals), objf_value(_objf_value),
      metrics(std::move(_metrics)), status(_status), message(_message) {}

  Optimum(const Optimum& other) noexcept : loss(other.loss), penalty(other.penalty), coefs(other.coefs),
                                           residuals(other.residuals), objf_value(other.objf_value),
                                           metrics(other.metrics ? new Metrics(*other.metrics) : nullptr),
                                           status(other.status), message(other.message) {}

  Optimum(Optimum&& other) = default;

  //! Move-assignable operator must be explicity defined.
  Optimum& operator=(Optimum&& other) {
    loss = std::move(other.loss);
    penalty = std::move(other.penalty);
    coefs = std::move(other.coefs);
    residuals = std::move(other.residuals);
    objf_value = other.objf_value;
    status = other.status;
    message = std::move(other.message);
    metrics = std::move(other.metrics);
    return *this;
  }

  //! A copy of the loss function for which this optimum was attained.
  LossFunction loss;
  //! A copy of the penalty function for which this optimum was attained.
  PenaltyFunction penalty;
  //! The coefficients at which the objective function attains its optimum.
  Coefficients coefs;
  //! Residuals associated with the loss and coefficients.
  ResidualType residuals;
  //! The value of the objective function at this optimum.
  double objf_value = std::numeric_limits<double>::max();
  //! Optional metrics associated with this optimum.
  MetricsPtr metrics;
  //! The status of the optimizer at the time this optimum was found.
  OptimumStatus status = OptimumStatus::kError;
  //! An optional status message of the optimizer at the time this optimum was found.
  std::string message;
};
}  // namespace optimum_internal

//! Wrapper around the information at an optimum point.
template <typename LossFunction, typename PenaltyFunction, typename Coefficients>
using Optimum = optimum_internal::Optimum<typename std::decay<LossFunction>::type,
                                          typename std::decay<PenaltyFunction>::type,
                                          typename std::decay<Coefficients>::type>;


template <typename L, typename P, typename C>
Optimum<L, P, C> MakeOptimum(
    const L& loss, const P& penalty, const C& coefs, const typename L::ResidualType& residuals, const double objf_value,
    optimum_internal::MetricsPtr metrics,
    const OptimumStatus status = OptimumStatus::kOk, const std::string& message = {}) noexcept {
  return Optimum<L, P, C>(loss, penalty, coefs, residuals, objf_value, std::move(metrics), status, message);
}

template <typename L, typename P, typename C>
Optimum<L, P, C> MakeOptimum(
    const L& loss, const P& penalty, const C& coefs, const typename L::ResidualType& residuals,
    optimum_internal::MetricsPtr metrics,
    const OptimumStatus status = OptimumStatus::kOk, const std::string& message = {}) noexcept {
  return Optimum<L, P, C>(loss, penalty, coefs, residuals, loss(residuals) + penalty(coefs),
                          std::move(metrics), status, message);
}

template <typename L, typename P, typename C>
Optimum<L, P, C> MakeOptimum(
    const L& loss, const P& penalty, const C& coefs,
    optimum_internal::MetricsPtr metrics,
    const OptimumStatus status = OptimumStatus::kOk, const std::string& message = {}) noexcept {
  const auto residuals = loss.Residuals(coefs);
  return Optimum<L, P, C>(loss, penalty, coefs, residuals, loss(residuals) + penalty(coefs),
                          std::move(metrics), status, message);
}

//! Create an Optimum from the given arguments.
template <typename L, typename P, typename C>
Optimum<L, P, C> MakeOptimum(
    const L& loss, const P& penalty, const C& coefs, const typename L::ResidualType& residuals, const double objf_value,
    const OptimumStatus status = OptimumStatus::kOk, const std::string& message = {}) noexcept {
  return Optimum<L, P, C>(loss, penalty, coefs, residuals, objf_value, nullptr,
                                                              status, message);
}

//! Create an Optimum from the given arguments.
template <typename L, typename P, typename C>
Optimum<L, P, C> MakeOptimum(
    const L& loss, const P& penalty, const C& coefs, const typename L::ResidualType& residuals,
    const OptimumStatus status = OptimumStatus::kOk, const std::string& message = {}) noexcept {
  return Optimum<L, P, C>(loss, penalty, coefs, residuals,
                                                              loss(residuals) + penalty(coefs), nullptr, status,
                                                              message);
}

//! Create an Optimum from the given arguments.
template <typename L, typename P, typename C>
Optimum<L, P, C> MakeOptimum(
    const L& loss, const P& penalty, const C& coefs,
    const OptimumStatus status = OptimumStatus::kOk, const std::string& message = {}) noexcept {
  const auto residuals = loss.Residuals(coefs);
  return Optimum<L, P, C>(loss, penalty, coefs, residuals, loss(residuals) + penalty(coefs), nullptr, status, message);
}
}  // namespace nsoptim

#endif  // NSOPTIM_OPTIMIZER_OPTIMUM_HPP_
