//
//  coordinate_descent.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2022-02-24.
//  Copyright © 2022 David Kepplinger. All rights reserved.
//

#ifndef NSOPTIM_OPTIMIZER_COORDINATE_DESCENT_HPP_
#define NSOPTIM_OPTIMIZER_COORDINATE_DESCENT_HPP_

#include <exception>
#include <forward_list>
#include <algorithm>
#include <limits>
#include <memory>

#include "../armadillo.hpp"
#include "../utilities.hpp"
#include "../container/regression_coefficients.hpp"
#include "../container/data.hpp"
#include "optimizer_base.hpp"
#include "optimum.hpp"
#include "../objective/ls_regression_loss.hpp"
#include "../objective/en_penalty.hpp"
#include "soft_threshold.hpp"
#include "../traits/traits.hpp"

namespace nsoptim {
//! Configuration options for the DAL algorithm.
struct CDConfiguration {
  //! Maximum number of iterations allowed.
  int max_it;
  //! Re-compute the residuals every `reset_iter` iterations to avoid drift.
  int reset_iter;
};

namespace coorddesc {
constexpr CDConfiguration kDefaultCDConfiguration = { 1000, 8 };

template<class Coefficients>
struct State {
  Coefficients coefs;
  arma::vec residuals;
};

} // namespace coorddesc

//! Compute the EN regression estimate using the LARS algorithm on the
//! augmented response vector and predictor matrix.
template<class LossFunction, class PenaltyFunction, class Coefficients>
class CoordinateDescentOptimizer :
    public Optimizer<LossFunction, PenaltyFunction, Coefficients> {
  using Base = Optimizer<LossFunction, PenaltyFunction, Coefficients>;
  using LossFunctionPtr = std::unique_ptr<LossFunction>;
  using PenaltyPtr = std::unique_ptr<PenaltyFunction>;
  using IsWeightedTag = typename traits::is_weighted<LossFunction>::type;
  using IsAdaptiveTag = typename traits::is_adaptive<PenaltyFunction>::type;
  using IsSparseTag = typename
    std::is_same<typename Coefficients::SlopeCoefficient, arma::sp_vec>::type;
  using EnThreshold = typename
    std::conditional<IsAdaptiveTag::value, arma::vec, double>::type;

  static_assert(traits::is_en_penalty<PenaltyFunction>::value,
                "PenaltyFunction must be an EN-type penalty.");
  static_assert(traits::is_ls_regression_loss<LossFunction>::value,
                "LossFunction must be an least-squares-type loss.");

 public:
  using Optimum = typename Base::Optimum;

    //! Ininitialize the optimizer without a loss or penalty function.
  CoordinateDescentOptimizer(
    const CDConfiguration& config = coorddesc::kDefaultCDConfiguration) noexcept
      : config_(config) {}

  //! Ininitialize the optimizer using the given (weighted) LS loss function
  //! and penalty function.
  //! @param loss a weighted LS loss function.
  //! @param penalty penalty function.
  CoordinateDescentOptimizer(const LossFunction& loss,
    const PenaltyFunction& penalty,
    const CDConfiguration& config = coorddesc::kDefaultCDConfiguration) noexcept
    : loss_(new LossFunction(loss)),
      penalty_(new PenaltyFunction(penalty)), config_(config) {}

  //! Default copy constructor.
  //!
  //! The copied optimizer will share the identical loss and penalty
  //! functions after construction.
  CoordinateDescentOptimizer(const CoordinateDescentOptimizer& other) noexcept
    : loss_(other.loss_? new LossFunction(*other.loss_) : nullptr),
      penalty_(other.penalty_ ? new PenaltyFunction(*other.penalty_) : nullptr),
      config_(other.config_),
      state_(other.state_),
      convergence_tolerance_(other.convergence_tolerance_) {}

  //! Default copy assignment.
  //!
  //! The copied optimizer will share the identical loss and penalty
  //! functions after construction.
  CoordinateDescentOptimizer& operator=(const CoordinateDescentOptimizer& other)
    = default;

  //! Default move constructor.
  CoordinateDescentOptimizer(CoordinateDescentOptimizer&& other)
    = default;

  //! Default move assignment operator.
  CoordinateDescentOptimizer& operator=(CoordinateDescentOptimizer&& other)
    = default;

  ~CoordinateDescentOptimizer() = default;

  void Reset() {
    loss_.reset();
    penalty_.reset();
    state_.residuals.reset();
  }

  //! Get the convergence tolerance for the CD algorithm.
  //!
  //! @return convergence tolerance.
  double convergence_tolerance() const noexcept {
    return convergence_tolerance_;
  }

  //! Set the convergence tolerance for the CD algorithm.
  //!
  //! @param convergence_tolerance convergene tolerance for the MM algorithm.
  void convergence_tolerance(double convergence_tolerance) noexcept {
    convergence_tolerance_ = convergence_tolerance;
  }

  //! Get the current loss function.
  LossFunction& loss() const {
    if (!loss_) {
      throw std::logic_error("no loss set");
    }
    return *loss_;
  }

  //! Set the new loss function.
  void loss(const LossFunction& loss) noexcept {
    loss_.reset(new LossFunction(loss));
    ls_stepsize_.reset();
  }

  PenaltyFunction& penalty() const {
    if (!penalty_) {
      throw std::logic_error("no penalty set");
    }
    return *penalty_;
  }

  void penalty(const PenaltyFunction& penalty) noexcept {
    penalty_.reset(new PenaltyFunction(penalty));
    ResetEnThreshold(IsAdaptiveTag{});
    ls_stepsize_.reset();
  }

  //! Find the minimum of the objective function, using the previous solution
  //! (or the 0-vector if no previous solution exists) as starting point.
  //!
  //! @return information about the optimum.
  Optimum Optimize() {
    return Optimize(config_.max_it);
  }

  //! Find the minimum of the objective function, using the given coefficients
  //! as starting point and at most ``max_it`` iterations.
  //!
  //! @param start where to start the optimization from.
  //! @return information about the optimum.
  Optimum Optimize(const Coefficients& start) {
    ResetState(start);
    return Optimize(config_.max_it);
  }

  //! Find the minimum of the objective function, using the given coefficients
  //! as starting point and at most ``max_it`` iterations.
  //!
  //! @param start where to start the optimization from.
  //! @param max_it maximum number of iterations.
  //! @return information about the optimum.
  Optimum Optimize(const Coefficients& start, const int max_it) {
    ResetState(start);
    return Optimize(max_it);
  }

  //! Find the minimum of the objective function, using the previous solution
  //! (or the 0-vector if no previous solution exists) as starting point and
  //! at most ``max_it`` iterations.
  //!
  //! @param max_it maximum number of iterations.
  //! @return information about the optimum.
  Optimum Optimize(const int max_it) {
    if (!loss_) {
      throw std::logic_error("no loss set");
    }
    if (!penalty_) {
      throw std::logic_error("no penalty set");
    }

    auto metrics = std::make_unique<Metrics>("cd-ls_en");

    if (state_.residuals.n_elem == 0) {
      ResetState(loss_->template ZeroCoefficients<Coefficients>());
    }

    UpdateEnThreshold(IsAdaptiveTag{});
    if (ls_stepsize_.n_elem == 0) {
      UpdateLsStepSize(IsWeightedTag{}, IsAdaptiveTag{});
    }

    int iter = 0;
    const auto& data = loss_->data();
    while (iter++ < max_it) {
      Metrics& iteration_metrics = metrics->CreateSubMetrics("cd_iteration");

      auto prev_coefs = state_.coefs;
      double total_change = 0;
      if (loss_->IncludeIntercept()) {
        state_.coefs.intercept = UpdateIntercept(IsWeightedTag{});
        const auto diff = prev_coefs.intercept - state_.coefs.intercept;
        state_.residuals += diff;
        total_change += std::abs(diff);
      }

      for (arma::uword j = 0; j < data.n_pred(); ++j) {
        // @TODO -- this is inefficient if we have a sparse vector!
        state_.coefs.beta[j] = UpdateSlope(j, IsWeightedTag{}, IsAdaptiveTag{});
        const auto diff = prev_coefs.beta[j] - state_.coefs.beta[j];
        if (diff != 0) {
          state_.residuals += diff * data.cx().col(j);
          total_change += std::abs(diff);
        }
      }

      iteration_metrics.AddMetric("iter", iter);
      iteration_metrics.AddMetric("change", total_change);

      if (total_change < data.n_pred() * convergence_tolerance_) {
        metrics->AddMetric("iter", iter);
        state_.residuals = loss_->Residuals(state_.coefs);
        return MakeOptimum(*loss_, *penalty_, state_.coefs, state_.residuals,
                           std::move(metrics));
      }

      // Re-compute the residuals after every few cycles to avoid any drifts
      if (iter > 0 && iter % config_.reset_iter == 0) {
        state_.residuals = loss_->Residuals(state_.coefs);
      }
    }

    metrics->AddMetric("iter", iter);
    state_.residuals = loss_->Residuals(state_.coefs);
    return MakeOptimum(*loss_, *penalty_, state_.coefs, state_.residuals,
                       std::move(metrics), OptimumStatus::kWarning,
                       "Coordinate descent did not converge.");
  }

 private:
  double UpdateIntercept (std::false_type /* is_weighted */) {
    return arma::accu(state_.residuals + state_.coefs.intercept);
  }

  double UpdateIntercept (std::true_type /* is_weighted */) {
    return arma::mean((state_.residuals + state_.coefs.intercept) %
                      arma::square(loss_->sqrt_weights()));
  }

  double UpdateSlope (arma::uword j, std::false_type /* is_weighted */,
                      std::false_type /* is_adaptive */) {
    const double dir = (state_.coefs.beta[j] != 0) ?
      arma::dot(loss_->data().cx().col(j), state_.residuals +
                  state_.coefs.beta[j] * loss_->data().cx().col(j)) :
      arma::dot(loss_->data().cx().col(j), state_.residuals);
    return SoftThreshold(dir, en_softthresh_) / ls_stepsize_[j];
  }

  double UpdateSlope (arma::uword j, std::true_type /* is_weighted */,
                      std::false_type /* is_adaptive */) {
    const double dir = (state_.coefs.beta[j] != 0) ?
      arma::dot(arma::square(loss_->sqrt_weights()) % loss_->data().cx().col(j),
        state_.residuals + state_.coefs.beta[j] * loss_->data().cx().col(j)) :
      arma::dot(arma::square(loss_->sqrt_weights()) % loss_->data().cx().col(j),
                state_.residuals);

    return SoftThreshold(dir, en_softthresh_ / loss_->mean_weight()) /
      ls_stepsize_[j];
  }

  double UpdateSlope (arma::uword j, std::false_type /* is_weighted */,
                      std::true_type /* is_adaptive */) {
    const double dir = (state_.coefs.beta[j] != 0) ?
      arma::dot(loss_->data().cx().col(j), state_.residuals +
                  state_.coefs.beta[j] * loss_->data().cx().col(j)) :
      arma::dot(loss_->data().cx().col(j), state_.residuals);
    return SoftThreshold(dir, en_softthresh_[j]) / ls_stepsize_[j];
  }

  double UpdateSlope (arma::uword j, std::true_type /* is_weighted */,
                      std::true_type /* is_adaptive */) {
    const double dir = (state_.coefs.beta[j] != 0) ?
      arma::dot(arma::square(loss_->sqrt_weights()) % loss_->data().cx().col(j),
                state_.residuals + state_.coefs.beta[j] *
                  loss_->data().cx().col(j)) :
      arma::dot(arma::square(loss_->sqrt_weights()) % loss_->data().cx().col(j),
                state_.residuals);

    return SoftThreshold(dir, en_softthresh_[j] / loss_->mean_weight()) /
      ls_stepsize_[j];
  }

  void ResetState (const Coefficients &coefs) {
    if (!loss_) {
      throw std::logic_error("no loss set");
    }
    if (!penalty_) {
      throw std::logic_error("no penalty set");
    }
    state_ = { coefs, loss_->Residuals(coefs) };
  }

  void ResetEnThreshold(std::false_type /* is_adaptive */) {}

  void ResetEnThreshold(std::true_type /* is_adaptive */) {
    en_softthresh_.reset();
  }

  void UpdateLsStepSize(std::false_type /* is_weighted */,
                        std::false_type /* is_adaptive */) {
    ls_stepsize_ = arma::sum(arma::square(loss_->data().cx())).t() +
      loss_->data().n_obs() * penalty_->lambda() * (1 - penalty_->alpha());
  }

  void UpdateLsStepSize(std::true_type /* is_weighted */,
                        std::false_type /* is_adaptive */) {
    const arma::vec ls = arma::sum(arma::square(loss_->data().cx().each_col() %
      loss_->sqrt_weights())).t();
    const double en = loss_->data().n_obs() * penalty_->lambda() *
      (1 - penalty_->alpha()) / loss_->mean_weight();
    ls_stepsize_ = ls + en;
  }

  void UpdateLsStepSize(std::false_type /* is_weighted */,
                        std::true_type /* is_adaptive */) {
    ls_stepsize_ = arma::sum(arma::square(loss_->data().cx())).t() +
      penalty_->loadings() * loss_->data().n_obs() * penalty_->lambda() *
        (1 - penalty_->alpha());
  }

  void UpdateLsStepSize(std::true_type /* is_weighted */,
                        std::true_type /* is_adaptive */) {
    const arma::vec ls = arma::sum(arma::square(loss_->data().cx().each_col() %
      loss_->sqrt_weights())).t();
    const arma::vec adaen = penalty_->loadings() * loss_->data().n_obs() *
      penalty_->lambda() * (1 - penalty_->alpha()) / loss_->mean_weight();
    ls_stepsize_ = ls + adaen;
  }

  void UpdateEnThreshold(std::false_type /* is_adaptive */) {
    en_softthresh_ = loss_->data().n_obs() * penalty_->lambda() *
      penalty_->alpha();
  }

  void UpdateEnThreshold(std::true_type /* is_adaptive */) {
    if (en_softthresh_.n_elem == 0) {
      en_softthresh_ = loss_->data().n_obs() * penalty_->loadings() *
        penalty_->lambda() * penalty_->alpha();
    }
  }

  LossFunctionPtr loss_;
  PenaltyPtr penalty_;
  CDConfiguration config_;
  arma::vec ls_stepsize_;
  EnThreshold en_stepsize_;
  EnThreshold en_softthresh_;
  coorddesc::State<Coefficients> state_;
  double convergence_tolerance_ = 1e-8;
};
} // namespace nsoptim

#endif // NSOPTIM_OPTIMIZER_COORDINATE_DESCENT_HPP_
