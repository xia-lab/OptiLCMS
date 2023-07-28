//
//  coordinate_descent.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2022-02-24.
//  Copyright © 2022 David Kepplinger. All rights reserved.
//

#ifndef PENSE_CD_PENSE_HPP_
#define PENSE_CD_PENSE_HPP_

#include <exception>
#include <forward_list>
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>

#include "nsoptim.hpp"
#include "s_loss.hpp"
#include "robust_scale_location.hpp"

namespace pense {
//! Configuration options for the DAL algorithm.
struct CDPenseConfiguration {
  //! Maximum number of iterations allowed.
  int max_it;
  //! Multiplier for adjusting the step size during line search.
  double linesearch_ss_multiplier;
  //! Number of step sizes to be considered for line search.
  int linesearch_ss_num;
  //! Re-compute the residuals every `reset_iter` iterations to avoid drift.
  int reset_iter;
};

namespace coorddesc {
constexpr CDPenseConfiguration kDefaultCDConfiguration = { 1000, 0.5, 10, 8 };

struct SurrogateGradient {
  const double gradient;
  double lipschitz_constant;
};

template<class Coefficients>
struct State {
  Coefficients coefs;
  arma::vec residuals;
  double mscale;
  double objf_loss;
  double objf_pen;
};

} // namespace coorddesc

//! Compute the EN regression estimate using the LARS algorithm on the
//! augmented response vector and predictor matrix.
template<class PenaltyFunction, class Coefficients>
class CDPense :
    public nsoptim::Optimizer<SLoss, PenaltyFunction, Coefficients> {
  using Base = nsoptim::Optimizer<SLoss, PenaltyFunction, Coefficients>;
  using LossFunctionPtr = std::unique_ptr<SLoss>;
  using PenaltyPtr = std::unique_ptr<PenaltyFunction>;
  using IsAdaptiveTag = typename nsoptim::traits::is_adaptive<PenaltyFunction>::type;
  using IsSparseTag = typename
    std::is_same<typename Coefficients::SlopeCoefficient, arma::sp_vec>::type;

  static_assert(nsoptim::traits::is_en_penalty<PenaltyFunction>::value,
                "PenaltyFunction must be an EN-type penalty.");

 public:
  using Optimum = typename Base::Optimum;

    //! Ininitialize the optimizer without a loss or penalty function.
  CDPense(
    const CDPenseConfiguration& config = coorddesc::kDefaultCDConfiguration) noexcept
      : config_(config) {}

  //! Ininitialize the optimizer using the given (weighted) LS loss function
  //! and penalty function.
  //! @param loss a weighted LS loss function.
  //! @param penalty penalty function.
  CDPense(const SLoss& loss,
    const PenaltyFunction& penalty,
    const CDPenseConfiguration& config = coorddesc::kDefaultCDConfiguration) noexcept
    : loss_(new SLoss(loss)),
      penalty_(new PenaltyFunction(penalty)), config_(config) {}

  //! Default copy constructor.
  //!
  //! The copied optimizer will share the identical loss and penalty
  //! functions after construction.
  CDPense(const CDPense& other) noexcept
    : loss_(other.loss_? new SLoss(*other.loss_) : nullptr),
      penalty_(other.penalty_ ? new PenaltyFunction(*other.penalty_) : nullptr),
      config_(other.config_),
      lipschitz_bounds_(other.lipschitz_bounds_),
      lipschitz_bound_intercept_(other.lipschitz_bound_intercept_),
      state_(other.state_),
      convergence_tolerance_(other.convergence_tolerance_) {}

  //! Default copy assignment.
  //!
  //! The copied optimizer will share the identical loss and penalty
  //! functions after construction.
  CDPense& operator=(const CDPense& other) = default;

  //! Default move constructor.
  CDPense(CDPense&& other) = default;

  //! Default move assignment operator.
  CDPense& operator=(CDPense&& other) = default;

  ~CDPense() = default;

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
  SLoss& loss() const {
    if (!loss_) {
      throw std::logic_error("no loss set");
    }
    return *loss_;
  }

  //! Set the new loss function.
  void loss(const SLoss& loss) noexcept {
    loss_ = std::make_unique<SLoss>(loss);
    lipschitz_bounds_.reset();
  }

  PenaltyFunction& penalty() const {
    if (!penalty_) {
      throw std::logic_error("no penalty set");
    }
    return *penalty_;
  }

  void penalty(const PenaltyFunction& penalty) noexcept {
    penalty_ = std::make_unique<PenaltyFunction>(penalty);
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

    auto metrics = std::make_unique<nsoptim::Metrics>("cd-pense");

    if (state_.residuals.n_elem == 0) {
      ResetState(loss_->template ZeroCoefficients<Coefficients>());
    }

    if (lipschitz_bounds_.n_elem == 0) {
      UpdateLipschitzBounds();
    }

    int iter = 0;
    const auto& data = loss_->data();

    while (iter++ < max_it) {
      double coef_change = 0;
      auto& iteration_metrics = metrics->CreateSubMetrics("cd_iteration");
      iteration_metrics.AddMetric("iter", iter);

      const double objf_before_iter = state_.objf_loss + state_.objf_pen;

      for (arma::uword j = 0; j < data.n_pred(); ++j) {
        // @TODO -- this iteration is inefficient if we have a sparse vector!
        auto& cycle_metrics = iteration_metrics.CreateSubMetrics("coordinate");
        auto gradlip = GradientAndSurrogateLipschitz(j);
        int total_mscale_iterations = 0;
        const double objf_pen_prev = state_.objf_pen - PenaltyContribution(state_.coefs.beta[j], j, IsAdaptiveTag{});
        double updated_coef = state_.coefs.beta[j];

        cycle_metrics.AddMetric("index", static_cast<int>(j));
        cycle_metrics.AddMetric("gradient", gradlip.gradient);
        cycle_metrics.AddMetric("lipschitz", lipschitz_bounds_[j]);
        cycle_metrics.AddMetric("lipschitz_surrogate", gradlip.lipschitz_constant);

        int ls_step = 0;
        bool improved = false;

        while (ls_step++ < config_.linesearch_ss_num) {
          const double try_coef = UpdateSlope(j, gradlip.lipschitz_constant, gradlip.gradient, IsAdaptiveTag{});

          if (std::abs(try_coef - state_.coefs.beta[j]) > kNumericZero) {
            state_.residuals += (updated_coef - try_coef) * data.cx().col(j);

            const auto eval_loss = loss_->EvaluateResiduals(state_.residuals);
            const double new_objf_pen = objf_pen_prev + PenaltyContribution(try_coef, j, IsAdaptiveTag{});
            const double new_objf = eval_loss.loss + new_objf_pen;

            total_mscale_iterations += loss_->mscale().LastIterations();

            if (new_objf < state_.objf_loss + state_.objf_pen + convergence_tolerance_) {
              // The objective function improved or did not change much. Stop here.
              coef_change += std::abs(state_.coefs.beta[j] - try_coef);

              state_.coefs.beta[j] = try_coef;
              state_.objf_loss = eval_loss.loss;
              state_.objf_pen = new_objf_pen;
              state_.mscale = eval_loss.scale;

              improved = true;
              cycle_metrics.AddMetric("ls_stepsize", gradlip.lipschitz_constant);
              break;
            } else if (gradlip.lipschitz_constant >= lipschitz_bounds_[j]) {
              // We are at the upper end of the step size range and haven't seen an improvement.
              // Stop here.
              updated_coef = try_coef;
              break;
            }

            updated_coef = try_coef;
            gradlip.lipschitz_constant /= config_.linesearch_ss_multiplier;
          } else {
            break;
          }
        }

        if (!improved) {
          if (std::abs(updated_coef - state_.coefs.beta[j]) > kNumericZero) {
            state_.residuals += (updated_coef - state_.coefs.beta[j]) * data.cx().col(j);
          }
          cycle_metrics.AddMetric("ls_stepsize", 0.);
        }

        cycle_metrics.AddMetric("ls_steps", ls_step);
        cycle_metrics.AddMetric("mscale_iterations", total_mscale_iterations);
      }

      // After updating the slope coefficients, update the intercept.
      if (loss_->IncludeIntercept()) {
        // If we don't do linesearch, also compute the lipschitz constant of the surrogate WLS loss.
        auto gradlip = GradientAndSurrogateLipschitz();
        int total_mscale_iterations = 0;
        double updated_coef = state_.coefs.intercept;
        iteration_metrics.AddMetric("gradient_int", gradlip.gradient);
        iteration_metrics.AddMetric("lipschitz_int", lipschitz_bound_intercept_);
        iteration_metrics.AddMetric("lipschitz_int_surrogate", gradlip.lipschitz_constant);

        int ls_step = 0;
        bool improved = false;

        while (ls_step++ < config_.linesearch_ss_num) {
          const double try_coef = state_.coefs.intercept - gradlip.gradient / gradlip.lipschitz_constant;
          state_.residuals += updated_coef - try_coef;
          const auto eval_loss = loss_->EvaluateResiduals(state_.residuals);
          total_mscale_iterations += loss_->mscale().LastIterations();

          if (eval_loss.loss < state_.objf_loss + convergence_tolerance_) {
            // The objective function improved or did not change much. Stop here.
            coef_change += std::abs(state_.coefs.intercept - try_coef);

            state_.coefs.intercept = try_coef;
            state_.objf_loss = eval_loss.loss;
            state_.mscale = eval_loss.scale;

            improved = true;

            iteration_metrics.AddMetric("ls_stepsize_int", gradlip.lipschitz_constant);
            iteration_metrics.AddMetric("ls_steps_int", ls_step);
            break;
          } else if (gradlip.lipschitz_constant > lipschitz_bound_intercept_) {
            // We are at the upper end of the step size range and haven't seen an improvement.
            // Stop here.
            updated_coef = try_coef;
            state_.residuals += try_coef - state_.coefs.intercept;
            // The coefficient value did not change.
            break;
          }

          updated_coef = try_coef;
          gradlip.lipschitz_constant /= config_.linesearch_ss_multiplier;
        }

        if (!improved) {
          state_.residuals += updated_coef - state_.coefs.intercept;
          iteration_metrics.AddMetric("ls_stepsize_int", 0.);
        }

        iteration_metrics.AddMetric("ls_steps_int", ls_step);
        iteration_metrics.AddMetric("mscale_iterations_int", total_mscale_iterations);
      }

      const double objf_change = (state_.objf_loss + state_.objf_pen) - objf_before_iter;
      iteration_metrics.AddMetric("change", objf_change);
      iteration_metrics.AddMetric("coef_change", coef_change);

      if (objf_change * objf_change < convergence_tolerance_ * convergence_tolerance_) {
        // The objective function value did not change. Algorithm converged.
        metrics->AddMetric("iter", iter);
        return nsoptim::MakeOptimum(*loss_, *penalty_, state_.coefs, state_.residuals,
                                    std::move(metrics));
      }

      // Re-compute the residuals after every few cycles to avoid any drifts
      if (iter > 0 && iter % config_.reset_iter == 0) {
        state_.residuals = loss_->Residuals(state_.coefs);
      }
    }

    metrics->AddMetric("iter", iter);
    state_.residuals = loss_->Residuals(state_.coefs);
    return nsoptim::MakeOptimum(*loss_, *penalty_, state_.coefs, state_.residuals,
                                std::move(metrics), nsoptim::OptimumStatus::kWarning,
                                "Coordinate descent did not converge.");
  }

 private:
  void UpdateLipschitzBounds() {
    const auto data = loss_->data();
    const auto ms = loss_->mscale();
    const double eff_n = data.n_obs() * (1. - ms.delta());
    const double separation = eff_n - std::floor(eff_n);
    const double mult = std::log(separation * (1 - separation)) / std::cbrt(eff_n);
    const double u1 = std::min(80., -40. * mult) / ms.rho().cc();
    const double u2 = std::min(50., 100. * mult * mult * mult * mult) / ms.rho().cc();
    lipschitz_bounds_ = u1 * u1 * arma::square(arma::sum(data.cx())).t();
    for (arma::uword j = 0; j < data.n_pred(); ++j) {
      lipschitz_bounds_[j] += u2 * state_.mscale * std::abs(arma::accu(data.cx().col(j) * data.cx().col(j).t()));
    }

    lipschitz_bound_intercept_ = (u1 * u1 + u2 * state_.mscale) * data.n_obs() * data.n_obs();
  }

  coorddesc::SurrogateGradient GradientAndSurrogateLipschitz() {
    const arma::vec wgt = loss_->mscale().rho().Weight(state_.residuals, state_.mscale);
    const double gradient = -state_.mscale * state_.mscale * arma::dot(wgt, state_.residuals) /
      arma::dot(wgt, arma::square(state_.residuals));
    const double lipschitz = 2 * arma::mean(wgt);
    return coorddesc::SurrogateGradient { gradient, lipschitz };
  }

  coorddesc::SurrogateGradient GradientAndSurrogateLipschitz(const arma::uword j) {
    auto&& xmat = loss_->data().cx();
    const arma::vec wgt = loss_->mscale().rho().Weight(state_.residuals, state_.mscale);
    const double gradient = -state_.mscale * state_.mscale * arma::dot(wgt % xmat.col(j), state_.residuals) /
      arma::dot(wgt, arma::square(state_.residuals));
    const double lipschitz = 2 * arma::mean(wgt % arma::square(xmat.col(j)));
    return coorddesc::SurrogateGradient { gradient, lipschitz };
  }

  double UpdateSlope (const arma::uword j, const double stepsize, const double gradient,
                      std::false_type /* is_adaptive */) {
    const double dir = stepsize * state_.coefs.beta[j] - gradient;
    return nsoptim::SoftThreshold(dir, penalty_->lambda() * penalty_->alpha()) /
      (stepsize + penalty_->lambda() * (1 - penalty_->alpha()));
  }

  double UpdateSlope (const arma::uword j, const double stepsize, const double gradient,
                      std::true_type /* is_adaptive */) {
    const double dir = stepsize * state_.coefs.beta[j] - gradient;
    const double penalty_level = penalty_->loadings()[j] * penalty_->lambda();
    return nsoptim::SoftThreshold(dir, penalty_level * penalty_->alpha()) /
      (stepsize + penalty_level * (1 - penalty_->alpha()));
  }

  void ResetState (const Coefficients &coefs) {
    if (!loss_) {
      throw std::logic_error("no loss set");
    }
    if (!penalty_) {
      throw std::logic_error("no penalty set");
    }
    state_ = { coefs, loss_->Residuals(coefs), 0, 0, penalty_->Evaluate(coefs) };
    auto loss_eval = loss_->EvaluateResiduals(state_.residuals);
    state_.mscale = loss_eval.scale;
    state_.objf_loss = loss_eval.loss;
  }

  double PenaltyContribution(const double value, const arma::uword, std::false_type /* is_adaptive */) {
    return penalty_->lambda() * (penalty_->alpha() * std::abs(value) +
      0.5 * (1 - penalty_->alpha()) * value * value);
  }

  double PenaltyContribution(const double value, const arma::uword j, std::true_type /* is_adaptive */) {
    return penalty_->loadings()[j] * penalty_->lambda() *
      (penalty_->alpha() * std::abs(value) +
        0.5 * (1 - penalty_->alpha()) * value * value);
  }

  LossFunctionPtr loss_;
  PenaltyPtr penalty_;
  CDPenseConfiguration config_;
  arma::vec lipschitz_bounds_;
  double lipschitz_bound_intercept_;
  coorddesc::State<Coefficients> state_;
  double convergence_tolerance_ = kDefaultConvergenceTolerance;
};
} // namespace pense

#endif // PENSE_CD_PENSE_HPP_
