//
//  optim_dal.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2018-11-30.
//  Copyright © 2018 David Kepplinger. All rights reserved.
//

#ifndef NSOPTIM_OPTIMIZER_DAL_HPP_
#define NSOPTIM_OPTIMIZER_DAL_HPP_

#include <exception>
#include <type_traits>
#include <memory>
#include <algorithm>

#include "../armadillo.hpp"
#include "../traits/traits.hpp"
#include "optimizer_base.hpp"
#include "optimum.hpp"
#include "../objective/ls_regression_loss.hpp"
#include "../container/metrics.hpp"
#include "../container/data.hpp"
#include "../container/regression_coefficients.hpp"
#include "dal_helper.hpp"
#include "soft_threshold.hpp"

namespace nsoptim {

//! Configuration options for the DAL algorithm.
struct DalEnConfiguration {
  //! Maximum number of iterations allowed.
  int max_it;
  //! Maximum number of inner iterations allowed.
  int max_inner_it;
  //! Conservative setting for the numerator when computing the initial value of the proximity parameters.
  double eta_start_numerator_conservative;
  //! Aggressive setting for the numerator when computing the initial value of the proximity parameters.
  double eta_start_numerator_aggressive;
  //! Maximum relative change in `lambda` that allows the use of the aggressive numerator.
  double lambda_relchange_aggressive;
  //! Multiplier to scale the proximity parameters at each outer iteration.
  double eta_multiplier;
};

namespace _optim_dal_internal {
constexpr DalEnConfiguration kDefaultDalEnConfiguration = {100, 100, 0.01, 1, 0.25, 2};
}  // namespace _optim_dal_internal

template <class LossFunction, class PenaltyFunction>
class DalEnOptimizer : public Optimizer<LossFunction, PenaltyFunction, RegressionCoefficients<arma::sp_vec>> {
  using Base = Optimizer<LossFunction, PenaltyFunction, RegressionCoefficients<arma::sp_vec>>;
  using LossFunctionPtr = std::unique_ptr<LossFunction>;
  using PenaltyFunctionPtr = std::unique_ptr<PenaltyFunction>;

  static_assert(traits::is_en_penalty<PenaltyFunction>::value, "PenaltyFunction must be an EN-type penalty.");
  static_assert(traits::is_ls_regression_loss<LossFunction>::value, "LossFunction must be an least-squares-type loss.");

 public:
  using Coefficients = typename Base::Coefficients;
  using Optimum = typename Base::Optimum;
  using Configuration = DalEnConfiguration;

  //! Initialize the DAL optimizer for EN penalties without initializing the loss/penalty functions.
  //!
  //! @param config configuration for the algorithm.
  explicit DalEnOptimizer(const Configuration& config = _optim_dal_internal::kDefaultDalEnConfiguration) noexcept
      : config_(config) {}

  //! Initialize the DAL optimizer for EN penalties.
  //!
  //! @param loss the loss function to optimize over.
  //! @param penalty the penalty function to optimizer over.
  //! @param config configuration for the algorithm.
  DalEnOptimizer(const LossFunction& loss, const PenaltyFunction& penalty,
                 const Configuration& config = _optim_dal_internal::kDefaultDalEnConfiguration) noexcept
    : config_(config), loss_(LossFunctionPtr(new LossFunction(loss))),
      penalty_(PenaltyFunctionPtr(new PenaltyFunction(penalty))), data_(loss_.get()) {}

  //! Copy constructor.
  //!
  //! This creates shallow copies of the loss and penalty functions and a deep copy of the hessian as well as possibly
  //! the weighted data.
  DalEnOptimizer(const DalEnOptimizer& other) noexcept :
    config_(other.config_),
    loss_(other.loss_ ? LossFunctionPtr(new LossFunction(*other.loss_)) : nullptr),
    penalty_(other.penalty_ ? PenaltyFunctionPtr(new PenaltyFunction(*other.penalty_)) : nullptr),
    coefs_(other.coefs_),
    data_(loss_.get()), eta_(other.eta_),
    convergence_tolerance_(other.convergence_tolerance_) {}

  DalEnOptimizer(DalEnOptimizer&& other) = default;
  DalEnOptimizer& operator=(DalEnOptimizer&& other) = default;

  //! Reset the optimizier. This compeletely purges the current *state*.
  void Reset() {
    eta_.nxlambda = -1;
    coefs_.Reset();
  }

  //! Access the loss function.
  //!
  //! @return the loss function currently in use by the optimizer.
  LossFunction& loss() const {
    if (!loss_) {
      throw std::logic_error("no loss set");
    }
    return *loss_;
  }

  //! Set the loss function.
  //!
  //! @param loss the new loss function to optimize over.
  void loss(const LossFunction& loss) noexcept {
    // Check if the dimensionality of the problem changed. If so, reset the coefficients.
    if (loss_ && (loss_->data().n_pred() != loss.data().n_pred())) {
      coefs_.Reset();
    }

    auto changes = data_.Update(loss);
    loss_.reset(new LossFunction(loss));

    if (changes.data_changed || changes.weights_changed > 1) {
      // If the data changed, the proximity parameters must be reset.
      eta_.nxlambda = -1;
    }
  }

  //! Access the penalty function.
  //!
  //! @return the penalty function currently in use by the optimizer.
  PenaltyFunction& penalty() const {
    if (!penalty_) {
      throw std::logic_error("no penalty set");
    }
    return *penalty_;
  }

  //! Set the penalty function.
  //!
  //! @param penalty the new penalty function to optimize over.
  void penalty(const PenaltyFunction& penalty) noexcept {
    penalty_.reset(new PenaltyFunction(penalty));
  }

  //! Get the convergence tolerance for the DAL algorithm.
  //!
  //! @return convergence tolerance.
  double convergence_tolerance() const noexcept {
    return convergence_tolerance_;
  }

  //! Set the convergence tolerance for the DAL algorithm.
  //!
  //! @param convergence_tolerance convergene tolerance for the MM algorithm.
  void convergence_tolerance(double convergence_tolerance) noexcept {
    convergence_tolerance_ = convergence_tolerance;
  }

  //! Find the minimum of the objective function, using the given coefficients as starting point.
  //!
  //! @param start where to start the optimization from.
  //! @return information about the optimum.
  Optimum Optimize(const Coefficients& start) {
    coefs_ = start;
    // If we start from a new place, reset the proximity parameters.
    eta_.nxlambda = -1;
    return Optimize(config_.max_it);
  }

  //! Find the minimum of the objective function, using the given coefficients as starting point
  //! and at most ``max_it`` iterations.
  //!
  //! @param start where to start the optimization from.
  //! @param max_it maximum number of iterations.
  //! @return information about the optimum.
  Optimum Optimize(const Coefficients& start, const int max_it) {
    coefs_ = start;
    // If we start from a new place, reset the proximity parameters.
    eta_.nxlambda = -1;
    return Optimize(max_it);
  }

  //! Find the minimum of the objective function, using the previous solution (or the 0-vector if no
  //! previous solution exists) as starting point.
  //!
  //! @return information about the optimum.
  Optimum Optimize() {
    return Optimize(config_.max_it);
  }

  //! Find the minimum of the objective function, using the previous solution (or the 0-vector if no
  //! previous solution exists) as starting point and at most ``max_it`` iterations.
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

    const double alpha = penalty_->alpha();
    const double nxlambda = data_->n_obs() * penalty_->lambda() / data_.mean_weight();
    const auto softthr_cutoff = SoftthresholdCutoff(nxlambda * alpha, IsAdaptiveTag{});
    const auto update_denom_mult = OneOverShiftPlusC(0., nxlambda * (1 - alpha), IsAdaptiveTag{});

    // If the coefficients are not yet initialized, reset them to the 0-vector.
    if (coefs_.beta.n_elem == 0) {
      coefs_ = loss_->template ZeroCoefficients<Coefficients>();
    }

    // Create new optional metrics
    auto metrics = std::make_unique<Metrics>("dal-ls_en");

    // Initialize the proximity parameters.
    double eta_start_numerator = config_.eta_start_numerator_conservative;
    if (eta_.nxlambda > 0) {
      // If lambda (over the mean observation weights) is close enough to the previous lambda,
      // we can start with larger proximity parameters and thus speed up convergence.
      const double lambda_rel_change = std::abs(1. - eta_.nxlambda / nxlambda);
      if (lambda_rel_change < config_.lambda_relchange_aggressive) {
        eta_start_numerator = config_.eta_start_numerator_aggressive;
      }
    }

    eta_.slope = std::min(eta_start_numerator / nxlambda, _optim_dal_internal::kMaximumEtaStart);
    eta_.intercept = eta_.slope;
    eta_.nxlambda = nxlambda;

    metrics->AddDetail("eta_start_numerator", eta_start_numerator);
    metrics->AddDetail("convergence_tol", convergence_tolerance_);

    // @TODO: maybe `phi_argmin` could be re-used from the previous run if nothing changed!
    arma::vec phi_argmin = data_->cy() - data_->cx() * coefs_.beta;

    if (loss_->IncludeIntercept()) {
      coefs_.intercept = InitializeIntercept(&phi_argmin, HasWeightsTag{});
    } else {
      coefs_.intercept = 0;
    }

    arma::vec dual_constraint_rhs = data_->cx().t() * phi_argmin;
    int outer_iterations = 0;
    double dual_fun_value = 0;
    double eq_constr_violation = 0;

    while (true) {
      Metrics& iteration_metrics = metrics->CreateSubMetrics("dal_iteration");
      // Check if relative duality gap (rdg) is below the threshold
      const double dual_fun_val_prev = dual_fun_value;
      arma::vec dual_vec = ComputeDualVector(phi_argmin, HasWeightsTag{});
      arma::vec xtr_dual = data_->cx().t() * dual_vec;
      if (alpha < 1) {
        SoftThreshold(softthr_cutoff, &xtr_dual);
        dual_fun_value = _optim_dal_internal::DualLoss(dual_vec, data_->cy()) +
          0.5 * arma::dot(linalg::ElementwiseProduct(update_denom_mult, xtr_dual), xtr_dual);
      } else {
        dual_vec *= DualVectorUpdate(nxlambda, xtr_dual, IsAdaptiveTag{});
        dual_fun_value = _optim_dal_internal::DualLoss(dual_vec, data_->cy());
      }

      // Ensure that the value of the dual function does not increase.
      if (outer_iterations > 0 && dual_fun_value > dual_fun_val_prev) {
        dual_fun_value = dual_fun_val_prev;
      }

      const arma::vec residuals = loss_->data().cy() - loss_->data().cx() * coefs_.beta - coefs_.intercept;
      const double orig_objf_value = loss_->Evaluate(residuals) + penalty_->Evaluate(coefs_);
      const double primal_fun_val = data_->n_obs() * orig_objf_value / data_.mean_weight();
      const double rel_duality_gap = (primal_fun_val + dual_fun_value) / primal_fun_val;

      // Record some metrics for the current iteration
      iteration_metrics.AddDetail("fval", primal_fun_val);
      iteration_metrics.AddDetail("dval", dual_fun_value);
      iteration_metrics.AddDetail("rdg", rel_duality_gap);
      iteration_metrics.AddDetail("eta(slope)", eta_.slope);
      iteration_metrics.AddDetail("eta(intercept)", eta_.intercept);

      // Check for convergence or reaching the maximum nr. of iterations.
      if (std::abs(rel_duality_gap) < convergence_tolerance_) {
        metrics->AddMetric("iter", outer_iterations);
        return MakeOptimum(*loss_, *penalty_, coefs_, residuals, orig_objf_value, std::move(metrics));
      } else if (++outer_iterations > max_it) {
        metrics->AddMetric("iter", outer_iterations);
        metrics->AddDetail("warning", "Algorithm did not converge");
        return MakeOptimum(*loss_, *penalty_, coefs_, residuals, orig_objf_value, std::move(metrics),
                           OptimumStatus::kWarning, "Algorithm did not convergence");
      } else if (rel_duality_gap < -convergence_tolerance_) {
        metrics->AddMetric("iter", outer_iterations);
        metrics->AddDetail("warning", "Relative duality gap is negative");
        return MakeOptimum(*loss_, *penalty_, coefs_, residuals, orig_objf_value, std::move(metrics),
                           OptimumStatus::kError, "Relative duality gap is negative");
      }

      // Minimize phi
      const int phi_iter = MinimizePhi(softthr_cutoff, nxlambda, alpha, &phi_argmin, &dual_constraint_rhs,
                                       &iteration_metrics);

      if (phi_iter == _optim_dal_internal::kMinimizePhiHessianSingular) {
        metrics->AddMetric("iter", outer_iterations);
        metrics->AddDetail("warning", "Hessian matrix is singular");
        return MakeOptimum(*loss_, *penalty_, coefs_, std::move(metrics),
                           OptimumStatus::kError, "Hessian matrix singular");
      } else if (phi_iter == _optim_dal_internal::kMinimizePhiGradientTooSteep) {
        metrics->AddMetric("iter", outer_iterations);
        metrics->AddDetail("warning", "Gradient is too steep");
        return MakeOptimum(*loss_, *penalty_, coefs_, std::move(metrics),
                           OptimumStatus::kError, "Gradient is too steep");
      }

      iteration_metrics.AddDetail("phi_iter", phi_iter);

      // Update `eta` parameters
      eta_.slope *= config_.eta_multiplier;
      // Update `eta (intercept)`
      const double new_eq_constr_violation = EqualityConstraintViolation(phi_argmin, HasWeightsTag{});
      if ((outer_iterations > 1) && (new_eq_constr_violation > convergence_tolerance_) &&
          (new_eq_constr_violation > 0.5 * eq_constr_violation)) {
        eta_.intercept *= 10 * config_.eta_multiplier;
      } else {
        eta_.intercept *= config_.eta_multiplier;
      }
      eq_constr_violation = new_eq_constr_violation;
    }
  }

 private:
  using HasWeightsTag = typename traits::is_weighted<LossFunction>::type;
  using IsAdaptiveTag = typename traits::is_adaptive<PenaltyFunction>::type;
  using PenaltyLevelType = typename std::conditional<traits::is_adaptive<PenaltyFunction>::value,
                                                     arma::vec, double>::type;

  //! Minimize the AL function `phi`.
  //!
  //! @param softthr_cutoff the cutoff value(s) for the soft-thresholding function.
  //! @param nxlambda the regularization parameter (scaled by `n` over the mean of the weights).
  //! @param alpha alpha of the penalty term.
  //! @param phi_argmin a pointer to the initial guess of the value of the minimizing argument. On output, the value
  //!                    of this vector is updated with the minimizing argument.
  //! @param dual_constraint_rhs a pointer to the value of the right-hand-side of the current duality constraint.
  //!                            On output, the value of this of this vector is updated with the duality constraint at
  //!                            `phi_argmin`.
  //! @param metrics a pointer to the metrics for the current iteration of DAL.
  //! @return the number of iterations required to find the minimum of `phi`. If < 0, an error occured.
  int MinimizePhi(PenaltyLevelType softthr_cutoff, const double nxlambda, const double alpha,
                  arma::vec * const phi_argmin, arma::vec * const dual_constraint_rhs, Metrics * const metrics) {
    const auto moreau_factor = OneOverShiftPlusC(1., nxlambda * eta_.slope * (1 - alpha), IsAdaptiveTag{});
    const Coefficients prev_coefs(coefs_);
    softthr_cutoff *= eta_.slope;

    int iter = 0;
    arma::vec step_dir;
    while (true) {
      Metrics& inner_metrics = metrics->CreateSubMetrics("phi_iteration");
      // Update coefficient values
      coefs_.beta = SoftThreshold(prev_coefs.beta, eta_.slope, *dual_constraint_rhs, softthr_cutoff);
      coefs_.intercept = prev_coefs.intercept + ComputeInterceptChange(*phi_argmin, HasWeightsTag{});

      // Evaluate the phi function and its gradient.
      const double phi_fun_value = _optim_dal_internal::DualLoss(*phi_argmin, data_->cy()) +
        0.5 * arma::dot(linalg::ElementwiseProduct(coefs_.beta, moreau_factor), coefs_.beta) / eta_.slope +
        0.5 * coefs_.intercept * coefs_.intercept / eta_.intercept;

      // Check if we reached the maximum number of inner iterations.
      if (++iter > config_.max_inner_it) {
        break;
      }

      // Determine gradient and step direction.
      const arma::vec gradient = EvaluatePhiGradient(*phi_argmin, moreau_factor, HasWeightsTag{});
      const double gradient_norm_sq = arma::dot(gradient, gradient);
      const double intercept_diff = prev_coefs.intercept - coefs_.intercept;
      const arma::sp_vec beta_diff = prev_coefs.beta - linalg::ElementwiseProduct(moreau_factor, coefs_.beta);
      const double conv_threshold = std::max(0.5 * convergence_tolerance_,
        arma::dot(beta_diff, beta_diff) / eta_.slope + (intercept_diff * intercept_diff) / eta_.intercept);

      inner_metrics.AddDetail("gradient_norm_sq", gradient_norm_sq);
      inner_metrics.AddDetail("threshold", conv_threshold);
      if (iter > 1 && gradient_norm_sq < conv_threshold) {
        break;
      }

      const bool hessian_inverted = PhiStepDir(coefs_.beta, gradient, moreau_factor, &step_dir);
      if (!hessian_inverted) {
        iter = _optim_dal_internal::kMinimizePhiHessianSingular;
        break;
      }

      // Backtracking line search for step size to update `a`
      const double linesearch_step = arma::dot(step_dir, gradient);
      const double intercept_step_base = prev_coefs.intercept + ComputeInterceptChange(*phi_argmin, HasWeightsTag{});
      const double intercept_step_decrease = ComputeInterceptChange(step_dir, HasWeightsTag{});
      const arma::vec dual_constraint_rhs_step_dir = data_->cx().t() * step_dir;
      double step_size = 1;
      int line_search_iter = 0;

      *phi_argmin -= step_dir;
      *dual_constraint_rhs -= dual_constraint_rhs_step_dir;

      while (++line_search_iter <= _optim_dal_internal::kLinesearchMaxSteps) {
        // Update coefficient values
        coefs_.beta = SoftThreshold(prev_coefs.beta, eta_.slope, *dual_constraint_rhs, softthr_cutoff);
        coefs_.intercept = intercept_step_base - step_size * intercept_step_decrease;

        // Evaluate the phi function and its gradient.
        const double phi_fun_value_step = _optim_dal_internal::DualLoss(*phi_argmin, data_->cy()) +
          0.5 * arma::dot(linalg::ElementwiseProduct(coefs_.beta, moreau_factor), coefs_.beta) / eta_.slope +
          0.5 * coefs_.intercept * coefs_.intercept / eta_.intercept;

        const double phi_fun_value_target = phi_fun_value -
          _optim_dal_internal::kLinesearchStepCondition * step_size * linesearch_step;
        if (phi_fun_value_step < phi_fun_value_target) {
          break;
        }

        const double prev_step_size = step_size;
        step_size *= _optim_dal_internal::kLinesearchStepsizeMultiplier;
        *phi_argmin += (prev_step_size - step_size) * step_dir;
        *dual_constraint_rhs += (prev_step_size - step_size) * dual_constraint_rhs_step_dir;
      }
      inner_metrics.AddDetail("step_size", step_size);
      inner_metrics.AddDetail("linesearch_steps", line_search_iter);

      if (line_search_iter > _optim_dal_internal::kLinesearchMaxSteps) {
        iter = _optim_dal_internal::kMinimizePhiGradientTooSteep;
        break;
      }
    }
    linalg::InplaceElementwiseProduct(moreau_factor, &coefs_.beta);
    return iter;
  }

  //! Compute the initial intercept and update the `a` vector accordingly for the weighted LS loss.
  double InitializeIntercept(arma::vec* phi_argmin, std::true_type) const {
    const double new_intercept = arma::mean(data_.sqrt_weights() % *phi_argmin);
    *phi_argmin -= data_.sqrt_weights() * new_intercept;
    return new_intercept;
  }

  //! Compute the initial intercept and update the `a` vector accordingly for the un-weighted LS loss.
  double InitializeIntercept(arma::vec* phi_argmin, std::false_type) const {
    const double new_intercept = arma::mean(*phi_argmin);
    *phi_argmin -= new_intercept;
    return new_intercept;
  }

  //! Compute the new intercept term for a weighted LS loss.
  double ComputeInterceptChange(const arma::vec& phi_argmin, std::true_type) const {
    if (loss_->IncludeIntercept()) {
      return eta_.intercept * arma::dot(data_.sqrt_weights(), phi_argmin);
    }
    return 0.;
  }

  //! Compute the new intercept term for an un-weighted LS loss.
  double ComputeInterceptChange(const arma::vec& phi_argmin, std::false_type) const {
    if (loss_->IncludeIntercept()) {
      return eta_.intercept * arma::accu(phi_argmin);
    }
    return 0.;
  }

  //! Compute the dual vector for the weighted LS loss.
  arma::vec ComputeDualVector(const arma::vec& phi_argmin, std::true_type) const {
    return phi_argmin - data_.sqrt_weights() * arma::mean(data_.sqrt_weights() % phi_argmin);
  }

  //! Compute the dual vector for the un-weighted LS loss.
  arma::vec ComputeDualVector(const arma::vec& phi_argmin, std::false_type) const {
    return phi_argmin - arma::mean(phi_argmin);
  }

  //! Compute the update to the dual vector for adaptive penalties.
  double DualVectorUpdate(const double nxlambda, const arma::vec& xtr_dual, std::true_type) {
    return nxlambda / arma::norm(xtr_dual / penalty_->loadings(), "inf");
  }

  //! Compute the update to the dual vector for non-adaptive penalties.
  double DualVectorUpdate(const double nxlambda, const arma::vec& xtr_dual, std::false_type) {
    return std::min(nxlambda / arma::norm(xtr_dual, "inf"), 1.);
  }

  //! Compute the norm of the dual vector for weighted LS.
  double EqualityConstraintViolation(const arma::vec& phi_argmin, std::true_type) {
    return arma::dot(phi_argmin, data_.sqrt_weights());
  }

  //! Compute the norm of the dual vector for un-weighted LS.
  double EqualityConstraintViolation(const arma::vec& phi_argmin, std::false_type) {
    return arma::accu(phi_argmin);
  }

  //! Evaluate the gradient of the `phi` function at (a, coefs) if the LS loss includes weights.
  arma::vec EvaluatePhiGradient(const arma::vec& phi_argmin, const PenaltyLevelType& moreau_factor,
                                std::true_type) const {
    if (loss_->IncludeIntercept()) {
      return phi_argmin - data_->cy() + data_->cx() * linalg::ElementwiseProduct(moreau_factor, coefs_.beta) +
        coefs_.intercept * data_.sqrt_weights();
    } else {
      return phi_argmin - data_->cy() + data_->cx() * linalg::ElementwiseProduct(moreau_factor, coefs_.beta);
    }
  }

  //! Evaluate the gradient of the `phi` function at (a, coefs) if the LS loss does not include weights.
  arma::vec EvaluatePhiGradient(const arma::vec& phi_argmin, const PenaltyLevelType& moreau_factor,
                                std::false_type) const {
    if (loss_->IncludeIntercept()) {
      return phi_argmin - data_->cy() + data_->cx() * linalg::ElementwiseProduct(moreau_factor, coefs_.beta) +
        coefs_.intercept;
    } else {
      return phi_argmin - data_->cy() + data_->cx() * linalg::ElementwiseProduct(moreau_factor, coefs_.beta);
    }
  }

  //! Get the soft-thresholding cutoff for adaptive penalties.
  arma::vec SoftthresholdCutoff(const double c, std::true_type) const noexcept {
    return c * penalty_->loadings();
  }

  //! Get the soft-thresholding cutoff for non-adaptive penalties.
  double SoftthresholdCutoff(const double c, std::false_type) const noexcept {
    return c;
  }

  //! Get the Moreau factor for adaptive penalties.
  arma::vec OneOverShiftPlusC(const double shift, const double c, std::true_type) const noexcept {
    return 1 / (shift + c * penalty_->loadings());
  }

  //! Get the Moreau factor for non-adaptive penalties.
  double OneOverShiftPlusC(const double shift, const double c, std::false_type) const noexcept {
    return 1 / (shift + c);
  }

  //! Find the appropriate step direction.
  //! @return boolean indicating if a solution was found or not.
  bool PhiStepDir(const arma::sp_vec& beta, const arma::vec& gradient, const PenaltyLevelType& moreau_factor,
                  arma::vec * const step_dir) const {
    if (beta.n_nonzero == 0) {
      // No predictors are active. We can compute the inverse of the Hessian and the step direction directly.
      if (loss_->IncludeIntercept()) {
        *step_dir = PhiStepDirNoPredictors(gradient, HasWeightsTag{});
      } else {
        // If there's no intercept and no active predictors, the Hessian is the identity matrix.
        *step_dir = gradient;
      }
      return true;
    }

    // Ensure that the row indicies of the non-zero elements are upated.
    beta.sync();
    const arma::uvec active(const_cast<arma::uword*>(beta.row_indices), beta.n_nonzero, false, true);

    arma::mat hessian = PhiHessian(active, moreau_factor, HasWeightsTag{});
    hessian.diag() += 1;  //< this is faster than adding a diagonal matrix.
    return arma::solve(*step_dir, hessian, gradient, arma::solve_opts::likely_sympd);
  }

  //! Compute the Hessian for a weighted LS loss and non-adaptive penalty.
  arma::mat PhiHessian(const arma::uvec& active_predictors, const double moreau_factor,
                       std::true_type) const {
    auto active_view = data_->cx().cols(active_predictors);
    if (loss_->IncludeIntercept()) {
      return eta_.slope * moreau_factor * active_view * active_view.t() + eta_.intercept * data_.sqrt_weights_outer();
    } else {
      return eta_.slope * moreau_factor * active_view * active_view.t();
    }
  }

  //! Compute the Hessian for an un-weighted LS loss and non-adaptive penalty.
  arma::mat PhiHessian(const arma::uvec& active_predictors, const double moreau_factor,
                       std::false_type) const {
    auto active_view = data_->cx().cols(active_predictors);
    if (loss_->IncludeIntercept()) {
      return eta_.slope * moreau_factor * active_view * active_view.t() + eta_.intercept;
    } else {
      return eta_.slope * moreau_factor * active_view * active_view.t();
    }
  }

  //! Compute the Hessian for a weighted LS loss and adaptive penalty.
  arma::mat PhiHessian(const arma::uvec& active_predictors, const arma::vec& moreau_factor,
                       std::true_type) const {
    const arma::mat active_view = data_->cx().cols(active_predictors);
    if (loss_->IncludeIntercept()) {
      return eta_.slope * (active_view.each_row() % moreau_factor.rows(active_predictors).t()) * active_view.t() +
        eta_.intercept * data_.sqrt_weights_outer();
    } else {
      return eta_.slope * (active_view.each_row() % moreau_factor.rows(active_predictors).t()) * active_view.t();
    }
  }

  //! Compute the Hessian for an un-weighted LS loss and adaptive penalty.
  arma::mat PhiHessian(const arma::uvec& active_predictors, const arma::vec& moreau_factor,
                       std::false_type) const {
    const arma::mat active_view = data_->cx().cols(active_predictors);
    if (loss_->IncludeIntercept()) {
      return eta_.slope * (active_view.each_row() % moreau_factor.rows(active_predictors).t()) * active_view.t() +
        eta_.intercept;
    } else {
      return eta_.slope * (active_view.each_row() % moreau_factor.rows(active_predictors).t()) * active_view.t();
    }
  }

  //! Find the appropriate step direction if no predictors are active and a weighted LS loss is used.
  arma::vec PhiStepDirNoPredictors(const arma::vec& gradient, std::true_type) const {
    // If there's an intercept and weights, the preconditioner is set to the inverse of the Hessian.
    const double sum_weights = data_->n_obs() * data_.mean_weight();
    arma::mat inv_hessian = data_.sqrt_weights_outer() * (-eta_.intercept / (1 + eta_.intercept * sum_weights));
    inv_hessian.diag() += 1.;
    // The preconditioner is now the actual inverse itself and can be used directly to compute the step direction.
    return inv_hessian * gradient;
  }

  //! Find the appropriate step direction if no predictors are active and a non-weighted LS loss is used.
  arma::vec PhiStepDirNoPredictors(const arma::vec& gradient, std::false_type) const noexcept {
    return gradient / (1 + eta_.intercept);
  }

  const Configuration config_;
  LossFunctionPtr loss_;
  PenaltyFunctionPtr penalty_;
  Coefficients coefs_;
  _optim_dal_internal::DataProxy<LossFunction> data_;
  _optim_dal_internal::ProximityParameters eta_;
  double convergence_tolerance_ = 1e-8;
};
}  // namespace nsoptim

#endif  // NSOPTIM_OPTIMIZER_DAL_HPP_
