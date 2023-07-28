//
//  admm.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2019-07-10.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifndef NSOPTIM_OPTIMIZER_ADMM_HPP_
#define NSOPTIM_OPTIMIZER_ADMM_HPP_

#include <string>
#include <exception>
#include <type_traits>
#include <algorithm>
#include <limits>

#include "../armadillo.hpp"
#include "../container/regression_coefficients.hpp"
#include "../container/data.hpp"
#include "optimizer_base.hpp"
#include "optimum.hpp"
#include "../objective/ls_regression_loss.hpp"
#include "../objective/en_penalty.hpp"
#include "../traits/traits.hpp"
#include "soft_threshold.hpp"
#include "linear_algebra_utilities.hpp"

namespace nsoptim {

//! Configuration options for the variable-stepsize ADMM algorithm.
struct AdmmVarStepConfiguration {
  int max_it;  //!< Maximum number of iterations allowed.
  double tau;  //!< Initial step size. If negative (the default), use the square L_2 norm of `x`.
  double tau_lower_mult;  //!< Lower bound for the step size, defined as a multiple of `tau`.
  double tau_adjustment_lower;  //!< Lower bound of the step-size adjustment factor.
  double tau_adjustment_upper;  //!< Upper bound of the step-size adjustment factor.
};

//! Configuration options for the linearized ADMM algorithm.
struct AdmmLinearConfiguration {
  int max_it;  //!< maximum number of iterations allowed.
  double accelerate;  //!< acceleration factor.
};

namespace admm_optimizer {
//! Default configuration for the variable-stepsize ADMM algorithm
constexpr AdmmVarStepConfiguration kDefaultVarStepConfig { 1000, -1, 0.01, 0.98, 0.999 };
//! Default configuration for the variable-stepsize ADMM algorithm
constexpr AdmmLinearConfiguration kDefaultLinConfig { 1000, 1. };

//! How often does the secondary convergence criterion (relative to the maximum number of iterations) need to be
//! fulfulled to stop the linearized ADMM early.
constexpr int kSecondCriterionMultiplier = 10;

// Data cache for inner-products that don't change (often)
struct DataCache {
  arma::mat xtx;
  arma::vec xty;
  arma::vec xtwgt;
  double chol_xtx_tau;
  arma::mat chol_xtx;
};

//! Check whether any of the predictors in `x` violates the KKT conditions for a EN-type problem.
//!
//! @param x matrix of predictor values.
//! @param residuals vector of residuals.
//! @param lambda the lambda value to use (overrides the lambda in the EN penalty)
//! @param penalty elastic net penalty object.
inline bool AnyViolateKKT(const arma::mat&, const arma::vec&, const double, const RidgePenalty&) {
  return true;
}

//! Check whether any of the predictors in `x` violates the KKT conditions for a EN-type problem.
//!
//! @param x matrix of predictor values.
//! @param residuals vector of residuals.
//! @param lambda the lambda value to use (overrides the lambda in the EN penalty)
//! @param penalty elastic net penalty object.
inline bool AnyViolateKKT(const arma::mat& x, const arma::vec& residuals, const double lambda,
                          const EnPenalty& penalty) {
  // const double lambda_1 = residuals.n_elem * penalty.alpha() * penalty.lambda();
  const double lambda_1 = penalty.alpha() * lambda * residuals.n_elem;
  for (arma::uword j = 0; j < x.n_cols; ++j) {
    // const double cutoff = lambda1 * penalty.loadings()[j];
    const double inner = arma::dot(x.col(j), residuals);
    if ((inner < -lambda_1) || (inner > lambda_1)) {
      return true;
    }
  }
  return false;
}

//! Check whether any of the predictors in `x` violates the KKT conditions for a adaptive EN-type problem.
//!
//! @param x matrix of predictor values.
//! @param residuals vector of residuals.
//! @param lambda the lambda value to use (overrides the lambda in the EN penalty)
//! @param penalty adaptive elastic net penalty object.
inline bool AnyViolateKKT(const arma::mat& x, const arma::vec& residuals, const double lambda,
                          const AdaptiveEnPenalty& penalty) {
  // const double lambda_1 = residuals.n_elem * penalty.alpha() * penalty.lambda();
  const double lambda_1 = penalty.alpha() * lambda * residuals.n_elem;
  for (arma::uword j = 0; j < x.n_cols; ++j) {
    const double cutoff = lambda_1 * penalty.loadings()[j];
    const double inner = arma::dot(x.col(j), residuals);
    if ((inner < -cutoff) || (inner > cutoff)) {
      return true;
    }
  }
  return false;
}

//! Check if all entries in the vector are zero.
//! For dense vectors, always return false.
inline constexpr bool AllZero(const arma::vec&) noexcept {
  return false;
}

//! Check if all entries in the vector are zero.
//! For sparse vectors, check the n_nonzero field.
inline constexpr bool AllZero(const arma::sp_vec& v) noexcept {
  return v.n_nonzero == 0;
}

}  // namespace admm_optimizer

//! Proximal operator for the unweighted LS loss.
class LsProximalOperator {
 public:
  using LossFunction = LsRegressionLoss;

  //! Initialize the proximal operator for the unweighted LS loss.
  //!
  //! Deng & Yin (2016) provide a good reference on how to choose the step size (`beta` in their paper).
  //!
  //! @param step_size Fixed step size. If negative (the default), compute a good step size based on the data.
  //! @see Deng, W. and Yin, W. (2016). On the global and linear convergence of the generalized alternating direction
  //!        method of multipliers. Journal of Scientific Computing, 66(3):889–916.
  explicit LsProximalOperator(const double step_size = -1) noexcept : config_step_size_(step_size) {}

  //! Set the loss function for the proximal operator.
  //!
  //! @param loss the LS-loss for optimization. The object retains only a reference to the loss, so it is
  //!             the user's responsibility to not use the object after the loss is removed!
  inline void loss(LsRegressionLoss* loss) noexcept {
    loss_ = loss;
  }

  //! Evaluate the proximal operator at `u`.
  //!
  //! @param u Where to evaluate the proximal operator.
  //! @param intercept Current value of the intercept.
  //! @param scale Scaling factor for the proximal operator.
  //! @param metrics Optional metrics object to collect metrics of the proximal operator.
  inline arma::vec operator()(const arma::vec& u, const double intercept, const double scale,
                              Metrics * const = nullptr) const {
    const double mult_fact = 1 / (1 + scale);
    if (loss_->IncludeIntercept()) {
      return mult_fact * u + mult_fact * scale * loss_->data().cy() + mult_fact * intercept;
    } else {
      return mult_fact * u + mult_fact * scale * loss_->data().cy();
    }
  }

  //! Evaluate the proximal operator at `u`.
  //!
  //! @param u Where to evaluate the proximal operator.
  //! @param prev Previous result of the proximal operator.
  //! @param intercept Current value of the intercept.
  //! @param scale Scaling factor for the proximal operator.
  //! @param metrics Optional metrics object to collect metrics of the proximal operator.
  inline arma::vec operator()(const arma::vec& u, const arma::vec&, const double intercept, const double scale,
                              Metrics * const = nullptr) const {
    return operator()(u, intercept, scale);
  }

  //! Get the proximal operator scaling (i.e., 1 over step size) optimized for the currently set loss.
  //!
  //! @return Loss-specific operator scaling.
  double OperatorScaling() const {
    if (config_step_size_ < 0) {
      return 1;
    }
    return 1 / config_step_size_;
  }

  //! Compute a the scaling constant of the penalty parameter for the currently set loss.
  //!
  //! @return Loss-specific penalty scaling.
  double PenaltyScaling() const {
    return loss_->data().n_obs();
  }

 private:
  double config_step_size_;
  LsRegressionLoss * loss_;
};

//! Proximal operator for the weighted LS loss.
class WeightedLsProximalOperator {
 public:
  using LossFunction = WeightedLsRegressionLoss;

  //! Initialize the proximal operator for the unweighted LS loss.
  //!
  //! Deng & Yin (2016) provide a good reference on how to choose the step size (`beta` in their paper).
  //!
  //! @param step_size Fixed step size. If negative (the default), compute a good step size based on the data.
  //! @see Deng, W. and Yin, W. (2016). On the global and linear convergence of the generalized alternating direction
  //!        method of multipliers. Journal of Scientific Computing, 66(3):889–916.
  explicit WeightedLsProximalOperator(const double step_size = -1) noexcept
    : automatic_step_size_(step_size < 0), op_scaling_(1 / step_size) {}

  //! Set the loss function for the proximal operator.
  //!
  //! @param loss the LS-loss for optimization. The object retains only a reference to the loss, so it is
  //!             the user's responsibility to not use the object after the loss is removed!
  inline void loss(WeightedLsRegressionLoss* loss) noexcept {
    loss_ = loss;

    // Manually count the number of non-zero weights and find the minimum/maximum non-zero weight
    if (automatic_step_size_) {
      nnz_weights_ = 0;
      double min_nnz_weight = std::numeric_limits<double>::max();
      double max_nnz_weight = 0;
      for (auto&& sqrt_wgt : loss_->sqrt_weights()) {
        if (sqrt_wgt > 0) {
          ++nnz_weights_;
          if (sqrt_wgt < min_nnz_weight) {
            min_nnz_weight = sqrt_wgt;
          }
          if (sqrt_wgt > max_nnz_weight) {
            max_nnz_weight = sqrt_wgt;
          }
        }
      }
      if (nnz_weights_ > 0) {
        op_scaling_ = 1 / (min_nnz_weight * max_nnz_weight);
      }
    }
  }

  //! Evaluate the proximal operator at `u`.
  //!
  //! @param u Where to evaluate the proximal operator.
  //! @param intercept Current value of the intercept.
  //! @param scale Scaling factor for the proximal operator.
  //! @param metrics Optional metrics object to collect metrics of the proximal operator.
  inline arma::vec operator()(const arma::vec& u, const double intercept, const double scale,
                              Metrics * const = nullptr) const {
    const arma::vec scaled_squared_weights = scale * arma::square(loss_->sqrt_weights());
    if (loss_->IncludeIntercept()) {
      return (u + scaled_squared_weights % loss_->data().cy() + intercept) / (1 + scaled_squared_weights);
    } else {
      return (u + scaled_squared_weights % loss_->data().cy()) / (1 + scaled_squared_weights);
    }
  }

  //! Evaluate the proximal operator at `u`.
  //!
  //! @param u Where to evaluate the proximal operator.
  //! @param prev Previous result of the proximal operator.
  //! @param intercept Current value of the intercept.
  //! @param scale Scaling factor for the proximal operator.
  //! @param metrics Optional metrics object to collect metrics of the proximal operator.
  inline arma::vec operator()(const arma::vec& u, const arma::vec&, const double intercept, const double scale,
                              Metrics * const = nullptr) const {
    return operator()(u, intercept, scale);
  }

  //! Get the proximal operator scaling (i.e., 1 over step size) optimized for the currently set loss.
  //!
  //! @return Loss-specific operator scaling.
  double OperatorScaling() const {
    return op_scaling_;
  }

  //! Compute a the scaling constant of the penalty parameter for the currently set loss.
  //!
  //! @return Loss-specific penalty scaling.
  double PenaltyScaling() const {
    return loss_->data().n_obs() / loss_->mean_weight();
  }

 private:
  bool automatic_step_size_;
  double op_scaling_;
  WeightedLsRegressionLoss * loss_;
  int nnz_weights_;
};

namespace admm_optimizer {
//! Type trait mapping the LsRegressionLoss to the LsProximalOperator, the WeightedLsRegressionLoss to the
//! WeightedLsProximalOperator,
//! and any other type to itself.
template <typename T>
using ProximalOperator = typename std::conditional<
  std::is_same<T, LsRegressionLoss>::value, LsProximalOperator,
    typename std::conditional<
      std::is_same<T, WeightedLsRegressionLoss>::value, WeightedLsProximalOperator, T>::type >::type;
}  // namespace admm_optimizer


//! Compute the EN regression estimate using the alternating direction method of multiplier (ADMM)
//! with linearization. This optimizer uses the given proximal operator class `ProxOp`.
template <typename ProxOp, typename PenaltyFunction, typename Coefficients>
class GenericLinearizedAdmmOptimizer : public Optimizer<typename ProxOp::LossFunction, PenaltyFunction, Coefficients> {
 public:
  using ProximalOperator = ProxOp;
  using LossFunction = typename ProxOp::LossFunction;

 private:
  using Base = Optimizer<LossFunction, PenaltyFunction, Coefficients>;
  using LossFunctionPtr = std::unique_ptr<LossFunction>;
  using PenaltyPtr = std::unique_ptr<PenaltyFunction>;
  using IsAdaptiveTag = typename traits::is_adaptive<PenaltyFunction>::type;

  static_assert(traits::is_en_penalty<PenaltyFunction>::value, "PenaltyFunction must be an EN-type penalty.");
  // ADMM state variables
  struct State {
    arma::vec fitted;
    arma::vec lagrangian;
  };

  // Helper-traits to identify constructor arguments.
  template<typename T>
  using IsConfiguration = std::is_same<T, AdmmLinearConfiguration>;
  template<typename T>
  using IsLossFunction = std::is_same<T, LossFunction>;
  template<typename T>
  using IsProximalOperator = std::is_same<T, ProximalOperator>;

 public:
  using Optimum = typename Base::Optimum;

  //! Initialize the ADMM optimizer without setting a loss or penalty function.
  GenericLinearizedAdmmOptimizer() noexcept
      : config_(admm_optimizer::kDefaultLinConfig), prox_(), loss_(nullptr), penalty_(nullptr) {}

  //! Initialize the ADMM optimizer with the given loss and penalty functions.
  //!
  //! @param loss the loss function object.
  //! @param penalty the penalty function object.
  GenericLinearizedAdmmOptimizer(const LossFunction& loss, const PenaltyFunction& penalty) noexcept
      : config_(admm_optimizer::kDefaultLinConfig), prox_(),
        loss_(new LossFunction(loss)), penalty_(new PenaltyFunction(penalty)) {}

  //! Initialize the ADMM optimizer with the given loss and penalty functions.
  //!
  //! @param loss the loss function object.
  //! @param penalty the penalty function object.
  //! @param prox proximal operator object.
  //! @param config optional ADMM configuration object.
  GenericLinearizedAdmmOptimizer(const LossFunction& loss, const PenaltyFunction& penalty, const ProximalOperator& prox,
                                 const AdmmLinearConfiguration& config = admm_optimizer::kDefaultLinConfig) noexcept
      : config_(config), prox_(prox), loss_(new LossFunction(loss)), penalty_(new PenaltyFunction(penalty)) {}

  //! Initialize the ADMM optimizer without setting a loss or penalty function.
  //!
  //! @param config ADMM configuration object.
  //! @param prox_args... arguments to the constructor of the proximal operator.
  template<typename C, typename... Args, typename = typename
           std::enable_if<IsConfiguration<C>::value, void>::type >
  explicit GenericLinearizedAdmmOptimizer(const C& config, Args&&... prox_args) noexcept
      : config_(config), prox_(std::forward<Args>(prox_args)...), loss_(nullptr), penalty_(nullptr) {}

  //! Initialize the ADMM optimizer without setting a loss or penalty function.
  //!
  //! @param config ADMM configuration object.
  //! @param prox_arg_1 first argument to constructor of the proximal operator.
  //! @param prox_args... further arguments to the constructor of the proximal operator.
  template<typename C, typename T1, typename... Args, typename = typename
           std::enable_if<IsConfiguration<C>::value, void>::type >
  GenericLinearizedAdmmOptimizer(const C& config, T1&& prox_arg_1, Args&&... prox_args) noexcept
      : config_(config), prox_(std::forward<T1>(prox_arg_1), std::forward<Args>(prox_args)...),
        loss_(nullptr), penalty_(nullptr) {}

  //! Initialize the ADMM optimizer without setting a loss or penalty function.
  //!
  //! @param config ADMM configuration object.
  //! @param prox_arg_1 first argument to constructor of the proximal operator.
  //! @param prox_arg_2 second argument to constructor of the proximal operator.
  //! @param prox_args... further arguments to the constructor of the proximal operator.
  template<typename C, typename T1, typename T2, typename... Args, typename = typename
           std::enable_if<IsConfiguration<C>::value, void>::type >
  GenericLinearizedAdmmOptimizer(const C& config, T1&& prox_arg_1, T2&& prox_arg_2, Args&&... prox_args) noexcept
      : config_(config),
        prox_(std::forward<T1>(prox_arg_1), std::forward<T2>(prox_arg_2), std::forward<Args>(prox_args)...),
        loss_(nullptr), penalty_(nullptr) {}

  //! Initialize the ADMM optimizer without setting a loss or penalty function.
  //!
  //! @param loss the loss function object.
  //! @param penalty the penalty function object.
  //! @param prox_arg_1 first argument to constructor of the proximal operator.
  //! @param prox_args... further arguments to the constructor of the proximal operator.
  template<typename L, typename P, typename T1, typename... Args, typename = typename
           std::enable_if<IsLossFunction<L>::value && !IsConfiguration<T1>::value, void>::type >
  GenericLinearizedAdmmOptimizer(const L& loss, const P& penalty, T1&& prox_arg_1, Args&&... prox_args) noexcept
      : config_(admm_optimizer::kDefaultLinConfig),
        prox_(std::forward<T1>(prox_arg_1), std::forward<Args>(prox_args)...),
        loss_(new LossFunction(loss)), penalty_(new PenaltyFunction(penalty)) {}

  //! Initialize the ADMM optimizer without setting a loss or penalty function.
  //!
  //! @param loss the loss function object.
  //! @param penalty the penalty function object.
  //! @param config ADMM configuration object.
  //! @param prox_args... further arguments to the constructor of the proximal operator.
  template<typename L, typename P, typename C, typename... Args, typename = typename
           std::enable_if<IsLossFunction<L>::value && IsConfiguration<C>::value, void>::type >
  GenericLinearizedAdmmOptimizer(const L& loss, const P& penalty, const C& config, Args&&... prox_args) noexcept
      : config_(config), prox_(std::forward<Args>(prox_args)...),
        loss_(new LossFunction(loss)), penalty_(new PenaltyFunction(penalty)) {}

  //! Default copy constructor.
  //!
  //! The copied optimizer will share the identical loss and penalty functions after construction.
  //! In case the loss or penalty function are mutated in any way, the change will affect both optimizers.
  //! If the loss/penalty function is changed on one of the optimizers (using the `loss()` or `penalty()` methods),
  //! the two optimizers will *not* share the new loss/penalty function.
  GenericLinearizedAdmmOptimizer(const GenericLinearizedAdmmOptimizer& other) noexcept
      : config_(other.config_),
        prox_(other.prox_),
        loss_(other.loss_? new LossFunction(*other.loss_) : nullptr),
        penalty_(other.penalty_ ? new PenaltyFunction(*other.penalty_) : nullptr),
        coefs_(other.coefs_),
        state_(other.state_),
        x_col_sum_(other.x_col_sum_),
        operator_scaling_g_(other.operator_scaling_g_),
        operator_scaling_f_(other.operator_scaling_f_),
        convergence_tolerance_(other.convergence_tolerance_) {}

  //! Default copy assignment.
  //!
  //! The copied optimizer will share the identical loss and penalty functions after construction.
  //! In case the loss or penalty function are mutated in any way, the change will affect both optimizers.
  //! If the loss/penalty function is changed on one of the optimizers (using the `loss()` or `penalty()` methods),
  //! the two optimizers will *not* share the new loss/penalty function.
  GenericLinearizedAdmmOptimizer& operator=(const GenericLinearizedAdmmOptimizer& other) = default;

  //! Default move constructor.
  GenericLinearizedAdmmOptimizer(GenericLinearizedAdmmOptimizer&& other) = default;

  //! Default move assignment operator.
  GenericLinearizedAdmmOptimizer& operator=(GenericLinearizedAdmmOptimizer&& other) = default;

  ~GenericLinearizedAdmmOptimizer() = default;

  void Reset() {
    state_.lagrangian.reset();
  }

  LossFunction& loss() const {
    if (!loss_) {
      throw std::logic_error("no loss set");
    }
    return *loss_;
  }

  void loss(const LossFunction& loss) noexcept {
    loss_.reset(new LossFunction(loss));
    prox_.loss(loss_.get());

    x_col_sum_ = arma::sum(loss_->data().cx()).t();
    const double norm_x = loss_->IncludeIntercept() ?
      arma::norm(arma::join_rows(arma::ones(loss_->data().n_obs(), 1), loss_->data().cx()), 2) :
      arma::norm(loss_->data().cx(), 2);

    operator_scaling_g_ = 1 / (norm_x * norm_x);  // this is `tau` in Deng & Yin (2016)
  }

  PenaltyFunction& penalty() const {
    if (!penalty_) {
      throw std::logic_error("no penalty set");
    }
    return *penalty_;
  }

  void penalty(const PenaltyFunction& penalty) noexcept {
    penalty_.reset(new PenaltyFunction(penalty));
  }

  //! Get the convergence tolerance for the ADMM algorithm.
  //!
  //! @return convergence tolerance.
  double convergence_tolerance() const noexcept {
    return convergence_tolerance_;
  }

  //! Set the convergence tolerance for the ADMM algorithm.
  //!
  //! @param convergence_tolerance convergene tolerance for the ADMM algorithm.
  void convergence_tolerance(double convergence_tolerance) noexcept {
    convergence_tolerance_ = convergence_tolerance;
  }

  //! Find the minimum of the objective function, using the previous solution (or the 0-vector if no
  //! previous solution exists) as starting point.
  //!
  //! @return information about the optimum.
  Optimum Optimize() {
    return Optimize(config_.max_it);
  }

  //! Find the minimum of the objective function, using the given coefficients as starting point.
  //!
  //! @param start where to start the optimization from.
  //! @return information about the optimum.
  Optimum Optimize(const Coefficients& start) {
    return Optimize(start, config_.max_it);
  }

  //! Find the minimum of the objective function, using the given coefficients as starting point
  //! and at most ``max_it`` iterations.
  //!
  //! @param start where to start the optimization from.
  //! @param max_it maximum number of iterations.
  //! @return information about the optimum.
  Optimum Optimize(const Coefficients& start, const int max_it) {
    if (!loss_) {
      throw std::logic_error("no loss set");
    }

    coefs_ = start;
    // Reset the lagrangian
    state_.lagrangian.reset();
    return Optimize(max_it);
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

    const PredictorResponseData& data = loss_->data();
    const bool include_intercept = loss_->IncludeIntercept();

    // Check if the coefficients are correct.
    if (coefs_.beta.n_elem != data.n_pred()) {
      coefs_.beta.zeros(data.n_pred());
      coefs_.intercept = 0;
    }

    auto metrics = std::make_unique<Metrics>("admm");
    metrics->AddDetail("type", "linearized");

    operator_scaling_f_ = prox_.OperatorScaling();  // this is (1/beta) in Deng & Yin (2016)
    const double scaled_lambda = penalty_->lambda() * prox_.PenaltyScaling();

    const auto en_cutoff = DetermineCutoff(scaled_lambda, IsAdaptiveTag{});
    const auto en_multiplier = DetermineEnMultiplier(scaled_lambda, IsAdaptiveTag{});

    double gap = 0;

    arma::vec fitted_step_1 = data.cx() * coefs_.beta;
    // Check if the lagrangian can be reused or if it needs to be reinitialized.
    if (state_.lagrangian.n_elem != data.n_obs()) {
      state_.lagrangian.zeros(data.n_obs());
      state_.fitted = prox_(fitted_step_1, coefs_.intercept, operator_scaling_f_);
    } else {
      // Continue to use the previous lagrangian
      state_.fitted = prox_(fitted_step_1 + state_.lagrangian * operator_scaling_f_, coefs_.intercept,
                            operator_scaling_f_);
    }

    fitted_step_1 -= state_.fitted;

    metrics->AddDetail("convergence_tolerance", convergence_tolerance_);
    metrics->AddDetail("op_scaling_g", operator_scaling_g_);
    metrics->AddDetail("op_scaling_f", operator_scaling_f_);

    int iter = 0;
    State prev_state;

    while (iter++ < max_it) {
      Metrics& iter_metrics = metrics->CreateSubMetrics("admm-iteration");
      prev_state = state_;

      if (include_intercept) {
        const double intercept = coefs_.intercept;

        coefs_.intercept -= operator_scaling_g_ * (data.n_obs() * coefs_.intercept +
          arma::dot(coefs_.beta, x_col_sum_) - arma::accu(state_.fitted - state_.lagrangian * operator_scaling_f_));

        // remember: fitted_step_1 is already fitted_step_1 - state_.fitted
        coefs_.beta = UpdateSlope(en_multiplier, SoftThreshold(coefs_.beta, -operator_scaling_g_,
          intercept * x_col_sum_ + data.cx().t() * (fitted_step_1 + operator_scaling_f_ * state_.lagrangian),
          en_cutoff));
      } else {
        // remember: fitted_step_1 is already fitted_step_1 - state_.fitted
        coefs_.beta = UpdateSlope(en_multiplier, SoftThreshold(coefs_.beta, -operator_scaling_g_,
          data.cx().t() * (fitted_step_1 + operator_scaling_f_ * state_.lagrangian), en_cutoff));
      }

      fitted_step_1 = data.cx() * coefs_.beta;

      state_.fitted = prox_(fitted_step_1 + state_.lagrangian * operator_scaling_f_, state_.fitted,
                            coefs_.intercept, operator_scaling_f_, &(iter_metrics.CreateSubMetrics("prox")));

      // Instead of `lagrangian -= accelerate * (fitted - fitted_step_1 - intercept) / operator_scaling_f_`, do
      fitted_step_1 -= state_.fitted;
      state_.lagrangian += (fitted_step_1 + coefs_.intercept) * config_.accelerate / operator_scaling_f_;

      const double fitted_diff = arma::accu(arma::square(state_.fitted - prev_state.fitted));
      const double lagrangian_diff = arma::accu(arma::square(state_.lagrangian - prev_state.lagrangian));

      gap = (fitted_diff + lagrangian_diff);

      iter_metrics.AddDetail("fitted_diff", fitted_diff);
      iter_metrics.AddDetail("lagrangian_diff", lagrangian_diff);
      iter_metrics.AddDetail("gap", gap);

      if (gap < convergence_tolerance_) {
        return FinalizeResult(iter, gap, state_.fitted, OptimumStatus::kOk, std::move(metrics));
      }
    }
    return FinalizeResult(--iter, gap, state_.fitted,
                          OptimumStatus::kWarning, "ADMM-algorithm did not converge.", std::move(metrics));
  }

 private:
  Optimum FinalizeResult(const int iter, const double gap, const arma::vec& fitted, const OptimumStatus status,
                         std::unique_ptr<Metrics> metrics) {
    // Update the active set.
    metrics->AddMetric("iter", iter);
    metrics->AddMetric("gap", gap);
    return MakeOptimum(*loss_, *penalty_, coefs_, loss_->data().cy() - fitted, std::move(metrics), status);
  }

  Optimum FinalizeResult(const int iter, const double gap, const arma::vec& fitted, const OptimumStatus status,
                         const std::string& message, std::unique_ptr<Metrics> metrics) {
    metrics->AddMetric("iter", iter);
    metrics->AddMetric("gap", gap);
    return MakeOptimum(*loss_, *penalty_, coefs_, loss_->data().cy() - fitted, std::move(metrics), status, message);
  }

  //! Determine the cutoff for the soft-threshold function for adaptive penalties
  //!
  arma::vec DetermineCutoff(const double scaled_lambda, std::true_type /* is_adaptive */) const noexcept {
    return penalty_->loadings() * DetermineCutoff(scaled_lambda, std::false_type{});
  }

  //! Determine the cutoff for the soft-threshold function for non-adaptive penalties
  //!
  //! @param scaled_lambda scaled adaptive EN penalty parameter
  double DetermineCutoff(const double scaled_lambda, std::false_type) const noexcept {
    return penalty_->alpha() * scaled_lambda * operator_scaling_g_ * operator_scaling_f_;
  }

  //! Determine the EN multiplier for adaptive penalties
  //!
  arma::vec DetermineEnMultiplier(const double scaled_lambda, std::true_type /* is_adaptive */) const noexcept {
    return 1 / (1 + penalty_->loadings() * scaled_lambda * (1 - penalty_->alpha()) *
      operator_scaling_g_ * operator_scaling_f_);
  }

  //! Determine the EN multiplier for non-adaptive penalties
  //!
  //! @param scaled_lambda scaled adaptive EN penalty parameter
  double DetermineEnMultiplier(const double scaled_lambda, std::false_type) const noexcept {
    return 1 / (1 + scaled_lambda * (1 - penalty_->alpha()) * operator_scaling_g_ * operator_scaling_f_);
  }

  //! Update the slope coefficients for adaptive penalties
  //!
  template<typename T>
  T UpdateSlope(const arma::vec& en_mult, const T& soft_thresh) const noexcept {
    return en_mult % soft_thresh;
  }

  //! Update the slope coefficients for non-adaptive penalties
  //!
  //! @param scaled_lambda scaled adaptive EN penalty parameter
  template<typename T>
  T UpdateSlope(const double en_mult, const T& soft_thresh) const noexcept {
    return en_mult * soft_thresh;
  }

  const AdmmLinearConfiguration config_;
  ProximalOperator prox_;
  LossFunctionPtr loss_;
  PenaltyPtr penalty_;
  Coefficients coefs_;
  State state_;
  arma::vec x_col_sum_;
  double operator_scaling_g_;  //< this is \tau in Deng (2016)
  double operator_scaling_f_;  //< this is 1 / \rho in Deng (2016)
  double convergence_tolerance_ = 1e-6;
};

//! Alias to make the GenericLinearizedAdmmOptimizer template more versatile. It either excepts a LS-type loss
//! function class, or a proximal operator.
template <typename LossFunction, typename PenaltyFunction, typename Coefficients>
using LinearizedAdmmOptimizer = GenericLinearizedAdmmOptimizer<admm_optimizer::ProximalOperator<LossFunction>,
                                                               PenaltyFunction, Coefficients>;

//! Compute the EN regression estimate using the alternating direction method of multiplier (ADMM)
//! with variable step-size.
template <typename LossFunction, typename PenaltyFunction, typename Coefficients>
class AdmmVarStepOptimizer : public Optimizer<LossFunction, PenaltyFunction, Coefficients> {
  using Base = Optimizer<LossFunction, PenaltyFunction, Coefficients>;
  using LossFunctionPtr = std::unique_ptr<LossFunction>;
  using PenaltyPtr = std::unique_ptr<PenaltyFunction>;
  using IsWeightedTag = typename traits::is_weighted<LossFunction>::type;
  using IsAdaptiveTag = typename traits::is_adaptive<PenaltyFunction>::type;
  using IsSparseTag = typename std::is_same<typename Coefficients::SlopeCoefficient, arma::sp_vec>::type;
  using Weights = typename std::conditional<IsWeightedTag::value, arma::vec, double>::type;
  // using Weights = typename std::conditional<IsWeightedTag::value, arma::vec, double>::type;

  static_assert(traits::is_en_penalty<PenaltyFunction>::value, "PenaltyFunction must be an EN-type penalty.");
  static_assert(traits::is_ls_regression_loss<LossFunction>::value, "LossFunction must be an least-squares-type loss.");

  // ADMM state structure
  struct State {
    double gap = -1;
    double tau = -1;
    double tau_lower;
    double gamma;
    arma::vec v;
    arma::vec l;
  };

 public:
  using Optimum = typename Base::Optimum;

  //! Ininitialize the optimizer using the given (weighted) LS loss function and penalty function.
  //!
  //! @param loss a weighted LS loss function.
  //! @param penalty penalty function.
  explicit AdmmVarStepOptimizer(const AdmmVarStepConfiguration& config = admm_optimizer::kDefaultVarStepConfig) noexcept
      : config_(config), loss_(nullptr), penalty_(nullptr), data_(nullptr) {}

  //! Ininitialize the optimizer using the given (weighted) LS loss function and penalty function.
  //!
  //! @param loss a weighted LS loss function.
  //! @param penalty penalty function.
  AdmmVarStepOptimizer(const LossFunction& loss, const PenaltyFunction& penalty,
                       const AdmmVarStepConfiguration& config = admm_optimizer::kDefaultVarStepConfig) noexcept
      : config_(config), loss_(new LossFunction(loss)), penalty_(new PenaltyFunction(penalty)), data_(nullptr) {}

  //! Default copy constructor.
  //!
  //! The copied optimizer will use the same loss and penalty functions after construction.
  AdmmVarStepOptimizer(const AdmmVarStepOptimizer& other) noexcept
      : config_(other.config_),
        loss_(other.loss_? new LossFunction(*other.loss_) : nullptr),
        penalty_(other.penalty_ ? new PenaltyFunction(*other.penalty_) : nullptr),
        coefs_(other.coefs_),
        state_(other.state_),
        weighted_data_(other.weighted_data_ ? new PredictorResponseData(*other.weighted_data_) : nullptr),
        data_(weighted_data_.get()),
        convergence_tolerance_(other.convergence_tolerance_) {}

  //! Default copy assignment.
  //!
  //! The copied optimizer will use the same loss and penalty functions after construction.
  AdmmVarStepOptimizer& operator=(const AdmmVarStepOptimizer& other) = default;

  //! Default move constructor.
  AdmmVarStepOptimizer(AdmmVarStepOptimizer&& other) = default;

  //! Default move assignment operator.
  AdmmVarStepOptimizer& operator=(AdmmVarStepOptimizer&& other) = default;

  ~AdmmVarStepOptimizer() = default;

  void Reset() {
    state_.gap = -1;
    data_ = nullptr;
  }

  LossFunction& loss() const {
    if (!loss_) {
      throw std::logic_error("no loss set");
    }
    return *loss_;
  }

  void loss(const LossFunction& loss) noexcept {
    loss_.reset(new LossFunction(loss));
    Reset();
  }

  PenaltyFunction& penalty() const {
    if (!penalty_) {
      throw std::logic_error("no penalty set");
    }
    return *penalty_;
  }

  void penalty(const PenaltyFunction& penalty) noexcept {
    penalty_.reset(new PenaltyFunction(penalty));
  }

  //! Get the convergence tolerance for the ADMM algorithm.
  //!
  //! @return convergence tolerance.
  double convergence_tolerance() const noexcept {
    return convergence_tolerance_;
  }

  //! Set the convergence tolerance for the ADMM algorithm.
  //!
  //! @param convergence_tolerance convergene tolerance for the ADMM algorithm.
  void convergence_tolerance(double convergence_tolerance) noexcept {
    convergence_tolerance_ = convergence_tolerance;
  }

  //! Find the minimum of the objective function, using the previous solution (or the 0-vector if no
  //! previous solution exists) as starting point.
  //!
  //! @return information about the optimum.
  Optimum Optimize() {
    return Optimize(config_.max_it);
  }

  //! Find the minimum of the objective function, using the given coefficients as starting point.
  //!
  //! @param start where to start the optimization from.
  //! @return information about the optimum.
  Optimum Optimize(const Coefficients& start) {
    return Optimize(start, config_.max_it);
  }

  //! Find the minimum of the objective function, using the given coefficients as starting point
  //! and at most ``max_it`` iterations.
  //!
  //! @param start where to start the optimization from.
  //! @param max_it maximum number of iterations.
  //! @return information about the optimum.
  Optimum Optimize(const Coefficients& start, const int max_it) {
    if (!loss_) {
      throw std::logic_error("no loss set");
    }

    if (!data_) {
      UpdateData(IsWeightedTag{});
    }

    coefs_ = start;
    state_.gap = -1;
    return Optimize(max_it);
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

    // Check if the data needs to be updated
    if (!data_) {
      UpdateData(IsWeightedTag{});
    }

    const bool include_intercept = loss_->IncludeIntercept();
    const double scaled_lambda = ScaledLambda(IsWeightedTag{});
    const bool check_empty = admm_optimizer::AllZero(coefs_.beta) || (coefs_.beta.n_elem != data_->n_pred());

    auto metrics = std::make_unique<Metrics>("admm");
    metrics->AddDetail("type", "var-stepsize");

    // Check if the coefficients are correct.
    if (coefs_.beta.n_elem != data_->n_pred()) {
      coefs_.beta.zeros(data_->n_pred());
      coefs_.intercept = 0;
    }

    if (include_intercept) {
      EmptyModelIntercept(IsWeightedTag{});
    }
    // Check if any of the predictors might be active.
    if (check_empty && !admm_optimizer::AnyViolateKKT(data_->cx(), EmptyModelResiduals(IsWeightedTag{}),
                                                      scaled_lambda / data_->n_obs(), *penalty_)) {
      // None of the predictors will be activated for the current penalty. Return the current coefficient value.
      return FinalizeResult(0, OptimumStatus::kOk, std::move(metrics));
    }

    // Compute the upper limit for tau (if we would do a restart)
    const double tau_upper = (config_.tau > 0) ? config_.tau : ((1.1 - penalty_->alpha()) * scaled_lambda);
    state_.tau = tau_upper;
    state_.gamma = config_.tau_adjustment_lower;
    state_.tau_lower = state_.tau * config_.tau_lower_mult;

    // Check if the state needs to be re-initialized
    if (state_.gap < 0) {
      state_.l.zeros(data_->n_pred());
    }

    // This is the convergence tolerance for the "standardized" residual.
    const double conv_tol = convergence_tolerance_ * penalty_->alpha();

    double tau_inv = 1 / state_.tau;
    auto en_cutoff = DetermineCutoff(scaled_lambda, IsAdaptiveTag{});
    double en_multiplier = state_.tau / (state_.tau + scaled_lambda * (1 - penalty_->alpha()));

    metrics->AddDetail("convergence_tolerance", conv_tol);
    metrics->AddDetail("tau-start", state_.tau);
    metrics->AddDetail("gamma-start", state_.gamma);

    if (state_.gap < 0 || cache_.chol_xtx_tau != state_.tau) {
      if (!UpdateCholesky()) {
        return FinalizeResult(0, OptimumStatus::kError, "Predictor matrix is singular. Cholesky decomposition failed.",
                              std::move(metrics));
      }
    }

    // Initialize the `v` vector
    if (state_.gap < 0) {
      const bool prox_ls_success = ProximalLs();
      if (!prox_ls_success) {
        return FinalizeResult(0, OptimumStatus::kError, "Predictor matrix is singular. Can not solve equations.",
                              std::move(metrics));
      }
    }

    int iter = 0;
    double intercept_change = 0;
    while (iter++ < max_it) {
      State prev_state = state_;
      coefs_.beta = en_multiplier * SoftThreshold(state_.v, tau_inv, state_.l, en_cutoff);
      if (include_intercept) {
        const double new_intercept = ComputeIntercept(IsWeightedTag{});
        intercept_change = coefs_.intercept - new_intercept;
        intercept_change *= intercept_change;
        coefs_.intercept = new_intercept;
      }

      const bool prox_ls_success = ProximalLs(state_.l);
      if (!prox_ls_success) {
        return FinalizeResult(iter, OptimumStatus::kError, "Singular predictor matrix.", std::move(metrics));
      }

      UpdateStateL(IsSparseTag{});

      const double diff_scaling = arma::norm(state_.v, 2) + state_.tau * arma::norm(state_.l, 2);
      state_.gap = arma::accu(state_.tau * state_.tau * arma::square(state_.v - prev_state.v) +
                              arma::square(state_.l - prev_state.l)) + intercept_change;

      if (state_.gap * diff_scaling * diff_scaling < conv_tol) {
        return FinalizeResult(iter, OptimumStatus::kOk, std::move(metrics));
      }

      if (iter > 1 && state_.gap > prev_state.gap * state_.gamma) {
        if (state_.tau > state_.tau_lower) {
          state_.tau = std::max(state_.tau_lower, state_.tau * state_.gamma);
        } else {
          state_.tau = tau_upper;
          state_.gamma = std::min(0.5 * state_.gamma + 0.5, config_.tau_adjustment_upper);
          state_.l.zeros(data_->n_pred());
          state_.v.zeros(data_->n_pred());
        }

        if (!UpdateCholesky()) {
          return FinalizeResult(iter, OptimumStatus::kError, "Predictor matrix is singular. "
                                "Cholesky decomposition failed.", std::move(metrics));
        }
        tau_inv = 1 / state_.tau;
        en_cutoff = DetermineCutoff(scaled_lambda, IsAdaptiveTag{});
        en_multiplier = state_.tau / (state_.tau + scaled_lambda * (1 - penalty_->alpha()));
      }
    }

    return FinalizeResult(--iter, OptimumStatus::kWarning, "ADMM algorithm did not converge.", std::move(metrics));
  }

 private:
  Optimum FinalizeResult(const int iter, const OptimumStatus status,
                         std::unique_ptr<Metrics> metrics) {
    metrics->AddMetric("iter", iter);
    metrics->AddMetric("gap", state_.gap);
    metrics->AddDetail("tau-end", state_.tau);
    metrics->AddDetail("gamma-end", state_.gamma);
    return MakeOptimum(*loss_, *penalty_, coefs_, std::move(metrics), status);
  }

  Optimum FinalizeResult(const int iter, const OptimumStatus status,
                         const std::string& message, std::unique_ptr<Metrics> metrics) {
    metrics->AddMetric("iter", iter);
    metrics->AddMetric("gap", state_.gap);
    metrics->AddDetail("tau-end", state_.tau);
    metrics->AddDetail("gamma-end", state_.gamma);
    return MakeOptimum(*loss_, *penalty_, coefs_, std::move(metrics), status, message);
  }

  //! Get the appropriately scaled lambda for a weighted LS loss
  inline double ScaledLambda(std::true_type) const noexcept {
    return data_->n_obs() * penalty_->lambda() / loss_->mean_weight();
  }
  //! Get the appropriately scaled lambda for an unweighted LS loss
  inline double ScaledLambda(std::false_type) const noexcept {
    return data_->n_obs() * penalty_->lambda();
  }

  //! Determine the cutoff for the soft-threshold function for adaptive penalties
  arma::vec DetermineCutoff(const double scaled_lambda, std::true_type) const {
    return penalty_->loadings() * (penalty_->alpha() * scaled_lambda  / state_.tau);
  }

  //! Determine the cutoff for the soft-threshold function for non-adaptive penalties
  double DetermineCutoff(const double scaled_lambda, std::false_type) const {
    return penalty_->alpha() * scaled_lambda / state_.tau;
  }

  //! Update the data for a weighted LS loss
  void UpdateData(std::true_type) {
    weighted_data_.reset(new PredictorResponseData(loss_->data().cx().each_col() % loss_->sqrt_weights(),
                                                   loss_->data().cy() % loss_->sqrt_weights()));
    data_ = weighted_data_.get();
    cache_.xtwgt = data_->cx().t() * loss_->sqrt_weights();
    UpdateCache();
  }

  //! Update the data for an un-weighted LS loss
  void UpdateData(std::false_type) {
    data_ = &(loss_->data());
    cache_.xtwgt = arma::sum(data_->cx(), 0).t();
    UpdateCache();
  }

  //! Update the cache of inner products
  void UpdateCache() {
    cache_.xty = data_->cx().t() * data_->cy();
    cache_.xtx = data_->cx().t() * data_->cx();
    UpdateCholesky();
  }

  //! Update the Choleskey decomposition of `X'X + tau I` for a new `tau`
  bool UpdateCholesky() {
    cache_.chol_xtx_tau = state_.tau;
    if (state_.tau > 0) {
      cache_.chol_xtx = cache_.xtx;
      cache_.chol_xtx.diag() += state_.tau;
      // Manually compute Cholesky decomposition to avoid unnecessary copying.
      return arma::auxlib::chol_simple(cache_.chol_xtx);
    }
    return true;
  }

  //! Apply the proximal operator to the vector `tau * beta + X'y - intercept w`
  bool ProximalLs() {
    state_.v = loss_->IncludeIntercept() ?
      arma::vec(state_.tau * coefs_.beta + cache_.xty - coefs_.intercept * cache_.xtwgt) :
      arma::vec(state_.tau * coefs_.beta + cache_.xty);
    return linalg::SolveChol(cache_.chol_xtx, &state_.v);
  }

  //! Apply the proximal operator to the vector `tau * beta + X'y - intercept w - l`
  bool ProximalLs(const arma::vec& l) {
    state_.v = loss_->IncludeIntercept() ?
      arma::vec(state_.tau * coefs_.beta + cache_.xty - coefs_.intercept * cache_.xtwgt - l) :
      arma::vec(state_.tau * coefs_.beta + cache_.xty - l);
    return linalg::SolveChol(cache_.chol_xtx, &state_.v);
  }

  //! Compute the intercept for weighted LS
  double ComputeIntercept(std::true_type) const {
    return arma::mean((data_->cy() - data_->cx() * coefs_.beta) % loss_->sqrt_weights());
  }

  //! Compute the intercept for unweighted LS
  double ComputeIntercept(std::false_type) const {
    return arma::mean(data_->cy() - data_->cx() * coefs_.beta);
  }

  //! Update the `l` vector in the state structure if using sparse coefficients.
  //! This is required because otherwise SAN/UBSAN complains about a stack-use-after-scope
  void UpdateStateL(std::true_type) {
    state_.l += state_.tau * (state_.v - arma::vec(coefs_.beta));
  }

  //! Update the `l` vector in the state structure if using dense coefficients.
  void UpdateStateL(std::false_type) {
    state_.l += state_.tau * (state_.v - coefs_.beta);
  }

  //! Compute the intercept in the empty model.
  void EmptyModelIntercept(std::true_type) {
    coefs_.intercept = arma::mean(data_->cy() % loss_->sqrt_weights());
  }

  //! Compute the intercept in the empty model.
  void EmptyModelIntercept(std::false_type) {
    coefs_.intercept = arma::mean(data_->cy());
  }

  //! Determine the residuals in the empty model when using weighted LS
  arma::vec EmptyModelResiduals(std::true_type) const {
    if (loss_->IncludeIntercept()) {
      return data_->cy() - loss_->sqrt_weights() * coefs_.intercept;
    } else {
      return data_->cy();
    }
  }

  //! Determine the residuals in the empty model when using unweighted LS
  arma::vec EmptyModelResiduals(std::false_type) {
    if (loss_->IncludeIntercept()) {
      return data_->cy() - coefs_.intercept;
    } else {
      return data_->cy();
    }
  }
  //! The ADMM algorithm uses a different criterion to determine convergence than other optimizers.
  //! To make the precision comparable to other optimizers, the convergence tolerance needs to be adjusted.
  // constexpr static kConvergenceToleranceAdjustment = 1e-2;
  const AdmmVarStepConfiguration config_;
  LossFunctionPtr loss_;
  PenaltyPtr penalty_;
  Coefficients coefs_;
  State state_;
  std::unique_ptr<PredictorResponseData> weighted_data_;
  PredictorResponseData const * data_;
  admm_optimizer::DataCache cache_;
  double convergence_tolerance_ = 1e-6;
};

}  // namespace nsoptim

#endif  // NSOPTIM_OPTIMIZER_ADMM_HPP_
