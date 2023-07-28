//
//  mm.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2019-01-02.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifndef NSOPTIM_OPTIMIZER_MM_HPP_
#define NSOPTIM_OPTIMIZER_MM_HPP_

#include <string>
#include <memory>
#include <type_traits>
#include <algorithm>

#include "../armadillo.hpp"
#include "../traits/traits.hpp"
#include "../container/metrics.hpp"
#include "optimizer_base.hpp"
#include "optimum.hpp"

namespace nsoptim {

//! Configuration options for the MM algorithm.
struct MMConfiguration {
  //! Type of tightening for inner optimization.
  enum class TighteningType {
    //! No tightening, i.e., always use the configured numeric tolerance for the inner optimization.
    //! This is automatically chosen if the inner optimizer does not support chaning the numeric tolerance.
    kNone = 0,

    //! Start with a large inner tolerance and at each iteration make the inner optimization more precise by
    //! reducing the tolerance level of the inner optimizer by a constant factor, up until the minimum inner
    //! tolerance level is reached.
    kExponential = 1,

    //! Start with a large inner tolerance and reduce by a constant factor as soon as parameter change in the outer
    //! optimization is less than the inner tolerance.
    kAdaptive = 2
  };

  //! Maximum number of iterations allowed.
  int max_it;
  //! Type of tightening for inner optimization.
  TighteningType tightening;
  //! Number of tightening steps if using adaptive thightening.
  int adaptive_tightening_steps;
};

namespace mm_optimizer {
//! Default configuration for the MM algorithm.
constexpr MMConfiguration kDefaultMMConfiguration = {500, MMConfiguration::TighteningType::kNone, 10};

template<typename Optimizer>
class InnerToleranceTightening {
  using IsIterativeAlgorithmTag = typename traits::is_iterative_algorithm<Optimizer>::type;

 public:
  InnerToleranceTightening() noexcept : inner_tol_(0) {}
  explicit InnerToleranceTightening(Optimizer* optimizer, const double inner_tol) noexcept
      : optimizer_(optimizer), inner_tol_(inner_tol) {}

  virtual ~InnerToleranceTightening() noexcept = default;

  //! Get the current convergence tolerance for the (inner) optimizer.
  //! If the optimizer does not support convergence tolerance, returns 0.
  //!
  //! @return the currently set inner convergence tolerance, or 0 if not supported.
  double current_tolerance() const noexcept {
    return current_tolerance(IsIterativeAlgorithmTag{});
  }

  //! Check if further tightening is possible.
  bool CanTightenFurther() const noexcept {
    return CanTightenFurther(IsIterativeAlgorithmTag{});
  }

  //! Tighten the inner convergence tolerance.
  //!
  //! @param outer_change the change in the objective function after the last iteration.
  virtual void Tighten(const double) noexcept {}

  //! Tighten the inner convergence tolerance more aggressively.
  virtual void FastTighten() noexcept {
    FullyTighten(IsIterativeAlgorithmTag{});
  }

  //! Make the inner convergence tolerance as small as possible.
  void FullyTighten() noexcept {
    FullyTighten(IsIterativeAlgorithmTag{});
  }

 protected:
  //! Set the new convergence tolerance for the inner optimizer, if supported.
  void PropagateUpdate(const double new_tolerance) noexcept {
    PropagateUpdate(new_tolerance, IsIterativeAlgorithmTag{});
  }

  //! Get the minimum (and desired) inner tolerance level.
  //!
  //! @return the inner tolerance level.
  double inner_tolerance() const noexcept {
    return inner_tol_;
  }

 private:
  bool CanTightenFurther(std::true_type) const noexcept {
    return inner_tol_ < optimizer_->convergence_tolerance();
  }

  constexpr bool CanTightenFurther(std::false_type) const noexcept {
    return false;
  }

  void FullyTighten(std::true_type) noexcept {
    optimizer_->convergence_tolerance(inner_tol_);
  }

  void FullyTighten(std::false_type) const noexcept {}

  void PropagateUpdate(const double new_tolerance, std::true_type) noexcept {
    optimizer_->convergence_tolerance(new_tolerance);
  }

  void PropagateUpdate(const double, std::false_type) const noexcept {}

  double current_tolerance(std::true_type) const noexcept {
    return optimizer_->convergence_tolerance();
  }

  constexpr double current_tolerance(std::false_type) const noexcept {
    return 0;
  }

  Optimizer * const optimizer_;
  const double inner_tol_;
};

//! At each iteration, make the inner optimization tighter, up until the minimum inner tolerance level
//! is reached.
//! The initial inner tolerance level is set the square root of the outer tolerance.
//! A tightening step is such that the minimum inner convergence tolerance is reached after half the maximum number of
//! iterations.
//! A fast tightening step is such that 10 fast tightening steps lead to the minimum inner tolerance level.
template<class Optimizer>
class ExponentialTightening : public InnerToleranceTightening<Optimizer> {
 public:
  ExponentialTightening(Optimizer* optimizer, const double outer_tolerance, const double inner_tolerance,
                        const int max_it) noexcept
      : InnerToleranceTightening<Optimizer>(optimizer, inner_tolerance),
        multiplier_(std::pow(inner_tolerance, 2. / max_it)), fast_multiplier_(std::pow(inner_tolerance, 0.1)) {
    // Start with the sqrt of the outer tolerance.
    this->PropagateUpdate(std::sqrt(outer_tolerance));
  }

  virtual ~ExponentialTightening() noexcept = default;

  void Tighten(const double outer_change) noexcept override {
    const double new_tol = this->current_tolerance() * multiplier_;
    this->PropagateUpdate(std::max(this->inner_tolerance(), std::min(new_tol, outer_change)));
  }

  void FastTighten() noexcept override {
    this->PropagateUpdate(this->current_tolerance() * fast_multiplier_);
  }

 private:
  const double multiplier_;
  const double fast_multiplier_;
};


//! Make the inner optimization tighther as soon as the change in the objective function is less
//! than the inner optimization tolerance.
//! Fast tightening takes two tightening steps at once.
template<class Optimizer>
class AdaptiveTightening : public InnerToleranceTightening<Optimizer> {
 public:
  explicit AdaptiveTightening(Optimizer* optimizer, const double outer_tolerance, const double inner_tolerance,
                              const int steps) noexcept
      : InnerToleranceTightening<Optimizer>(optimizer, inner_tolerance),
        multiplier_(std::pow(inner_tolerance / std::sqrt(outer_tolerance), 1. / steps)),
        min_inner_tolerance_(inner_tolerance * 0.1) {
    // Start with the sqrt of the outer tolerance.
    this->PropagateUpdate(std::sqrt(outer_tolerance));
  }
  virtual ~AdaptiveTightening() noexcept = default;

  void Tighten(const double outer_change) noexcept override {
    const double tol = this->current_tolerance();
    if (outer_change < tol) {
      this->PropagateUpdate(std::max(tol * multiplier_, min_inner_tolerance_));
    }
  }

  void FastTighten() noexcept override {
    this->PropagateUpdate(std::max(this->current_tolerance() * multiplier_ * multiplier_, min_inner_tolerance_));
  }

 private:
  const double multiplier_;
  const double min_inner_tolerance_;
};

}  // namespace mm_optimizer

//! Compute the minimum of a non-convex objective function (the loss and/or the penalty can be non-convex) using
//! the Minimization by Majorization (MM) algorithm.
template <typename LossFunction, typename PenaltyFunction, typename InnerOptimizerType, typename Coefficients>
class MMOptimizer : public Optimizer<LossFunction, PenaltyFunction, Coefficients> {
  using LossFunctionPtr = std::unique_ptr<LossFunction>;
  using PenaltyFunctionPtr = std::unique_ptr<PenaltyFunction>;
  using IsIterativeAlgorithmTag = typename traits::is_iterative_algorithm<InnerOptimizerType>::type;
  using ResidType = typename LossFunction::ResidualType;

  static_assert(traits::has_convex_surrogate<LossFunction, Coefficients>::value,
                "LossFunction does not provide a convex surrogate.");
  static_assert(traits::has_convex_surrogate<PenaltyFunction, Coefficients>::value,
                "PenaltyFunction does not provide a convex surrogate.");
  static_assert(traits::can_optimize<InnerOptimizerType, Coefficients>(),
                "InnerOptimizerType can not optimize the requested coefficient type.");

 public:
  using InnerOptimizer = InnerOptimizerType;
  using Optimum = typename Optimizer<LossFunction, PenaltyFunction, Coefficients>::Optimum;

  //! Ininitialize the MM algorithm without loss or penalty function.
  //!
  //! @param config configuration for the MM optimizer.
  explicit MMOptimizer(const MMConfiguration& config = mm_optimizer::kDefaultMMConfiguration) noexcept
      : config_(config), optimizer_(),
        inner_convergence_tolerance_(InnerConvergenceTolerance(IsIterativeAlgorithmTag{})) {}

  //! Ininitialize the MM algorithm using the given inner optimizer,
  //!
  //! @param optimizer optimizer to use for the inner optimization.
  //! @param config configuration for the MM optimizer.
  MMOptimizer(const InnerOptimizerType& optimizer,
              const MMConfiguration& config = mm_optimizer::kDefaultMMConfiguration) noexcept
      : config_(config), optimizer_(optimizer),
        inner_convergence_tolerance_(InnerConvergenceTolerance(IsIterativeAlgorithmTag{})) {}

  //! Ininitialize the MM algorithm using the given loss function and the penalty function.
  //!
  //! @param loss a loss function.
  //! @param penalty a penalty function.
  //! @param config configuration for the MM optimizer.
  MMOptimizer(const LossFunction& loss, const PenaltyFunction& penalty,
              const MMConfiguration& config = mm_optimizer::kDefaultMMConfiguration) noexcept
      : config_(config), loss_(new LossFunction(loss)), penalty_(new PenaltyFunction(penalty)),
        optimizer_(), inner_convergence_tolerance_(InnerConvergenceTolerance(IsIterativeAlgorithmTag{})) {}

  //! Ininitialize the MM algorithm using the given loss function, the penalty function, and inner optimizer.
  //!
  //! @param loss a loss function.
  //! @param penalty a penalty function.
  //! @param optimizer optimizer to use for the inner optimization.
  //! @param config configuration for the MM optimizer.
  MMOptimizer(const LossFunction& loss, const PenaltyFunction& penalty, const InnerOptimizerType& optimizer,
              const MMConfiguration& config = mm_optimizer::kDefaultMMConfiguration) noexcept
      : config_(config), loss_(new LossFunction(loss)), penalty_(new PenaltyFunction(penalty)),
        optimizer_(optimizer), inner_convergence_tolerance_(InnerConvergenceTolerance(IsIterativeAlgorithmTag{})) {}

  //! Default copy constructor.
  MMOptimizer(const MMOptimizer& other)
    : config_(other.config_),
      loss_(other.loss_ ? new LossFunction(*other.loss_) : nullptr),
      penalty_(other.penalty_ ? new PenaltyFunction(*other.penalty_) : nullptr),
      optimizer_(other.optimizer_), coefs_(other.coefs_), convergence_tolerance_(other.convergence_tolerance_),
      inner_convergence_tolerance_(other.inner_convergence_tolerance_) {}

  //! Default move constructor.
  MMOptimizer(MMOptimizer&& other) = default;

  //! Default move assignment operator.
  MMOptimizer& operator=(MMOptimizer&& other) = default;

  ~MMOptimizer() = default;

  //! Reset the optimizier. This compeletely purges the current *state*.
  void Reset() {
    coefs_.Reset();
    optimizer_.Reset();
  }

  //! Access the loss function.
  //!
  //! @return the loss function currently in use by the MM algorithm.
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
    loss_.reset(new LossFunction(loss));
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

  //! Get the convergence tolerance for the MM algorithm.
  //!
  //! @return convergence tolerance.
  double convergence_tolerance() const noexcept {
    return convergence_tolerance_;
  }

  //! Set the convergence tolerance for the MM algorithm.
  //!
  //! @param convergence_tolerance convergene tolerance for the MM algorithm.
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
    coefs_ = start;
    optimizer_.Reset();
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
    optimizer_.Reset();
    return Optimize(max_it);
  }

  //! Find the minimum of the objective function, using the previous solution (or the 0-vector if no
  //! previous solution exists) as starting point and at most ``max_it`` iterations.
  //!
  //! @param max_it maximum number of iterations.
  //! @return information about the optimum.
  Optimum Optimize(const int max_it) {
    using LossHasDifferenceOp = typename traits::has_difference_op<LossFunction, Coefficients>::type;
    using PenaltyHasDifferenceOp = typename traits::has_difference_op<LossFunction, Coefficients>::type;

    if (!loss_) {
      throw std::logic_error("no loss set");
    }

    if (!penalty_) {
      throw std::logic_error("no penalty set");
    }

    auto metrics = std::make_unique<Metrics>("mm-algorithm");
    // If the coefficients are not yet initialized, reset them to the 0-vector.
    if (coefs_.beta.n_elem == 0) {
      coefs_ = loss_->template ZeroCoefficients<Coefficients>();
      // Ensure that the optimizer also looses it's state.
      optimizer_.Reset();
    }

    // Set the convex surrogates for the internal optimizer.
    auto residuals = loss_->Residuals(coefs_);
    optimizer_.loss(loss_->GetConvexSurrogate(residuals));
    optimizer_.penalty(penalty_->GetConvexSurrogate(coefs_));

    std::unique_ptr<mm_optimizer::InnerToleranceTightening<InnerOptimizerType>> tightener;

    switch (config_.tightening) {
      case MMConfiguration::TighteningType::kExponential:
        tightener.reset(new mm_optimizer::ExponentialTightening<InnerOptimizerType>(
          &optimizer_, convergence_tolerance_, inner_convergence_tolerance_, config_.max_it));
        break;
      case MMConfiguration::TighteningType::kAdaptive:
        tightener.reset(new mm_optimizer::AdaptiveTightening<InnerOptimizerType>(
          &optimizer_, convergence_tolerance_, inner_convergence_tolerance_, config_.adaptive_tightening_steps));
        break;
      case MMConfiguration::TighteningType::kNone:
      default:
        tightener.reset(new mm_optimizer::InnerToleranceTightening<InnerOptimizerType>(
          &optimizer_, inner_convergence_tolerance_));
    }

    // Evaluate the loss and penalty functions after computing the convex surrogate.
    // Often, the results from the convex surrogate can be used to quickly evaluate the loss/penalty function.
    double objf_value = loss_->Evaluate(residuals) + penalty_->Evaluate(coefs_);
    double rel_difference = 0;
    bool final_iterations = false;
    int iter = 0;

    while (iter++ < max_it) {
      // UpdateInnerConvergenceTolerance(rel_tol, IsIterativeAlgorithmTag{});
      auto&& iter_metrics = metrics->CreateSubMetrics("mm_iteration");

      // Compute the minimizer of the convex surrogates.
      auto optimum = optimizer_.Optimize();
      if (optimum.metrics) {
        iter_metrics.AddSubMetrics(std::move(*optimum.metrics));
        optimum.metrics.reset();
      }

      // Check for any problems.
      if (optimum.status == OptimumStatus::kError) {
        metrics->AddDetail("final_rel_difference", rel_difference);
        metrics->AddDetail("final_innner_tol", tightener->current_tolerance());
        metrics->AddMetric("iter", iter);
        return MakeOptimum(*loss_, *penalty_, coefs_, std::move(metrics), OptimumStatus::kError,
                           std::string("MM-iteration failed: ").append(optimum.message));
      }

      // Check for convergence.
      rel_difference = Difference(coefs_, residuals, optimum.coefs, optimum.residuals,
                                  LossHasDifferenceOp{}, PenaltyHasDifferenceOp{});

      iter_metrics.AddDetail("iter", iter);
      iter_metrics.AddDetail("rel_difference", rel_difference);
      // iter_metrics.AddDetail("coef_norm_diff", arma::norm(coefs_.beta - optimum.coefs.beta, 2));
      // iter_metrics.AddDetail("resid_norm_diff", arma::norm(residuals - optimum.residuals, 2));
      iter_metrics.AddDetail("objf_value", optimum.objf_value);
      iter_metrics.AddDetail("inner_tol", tightener->current_tolerance());

      // The solutions are good enough. If we can make the inner optimizer tighter, add a few iteration and then
      // return the obtained optimum.
      if (rel_difference < convergence_tolerance_) {
        if (final_iterations || !tightener->CanTightenFurther()) {
          coefs_ = std::move(optimum.coefs);
          metrics->AddMetric("iter", iter);
          metrics->AddDetail("final_rel_difference", rel_difference);
          metrics->AddDetail("final_innner_tol", tightener->current_tolerance());
          // The value of the convex surrogate should be close enough (actually it should be equal) to the
          // true objective. Therefore, there is no need to re-evaluate the loss.
          return MakeOptimum(*loss_, *penalty_, coefs_, optimum.residuals, optimum.objf_value, std::move(metrics));
        } else {
          final_iterations = true;
          tightener->FullyTighten();
        }
      }

      // Check if the value of the convex surrogate decreased. This can only fail if the inner optimizer is using an
      // iterative scheme and the relative tolerance is too large.
      if (IsIterativeAlgorithmTag::value) {
        if (objf_value > 0 && (objf_value - optimum.objf_value) < -convergence_tolerance_) {
          // The value of the objective function increased substantially.
          // Decrease the inner convergence tolerance considerably and continue iterating the current surrogate,
          // i.e., don't update the coefficients or the surrogate.
          // If the inner convergence tolerance is already very small, return the result as-is.
          // (Any step increases the objective function, hence we are at a minimum.)
          if (!tightener->CanTightenFurther()) {
            metrics->AddMetric("iter", iter);
            metrics->AddDetail("final_rel_difference", rel_difference);
            metrics->AddDetail("final_innner_tol", tightener->current_tolerance());
            return MakeOptimum(*loss_, *penalty_, coefs_, residuals, objf_value, std::move(metrics),
                               OptimumStatus::kOk);
          }
          iter_metrics.AddDetail("tighten_faster", "yes");
          tightener->FastTighten();
          continue;
        }
      }

      // Continue iterations.
      coefs_ = std::move(optimum.coefs);
      residuals = std::move(optimum.residuals);

      // Make inner iteration more precise, i.e., decrease the relative tolerance.
      tightener->Tighten(rel_difference);

      // Update the convex surrogates for the internal optimizer.
      try {
        optimizer_.loss(loss_->GetConvexSurrogate(residuals));
      } catch(...) {
        metrics->AddMetric("iter", iter);
        metrics->AddDetail("final_rel_difference", rel_difference);
        metrics->AddDetail("final_innner_tol", tightener->current_tolerance());
        return MakeOptimum(*loss_, *penalty_, coefs_, residuals, std::move(metrics), OptimumStatus::kWarning,
                           "MM-algorithm did not converge");
      }
      optimizer_.penalty(penalty_->GetConvexSurrogate(residuals));

      // Retain value of the objective function to check for improvement.
      objf_value = loss_->Evaluate(residuals) + penalty_->Evaluate(coefs_);
    }

    metrics->AddMetric("iter", iter);
    metrics->AddDetail("final_rel_difference", rel_difference);
    metrics->AddDetail("final_innner_tol", tightener->current_tolerance());
    return MakeOptimum(*loss_, *penalty_, coefs_, residuals, std::move(metrics), OptimumStatus::kWarning,
                       "MM-algorithm did not converge");
  }

 private:
  //! Get the difference between two coefficients if the loss and penalty functions both have a `Difference` operator.
  double Difference(const Coefficients& c1, const ResidType& r1, const Coefficients& c2, const ResidType&,
                    std::true_type, std::true_type) const {
    return loss_->Difference(c1, c2) + penalty_->Difference(c1, c2);
  }
  //! Get the difference between two coefficients if only the loss has a `Difference` operator.
  double Difference(const Coefficients& c1, const ResidType& r1, const Coefficients& c2, const ResidType&,
                    std::true_type, std::false_type) const {
    return loss_->Difference(c1, c2) + std::abs(penalty_->Evaluate(c1) - penalty_->Evaluate(c2));
  }
  //! Get the difference between two coefficients if only the penalty has a `Difference` operator.
  double Difference(const Coefficients& c1, const ResidType& r1, const Coefficients& c2, const ResidType& r2,
                    std::false_type, std::true_type) const {
    return std::abs(loss_->Evaluate(r1) - loss_->Evaluate(r2)) + penalty_->Difference(c1, c2);
  }
  //! Get the difference between two coefficients if neither the loss nor the penalty has a `Difference` operator.
  double Difference(const Coefficients& c1, const ResidType& r1, const Coefficients& c2, const ResidType& r2,
                    std::false_type, std::false_type) const {
    return std::abs(loss_->Evaluate(r1) + penalty_->Evaluate(c1) - (loss_->Evaluate(r2) + penalty_->Evaluate(c2)));
  }

  double InnerConvergenceTolerance(std::true_type) const noexcept {
    return 0.5 * optimizer_.convergence_tolerance();
  }

  constexpr double InnerConvergenceTolerance(std::false_type) const noexcept {
    return 0;
  }

  // //! Update the convergence tolerance for inner optimizers which are iteratively computing the optimum.
  // void UpdateInnerConvergenceTolerance(const double rel_tol, std::true_type) noexcept {
  //   optimizer_.convergence_tolerance(rel_tol);
  // }

  // //! Update the convergence tolerance for inner optimizers that do not support that.
  // void UpdateInnerConvergenceTolerance(const double, std::false_type) noexcept {}

  const MMConfiguration config_;
  LossFunctionPtr loss_;
  PenaltyFunctionPtr penalty_;
  InnerOptimizerType optimizer_;
  Coefficients coefs_;
  double convergence_tolerance_ = 1e-8;
  double inner_convergence_tolerance_ = 0;
};
}  // namespace nsoptim

#endif  // NSOPTIM_OPTIMIZER_MM_HPP_
