//
//  auglars.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2019-11-14.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifndef NSOPTIM_OPTIMIZER_AUGLARS_HPP_
#define NSOPTIM_OPTIMIZER_AUGLARS_HPP_

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
#include "../traits/traits.hpp"
#include "linear_algebra_utilities.hpp"

namespace nsoptim {

namespace auglars {
//! Proxy class to get the slope coefficients from the LarsPath.
class BetaProxy {
 public:
  //! Create a proxy to an all-0 slope coefficient.
  //!
  //! @param size number of slope coefficients.
  explicit BetaProxy(const arma::uword size) noexcept : size_(size) {}

  //! Create a proxy to the slope coefficients.
  //!
  //! @param size number of slope coefficients.
  //! @param values values of the slope coefficients. Must have the same number of elements as *indices*.
  //! @param indices the indicies of the slope coefficients in *values*.
  BetaProxy(const arma::uword size, double const * values, const arma::uvec& indices) noexcept
    : size_(size), indices_(indices), values_(new double[indices_.n_elem]) {
    std::copy(values, values + indices_.n_elem, values_.get());
  }

  //! Create a proxy to the slope coefficients.
  //!
  //! @param size number of slope coefficients.
  //! @param values values of the slope coefficients. Must have the same number of elements as *indices*.
  //! @param indices the indicies of the slope coefficients in *values*.
  BetaProxy(const arma::uword size, double const * values, arma::uvec&& indices) noexcept
    : size_(size), indices_(std::move(indices)), values_(new double[indices_.n_elem]) {
    std::copy(values, values + indices_.n_elem, values_.get());
  }

  BetaProxy(const BetaProxy& other) = delete;
  BetaProxy& operator=(const BetaProxy& other) = delete;

  BetaProxy(BetaProxy&& other) = default;
  BetaProxy& operator=(BetaProxy&& other) = default;

  //! Get a dense vector representation of the slope coefficients.
  arma::vec beta(std::false_type) {
    arma::vec beta(size_, arma::fill::zeros);
    double * value_ptr = values_.get();
    for (arma::uword i = 0; i < indices_.n_elem; ++i) {
      beta[indices_[i]] = *value_ptr++;
    }
    return beta;
  }

  //! Get a sparse vector representation of the slope coefficients.
  arma::sp_vec beta(std::true_type) {
    arma::sp_vec beta(size_);
    double * value_ptr = values_.get();
    for (arma::uword i = 0; i < indices_.n_elem; ++i) {
      beta[indices_[i]] = *value_ptr++;
    }
    return beta;
  }

 private:
  arma::uword size_;
  arma::uvec indices_;
  std::unique_ptr<double[]> values_;
};

//! A temporary beta proxy not to be used except as an rvalue.
class TemporaryBetaProxy {
 public:
  TemporaryBetaProxy(const TemporaryBetaProxy&) = delete;
  TemporaryBetaProxy& operator=(const TemporaryBetaProxy&) = delete;

  operator BetaProxy() {
    return BetaProxy(size_, values_, std::move(indices_));
  }

  //! Get a dense vector representation of the slope coefficients.
  arma::vec beta(std::false_type) {
    arma::vec beta(size_, arma::fill::zeros);
    double const * value_ptr = values_;
    for (arma::uword i = 0; i < indices_.n_elem; ++i) {
      beta[indices_[i]] = *value_ptr++;
    }
    return beta;
  }

  //! Get a sparse vector representation of the slope coefficients.
  arma::sp_vec beta(std::true_type) {
    arma::sp_vec beta(size_);
    double const * value_ptr = values_;
    for (arma::uword i = 0; i < indices_.n_elem; ++i) {
      beta[indices_[i]] = *value_ptr++;
    }
    return beta;
  }

 private:
  TemporaryBetaProxy(const arma::uword size, double const * values, const arma::uvec& indices) noexcept
    : size_(size), indices_(indices), values_(values) {}

  TemporaryBetaProxy(TemporaryBetaProxy&&) = default;
  TemporaryBetaProxy& operator=(TemporaryBetaProxy&&) = default;

  arma::uword size_;
  arma::uvec indices_;
  double const * values_;

  friend class LarsPath;
};

//! The LARS regularization path.
class LarsPath {
 public:
  LarsPath(const arma::mat& gram, const arma::vec& _cor_y_, const arma::uword max_active) noexcept
    : chol_(gram, max_active), cor_y_(_cor_y_), max_cor_(arma::norm(cor_y_, "inf")),
      cor_signs_(max_active), active_beta_(new double[max_active]), max_active_(max_active),
      remaining_usable_vars_(gram.n_cols), dropped_(false) {
      auto insert_it = inactive_.before_begin();
      for (arma::uword ind = 0; ind < gram.n_cols; ++ind) {
        insert_it = inactive_.insert_after(insert_it, ind);
      }
    }

  LarsPath(const LarsPath& other) noexcept
    : chol_(other.chol_), cor_y_(other.cor_y_), max_cor_(other.max_cor_),
      cor_signs_(other.cor_signs_), inactive_(other.inactive_), active_beta_(new double[other.max_active_]),
      max_active_(other.max_active_), remaining_usable_vars_(other.remaining_usable_vars_), dropped_(other.dropped_) {
    std::copy(other.active_beta_.get(), other.active_beta_.get() + max_active_, active_beta_.get());
  }

  LarsPath& operator=(const LarsPath& other) noexcept {
    chol_ = other.chol_;
    cor_y_ = other.cor_y_;
    max_cor_ = other.max_cor_;
    remaining_usable_vars_ = other.remaining_usable_vars_;
    cor_signs_ = other.cor_signs_;
    active_beta_.reset(new double[other.max_active_]);
    max_active_ = other.max_active_;
    std::copy(other.active_beta_.get(), other.active_beta_.get() + max_active_, active_beta_.get());
    inactive_ = other.inactive_;
    dropped_ = other.dropped_;
    return *this;
  }

  LarsPath(LarsPath&&) = default;
  LarsPath& operator=(LarsPath&&) = default;

  //! Update the diagonal of the Gram matrix.
  //! No other manipulations are performed, hence, the path is invalid afterwards and needs to be reset with `Reset`.
  //!
  //! @param add_gram_diagonal add this number to each diagonal element of the Gram matrix.
  void UpdateGram(const double add_gram_diagonal) noexcept {
    chol_.UpdateMatrixDiagonal(add_gram_diagonal);
  }

  //! Update the diagonal of the Gram matrix.
  //! No other manipulations are performed, hence, the path is invalid afterwards and needs to be reset with `Reset`.
  //!
  //! @param add_gram_diagonal add this vector to the diagonal of the Gram matrix.
  void UpdateGram(const arma::vec& add_gram_diagonal) {
    chol_.UpdateMatrixDiagonal(add_gram_diagonal);
  }

  //! Reset the lars path to the beginning, i.e., no active variables, but use the same Gram matrix.
  //!
  //! @param cor_y the new correlation between X and y.
  void Reset(const arma::vec& cor_y) noexcept {
    chol_.Reset();
    cor_y_ = cor_y;
    max_cor_ = arma::norm(cor_y_, "inf");
    dropped_ = false;
    inactive_.clear();
    remaining_usable_vars_ = chol_.matrix().n_cols;
    auto insert_it = inactive_.before_begin();
    for (arma::uword ind = 0; ind < remaining_usable_vars_; ++ind) {
      insert_it = inactive_.insert_after(insert_it, ind);
    }
  }

  //! Take next step along the LARS path.
  void Next() {
    // If no predictors were dropped previously, find next predictor(s) to activate
    // Compute step direction
    // Find correlations between direction and active/inactive variables
    // Compute step size
    // Check if any predictors are dropped along the step.
    // if so, only take step up until the predictor(s) are dropped
    if (!dropped_) {
      ActivateNext();
    }

    // Compute the step direction.
    arma::vec step_dir = cor_signs_.head(chol_.active_size());
    chol_.Solve(&step_dir);
    // Correlation of active variable with equiangular vector:
    const double cor_active_equi = 1. / std::sqrt(arma::dot(step_dir, cor_signs_.head(chol_.active_size())));
    step_dir *= cor_active_equi;

    // Compute correlations between the equiangular vector and all other variables
    const arma::vec cor_equi_v = chol_.matrix().cols(chol_.active()) * step_dir;

    // Compute step size in the direction of the equiangular vector.
    double step = FindStepSize(cor_active_equi, cor_equi_v);

    // Check if any predictors need to be dropped before taking the step!
    dropped_ = DropAlong(&step_dir, &step);

    // Compute new coefficients
    arma::vec(active_beta_.get(), chol_.active_size(), false, true) += step * step_dir;
    cor_y_ -= step * cor_equi_v;

    // Find new maximum correlation
    max_cor_ = std::numeric_limits<double>::epsilon();
    for (const arma::uword inactive_pred : inactive_) {
      const double tmp = std::abs(cor_y_[inactive_pred]);
      if (max_cor_ < tmp) {
        max_cor_ = tmp;
      }
    }
  }

  TemporaryBetaProxy CurrentSlope() const noexcept {
    return TemporaryBetaProxy(chol_.matrix().n_cols, active_beta_.get(), chol_.active());
  }

  arma::uword max_active() const noexcept {
    return max_active_;
  }

  arma::uword active_size() const noexcept {
    return chol_.active_size();
  }

  double max_cor() const noexcept {
    return max_cor_;
  }

  //! Get the correlation between the predictor and the response.
  double cor_y(const arma::uword predictor) const noexcept {
    return cor_y_[predictor];
  }

  //! Get the diagonal element of the gram matrix associated with *predictor*.
  double GramDiagonal(const arma::uword predictor) const noexcept {
    return chol_.matrix().at(predictor, predictor);
  }

 private:
  //! Activate the next inactive variable(s) with largest correlation with the response.
  void ActivateNext() {
    auto inactive_it = inactive_.begin();
    auto inactive_drop_it = inactive_.before_begin();

    while (inactive_it != inactive_.end()) {
      const arma::uword inactive_pred = *inactive_it++;
      if (max_cor_ <= std::abs(cor_y_[inactive_pred]) + std::numeric_limits<double>::epsilon()) {
        // This currently inactive predictor has maximum correlation with the response. Add it.
        // First add it to the decomposition.
        const bool not_singular = chol_.Add(inactive_pred);
        if (not_singular) {
          // No singularity. Add the variable to the active set.
          const arma::uword new_active_index = chol_.active_size() - 1;
          cor_signs_[new_active_index] = cor_y_[inactive_pred] < 0 ? -1. : 1.;
          active_beta_[new_active_index] = 0;
        } else {
          // Singularity detected! Drop variable for good.
          --remaining_usable_vars_;
          if (remaining_usable_vars_ < max_active_) {
            --max_active_;
          }
        }
        // Remove element from inactive set. Ensure that `inactive_it` is incremented before erasing it!
        inactive_.erase_after(inactive_drop_it);
      } else {
        ++inactive_drop_it;
      }
    }
  }

  //! Compute step size in the direction of the equiangular vector.
  //!
  //! @param cor_active_equi correlation of the newly active predictor with the equiangular vector
  //! @param cor_with_equi correlations of all predictors with the equiangular vector
  //! @return step size
  double FindStepSize(const double cor_active_equi, const arma::vec& cor_with_equi) {
    double step = max_cor_ / cor_active_equi;

    if (chol_.active_size() == max_active_) {
      return step;
    }

    for (const arma::uword inactive_pred : inactive_) {
      const double inactive_cor_with_y = cor_y_[inactive_pred];
      const double inactive_cor_with_equiv = cor_with_equi[inactive_pred];

      double possible_step = (max_cor_ - inactive_cor_with_y) / (cor_active_equi - inactive_cor_with_equiv);
      if (possible_step > std::numeric_limits<double>::epsilon() && possible_step < step) {
        step = possible_step;
      }
      possible_step = (max_cor_ + inactive_cor_with_y) / (cor_active_equi + inactive_cor_with_equiv);
      if (possible_step > std::numeric_limits<double>::epsilon() && possible_step < step) {
        step = possible_step;
      }
    }
    return step;
  }

  //! Drop predictors along the given direction.
  //!
  //! @param direction direction of the next step for all currently active predictors. On output, it will only have
  //!                  as many elements for the retained predictors (if any are dropped).
  //! @param step the desired step size. On output will be updated with the step size until drops happen.
  //! @return true if variables are dropped and the step size & direciton were adjusted. False otherwise.
  bool DropAlong(arma::vec* direction, double* step) {
    std::forward_list<arma::uword> drops;
    for (arma::uword i = 0; i < chol_.active_size(); ++i) {
      // Find smallest *positive* step size which results in a change of signs and is smaller than *step*.
      const double sign_change_step = -active_beta_[i] / direction->at(i);
      if (sign_change_step > 0 && sign_change_step < *step) {
        if (*step > sign_change_step + std::numeric_limits<double>::epsilon()) {
          // Previous smallest positive step is "much" larger than this step. Reset drops.
          drops.clear();
        }
        *step = sign_change_step;
        drops.push_front(i);
      }
    }

    // Drop variables if necessary. drops is in descending order!
    if (!drops.empty()) {
      arma::uword prev_active_size = chol_.active_size();
      const arma::uvec prev_active = chol_.active();
      chol_.Drop(drops.begin(), drops.end());

      // Add dropped variables back to the inactive set and remove traces from the dropped predictor.
      // This implementation requires that the drops are in descending order! It is efficient only if the number of
      // drops is small (which is usually the case)!
      for (const arma::uword drop_ind : drops) {
        std::copy(&active_beta_[drop_ind + 1], &active_beta_[prev_active_size], &active_beta_[drop_ind]);
        inactive_.push_front(prev_active[drop_ind]);
        direction->shed_row(drop_ind);
        // Don't use `shed_row` on the cor_signs_ vector because we don't want a size change!
        std::copy(cor_signs_.memptr() + drop_ind + 1, cor_signs_.memptr() + prev_active_size,
                  cor_signs_.memptr() + drop_ind);
        --prev_active_size;
      }
      return true;
    }
    return false;
  }

  linalg::Cholesky chol_;
  arma::vec cor_y_;
  double max_cor_ = 0;
  arma::vec cor_signs_;
  std::forward_list<arma::uword> inactive_;
  std::unique_ptr<double[]> active_beta_;
  arma::uword max_active_;
  arma::uword remaining_usable_vars_;
  bool dropped_;
};
}  // namespace auglars

//! Compute the EN regression estimate using the LARS algorithm on the augmented response vector and predictor
//! matrix.
template<class LossFunction, class PenaltyFunction, class Coefficients>
class AugmentedLarsOptimizer : public Optimizer<LossFunction, PenaltyFunction, Coefficients> {
  using Base = Optimizer<LossFunction, PenaltyFunction, Coefficients>;
  using LossFunctionPtr = std::unique_ptr<LossFunction>;
  using PenaltyPtr = std::unique_ptr<PenaltyFunction>;
  using IsWeightedTag = typename traits::is_weighted<LossFunction>::type;
  using IsAdaptiveTag = typename traits::is_adaptive<PenaltyFunction>::type;
  using IsSparseTag = typename std::is_same<typename Coefficients::SlopeCoefficient, arma::sp_vec>::type;
  using Weights = typename std::conditional<IsWeightedTag::value, arma::vec, double>::type;

  static_assert(traits::is_en_penalty<PenaltyFunction>::value, "PenaltyFunction must be an EN-type penalty.");
  static_assert(traits::is_ls_regression_loss<LossFunction>::value, "LossFunction must be an least-squares-type loss.");

 public:
  using Optimum = typename Base::Optimum;

  //! Ininitialize the optimizer without a loss or penalty function.
  AugmentedLarsOptimizer() noexcept {}

  //! Ininitialize the optimizer using the given (weighted) LS loss function and penalty function.
  //!
  //! @param loss a weighted LS loss function.
  //! @param penalty penalty function.
  AugmentedLarsOptimizer(const LossFunction& loss, const PenaltyFunction& penalty) noexcept
    : loss_(new LossFunction(loss)), penalty_(new PenaltyFunction(penalty)) {}

  //! Default copy constructor.
  //!
  //! The copied optimizer will share the identical loss and penalty functions after construction.
  AugmentedLarsOptimizer(const AugmentedLarsOptimizer& other) noexcept
    : loss_(other.loss_? new LossFunction(*other.loss_) : nullptr),
      penalty_(other.penalty_ ? new PenaltyFunction(*other.penalty_) : nullptr),
      path_(other.path_ ? new auglars::LarsPath(*other.path_) : nullptr), mean_x_(other.mean_x_),
      mean_y_(other.mean_y_) {}

  //! Default copy assignment.
  //!
  //! The copied optimizer will share the identical loss and penalty functions after construction.
  AugmentedLarsOptimizer& operator=(const AugmentedLarsOptimizer& other) = default;

  //! Default move constructor.
  AugmentedLarsOptimizer(AugmentedLarsOptimizer&& other) = default;

  //! Default move assignment operator.
  AugmentedLarsOptimizer& operator=(AugmentedLarsOptimizer&& other) = default;

  ~AugmentedLarsOptimizer() = default;

  void Reset() {
    loss_.reset();
    penalty_.reset();
    path_.reset();
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
    path_.reset();
    loss_.reset(new LossFunction(loss));
  }

  PenaltyFunction& penalty() const {
    if (!penalty_) {
      throw std::logic_error("no penalty set");
    }
    return *penalty_;
  }

  void penalty(const PenaltyFunction& penalty) noexcept {
    if (penalty_ && loss_ && path_) {
      path_->UpdateGram(LambdaRidge(penalty, IsWeightedTag{}, IsAdaptiveTag{}) -
                          LambdaRidge(*penalty_, IsWeightedTag{}, IsAdaptiveTag{}));
    }
    penalty_.reset(new PenaltyFunction(penalty));
  }

  //! Find the minimum of the objective function, using the previous solution (or the 0-vector if no
  //! previous solution exists) as starting point and at most ``max_it`` iterations.
  //!
  //! @param max_it maximum number of iterations.
  //! @return information about the optimum.
  Optimum Optimize(const int /* max_it */) {
    return Optimize();
  }

  //! Find the minimum of the objective function, using the given coefficients as starting point
  //! and at most ``max_it`` iterations.
  //!
  //! @param start where to start the optimization from.
  //! @return information about the optimum.
  Optimum Optimize(const Coefficients& /* start */) {
    return Optimize();
  }

  //! Find the minimum of the objective function, using the given coefficients as starting point
  //! and at most ``max_it`` iterations.
  //!
  //! @param start where to start the optimization from.
  //! @param max_it maximum number of iterations.
  //! @return information about the optimum.
  Optimum Optimize(const Coefficients& /* start */, const int /* max_it */) {
    return Optimize();
  }

  //! Find the minimum of the objective function, using the previous solution (or the 0-vector if no
  //! previous solution exists) as starting point.
  //!
  //! @return information about the optimum.
  Optimum Optimize() {
    if (!loss_) {
      throw std::logic_error("no loss set");
    }
    if (!penalty_) {
      throw std::logic_error("no penalty set");
    }

    InitializeLarsPath(IsWeightedTag{}, IsAdaptiveTag{});

    auto&& data = loss_->data();
    if (data.n_pred() == 1) {
      const Coefficients coefs = OptimizeSinglePredictor(IsWeightedTag{}, IsAdaptiveTag{});
      return MakeOptimum(*loss_, *penalty_, coefs);
    }

    const double lambda_lasso = LambdaLasso(*penalty_, IsWeightedTag{});
    double prev_max_cor = path_->max_cor();
    auglars::BetaProxy prev_beta(data.n_pred());

    // Walk along the LARS path until all predictors are added or the LASSO lambda is passed.
    while (path_->active_size() < path_->max_active() && path_->max_cor() > lambda_lasso &&
           path_->max_cor() <= prev_max_cor + std::numeric_limits<double>::epsilon()) {
      prev_beta = path_->CurrentSlope();
      prev_max_cor = path_->max_cor();
      path_->Next();
    }

    // Either the maximum number of predictors are added, or the maximum correlation is below the desired LASSO lambda.
    Coefficients coefs(path_->CurrentSlope().beta(IsSparseTag{}));
    if (path_->active_size() > 0) {
      if (path_->max_cor() < lambda_lasso) {
        // The requested LASSO lambda is somewhere between the previous step and the current.
        // This means the solution is the linear interpolation between the previous and the current slope.
        // If the current coefficients are the LS solution (or as close as it gets),
        // the interpolation is a bit different.
        const double mixing = (path_->active_size() == path_->max_active()) ?
          (lambda_lasso / prev_max_cor) : ((path_->max_cor() - lambda_lasso) / (path_->max_cor() - prev_max_cor));
        coefs.beta = mixing * prev_beta.beta(IsSparseTag{}) + (1 - mixing) * coefs.beta;
      }
    }

    const arma::vec residuals = FinalizeCoefficients(&coefs, IsWeightedTag{}, IsAdaptiveTag{});
    if (path_->max_cor() > prev_max_cor + std::numeric_limits<double>::epsilon()) {
      return MakeOptimum(*loss_, *penalty_, coefs, residuals,
                         OptimumStatus::kWarning, "Penalization level below numerical precision.");
    }
    return MakeOptimum(*loss_, *penalty_, coefs, residuals);
  }

 private:
  arma::vec FinalizeCoefficients(Coefficients* coefs, std::true_type /* is_weighted */,
                                 std::true_type /* is_adaptive */) const {
    auto&& data = loss_->data();
    coefs->beta /= penalty_->loadings();

    const arma::vec slope_prod = data.cx() * coefs->beta;
    coefs->intercept = loss_->IncludeIntercept() ?
      mean_y_ - arma::dot(arma::square(loss_->sqrt_weights()), slope_prod) / data.n_obs() : 0;
    return data.cy() - slope_prod - coefs->intercept;
  }

  arma::vec FinalizeCoefficients(Coefficients* coefs, std::true_type /* is_weighted */,
                                 std::false_type /* is_adaptive */) {
    auto&& data = loss_->data();
    const arma::vec slope_prod = data.cx() * coefs->beta;
    coefs->intercept = loss_->IncludeIntercept() ?
      mean_y_ - arma::dot(arma::square(loss_->sqrt_weights()), slope_prod) / data.n_obs() : 0;
    return data.cy() - slope_prod - coefs->intercept;
  }

  arma::vec FinalizeCoefficients(Coefficients* coefs, std::false_type /* is_weighted */,
                                 std::true_type /* is_adaptive */) {
    auto&& data = loss_->data();
    coefs->intercept = loss_->IncludeIntercept() ? mean_y_ - arma::as_scalar(mean_x_ * coefs->beta) : 0;

    coefs->beta /= penalty_->loadings();
    return data.cy() - data.cx() * coefs->beta - coefs->intercept;
  }

  arma::vec FinalizeCoefficients(Coefficients* coefs, std::false_type /* is_weighted */,
                                 std::false_type /* is_adaptive */) {
    auto&& data = loss_->data();
    coefs->intercept = loss_->IncludeIntercept() ? mean_y_ - arma::as_scalar(mean_x_ * coefs->beta) : 0;
    return data.cy() - data.cx() * coefs->beta - coefs->intercept;
  }

  //! Special case of only a single predictor with sufficiently large norm.
  //! Compute coefficients of the least squares solution.
  Coefficients OptimizeSinglePredictor(std::true_type /* is_weighted */, std::true_type /* is_adaptive */) {
    const double lambda_lasso = LambdaLasso(*penalty_, IsWeightedTag{});
    Coefficients coefs(mean_y_, typename Coefficients::SlopeCoefficient(1));

    // Compute the LASSO coefficients through linear interpolation
    if (lambda_lasso < path_->cor_y(0)) {
      auto&& data = loss_->data();
      // Back-transform scaling from penalty loadings *before* the intercept is computed.
      coefs.beta[0] = (path_->cor_y(0) - lambda_lasso) / (path_->GramDiagonal(0) * penalty_->loadings()[0]);
      // sqrt_weights^2 sum to n!
      coefs.intercept -= coefs.beta[0] * arma::dot(arma::square(loss_->sqrt_weights()), data.cx().col(0)) /
        data.n_obs();
    } else {
      coefs.beta.zeros();
    }
    return coefs;
  }

  //! Special case of only a single predictor with sufficiently large norm.
  //! Compute coefficients of the least squares solution.
  Coefficients OptimizeSinglePredictor(std::true_type /* is_weighted */, std::false_type /* is_adaptive */) {
    const double lambda_lasso = LambdaLasso(*penalty_, IsWeightedTag{});
    Coefficients coefs(mean_y_, typename Coefficients::SlopeCoefficient(1));

    // Compute the LASSO coefficients through linear interpolation
    if (lambda_lasso < path_->cor_y(0)) {
      auto&& data = loss_->data();
      coefs.beta[0] = (path_->cor_y(0) - lambda_lasso) / path_->GramDiagonal(0);
      // sqrt_weights^2 sum to n!
      coefs.intercept -= coefs.beta[0] * arma::dot(arma::square(loss_->sqrt_weights()), data.cx().col(0)) /
        data.n_obs();
    } else {
      coefs.beta.zeros();
    }
    return coefs;
  }

  //! Special case of only a single predictor with sufficiently large norm.
  //! Compute coefficients of the least squares solution.
  Coefficients OptimizeSinglePredictor(std::false_type /* is_weighted */, std::true_type /* is_adaptive */) {
    const double lambda_lasso = LambdaLasso(*penalty_, IsWeightedTag{});
    Coefficients coefs(mean_y_, typename Coefficients::SlopeCoefficient(1));

    // Compute the LASSO coefficients through linear interpolation
    if (lambda_lasso < path_->cor_y(0)) {
      coefs.beta[0] = (path_->cor_y(0) - lambda_lasso) / path_->GramDiagonal(0);
      coefs.intercept -= mean_x_[0] * coefs.beta[0];
      // Back-transform scaling from penalty loadings *after* the intercept is computed.
      coefs.beta[0] /= penalty_->loadings()[0];
    } else {
      coefs.beta.zeros();
    }
    return coefs;
  }

  //! Special case of only a single predictor with sufficiently large norm.
  //! Compute coefficients of the least squares solution.
  Coefficients OptimizeSinglePredictor(std::false_type /* is_weighted */, std::false_type /* is_adaptive */) {
    const double lambda_lasso = LambdaLasso(*penalty_, IsWeightedTag{});
    Coefficients coefs(mean_y_, typename Coefficients::SlopeCoefficient(1) );

    // Compute the LASSO coefficients through linear interpolation
    if (lambda_lasso < path_->cor_y(0)) {
      coefs.beta[0] = (path_->cor_y(0) - lambda_lasso) / path_->GramDiagonal(0);
      coefs.intercept -= mean_x_[0] * coefs.beta[0];
    } else {
      coefs.beta.zeros();
    }
    return coefs;
  }

  void InitializeLarsPath(std::true_type /* is_weighted */, std::true_type /* is_adaptive */) {
    auto&& data = loss_->data();
    if (!path_) {
      // Compute the gram matrix for the new data.
      const arma::uword max_active = (penalty_->alpha() < 1) ? data.n_pred() : std::min(data.n_obs(), data.n_pred());
      // TODO(me): does it matter that the weights are normalized??
      const arma::vec weighted_y = data.cy() % loss_->sqrt_weights();
      arma::mat weighted_x = data.cx().each_row() / penalty_->loadings().t();
      // sqrt_weights^2 sum to n!
      const double inv_n = 1. / data.n_obs();

      if (loss_->IncludeIntercept()) {
        mean_x_ = arma::mean(weighted_x);
        mean_y_ = arma::dot(weighted_y, loss_->sqrt_weights()) * inv_n;
        weighted_x.each_row() -= mean_x_;
        // Weight the predictors ...
        weighted_x.each_col() %= loss_->sqrt_weights();
        // ... and make the design orthogonal to y - intercept
        weighted_x -= inv_n * loss_->sqrt_weights() * loss_->sqrt_weights().t() * weighted_x;
      } else {
        mean_x_.reset();
        mean_y_ = 0;
        // Only weight the predictors.
        weighted_x.each_col() %= loss_->sqrt_weights();
      }

      path_.reset(new auglars::LarsPath(weighted_x.t() * weighted_x, weighted_x.t() * weighted_y, max_active));
      path_->UpdateGram(LambdaRidge(*penalty_, IsWeightedTag{}, IsAdaptiveTag{}));
    } else {
      // Re-use the gram matrix from the current path.
      arma::vec weighted_y = data.cy() % arma::square(loss_->sqrt_weights());
      if (loss_->IncludeIntercept()) {
        const double weighted_y_mean = arma::mean(weighted_y);
        weighted_y -= weighted_y_mean * arma::square(loss_->sqrt_weights());
      }
      const arma::vec cor_y = data.cx().t() * weighted_y;
      path_->Reset(cor_y / penalty_->loadings());
    }
  }

  void InitializeLarsPath(std::true_type /* is_weighted */, std::false_type /* is_adaptive */) {
    auto&& data = loss_->data();
    if (!path_) {
      // Compute the gram matrix for the new data.
      const arma::uword max_active = (penalty_->alpha() < 1) ? data.n_pred() : std::min(data.n_obs(), data.n_pred());
      const arma::vec weighted_y = data.cy() % loss_->sqrt_weights();
      // sqrt_weights^2 sum to n!
      const double inv_n = 1. / data.n_obs();

      arma::mat weighted_x;
      if (loss_->IncludeIntercept()) {
        mean_x_ = arma::mean(data.cx());
        mean_y_ = arma::dot(weighted_y, loss_->sqrt_weights()) * inv_n;
        weighted_x = data.cx().each_row() - mean_x_;
        // Weight the predictors...
        weighted_x.each_col() %= loss_->sqrt_weights();
        // ... and make the design orthogonal to y - intercept
        weighted_x -= inv_n * loss_->sqrt_weights() * loss_->sqrt_weights().t() * weighted_x;
      } else {
        mean_x_.reset();
        mean_y_ = 0;
        // Only weight the predictors.
        weighted_x = data.cx().each_col() % loss_->sqrt_weights();
      }

      path_.reset(new auglars::LarsPath(weighted_x.t() * weighted_x, weighted_x.t() * weighted_y, max_active));
      path_->UpdateGram(LambdaRidge(*penalty_, IsWeightedTag{}, IsAdaptiveTag{}));
    } else {
      // Re-use the gram matrix from the current path.
      arma::vec weighted_y = data.cy() % arma::square(loss_->sqrt_weights());
      if (loss_->IncludeIntercept()) {
        const double weighted_y_mean = arma::mean(weighted_y);
        weighted_y -= weighted_y_mean * arma::square(loss_->sqrt_weights());
      }
      const arma::vec cor_y = data.cx().t() * weighted_y;
      path_->Reset(cor_y);
    }
  }

  void InitializeLarsPath(std::false_type /* is_weighted */, std::true_type /* is_adaptive */) {
    auto&& data = loss_->data();
    if (!path_) {
      const arma::uword max_active = (penalty_->alpha() < 1) ? data.n_pred() : std::min(data.n_obs(), data.n_pred());
      arma::mat weighted_x = data.cx().each_row() / penalty_->loadings().t();

      if (loss_->IncludeIntercept()) {
        mean_x_ = arma::mean(weighted_x);
        mean_y_ = arma::mean(data.cy());
        weighted_x.each_row() -= mean_x_;
      } else {
        mean_x_.reset();
        mean_y_ = 0;
      }

      path_.reset(new auglars::LarsPath(weighted_x.t() * weighted_x, weighted_x.t() * data.cy(), max_active));
      path_->UpdateGram(LambdaRidge(*penalty_, IsWeightedTag{}, IsAdaptiveTag{}));
    } else {
      // Re-use the gram matrix from the current path.
      const arma::vec cor_y = (data.cx().t() * data.cy()) / penalty_->loadings();
      if (loss_->IncludeIntercept()) {
        path_->Reset(cor_y - data.n_obs() * mean_y_ * mean_x_.t());
      } else {
        path_->Reset(cor_y);
      }
    }
  }

  void InitializeLarsPath(std::false_type /* is_weighted */, std::false_type /* is_adaptive */) {
    auto&& data = loss_->data();
    if (!path_) {
      // Compute the gram matrix for the new data.
      const arma::uword max_active = (penalty_->alpha() < 1) ? data.n_pred() : std::min(data.n_obs(), data.n_pred());

      if (loss_->IncludeIntercept()) {
        mean_x_ = arma::mean(data.cx());
        mean_y_ = arma::mean(data.cy());
        const arma::mat centered_x = data.cx().each_row() - mean_x_;
        path_.reset(new auglars::LarsPath(centered_x.t() * centered_x, centered_x.t() * data.cy(), max_active));
      } else {
        mean_x_.reset();
        mean_y_ = 0;
        path_.reset(new auglars::LarsPath(data.cx().t() * data.cx(), data.cx().t() * data.cy(), max_active));
      }

      path_->UpdateGram(LambdaRidge(*penalty_, IsWeightedTag{}, IsAdaptiveTag{}));
    } else {
      // Re-use the gram matrix from the current path.
      const arma::vec cor_y = data.cx().t() * data.cy();
      if (loss_->IncludeIntercept()) {
        path_->Reset(cor_y - data.n_obs() * mean_y_ * mean_x_.t());
      } else {
        path_->Reset(cor_y);
      }
    }
  }

  //! Compute the "Ridge" part of the penalty.
  //! Note: The formulation of the LS-loss in the *nsoptim* package includes the factor ``1/(2n)`` in front.
  //! To make that compatible with the LARS formulation used by this optimizer, the Ridge penalty term must be
  //! multiplied by ``2n``.
  inline double LambdaRidge(const PenaltyFunction& penalty, std::false_type /* is_weighted */,
                            std::false_type /* is_adaptive */) const noexcept {
    return loss_->data().n_obs() * (1 - penalty.alpha()) * penalty.lambda();
  }

  //! Compute the "Ridge" part of the penalty in case of unweighted observations and penalty loadings.
  //! See notes for the implementation for unweighted observations.
  inline double LambdaRidge(const PenaltyFunction& penalty, std::true_type /* is_weighted */,
                            std::false_type /* is_adaptive */) const noexcept {
    return LambdaRidge(penalty, std::false_type{}, std::false_type{}) / loss_->mean_weight();
  }

  //! Compute the "Ridge" part of the penalty.
  //! See notes for the implementation for unweighted observations.
  inline arma::vec LambdaRidge(const PenaltyFunction& penalty, std::false_type /* is_weighted */,
                               std::true_type /* is_adaptive */) const noexcept {
    return LambdaRidge(penalty, std::false_type{}, std::false_type{}) / penalty_->loadings();
  }

  //! Compute the "Ridge" part of the penalty in case of weighted observations and penalty loadings.
  //! See notes for the implementation for unweighted observations.
  inline arma::vec LambdaRidge(const PenaltyFunction& penalty, std::true_type /* is_weighted */,
                            std::true_type /* is_adaptive */) const noexcept {
    return LambdaRidge(penalty, std::true_type{}, std::false_type{}) / penalty_->loadings();
  }

  //! Compute the "Lasso" part of the penalty in case of weighted observations.
  //! Note: The formulation of the LS-loss in the *nsoptim* package includes the factor ``1/(2n)`` in front.
  //! To make that compatible with the LARS formulation used by this optimizer, the Lasso penalty term must be
  //! multiplied by ``n``.
  inline double LambdaLasso(const PenaltyFunction& penalty, std::false_type /* is_weighted */) noexcept {
    return loss_->data().n_obs() * penalty.alpha() * penalty.lambda();
  }

  //! Compute the "Lasso" part of the penalty in case of weighted observations.
  //! See notes for the implementation for unweighted observations.
  inline double LambdaLasso(const PenaltyFunction& penalty, std::true_type /* is_weighted */) noexcept {
    return LambdaLasso(penalty, std::false_type{}) / loss_->mean_weight();
  }

  LossFunctionPtr loss_;
  PenaltyPtr penalty_;
  std::unique_ptr<auglars::LarsPath> path_;

  arma::rowvec mean_x_;
  double mean_y_;
};

//! Specialization of the LARS algorithm for Ridge penalty.
//! The Ridge regression estimate is computed using standard linear algebra on the augmented response vector and
//! predictor matrix.
template<class LossFunction>
class AugmentedLarsOptimizer<LossFunction, RidgePenalty, RegressionCoefficients<arma::vec>>
  : public Optimizer<LossFunction, RidgePenalty, RegressionCoefficients<arma::vec>> {
  using AugmentedRidgeOptimizer = AugmentedLarsOptimizer<LossFunction, RidgePenalty, RegressionCoefficients<arma::vec>>;
  using Base = Optimizer<LossFunction, RidgePenalty, RegressionCoefficients<arma::vec>>;
  using LossFunctionPtr = std::unique_ptr<LossFunction>;
  using RidgePenaltyPtr = std::unique_ptr<RidgePenalty>;
  static_assert(traits::is_ls_regression_loss<LossFunction>::value, "LossFunction must be a LS-type loss.");

 public:
  using Coefficients = typename Base::Coefficients;
  using Optimum = typename Base::Optimum;

  //! Ininitialize the optimizer using the given (weighted) LS loss function and the Ridge penalty.
  //!
  //! @param loss a weighted LS loss function.
  //! @param penalty the Ridge penalty.
  AugmentedLarsOptimizer<LossFunction, RidgePenalty, RegressionCoefficients<arma::vec>>() noexcept
    : previous_data_id_(ObjectId::null()), loss_(nullptr), penalty_(nullptr) {}

  //! Ininitialize the optimizer using the given (weighted) LS loss function and the Ridge penalty.
  //!
  //! @param loss a weighted LS loss function.
  //! @param penalty the Ridge penalty.
  AugmentedLarsOptimizer<LossFunction, RidgePenalty, RegressionCoefficients<arma::vec>>(
    const LossFunction& loss, const RidgePenalty& penalty) noexcept
      : previous_data_id_(ObjectId::null()), loss_(LossFunctionPtr(new LossFunction(loss))),
        penalty_(RidgePenaltyPtr(new RidgePenalty(penalty))) {}

  //! Default copy constructor.
  //!
  //! The copied optimizer will share the identical loss and penalty functions after construction.
  //! In case the loss or penalty function are mutated in any way, the change will affect both optimizers.
  //! If the loss/penalty function is changed on one of the optimizers (using the `loss()` or `penalty()` methods),
  //! the two optimizers will *not* share the new loss/penalty function.
  AugmentedLarsOptimizer<LossFunction, RidgePenalty, RegressionCoefficients<arma::vec>>(
    const AugmentedRidgeOptimizer& other) noexcept
      : previous_data_id_(other.previous_data_id_),
        loss_(other.loss_? LossFunctionPtr(new LossFunction(*other.loss_)) : nullptr),
        penalty_(other.penalty_ ? RidgePenaltyPtr(new RidgePenalty(*other.penalty_)) : nullptr),
        weighted_xy_cov_(other.weighted_xy_cov_), weighted_gram_(other.weighted_gram_),
        centered_x_(other.centered_x_), centered_y_(other.centered_y_) {}

  //! Default copy assignment.
  //!
  //! The copied optimizer will share the identical loss and penalty functions after construction.
  //! In case the loss or penalty function are mutated in any way, the change will affect both optimizers.
  //! If the loss/penalty function is changed on one of the optimizers (using the `loss()` or `penalty()` methods),
  //! the two optimizers will *not* share the new loss/penalty function.
  AugmentedRidgeOptimizer& operator=(const AugmentedRidgeOptimizer& other) = default;

  //! Default move constructor.
  AugmentedLarsOptimizer<LossFunction, RidgePenalty, RegressionCoefficients<arma::vec>>(
    AugmentedRidgeOptimizer&& other) = default;

  //! Default move assignment operator.
  AugmentedRidgeOptimizer& operator=(AugmentedRidgeOptimizer&& other) = default;

  ~AugmentedLarsOptimizer<LossFunction, RidgePenalty, RegressionCoefficients<arma::vec>>() = default;

  void Reset() {}

  LossFunction& loss() const {
    if (!loss_) {
      throw std::logic_error("no loss set");
    }
    return *loss_;
  }

  void loss(const LossFunction& loss) noexcept {
    loss_.reset(new LossFunction(loss));
  }

  RidgePenalty& penalty() const {
    if (!penalty_) {
      throw std::logic_error("no penalty set");
    }
    return *penalty_;
  }

  void penalty(const RidgePenalty& penalty) noexcept {
    penalty_.reset(new RidgePenalty(penalty));
  }

  AugmentedRidgeOptimizer::Optimum Optimize() {
    using IsWeightedTag = typename traits::is_weighted<LossFunction>::type;

    if (!loss_) {
      throw std::logic_error("no loss set");
    }

    if (!penalty_) {
      throw std::logic_error("no penalty set");
    }

    UpdateData(IsWeightedTag{});

    RegressionCoefficients<arma::vec> coefs;
    bool success;
    if (penalty_->lambda() > 0) {
      arma::mat gram = weighted_gram_;
      gram.diag() += DiagonalUpdate(IsWeightedTag{});
      success = arma::solve(coefs.beta, gram, weighted_xy_cov_, arma::solve_opts::likely_sympd);
    } else {
      success = arma::solve(coefs.beta, weighted_gram_, weighted_xy_cov_, arma::solve_opts::likely_sympd);
    }

    // Compute the intercept and residuals
    arma::vec residuals = loss_->data().cy() - loss_->data().cx() * coefs.beta;
    if (loss_->IncludeIntercept()) {
      coefs.intercept = ComputeIntercept(coefs, residuals, IsWeightedTag{});
      residuals -= coefs.intercept;
    }

    return MakeOptimum(*loss_, *penalty_, coefs, residuals,
                       success ? OptimumStatus::kOk : OptimumStatus::kError,
                       success ? "" : "Could not solve system of linear equations");
  }

  AugmentedRidgeOptimizer::Optimum Optimize(const RegressionCoefficients<arma::vec>&) {
    return Optimize();
  }

 private:
  void UpdateData(std::true_type) {
    const PredictorResponseData& data = loss_->data();
    if (loss_->IncludeIntercept()) {
      UpdateCenteredData();

      // The observations have different weights. Make the predictor matrix orthogonal to the centered response.
      // arma::mat weighted_x = centered_x_.each_col() % loss_->sqrt_weights();
      // weighted_x -= loss_->sqrt_weights() * loss_->sqrt_weights().t() * weighted_x / data.n_obs();

      // NOTE: Performing the computation manually saves a few operations
      arma::mat weighted_x(arma::size(centered_x_));
      auto weighted_x_it = weighted_x.begin();
      auto centered_x_it = centered_x_.begin();
      const auto sqrt_weights_end = loss_->sqrt_weights().end();
      const arma::uword ncol = centered_x_.n_cols;
      for (arma::uword col = 0; col < ncol; ++col) {
        double center = 0;
        auto centered_x_mean_it = centered_x_it;
        auto sqrt_weights_it = loss_->sqrt_weights().begin();
        while (sqrt_weights_it != sqrt_weights_end) {
          center += (*centered_x_mean_it) * (*sqrt_weights_it) * (*sqrt_weights_it);
          ++centered_x_mean_it;
          ++sqrt_weights_it;
        }
        center /= centered_x_.n_rows;
        sqrt_weights_it = loss_->sqrt_weights().begin();
        while (sqrt_weights_it != sqrt_weights_end) {
          *weighted_x_it = (*sqrt_weights_it) * (*centered_x_it - center);
          ++weighted_x_it;
          ++centered_x_it;
          ++sqrt_weights_it;
        }
      }

      weighted_gram_ = weighted_x.t() * weighted_x;
      // Apply weights to the response
      weighted_xy_cov_ = weighted_x.t() * (centered_y_ % loss_->sqrt_weights());
    } else {
      weighted_gram_ = data.cx().each_col() % loss_->sqrt_weights();
      weighted_xy_cov_ = weighted_gram_.t() * (data.cy() % loss_->sqrt_weights());
      weighted_gram_ = weighted_gram_.t() * weighted_gram_;
    }
  }

  void UpdateData(std::false_type) {
    const PredictorResponseData& data = loss_->data();
    if (loss_->IncludeIntercept()) {
      UpdateCenteredData();
      // The weights are all identical. Center the predictors and the response.
      weighted_gram_ = centered_x_.t() * centered_x_;
      weighted_xy_cov_ = data.cx().t() * centered_y_;
    } else {
      weighted_gram_ = data.cx().t() * data.cx();
      weighted_xy_cov_ = data.cx().t() * data.cy();
    }
  }

  // Update the local copy of the centered data. Only perform the update if the data is from a different address.
  void UpdateCenteredData() {
    const PredictorResponseData& data = loss_->data();
    if (previous_data_id_ != data.id()) {
      centered_x_ = data.cx().each_row() - arma::mean(data.cx(), 0);
      centered_y_ = data.cy() - arma::mean(data.cy());
      previous_data_id_ = data.id();
    }
  }

  double DiagonalUpdate(std::true_type) const {
    const PredictorResponseData& data = loss_->data();
    return penalty_->lambda() * (data.n_obs() / loss_->mean_weight());
  }

  double DiagonalUpdate(std::false_type) const {
    const PredictorResponseData& data = loss_->data();
    return data.n_obs() * penalty_->lambda();
  }

  double ComputeIntercept(const RegressionCoefficients<arma::vec>& coefs, const arma::vec& residuals,
                          std::true_type) const {
    return arma::mean(loss_->sqrt_weights() % loss_->sqrt_weights() % residuals);
  }

  double ComputeIntercept(const RegressionCoefficients<arma::vec>& coefs, const arma::vec& residuals,
                          std::false_type) const {
    return arma::mean(residuals);
  }

  ObjectId previous_data_id_;
  LossFunctionPtr loss_;
  RidgePenaltyPtr penalty_;
  arma::vec weighted_xy_cov_;
  arma::mat weighted_gram_;
  arma::mat centered_x_;
  arma::vec centered_y_;
};

}  // namespace nsoptim

#endif  // NSOPTIM_OPTIMIZER_AUGLARS_HPP_
