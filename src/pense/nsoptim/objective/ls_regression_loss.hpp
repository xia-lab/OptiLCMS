//
//  ls_regression_loss.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2018-11-30.
//  Copyright © 2018 David Kepplinger. All rights reserved.
//

#ifndef NSOPTIM_OBJECTIVE_LS_REGRESSION_LOSS_HPP_
#define NSOPTIM_OBJECTIVE_LS_REGRESSION_LOSS_HPP_

#include <algorithm>

#include "../armadillo.hpp"
#include "../utilities.hpp"
#include "convex.hpp"
#include "../container/data.hpp"
#include "../container/regression_coefficients.hpp"
#include "loss.hpp"

namespace nsoptim {
namespace ls_regression_loss {

//! Compute the minimum of the 1- and infinity norm of the weighted matrix :math:`W X`, where :math:`W` is a diagonal
//! matrix with entries `sqrt_weights`.
inline double TwoNormUpper(const arma::mat& x, const arma::vec& sqrt_weights) {
  double norm_1 = 0, norm_inf = 0;

  // Compute the inf-norm of the weighted matrix
  for (arma::uword i = 0; i < x.n_rows; ++i) {
    const double tmp = sqrt_weights[i] * arma::norm(x.row(i), 1);
    if (tmp > norm_inf) {
      norm_inf = tmp;
    }
  }

  // Compute the 1-norm of the weighted matrix
  for (arma::uword j = 0; j < x.n_cols; ++j) {
    const double tmp = arma::norm(sqrt_weights % x.col(j), 1);
    if (tmp > norm_1) {
      norm_1 = tmp;
    }
  }
  if (norm_1 < norm_inf) {
    return norm_1;
  }
  return norm_inf;
}

}  // namespace ls_regression_loss

//! The weighted least-squares loss for regression.
class WeightedLsRegressionLoss : public LossFunction<PredictorResponseData>,
                                 public ConvexFunction<WeightedLsRegressionLoss> {
 public:
  //! Type of the gradient if evaluated with coefficients of any type.
  template<typename>
  using GradientType = RegressionCoefficients<arma::vec>;
  using ResidualType = arma::vec;
  //! Type of the weights vector.
  using WeightsType = arma::vec;

  //! Tag this loss function as an LS-type loss.
  struct is_ls_regression_loss_tag {};

  //! Initialize a weighted LS regression loss.
  //!
  //! @param data A shared pointer to the constant predictor-response data.
  //! @param weighs A shared pointer to the constant vector of observataion weights.
  //!               Must be the same length as the number of observations in *data*.
  //! @param include_intercept Include an intercept term in the loss? Default is yes.
  WeightedLsRegressionLoss(std::shared_ptr<const PredictorResponseData> data, std::shared_ptr<const arma::vec> weights,
                 const bool include_intercept = true) noexcept
      : include_intercept_(include_intercept), data_(data), mean_weight_(arma::mean(*weights)),
        sqrt_weights_(std::make_shared<const arma::vec>(arma::sqrt(*weights / mean_weight_))),
        weighted_pred_norm_(-1) {}

  //! Initialize a weighted LS regression loss.
  //!
  //! @param data A shared pointer to the constant predictor-response data.
  //! @param weighs A vector of observation weights. Must be the same length as the number of observations in *data*.
  //! @param include_intercept Include an intercept term in the loss? Default is yes.
  WeightedLsRegressionLoss(std::shared_ptr<const PredictorResponseData> data, const arma::vec& weights,
                 const bool include_intercept = true) noexcept
      : WeightedLsRegressionLoss(data, std::make_shared<const arma::vec>(weights), include_intercept) {}

  //! Default copy constructor.
  //! The new object will **share pointers to the data and the weights**!
  WeightedLsRegressionLoss(const WeightedLsRegressionLoss& other) = default;
  //! Default move constructor.
  //! The new object will **share pointers to the data and the weights**!
  WeightedLsRegressionLoss(WeightedLsRegressionLoss&& other) = default;
  //! Default copy assignment operator.
  //! This object will **share pointers to the data and the weights** with *other*.
  WeightedLsRegressionLoss& operator=(const WeightedLsRegressionLoss& other) = default;
  //! Default move assignment operator.
  //! This object will inherit the **shared pointers to the data and the weights**.
  WeightedLsRegressionLoss& operator=(WeightedLsRegressionLoss&& other) = default;

  ~WeightedLsRegressionLoss() = default;

  //! Get the zero coefficients for this loss type.
  //!
  //! @return Zero coefficients of type `T`.
  template<typename T, typename = typename std::enable_if<
    std::is_same<T, RegressionCoefficients<typename T::SlopeCoefficient>>::value, void>::type>
  T ZeroCoefficients() const noexcept {
    return RegressionCoefficients<typename T::SlopeCoefficient>(data_->n_pred());
  }

  //! Evaluate the weighted LS loss function.
  //!
  //! @param where Point where to evaluate the loss function.
  //! @return loss evaluated at `where`.
  template<typename VectorType>
  double operator()(const RegressionCoefficients<VectorType>& where) const {
    return Evaluate(where);
  }

  //! Evaluate the weighted LS loss function.
  //!
  //! @param residuals residuals for evaluating the loss
  //! @return loss evaluated with the given residuals.
  double operator()(const ResidualType& residuals) const {
    return Evaluate(residuals);
  }

  //! Evaluate the weighted LS loss function.
  //!
  //! @param where Point where to evaluate the loss function.
  //! @return loss evaluated at `where`.
  template<typename VectorType>
  double Evaluate(const RegressionCoefficients<VectorType>& where) const {
    if (include_intercept_) {
      return Evaluate(data_->cy() - data_->cx() * where.beta - where.intercept);
    }
    return Evaluate(data_->cy() - data_->cx() * where.beta);
  }

  //! Evaluate the weighted LS loss function.
  //!
  //! @param residuals residuals for evaluating the loss
  //! @return loss evaluated with the given residuals.
  double Evaluate(const arma::vec& residuals) const {
    return 0.5 * mean_weight_ * arma::mean(arma::square(residuals % *sqrt_weights_));
  }

  //! Get the residuals for the LS loss function.
  //!
  //! @param where Point where to compute the residuals.
  //! @return residuals at `where`.
  template<typename VectorType>
  ResidualType Residuals(const RegressionCoefficients<VectorType>& where) const {
    if (include_intercept_) {
      return data_->cy() - data_->cx() * where.beta - where.intercept;
    }
    return data_->cy() - data_->cx() * where.beta;
  }

  //! Get the difference between two sets of regression coefficients.
  //!
  //! For the weighted LS loss, the difference is an approximation to the 2-norm of the matrix-vector product
  //!
  //! @param x Regression coefficients.
  //! @param y Regression coefficients.
  //! @return the difference between `x` and `y`.
  //!
  //! @verbatim embed:rst:leading-slashes
  ///
  /// .. math::
  ///
  ///   \| W X' (x_\beta - y_\beta) + w ( x_{\beta_0} - y_{\beta_0} )\|_2 <
  ///   | x_{\beta_0} - y_{\beta_0} | \sqrt{\sum_{i=1}^n w_i} + \| W X \| \| x_\beta - y_\beta \|_2
  ///
  /// @endverbatim
  template<typename T>
  double Difference(const RegressionCoefficients<T>& x, const RegressionCoefficients<T>& y) const {
    const double weighted_pred_norm = (weighted_pred_norm_ < 0) ?
      ls_regression_loss::TwoNormUpper(data_->cx(), *sqrt_weights_) : weighted_pred_norm_;
    return std::sqrt(sqrt_weights_->n_elem * mean_weight_) * std::abs(x.intercept - y.intercept) +
      weighted_pred_norm * arma::norm(x.beta - y.beta, 2);
  }

  //! Get the difference between two sets of regression coefficients.
  //!
  //! For the weighted LS loss, the difference is an approximation to the 2-norm of the matrix-vector product
  //!
  //! @param x Regression coefficients.
  //! @param y Regression coefficients.
  //! @return the difference between `x` and `y`.
  //!
  //! @verbatim embed:rst:leading-slashes
  ///
  /// .. math::
  ///
  ///   \| W X' (x_\beta - y_\beta) + w ( x_{\beta_0} - y_{\beta_0} )\|_2 <
  ///   | x_{\beta_0} - y_{\beta_0} | \sqrt{\sum_{i=1}^n w_i} + \| W X \| \| x_\beta - y_\beta \|_2
  ///
  /// @endverbatim
  template<typename T>
  double Difference(const RegressionCoefficients<T>& x, const RegressionCoefficients<T>& y) {
    if (weighted_pred_norm_ < 0) {
      weighted_pred_norm_ = ls_regression_loss::TwoNormUpper(data_->cx(), *sqrt_weights_);
    }
    return std::sqrt(sqrt_weights_->n_elem * mean_weight_) * std::abs(x.intercept - y.intercept) +
      weighted_pred_norm_ * arma::norm(x.beta - y.beta, 2);
  }

  //! Evaluate the gradient of the weighted LS loss at the given coefficient value.
  //!
  //! @param coefs Point where the gradient should be evaluated.
  template<typename T>
  GradientType<T> Gradient(const RegressionCoefficients<T>& coefs) const {
    if (include_intercept_) {
      const arma::vec neg_weighted_residuals =  mean_weight_ * arma::square(*sqrt_weights_) %
        (data_->cx() * coefs.beta + coefs.intercept - data_->cy());
      return GradientType<T>(arma::mean(neg_weighted_residuals),
                             arma::mean(data_->cx().each_col() % neg_weighted_residuals, 0));
    }

    const arma::vec neg_weighted_residuals =  mean_weight_ * arma::square(*sqrt_weights_) %
      (data_->cx() * coefs.beta - data_->cy());
    return GradientType<T>(0, -arma::mean(data_->cx().each_col() % neg_weighted_residuals, 0));
  }

  //! Access the data used by this weighted LS loss function.
  //!
  //! @return a constant reference to this loss' data.
  const PredictorResponseData& data() const noexcept {
    return *data_;
  }

  //! Access the un-normalized weights used by this weighted LS loss function.
  //!
  //! @return a vector of weights.
  WeightsType weights() const noexcept {
    return mean_weight_ * arma::square(*sqrt_weights_);
  }

  //! Access the normalized sqrt-weights used by this weighted LS loss function.
  //! Normalized means that `sqrt_weights[i] = sqrt(weights[i] / mean(weights))`.
  //!
  //! @return a constant reference to this loss' sqrt-weights.
  const WeightsType& sqrt_weights() const noexcept {
    return *sqrt_weights_;
  }

  //! Get the average of the weights, i.e., the normalizing constant used for `sqrt_weights()`.
  //!
  //! @return average of the raw weights vector.
  double mean_weight() const noexcept {
    return mean_weight_;
  }

  //! Check if the intercept term is included in this weighted LS loss function.
  //!
  //! @return `true` if the intercept term is included, `false` otherwise.
  bool IncludeIntercept() const noexcept {
    return include_intercept_;
  }

 private:
  bool include_intercept_;
  std::shared_ptr<const PredictorResponseData> data_;
  double mean_weight_;
  std::shared_ptr<const arma::vec> sqrt_weights_;
  double weighted_pred_norm_;
};


//! A regression loss function implementing the un-weighted least-squares loss defined as
//!  1/(2n) * sum_{i = 1}^n (y_i - \hat{\mu} - x_i' . \beta)^2
class LsRegressionLoss : public LossFunction<PredictorResponseData>,
               public ConvexFunction<LsRegressionLoss> {
 public:
  template<typename>
  using GradientType = RegressionCoefficients<arma::vec>;
  using ResidualType = arma::vec;

  //! Tag this loss function as an LS-type loss.
  struct is_ls_regression_loss_tag {};

  explicit LsRegressionLoss(std::shared_ptr<const PredictorResponseData> data, const bool include_intercept = true)
      : include_intercept_(include_intercept), data_(data), pred_norm_(-1) {}

  LsRegressionLoss(const LsRegressionLoss& other) = default;
  LsRegressionLoss(LsRegressionLoss&& other) = default;
  LsRegressionLoss& operator=(const LsRegressionLoss& other) = default;
  LsRegressionLoss& operator=(LsRegressionLoss&& other) = default;
  ~LsRegressionLoss() = default;

  //! Cast the un-weighted LS loss as weighted LS loss.
  operator WeightedLsRegressionLoss() const {
    return WeightedLsRegressionLoss(data_, arma::ones(data_->n_obs()), include_intercept_);
  }

  //! Get the zero coefficients for this loss type.
  //!
  //! @return zero coefficients.
  template<typename T, typename = typename std::enable_if<
    std::is_same<T, RegressionCoefficients<typename T::SlopeCoefficient>>::value, void>::type>
  T ZeroCoefficients() const noexcept {
    return RegressionCoefficients<typename T::SlopeCoefficient>(data_->n_pred());
  }

  //! Evaluate the LS loss function.
  //!
  //! @param where point where to evaluate the loss function.
  //! @return loss evaluated at `where`.
  template<typename T>
  double operator()(const RegressionCoefficients<T>& where) const {
    return Evaluate(where);
  }

  //! Evaluate the LS loss function.
  //!
  //! @param residuals residuals for evaluating the loss
  //! @return loss evaluated with the given residuals.
  double operator()(const arma::vec& residuals) const {
    return Evaluate(residuals);
  }

  //! Evaluate the LS loss function.
  //!
  //! @param where point where to evaluate the loss function.
  //! @return loss evaluated at `where`.
  template<typename T>
  double Evaluate(const RegressionCoefficients<T>& where) const {
    if (include_intercept_) {
      return Evaluate(data_->cy() - data_->cx() * where.beta);
    }
    return Evaluate(data_->cy() - data_->cx() * where.beta - where.intercept);
  }

  //! Evaluate the LS loss function.
  //!
  //! @param residuals residuals for evaluating the loss
  //! @return loss evaluated with the given residuals.
  double Evaluate(const arma::vec& residuals) const {
    return 0.5 * arma::mean(arma::square(residuals));
  }

  //! Get the residuals for the LS loss function.
  //!
  //! @param where Point where to compute the residuals.
  //! @return residuals at `where`.
  template<typename VectorType>
  ResidualType Residuals(const RegressionCoefficients<VectorType>& where) const {
    if (include_intercept_) {
      return data_->cy() - data_->cx() * where.beta - where.intercept;
    }
    return data_->cy() - data_->cx() * where.beta;
  }

  //! Get the difference between two sets of regression coefficients.
  //!
  //! For the LS loss, the relative difference is an approximation to the 2-norm of the matrix-vector product
  //! ||X . (beta1 - beta2) + ( mu1 - mu2 )||_2 < |mu1 - mu2| sqrt(n) + ||X|| ||beta1 - beta2||_2
  //!
  //! @param x a set of regression coefficients.
  //! @param y the other set of regression coefficients.
  //! @return the relative difference between `x` and `y`.
  template<typename T>
  double Difference(const RegressionCoefficients<T>& x, const RegressionCoefficients<T>& y) const {
    const double pred_norm = (pred_norm_ < 0) ? std::min(arma::norm(data_->cx(), "inf"), arma::norm(data_->cx(), 1)) :
                                                pred_norm_;
    return std::sqrt(data_->n_obs()) * std::abs(x.intercept - y.intercept) +
      pred_norm * arma::norm(x.beta - y.beta, 2);
  }

  //! Get the difference between two sets of regression coefficients.
  //!
  //! For the LS loss, the relative difference is an approximation to the 2-norm of the matrix-vector product
  //! ||X . (beta1 - beta2) + ( mu1 - mu2 )||_2 < |mu1 - mu2| sqrt(n) + ||X|| ||beta1 - beta2||_2
  //!
  //! @param x a set of regression coefficients.
  //! @param y the other set of regression coefficients.
  //! @return the relative difference between `x` and `y`.
  template<typename T>
  double Difference(const RegressionCoefficients<T>& x, const RegressionCoefficients<T>& y) {
    if (pred_norm_ < 0) {
      pred_norm_ = std::min(arma::norm(data_->cx(), "inf"), arma::norm(data_->cx(), 1));
    }
    return std::sqrt(data_->n_obs()) * std::abs(x.intercept - y.intercept) +
      pred_norm_ * arma::norm(x.beta - y.beta, 2);
  }

  //! Evaluate the gradient of the weighted LS loss at the given coefficient value.
  //!
  //! The gradient of the LS loss is given by
  //! intercept = -1/n * sum_{i = 1}^n w_i (y_i - \hat{\mu} - x_i' . \beta)
  //! beta = -1/n sum_{i = 1}^n w_i (y_i - \hat{\mu} - x_i' . \beta) x_i
  //!
  //! @param coefs the point where the gradient should be evaluated.
  //! @param gradient_intercept a pointer where the gradient of the intercept should be stored at.
  //! @param gradient_beta a pointer where the gradient of the slope should be stored at.
  template<typename T>
  GradientType<T> Gradient(const RegressionCoefficients<T>& where) const {
    const arma::vec neg_residuals = data_->cx() * where.beta - data_->cy();
    if (include_intercept_) {
      return GradientType<T>(arma::mean(neg_residuals) + where.intercept,
                             arma::mean(data_->cx().each_col() % (neg_residuals + where.intercept), 0));
    }
    return GradientType<T>(0, arma::mean(data_->cx().each_col() % neg_residuals, 0));
  }

  //! Access the data used by this weighted LS loss function.
  //!
  //! @return a constant reference to this loss' data.
  const PredictorResponseData& data() const noexcept {
    return *data_;
  }

  //! Check if the intercept term is included in this weighted LS loss function.
  //!
  //! @return `true` if the intercept term is included, `false` otherwise.
  bool IncludeIntercept() const noexcept {
    return include_intercept_;
  }

 private:
  bool include_intercept_;
  std::shared_ptr<const PredictorResponseData> data_;
  double pred_norm_;
};
}  // namespace nsoptim

#endif  // NSOPTIM_OBJECTIVE_LS_REGRESSION_LOSS_HPP_
