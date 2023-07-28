//
//  s_loss.hpp
//  pense
//
//  Created by David Kepplinger on 2019-01-30.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifndef S_LOSS_HPP_
#define S_LOSS_HPP_

#include <memory>
#include <algorithm>
#include "nsoptim.hpp"

#include "constants.hpp"
#include "robust_scale_location.hpp"

namespace pense {
//! A regression loss function implementing the robust S-loss defined as
//!   1/2 $sigma_M(\mu, \beta)^2$
//! where $sigma_M$ is the robust M-scale of the residuals, defined implicitly by
//!   $(1/n) sum_{i=1}^n \rho((y_i - \mu - x_i' . \beta) / sigma_M(\mu, \beta)) = \delta
class SLoss : public nsoptim::LossFunction<nsoptim::PredictorResponseData> {
  //! Alias for the shared pointer to the regression data.
  using ConstDataPtr = std::shared_ptr<const nsoptim::PredictorResponseData>;

 public:
  using ConvexSurrogateType = nsoptim::WeightedLsRegressionLoss;
  using ResidualType = arma::vec;

  struct ExtendedEvaluation {
    double loss;
    double scale;
  };

  SLoss(ConstDataPtr data, const Mscale<RhoBisquare>& mscale, const bool include_intercept = true) noexcept
    : include_intercept_(include_intercept), data_(data), mscale_(mscale),
      pred_norm_(std::min(arma::norm(data->cx(), "inf"), arma::norm(data->cx(), 1))) {}

  SLoss(const SLoss& other) = default;
  SLoss& operator=(const SLoss& other) = delete;
  SLoss(SLoss&& other) = default;
  SLoss& operator=(SLoss&& other) = default;
  ~SLoss() = default;

  //! Get a copy of this loss with the data replaced by the given data.
  //!
  //! @param data a new data set.
  SLoss ReplaceData(ConstDataPtr data) const {
    return SLoss(*this, data);
  }

  //! Get the zero coefficients for this loss type.
  //!
  //! @return zero coefficients.
  template<typename T, typename = typename std::enable_if<
    std::is_same<T, nsoptim::RegressionCoefficients<typename T::SlopeCoefficient>>::value, void>::type>
  T ZeroCoefficients() const noexcept {
    return nsoptim::RegressionCoefficients<typename T::SlopeCoefficient>(data_->n_pred());
  }

  //! Evaluate the S loss function.
  //!
  //! @param where point where to evaluate the loss function.
  //! @return loss evaluated at `where`.
  template<typename T>
  double operator()(const nsoptim::RegressionCoefficients<T>& where) {
    return Evaluate(where);
  }

  //! Evaluate the S loss function.
  //!
  //! @param where point where to evaluate the loss function.
  //! @return loss evaluated at `where`.
  template<typename T>
  double operator()(const nsoptim::RegressionCoefficients<T>& where) const {
    return Evaluate(where);
  }

  //! Evaluate the S loss function.
  //!
  //! @param where point where to evaluate the loss function.
  //! @return loss evaluated at `where`.
  template<typename T>
  double Evaluate(const nsoptim::RegressionCoefficients<T>& where) {
    const double scale = mscale_(Residuals(where));
    return 0.5 * scale * scale;
  }

  //! Evaluate the S loss function.
  //!
  //! @param where point where to evaluate the loss function.
  //! @return loss evaluated at `where`.
  template<typename T>
  double Evaluate(const nsoptim::RegressionCoefficients<T>& where) const {
    const double scale = mscale_(Residuals(where));
    return 0.5 * scale * scale;
  }

  //! Evaluate the S loss function.
  //!
  //! @param residuals residuals at the point where to evaluate the loss function.
  //! @return loss evaluated at `residuals`.
  double operator()(const ResidualType& residuals) const {
    return Evaluate(residuals);
  }

  //! Evaluate the S loss function.
  //!
  //! @param residuals residuals at the point where to evaluate the loss function.
  //! @return loss evaluated at `residuals`.
  double operator()(const ResidualType& residuals) {
    return Evaluate(residuals);
  }

  //! Evaluate the S loss function.
  //!
  //! @param residuals residuals at the point where to evaluate the loss function.
  //! @return loss evaluated at `residuals`.
  double Evaluate(const ResidualType& residuals) {
    const double scale = mscale_(residuals);
    return 0.5 * scale * scale;
  }

  //! Evaluate the S loss function.
  //!
  //! @param residuals residuals at the point where to evaluate the loss function.
  //! @return loss evaluated at `residuals`.
  double Evaluate(const ResidualType& residuals) const {
    const double scale = mscale_(residuals);
    return 0.5 * scale * scale;
  }

  //! Evaluate the S loss function at the given *residuals*.
  //!
  //! @param residuals evaluate the S-loss function at the given residuals.
  //! @return the loss and the scale evaluated at the given residuals.
  ExtendedEvaluation EvaluateResiduals(const ResidualType& residuals) const {
    ExtendedEvaluation result;
    result.scale = mscale_(residuals);
    result.loss = 0.5 * result.scale * result.scale;
    return result;
  }

  //! Evaluate the S loss function at the given *residuals*.
  //!
  //! @param residuals evaluate the S-loss function at the given residuals.
  //! @return the loss and the scale evaluated at the given residuals.
  ExtendedEvaluation EvaluateResiduals(const ResidualType& residuals) {
    ExtendedEvaluation result;
    result.scale = mscale_(residuals);
    result.loss = 0.5 * result.scale * result.scale;
    return result;
  }

  template<typename Coefficients>
  arma::vec Residuals(const Coefficients& where) const {
    return data_->cy() - data_->cx() * where.beta - where.intercept;
  }

  //! Get the weights for the surrogate LS-loss at the given residuals.
  //!
  //! @param residuals residuals where the surrogate weights are computed.
  //! @return a vector of weights, the same length as `residuals`.
  arma::vec SurrogateWeights(const ResidualType& residuals) {
    const double scale = mscale_(residuals);

    // Check if the scale is 0.
    if (scale < kNumericZero) {
      throw ZeroWeightsException();
    }

    arma::vec weights = mscale_.rho().Weight(residuals, scale);
    const double denominator = arma::dot(weights, arma::square(residuals));

    return weights * residuals.n_elem * scale * scale / denominator;
  }

  //! Get the weights for the surrogate LS-loss at the given location.
  //!
  //! @param where location where the surrogate weights are computed.
  //! @return a vector of weights.
  template<typename T>
  arma::vec SurrogateWeights(const nsoptim::RegressionCoefficients<T>& where) {
    return SurrogateWeights(Residuals(where));
  }

  //! Get the convex surrogate for the S loss function
  //!
  //! @param where where the convex surrogate should be constructed.
  //! @return the weighted-ls loss surrogate loss function
  template<typename T>
  ConvexSurrogateType GetConvexSurrogate(const nsoptim::RegressionCoefficients<T>& where) {
    return ConvexSurrogateType(data_, SurrogateWeights(where), include_intercept_);
  }

  //! Get the convex surrogate for the S loss function
  //!
  //! @param residuals residuals of the solution where the convex surrogate should be constructed.
  //! @return the weighted-ls loss surrogate loss function
  ConvexSurrogateType GetConvexSurrogate(const ResidualType& residuals) {
    return ConvexSurrogateType(data_, SurrogateWeights(residuals), include_intercept_);
  }

  //! Clone the S loss function. The returned object does not share anything with this loss function.
  //!
  //! @return a deep copy of this S loss function.
  SLoss Clone() const {
    return SLoss(*this);
  }

  //! Access the data used by this S loss function.
  //!
  //! @return a constant reference to this loss' data.
  const nsoptim::PredictorResponseData& data() const noexcept {
    return *data_;
  }

  //! Access the data used by this S loss function.
  //!
  //! @return a shared pointer to this loss' data.
  std::shared_ptr<const nsoptim::PredictorResponseData> SharedData() const noexcept {
    return data_;
  }

  //! Check if the intercept term is included in this S loss function.
  //!
  //! @return `true` if the intercept term is included, `false` otherwise.
  bool IncludeIntercept() const noexcept {
    return include_intercept_;
  }

  //! Get the M-scale function used by this S-loss.
  //!
  //! @return a constant reference to the M-scale function.
  const Mscale<RhoBisquare>& mscale() const noexcept {
    return mscale_;
  }

  //! Get the M-scale function used by this S-loss.
  //!
  //! @return a constant reference to the M-scale function.
  Mscale<RhoBisquare>& mscale() noexcept {
    return mscale_;
  }

 private:
  SLoss(const SLoss& other, ConstDataPtr data) :
      include_intercept_(other.include_intercept_), data_(data), mscale_(other.mscale_),
      pred_norm_(std::min(arma::norm(data->cx(), "inf"), arma::norm(data->cx(), 1))) {}

  bool include_intercept_;
  ConstDataPtr data_;
  Mscale<RhoBisquare> mscale_;
  double pred_norm_;
};
}  // namespace pense

#endif  // S_LOSS_HPP_
