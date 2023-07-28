//
//  m_loss.hpp
//  pense
//
//  Created by David Kepplinger on 2020-06-08.
//  Copyright © 2020 David Kepplinger. All rights reserved.
//

#ifndef M_LOSS_HPP_
#define M_LOSS_HPP_

#include <memory>
#include <algorithm>
#include "nsoptim.hpp"

#include "constants.hpp"
#include "rho.hpp"

namespace pense {
//! A regression loss function implementing the M-loss defined as
//!   (1/n) sum_{i=1}^n \rho((y_i - \mu - x_i' . \beta) / scale)
template<class RhoFunction>
class MLoss : public nsoptim::LossFunction<nsoptim::PredictorResponseData> {
  //! Alias for the shared pointer to the regression data.
  using ConstDataPtr = std::shared_ptr<const nsoptim::PredictorResponseData>;

 public:
  using ConvexSurrogateType = nsoptim::WeightedLsRegressionLoss;
  using ResidualType = arma::vec;

  MLoss(ConstDataPtr data, const RhoFunction& rho, const double scale, const bool include_intercept = true) noexcept
    : include_intercept_(include_intercept), data_(data), rho_(rho), scale_(scale),
      pred_norm_(std::min(arma::norm(data->cx(), "inf"), arma::norm(data->cx(), 1))) {}

  MLoss(const MLoss& other) = default;
  MLoss& operator=(const MLoss& other) = delete;
  MLoss(MLoss&& other) = default;
  MLoss& operator=(MLoss&& other) = default;
  ~MLoss() = default;

  //! Get a copy of this loss with the data replaced by the given data.
  //!
  //! @param data a new data set.
  MLoss ReplaceData(ConstDataPtr data) const {
    return MLoss(*this, data);
  }

  //! Get the zero coefficients for this loss type.
  //!
  //! @return zero coefficients.
  template<typename T, typename = typename std::enable_if<
    std::is_same<T, nsoptim::RegressionCoefficients<typename T::SlopeCoefficient>>::value, void>::type>
  T ZeroCoefficients() const noexcept {
    return nsoptim::RegressionCoefficients<typename T::SlopeCoefficient>(data_->n_pred());
  }

  //! Evaluate the M loss function.
  //!
  //! @param residuals residuals of the point where to evaluate the loss function.
  //! @return loss evaluated at `residuals`.
  double operator()(const ResidualType& residuals) const {
    return Evaluate(residuals);
  }

  //! Evaluate the M loss function.
  //!
  //! @param residuals residuals of the point where to evaluate the loss function.
  //! @return loss evaluated at `where`.
  double Evaluate(const ResidualType& residuals) const {
    return arma::mean(rho_(residuals, scale_));
  }

  //! Evaluate the M loss function.
  //!
  //! @param where point where to evaluate the loss function.
  //! @return loss evaluated at `where`.
  template<typename T>
  double operator()(const nsoptim::RegressionCoefficients<T>& where) const {
    return Evaluate(where);
  }

  //! Evaluate the M loss function.
  //!
  //! @param where point where to evaluate the loss function.
  //! @return loss evaluated at `where`.
  template<typename T>
  double Evaluate(const nsoptim::RegressionCoefficients<T>& where) const {
    return arma::mean(rho_(Residuals(where), scale_));
  }

  //! Get the residuals for the LS loss function.
  //!
  //! @param where Point where to compute the residuals.
  //! @return residuals at `where`.
  template<typename VectorType>
  ResidualType Residuals(const nsoptim::RegressionCoefficients<VectorType>& where) const {
    if (include_intercept_) {
      return data_->cy() - data_->cx() * where.beta - where.intercept;
    }
    return data_->cy() - data_->cx() * where.beta;
  }

  //! Get the weights for the surrogate LS-loss at the given residuals.
  //!
  //! @param residuals residuals where the surrogate weights are computed.
  //! @return a vector of weights, the same length as `residuals`.
  arma::vec SurrogateWeights(const ResidualType& residuals) const {
    return rho_.Weight(residuals, scale_) / (scale_ * scale_);
  }

  //! Get the weights for the surrogate LS-loss at the given location.
  //!
  //! @param where location where the surrogate weights are computed.
  //! @return a vector of weights.
  template<typename T>
  arma::vec SurrogateWeights(const nsoptim::RegressionCoefficients<T>& where) const {
    return SurrogateWeights(Residuals(where));
  }

  //! Get the convex surrogate for the M loss function
  //!
  //! @param where where the convex surrogate should be constructed.
  //! @return the weighted-ls loss surrogate loss function
  template<typename T>
  ConvexSurrogateType GetConvexSurrogate(const nsoptim::RegressionCoefficients<T>& where) const {
    return ConvexSurrogateType(data_, SurrogateWeights(where), include_intercept_);
  }

  //! Get the convex surrogate for the M loss function
  //!
  //! @param residuals residuals of the point where the convex surrogate should be constructed.
  //! @return the weighted-ls loss surrogate loss function
  ConvexSurrogateType GetConvexSurrogate(const ResidualType& residuals) const {
    return ConvexSurrogateType(data_, SurrogateWeights(residuals), include_intercept_);
  }

  //! Clone the M loss function. The returned object does not share anything with this loss function.
  //!
  //! @return a deep copy of this M loss function.
  MLoss Clone() const {
    return MLoss(*this);
  }

  //! Access the data used by this M loss function.
  //!
  //! @return a constant reference to this loss' data.
  const nsoptim::PredictorResponseData& data() const noexcept {
    return *data_;
  }

  //! Access the data used by this M loss function.
  //!
  //! @return a shared pointer to this loss' data.
  std::shared_ptr<const nsoptim::PredictorResponseData> SharedData() const noexcept {
    return data_;
  }

  //! Check if the intercept term is included in this M loss function.
  //!
  //! @return `true` if the intercept term is included, `false` otherwise.
  bool IncludeIntercept() const noexcept {
    return include_intercept_;
  }

 private:
  MLoss(const MLoss& other, ConstDataPtr data) :
      include_intercept_(other.include_intercept_), data_(data), rho_(other.rho_), scale_(other.scale_),
      pred_norm_(std::min(arma::norm(data->cx(), "inf"), arma::norm(data->cx(), 1))) {}

  bool include_intercept_;
  ConstDataPtr data_;
  RhoFunction rho_;
  double scale_;
  double pred_norm_;
};
}  // namespace pense

#endif  // M_LOSS_HPP_
