//
//  en_penalty.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2018-11-30.
//  Copyright © 2018 David Kepplinger. All rights reserved.
//

#ifndef NSOPTIM_OBJECTIVE_EN_PENALTY_HPP_
#define NSOPTIM_OBJECTIVE_EN_PENALTY_HPP_

#include <memory>

#include "../armadillo.hpp"
#include "../container/regression_coefficients.hpp"
#include "penalty.hpp"
#include "adaptive_en_penalty.hpp"

namespace nsoptim {
//! The EN penalty function with hyper-parameters *alpha* and *lambda*.
class EnPenalty : public PenaltyFunction, public ConvexFunction<EnPenalty> {
 public:
  //! Declare this penalty function as an EN penalty.
  struct is_en_penalty_tag {};

  //! Type of the gradient if evaluated with coefficients of type ``T``.
  template<typename T>
  using GradientType = typename T::SlopeCoefficient;

  //! Initialize an EN penalty, setting *alpha* to 1 and *lambda* to 0.
  EnPenalty() noexcept : alpha_(1), lambda_(0) {}

  //! Initialize an EN penalty.
  //!
  //! @param alpha Value for the *alpha* hyper-parameter.
  //! @param lambda Value for the *lambda* hyper-parameter.
  EnPenalty(const double alpha, const double lambda) noexcept : alpha_(alpha), lambda_(lambda) {}

  //! Default copy constructor.
  EnPenalty(const EnPenalty& other) = default;
  //! Default copy assignment operator.
  EnPenalty& operator=(const EnPenalty& other) = default;
  //! Default move operator.
  EnPenalty(EnPenalty&& other) = default;
  //! Default move assignment.
  EnPenalty& operator=(EnPenalty&& other) = default;

  ~EnPenalty() = default;

  //! Cast the EN penalty to an adaptive EN penalty (with penalty loadings all equal to 1).
  operator AdaptiveEnPenalty() const {
    return AdaptiveEnPenalty(std::make_shared<const arma::vec>(), alpha_, lambda_);
  }

  //! Set the *alpha* hyper-parameter.
  //!
  //! @param alpha New *alpha* value.
  void alpha(const double alpha) noexcept {
    alpha_ = alpha;
  }

  //! Get the value of the *alpha* hyper-parameter.
  double alpha() const noexcept {
    return alpha_;
  }

  //! Set the *lambda* hyper-parameter.
  //!
  //! @param lambda New *lambda* value.
  void lambda(const double lambda) noexcept {
    lambda_ = lambda;
  }

  //! Get the value of the *lambda* hyper-parameter.
  double lambda() const noexcept {
    return lambda_;
  }

  //! Evaluate the elastic net penalty at *where*.
  //!
  //! @param where Point where to evaluate the penalty function.
  //! @return penalty evaluated at *where*.
  template<typename T>
  double operator()(const RegressionCoefficients<T>& where) const {
    return Evaluate(where);
  }

  //! Evaluate the elastic net penalty at *where*.
  //!
  //! @param where point where to evaluate the penalty function.
  //! @return penalty evaluated at *where*.
  template<typename T>
  double Evaluate(const RegressionCoefficients<T>& where) const {
    return lambda_ * (alpha_ * arma::norm(where.beta, 1) + 0.5 * (1. - alpha_) * arma::dot(where.beta, where.beta));
  }

  //! Evaluate the subgradient of the EN penalty at the given coefficient value.
  //!
  //! Elements of the slope that are 0 will be set to 0 in the gradient.
  //!
  //! @param where Coefficients where the subgradient should be evaluated.
  template<typename T>
  T Gradient(const RegressionCoefficients<T>& where) const {
    // Compute the gradient of the non-zero elements. The 0-elements are also 0 in the gradient!
    return lambda_ * (alpha_ * arma::sign(where.beta) + (1 - alpha_) * where.beta);
  }

 private:
  double alpha_;
  double lambda_;
};

//! The LASSO penalty function with hyper-parameter *lambda*.
class LassoPenalty : public PenaltyFunction, public ConvexFunction<LassoPenalty> {
 public:
  //! Declare this penalty function as an EN penalty.
  struct is_en_penalty_tag {};

  //! Type of the gradient if evaluated with coefficients of type ``T``.
  template<typename T>
  using GradientType = typename T::SlopeCoefficient;

  //! Initialize a LASSO penalty with *lambda* set to 0.
  LassoPenalty() noexcept : lambda_(0) {}

  //! Initialize a LASSO penalty.
  //!
  //! @param lambda Value for the *lambda* hyper-parameter.
  explicit LassoPenalty(const double lambda) noexcept : lambda_(lambda) {}

  //! Default copy constructor.
  LassoPenalty(const LassoPenalty& other) = default;

  //! Default copy assignment.
  LassoPenalty& operator=(const LassoPenalty& other) = default;

  //! Default move constructor.
  LassoPenalty(LassoPenalty&& other) = default;

  //! Default move assignment.
  LassoPenalty& operator=(LassoPenalty&& other) = default;

  ~LassoPenalty() = default;

  //! Get the value of the *alpha* hyper-parameter. This is always 1!
  double alpha() const noexcept {
    return 1.;
  }

  //! Set the *lambda* hyper-parameter.
  //!
  //! @param lambda New *lambda* value.
  void lambda(const double lambda) noexcept {
    lambda_ = lambda;
  }

  //! Get the value of the *lambda* hyper-parameter.
  double lambda() const noexcept {
    return lambda_;
  }

  //! Cast the LASSO penalty to an adaptive EN penalty, with penalty loadings all equal to 1 and *alpha* set to 1.
  operator AdaptiveEnPenalty() const {
    return AdaptiveEnPenalty(std::make_shared<const arma::vec>(), 1., lambda_);
  }

  //! Cast the LASSO penalty to an adaptive LASSO penalty, with penalty loadings all equal to 1.
  operator AdaptiveLassoPenalty() const {
    return AdaptiveLassoPenalty(std::make_shared<const arma::vec>(), lambda_);
  }

  //! Cast the LASSO penalty to an EN penalty, with *alpha* set to 1.
  operator EnPenalty() const {
    return EnPenalty(1, lambda_);
  }

  //! Evaluate the LASSO penalty at *where*.
  //!
  //! @param where Point where to evaluate the penalty function.
  //! @return Penalty evaluated at *where*.
  template<typename T>
  double operator()(const RegressionCoefficients<T>& where) const {
    return Evaluate(where);
  }

  //! Evaluate the LASSO penalty at *where*.
  //!
  //! @param where Point where to evaluate the penalty function.
  //! @return Penalty evaluated at *where*.
  template<typename T>
  double Evaluate(const RegressionCoefficients<T>& where) const {
    return lambda_ * arma::norm(where.beta, 1);
  }

  //! Evaluate the subgradient of the LASSO penalty at the given coefficient value.
  //!
  //! Elements of the slope that are 0 will be set to 0 in the gradient.
  //!
  //! @param where Coefficients where the subgradient should be evaluated.
  template<typename T>
  T Gradient(const RegressionCoefficients<T>& where) const {
    // Compute the gradient of the non-zero elements. The 0-elements are also 0 in the gradient!
    return lambda_ * arma::sign(where.beta);
  }

 private:
  double lambda_;
};

//! The LASSO penalty function with hyper-parameter *lambda*.
class RidgePenalty : public PenaltyFunction, public ConvexFunction<RidgePenalty> {
 public:
  //! Declare this penalty function as an EN penalty.
  struct is_en_penalty_tag {};

  //! Type of the gradient if evaluated with coefficients of type ``T``.
  template<typename T>
  using GradientType = typename T::SlopeCoefficient;

  //! Initialize a Ridge penalty with *lambda* set to 0.
  RidgePenalty() noexcept : lambda_(0) {}

  //! Initialize a Ridge penalty.
  //!
  //! @param lambda Value for the *lambda* hyper-parameter.
  explicit RidgePenalty(const double lambda) noexcept : lambda_(lambda) {}

  //! Default copy constructor.
  RidgePenalty(const RidgePenalty& other) = default;
  //! Default copy assignment.
  RidgePenalty& operator=(const RidgePenalty& other) = default;
  //! Default move constructor.
  RidgePenalty(RidgePenalty&& other) = default;
  //! Default move assignment.
  RidgePenalty& operator=(RidgePenalty&& other) = default;

  ~RidgePenalty() = default;

  //! Get the value of the *alpha* hyper-parameter. This is always 0s!
  double alpha() const noexcept {
    return 0.;
  }

  //! Set the *lambda* hyper-parameter.
  //!
  //! @param lambda New *lambda* value.
  void lambda(const double lambda) noexcept {
    lambda_ = lambda;
  }

  //! Get the value of the *lambda* hyper-parameter.
  double lambda() const noexcept {
    return lambda_;
  }

  //! Cast the Ridge penalty to an adaptive EN penalty, with penalty loadings all equal to 1 and *alpha* set to 0.
  operator AdaptiveEnPenalty() const {
    return AdaptiveEnPenalty(std::make_shared<const arma::vec>(), 0, lambda_);
  }

  //! Cast the Ridge penalty to an EN penalty, with *alpha* set to 0.
  operator EnPenalty() const {
    return EnPenalty(0., lambda_);
  }

  //! Evaluate the Ridge penalty at *where*.
  //!
  //! @param where Point where to evaluate the penalty function.
  //! @return Penalty evaluated at *where*.
  template<typename T>
  double operator()(const RegressionCoefficients<T>& where) const {
    return Evaluate(where);
  }

  //! Evaluate the Ridge penalty at *where*.
  //!
  //! @param where Point where to evaluate the penalty function.
  //! @return Penalty evaluated at *where*.
  template<typename T>
  double Evaluate(const RegressionCoefficients<T>& where) const {
    return 0.5 * lambda_ * arma::dot(where.beta, where.beta);
  }

  //! Evaluate the gradient of the Ridge penalty at the given coefficient value.
  //!
  //! @param where Coefficients where the gradient should be evaluated.
  template<typename T>
  T Gradient(const RegressionCoefficients<T>& where) const {
    // Compute the gradient of the non-zero elements. The 0-elements are also 0 in the gradient!
    return lambda_ * where.beta;
  }

 private:
  double lambda_;
};
}  // namespace nsoptim

#endif  // NSOPTIM_OBJECTIVE_EN_PENALTY_HPP_
