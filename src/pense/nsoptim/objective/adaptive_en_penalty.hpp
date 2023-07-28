//
//  adaptive_en_penalty.hpp
//  nsoptim
//
//  Created by David Kepplinger on 2018-11-30.
//  Copyright © 2018 David Kepplinger. All rights reserved.
//

#ifndef NSOPTIM_OBJECTIVE_ADAPTIVE_EN_PENALTY_HPP_
#define NSOPTIM_OBJECTIVE_ADAPTIVE_EN_PENALTY_HPP_

#include <memory>

#include "../armadillo.hpp"
#include "convex.hpp"
#include "../container/regression_coefficients.hpp"
#include "penalty.hpp"


namespace nsoptim {
//! The adaptive elastic net penalty function with hyper-parameters *alpha*, *lambda*, and
//! non-negative penalty loadings.
class AdaptiveEnPenalty : public PenaltyFunction, public ConvexFunction<AdaptiveEnPenalty> {
 public:
  //! Declare this penalty function as an EN penalty.
  struct is_en_penalty_tag {};

  //! Type of the gradient if evaluated with coefficients of ``T``.
  template<typename T>
  using GradientType = typename T::SlopeCoefficient;

  //! Initialize an adaptive EN penalty.
  //!
  //! @param loadings A shared pointer to a constant vector of penalty loadings.
  //! @param alpha Value of the *alpha* hyper-parameter.
  //! @param lambda Value of the *lambda* hyper-parameter. Default is 0.
  AdaptiveEnPenalty(std::shared_ptr<const arma::vec> loadings, const double alpha, const double lambda = 0) noexcept
      : loadings_(loadings), alpha_(alpha), lambda_(lambda) {}

  //! Default copy constructor.
  //! The copied penalty will **share the pointer to the penalty loadings**!
  AdaptiveEnPenalty(const AdaptiveEnPenalty& other) = default;
  //! Default copy assignment operator.
  //! The copied penalty will **share the pointer to the penalty loadings**!
  AdaptiveEnPenalty& operator=(const AdaptiveEnPenalty& other) = default;
  //! Default move operator.
  //! The copied penalty will **share the pointer to the penalty loadings**!
  AdaptiveEnPenalty(AdaptiveEnPenalty&& other) = default;
  //! Default move assignment.
  //! The copied penalty will **share the pointer to the penalty loadings**!
  AdaptiveEnPenalty& operator=(AdaptiveEnPenalty&& other) = default;

  ~AdaptiveEnPenalty() = default;

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

  //! Get a constant reference to the vector of penalty loadings.
  const arma::vec& loadings() const noexcept {
    return *loadings_;
  }

  //! Evaluate the adaptive EN penalty at `where`.
  //!
  //! @param where Point where to evaluate the penalty function.
  //! @return penalty evaluated at *where*.
  template<typename T>
  double operator()(const RegressionCoefficients<T>& where) const {
    return Evaluate(where);
  }

  //! Evaluate the adaptive EN penalty at `where`.
  //!
  //! @param where Point where to evaluate the penalty function.
  //! @return penalty evaluated at `where`.
  template<typename VectorType>
  double Evaluate(const RegressionCoefficients<VectorType>& where) const {
    if (loadings_->n_elem > 0) {
      return lambda_ * (alpha_ * arma::accu(*loadings_ % arma::abs(where.beta)) +
        0.5 * (1 - alpha_) * arma::dot(*loadings_ % where.beta, where.beta));
    }
    return lambda_ * (alpha_ * arma::norm(where.beta, 1) + 0.5 * (1 - alpha_) * arma::dot(where.beta, where.beta));
  }

  //! Evaluate the subgradient of the adaptive EN penalty at the given coefficient value.
  //!
  //! Elements of the slope that are 0 will be set to 0 in the gradient.
  //!
  //! @param where Coefficients where the subgradient should be evaluated.
  template<typename T>
  T Gradient(const RegressionCoefficients<T>& where) const {
    // The gradient is computed only for the non-zero coefficients. The other elements are set to 0.
    if (loadings_->n_elem > 0) {
      return lambda_ * (alpha_ * (*loadings_ % arma::sign(where.beta)) + (1 - alpha_) * (*loadings_) % where.beta);
    }
    return lambda_ * (alpha_ * arma::sign(where.beta) + (1 - alpha_) * where.beta);
  }

 private:
  std::shared_ptr<const arma::vec> loadings_;
  double alpha_;
  double lambda_;
};

//! The adaptive lasso penalty with hyper-parmaeter *lambda* and non-negative penalty loadings.
class AdaptiveLassoPenalty : public PenaltyFunction, public ConvexFunction<AdaptiveLassoPenalty> {
 public:
  //! Declare this penalty function as an EN penalty.
  struct is_en_penalty_tag {};

  //! Type of the gradient if evaluated with coefficients of type ``T``.
  template<typename T>
  using GradientType = typename T::SlopeCoefficient;

  //! Initialize an adaptive LASSO penalty.
  //!
  //! @param loadings A shared pointer to a constant vector of penalty loadings.
  //! @param lambda Value of the *lambda* hyper-parameter. Default is 0.
  explicit AdaptiveLassoPenalty(std::shared_ptr<const arma::vec> loadings, const double lambda = 0) noexcept
      : loadings_(loadings), lambda_(lambda) {}

  //! Default copy constructor.
  //! The copied penalty will **share the pointer to the penalty loadings**!
  AdaptiveLassoPenalty(const AdaptiveLassoPenalty& other) = default;
  //! Default copy assignment operator.
  //! The copied penalty will **share the pointer to the penalty loadings**!
  AdaptiveLassoPenalty& operator=(const AdaptiveLassoPenalty& other) = default;
  //! Default move operator.
  //! The copied penalty will **share the pointer to the penalty loadings**!
  AdaptiveLassoPenalty(AdaptiveLassoPenalty&& other) = default;
  //! Default move assignment.
  //! The copied penalty will **share the pointer to the penalty loadings**!
  AdaptiveLassoPenalty& operator=(AdaptiveLassoPenalty&& other) = default;

  ~AdaptiveLassoPenalty() = default;

  //! Set the *alpha* hyper-parameter. For adaptive LASSO penalties, this has no effect.
  void alpha(const double) const noexcept {}

  //! Get the value of the *alpha* hyper-parameter, which is always 1.
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

  //! Get a constant reference to the vector of penalty loadings.
  const arma::vec& loadings() const noexcept {
    return *loadings_;
  }

  //! Cast the adaptive LASSO penalty to an adaptive EN penalty, with *alpha* set to 1.
  operator AdaptiveEnPenalty() const {
    return AdaptiveEnPenalty(loadings_, 0., lambda_);
  }

  //! Evaluate the adaptive LASSO penalty at `where`.
  //!
  //! @param where Point where to evaluate the penalty function.
  //! @return penalty evaluated at *where*.
  template<typename T>
  double operator()(const RegressionCoefficients<T>& where) const {
    return Evaluate(where);
  }

  //! Evaluate the adaptive LASSO penalty at `where`.
  //!
  //! @param where Point where to evaluate the penalty function.
  //! @return penalty evaluated at *where*.
  template<typename T>
  double Evaluate(const RegressionCoefficients<T>& where) const {
    if (loadings_->n_elem > 0) {
      return lambda_ * arma::accu(*loadings_ % arma::abs(where.beta));
    }
    return lambda_ * arma::norm(where.beta, 1);
  }

  //! Evaluate the subgradient of the adaptive LASSO penalty at the given coefficient value.
  //!
  //! Elements of the slope that are 0 will be set to 0 in the gradient.
  //!
  //! @param where Coefficients where the subgradient should be evaluated.
  template<typename T>
  T Gradient(const RegressionCoefficients<T>& where) const {
    // The gradient is computed only for the non-zero coefficients. The other elements are set to 0.
    if (loadings_->n_elem > 0) {
      return lambda_ * (*loadings_ % arma::sign(where.beta));
    }
    return lambda_ * arma::sign(where.beta);
  }

 private:
  std::shared_ptr<const arma::vec> loadings_;
  double lambda_;
};
}  // namespace nsoptim

#endif  // NSOPTIM_OBJECTIVE_ADAPTIVE_EN_PENALTY_HPP_
