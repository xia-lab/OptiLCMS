//
//  rho.hpp
//  pense
//
//  Created by David Kepplinger on 2019-01-30.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifndef RHO_HPP_
#define RHO_HPP_

#include "nsoptim.hpp"

namespace pense {
//! Implementation of Huber's unbounded rho function defined by
//! rho(x) = 0.5 x^2 * [x < cc] + cc * (abs(x) - cc / 2) [x >= cc]
class RhoHuber {
 public:
  //! Create the Huber rho-function with cutoff `cc`
  //!
  //! @param cc cutoff value (outside the rho-function is linear)
  explicit RhoHuber(const double cc) noexcept : cc_(cc) {}

  RhoHuber(const RhoHuber&) = default;
  RhoHuber& operator=(const RhoHuber&) = default;
  RhoHuber(RhoHuber&&) = default;
  RhoHuber& operator=(RhoHuber&&) = default;

  ~RhoHuber() {}

  //! Get the value of the threshold paramter.
  //!
  //! @return threshold paramter value.
  double cc() const noexcept {
    return cc_;
  }

  //! Get the value of the rho function evaluated at x/scale.
  //! This is equivalent to `robustbase::Mpsi(x / scale, cc, 'huber', deriv = -1)`.
  double    operator()(double x, const double scale) const noexcept;
  arma::vec operator()(const arma::vec& x, const double scale) const noexcept;
  void      operator()(const arma::vec& x, const double scale, arma::vec* out) const noexcept;
  double    Sum(const arma::vec& x, const double scale) const noexcept;

  //! Get the derivative of the rho function evaluated at x/scale.
  //! Note that Derivative(x/scale, 1) == Derivative(x, scale) * scale
  //! This is equivalent to `robustbase::Mpsi(x / scale, cc, 'huber', deriv = 0)`.
  double    Derivative(double x, const double scale) const noexcept;
  arma::vec Derivative(const arma::vec& x, const double scale) const noexcept;
  void      Derivative(const arma::vec& x, const double scale, arma::vec* out) const noexcept;

  //! Get the second derivative of the rho function evaluated at x/scale.
  //! Note that SecondDerivative(x/scale, 1) == SecondDerivative(x, scale)
  //! This is equivalent to `robustbase::Mpsi(x / scale, cc, 'huber', deriv = 1)`.
  double    SecondDerivative(double x, const double scale) const noexcept;
  arma::vec SecondDerivative(const arma::vec& x, const double scale) const noexcept;
  void      SecondDerivative(const arma::vec& x, const double scale, arma::vec* out) const noexcept;

  //! Get the derivative of the rho function evaluated at (x/scale), divided by (x/scale), i.e.,
  //! Weight(x, scale) = Derivative(x/scale) / (x/scale).
  //! Note that Weight(x/scale, 1) == Weight(x, scale)
  //! This is equivalent to `robustbase::Mwgt(x / scale, cc, 'huber')`.
  double    Weight(double x, const double scale) const noexcept;
  arma::vec Weight(const arma::vec& x, const double scale) const noexcept;
  void      Weight(const arma::vec& x, const double scale, arma::vec* out) const noexcept;

 private:
  double cc_;
};

//! Implementation of Tukey's bisquare rho function defined by
//! rho(x) = min(1, 1 - (1 - x^2/cc^2)^3)
//! All of the virtual functions are declared final and can not be overwritten again!
class RhoBisquare {
 public:
  //! Create the bisquare-rho function with threshold `cc`
  //!
  //! @param cc threshold paramter value.
  explicit RhoBisquare(const double cc) noexcept : cc_(cc) {}

  RhoBisquare(const RhoBisquare&) = default;
  RhoBisquare& operator=(const RhoBisquare&) = default;
  RhoBisquare(RhoBisquare&&) = default;
  RhoBisquare& operator=(RhoBisquare&&) = default;

  ~RhoBisquare() {}

  //! Get the value of the threshold paramter.
  //!
  //! @return threshold paramter value.
  double cc() const noexcept {
    return cc_;
  }

  //! Get the value of the unstandardized rho function evaluated at x/scale.
  //! This is equivalent to `robustbase::Mpsi(x / scale, cc, 'bisquare', deriv = -1)`.
  double    operator()(double x, const double scale) const noexcept;
  arma::vec operator()(const arma::vec& x, const double scale) const noexcept;
  void      operator()(const arma::vec& x, const double scale,
                       arma::vec* out) const noexcept;
  double    Sum(const arma::vec& x, const double scale) const noexcept;

  //! Get the value of the *standardized* rho function evaluated at x/scale.
  //! Standardized means that rho(inf) = 1.
  //! This is equivalent to `robustbase::Mchi(x / scale, cc, 'bisquare', deriv = 0)`.
  double    EvaluateStd(double x, const double scale) const noexcept;
  arma::vec EvaluateStd(const arma::vec& x, const double scale) const noexcept;
  void      EvaluateStd(const arma::vec& x, const double scale,
                       arma::vec* out) const noexcept;
  double    SumStd(const arma::vec& x, const double scale) const noexcept;

  //! Get the derivative of the *unstandardized* rho function evaluated at x/scale.
  //! Note that Derivative(x/scale, 1) == Derivative(x, scale) * scale
  //! This is equivalent to `robustbase::Mpsi(x / scale, cc, 'bisquare', deriv = 0)`
  double    Derivative(double x, const double scale) const noexcept;
  arma::vec Derivative(const arma::vec& x, const double scale) const noexcept;
  void      Derivative(const arma::vec& x, const double scale, arma::vec* out) const noexcept;

  //! Get the derivative of the *standardized* rho function evaluated at x/scale.
  //! Note that DerivativeStd(x/scale, 1) == DerivativeStd(x, scale) / scale.
  //! This is equivalent to `robustbase::Mchi(x / scale, cc, 'bisquare', deriv = 1)`.
  double    DerivativeStd(double x, const double scale) const noexcept;
  arma::vec DerivativeStd(const arma::vec& x, const double scale) const noexcept;
  void      DerivativeStd(const arma::vec& x, const double scale, arma::vec* out) const noexcept;

  //! Get the second derivative of the *unstandardized* rho function evaluated at x/scale.
  //! Note that SecondDerivative(x/scale, 1) == SecondDerivative(x, scale)
  //! This is equivalent to `robustbase::Mpsi(x / scale, cc, 'bisquare', deriv = 1)`.
  double    SecondDerivative(double x, const double scale) const noexcept;
  arma::vec SecondDerivative(const arma::vec& x, const double scale) const noexcept;
  void      SecondDerivative(const arma::vec& x, const double scale, arma::vec* out) const noexcept;

  //! Get the second derivative of the *standardized* rho function evaluated at x/scale.
  //! Note that SecondDerivativeStd(x/scale, 1) == SecondDerivativeStd(x, scale) * scale^2
  //! This is equivalent to `robustbase::Mchi(x / scale, cc, 'bisquare', deriv = 2)`.
  double    SecondDerivativeStd(double x, const double scale) const noexcept;
  arma::vec SecondDerivativeStd(const arma::vec& x, const double scale) const noexcept;
  void      SecondDerivativeStd(const arma::vec& x, const double scale, arma::vec* out) const noexcept;

  //! Get the derivative of the *unstandardized* rho function evaluated at (x/scale), divided by (x/scale), i.e.,
  //! Weight(x, scale) = Derivative(x/scale) / (x/scale).
  //! Note that Weight(x/scale, 1) == Weight(x, scale)
  //! This is equivalent to `robustbase::Mwgt(x / scale, cc, 'bisquare')`.
  double    Weight(double x, const double scale) const noexcept;
  arma::vec Weight(const arma::vec& x, const double scale) const noexcept;
  void      Weight(const arma::vec& x, const double scale, arma::vec* out) const noexcept;

  //! Get the derivative of the *standardized* rho function evaluated at (x/scale), divided by (x/scale), i.e.,
  //! WeightStd(x, scale) = DerivativeStd(x/scale) / (x/scale).
  //! Note that WeightStd(x/scale, 1) == WeightStd(x, scale)
  //! This is equivalent to
  //!   `robustbase::Mwgt(x / scale, cc, 'bisquare') / robustbase::MrhoInf(cc, 'bisquare')`.
  double    WeightStd(double x, const double scale) const noexcept;
  arma::vec WeightStd(const arma::vec& x, const double scale) const noexcept;
  void      WeightStd(const arma::vec& x, const double scale, arma::vec* out) const noexcept;

  //! Get the upper bound of the rho function, i.e., the limiting value of `operator()(x)` for x to infinity.
  double UpperBound() const noexcept;

 private:
  double cc_;
};
}  // namespace pense
#endif  // RHO_HPP_
