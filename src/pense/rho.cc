//
//  rho.cc
//  pense
//
//  Created by David Kepplinger on 2019-01-30.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#include <cmath>
#include "nsoptim.hpp"

#include "rho.hpp"

using arma::vec;

namespace {
//! Actual implementations of the bisquare function.
//! Note that the derivatives are for the *unstandardized* bisquare function, while BisquareFunctionValueStd
//! gives the *standardized* value, i.e., BisquareFunctionValueStd(t, cc_scaled) == 1 for all |t| > cc_scaled!
double BisquareFunctionValueStd(double x, const double cc_scaled) noexcept;
double BisquareDerivativeValue(double x, const double cc_scaled) noexcept;
double BisquareSecondDerivativeValue(double x, const double cc_scaled) noexcept;
double BisquareWeightValue(double x, const double cc_scaled) noexcept;

//! Actual implementations of the Huber function.
double HuberFunctionValue(double x, const double scale, const double cc) noexcept;
double HuberDerivativeValue(double x, const double scale_sq, const double cc_inv_scaled) noexcept;
double HuberSecondDerivativeValue(double x, const double cc_scaled) noexcept;
double HuberWeightValue(double x, const double cc_scaled) noexcept;
}  // namespace

namespace pense {
// ================================== Huber's Rho Function ========================================================== //
double RhoHuber::operator()(double x, const double scale) const noexcept {
  return HuberFunctionValue(x, scale, cc_);
}

vec RhoHuber::operator()(const vec& x, const double scale) const noexcept {
  vec out;
  this->operator()(x, scale, &out);
  return out;
}

void RhoHuber::operator()(const vec& x, const double scale, vec* out) const noexcept {
  out->copy_size(x);
  auto read_it = x.cbegin();
  for (auto write_it = out->begin(), end = out->end(); write_it != end; ++write_it, ++read_it) {
    *write_it = HuberFunctionValue(*read_it, scale, cc_);
  }
}

double RhoHuber::Sum(const vec& x, const double scale) const noexcept {
  double tmp = 0.;
  for (auto read_it = x.cbegin(); read_it != x.cend(); ++read_it) {
    tmp += HuberFunctionValue(*read_it, scale, cc_);
  }
  return tmp;
}

double RhoHuber::Derivative(double x, const double scale) const noexcept {
  return HuberDerivativeValue(x, scale * scale, cc_ / scale);
}

vec RhoHuber::Derivative(const vec& x, const double scale) const noexcept {
  vec out;
  this->Derivative(x, scale, &out);
  return out;
}

void RhoHuber::Derivative(const vec& x, const double scale, vec* out) const
    noexcept {
  auto read_it = x.cbegin();
  const double scale_sq = scale * scale;
  const double cc_inv_scale = cc_ / scale;
  out->copy_size(x);
  for (auto write_it = out->begin(), end = out->end(); write_it != end; ++write_it, ++read_it) {
    *write_it = HuberDerivativeValue(*read_it, scale_sq, cc_inv_scale);
  }
}

double RhoHuber::Weight(double x, const double scale) const noexcept {
  return HuberWeightValue(x, scale * cc_);
}

vec RhoHuber::Weight(const vec& x, const double scale) const noexcept {
  vec out;
  this->Weight(x, scale, &out);
  return out;
}

void RhoHuber::Weight(const vec& x, const double scale, vec* out) const noexcept {
  const double cc_scaled = cc_ * scale;
  auto read_it = x.cbegin();
  out->copy_size(x);
  for (auto write_it = out->begin(), end = out->end(); write_it != end; ++write_it, ++read_it) {
    *write_it = HuberWeightValue(*read_it, cc_scaled);
  }
}

double RhoHuber::SecondDerivative(double x, const double scale) const noexcept {
  return HuberSecondDerivativeValue(x, cc_ * scale);
}

vec RhoHuber::SecondDerivative(const vec& x, const double scale) const noexcept {
  vec out;
  this->SecondDerivative(x, scale, &out);
  return out;
}

void RhoHuber::SecondDerivative(const vec& x, const double scale, vec* out) const noexcept {
  const double cc_scaled = cc_ * scale;
  auto read_it = x.cbegin();
  out->copy_size(x);
  for (auto write_it = out->begin(), end = out->end(); write_it != end; ++write_it, ++read_it) {
    *write_it = HuberSecondDerivativeValue(*read_it, cc_scaled);
  }
}

// ============================= Tukey's Bisquare Function ========================================================== //

double RhoBisquare::UpperBound() const noexcept {
  return cc_ * cc_ / 6.;
}

double RhoBisquare::operator()(double x, const double scale) const noexcept {
  return UpperBound() * BisquareFunctionValueStd(x, cc_ * scale);
}

vec RhoBisquare::operator()(const vec& x, const double scale) const noexcept {
  vec out;
  this->operator()(x, scale, &out);
  return out;
}

void RhoBisquare::operator()(const vec& x, const double scale, vec* out) const noexcept {
  out->copy_size(x);
  const double cc_scaled = cc_ * scale;
  const double rho_inf = UpperBound();
  auto read_it = x.cbegin();
  for (auto write_it = out->begin(), end = out->end(); write_it != end; ++write_it, ++read_it) {
    *write_it = rho_inf * BisquareFunctionValueStd(*read_it, cc_scaled);
  }
}

double RhoBisquare::Sum(const vec& x, const double scale) const noexcept {
  double tmp = 0.;
  const double cc_scaled = cc_ * scale;
  for (auto read_it = x.cbegin(); read_it != x.cend(); ++read_it) {
    tmp += BisquareFunctionValueStd(*read_it, cc_scaled);
  }
  return UpperBound() * tmp;
}

double RhoBisquare::EvaluateStd(double x, const double scale) const noexcept {
  return BisquareFunctionValueStd(x, cc_ * scale);
}

vec RhoBisquare::EvaluateStd(const vec& x, const double scale) const noexcept {
  vec out;
  this->operator()(x, scale, &out);
  return out;
}

void RhoBisquare::EvaluateStd(const vec& x, const double scale, vec* out) const noexcept {
  out->copy_size(x);
  const double cc_scaled = cc_ * scale;
  auto read_it = x.cbegin();
  for (auto write_it = out->begin(), end = out->end(); write_it != end; ++write_it, ++read_it) {
    *write_it = BisquareFunctionValueStd(*read_it, cc_scaled);
  }
}

double RhoBisquare::SumStd(const vec& x, const double scale) const noexcept {
  double tmp = 0.;
  const double cc_scaled = cc_ * scale;
  for (auto read_it = x.cbegin(); read_it != x.cend(); ++read_it) {
    tmp += BisquareFunctionValueStd(*read_it, cc_scaled);
  }
  return tmp;
}

double RhoBisquare::Derivative(double x, const double scale) const noexcept {
  return BisquareDerivativeValue(x, scale * cc_);
}

vec RhoBisquare::Derivative(const vec& x, const double scale) const noexcept {
  vec out;
  this->Derivative(x, scale, &out);
  return out;
}

void RhoBisquare::Derivative(const vec& x, const double scale, vec* out) const noexcept {
  const double cc_scaled = cc_ * scale;
  auto read_it = x.cbegin();
  out->copy_size(x);
  for (auto write_it = out->begin(), end = out->end(); write_it != end; ++write_it, ++read_it) {
    *write_it = BisquareDerivativeValue(*read_it, cc_scaled);
  }
}

double RhoBisquare::DerivativeStd(double x, const double scale) const noexcept {
  return BisquareDerivativeValue(x, scale * cc_) / UpperBound();
}

arma::vec RhoBisquare::DerivativeStd(const arma::vec& x, const double scale) const noexcept {
  vec out;
  this->DerivativeStd(x, scale, &out);
  return out;
}

void RhoBisquare::DerivativeStd(const arma::vec& x, const double scale, arma::vec* out) const noexcept {
  const double cc_scaled = cc_ * scale;
  const double rho_inf = UpperBound();
  auto read_it = x.cbegin();
  out->copy_size(x);
  for (auto write_it = out->begin(), end = out->end(); write_it != end; ++write_it, ++read_it) {
    *write_it = BisquareDerivativeValue(*read_it, cc_scaled) / rho_inf;
  }
}

double RhoBisquare::SecondDerivative(double x, const double scale) const noexcept {
  return BisquareSecondDerivativeValue(x, cc_ * scale);
}

vec RhoBisquare::SecondDerivative(const vec& x, const double scale) const noexcept {
  vec out;
  this->SecondDerivative(x, scale, &out);
  return out;
}

void RhoBisquare::SecondDerivative(const vec& x, const double scale, vec* out) const noexcept {
  const double cc_scaled = cc_ * scale;
  auto read_it = x.cbegin();
  out->copy_size(x);
  for (auto write_it = out->begin(), end = out->end(); write_it != end; ++write_it, ++read_it) {
    *write_it = BisquareSecondDerivativeValue(*read_it, cc_scaled);
  }
}

double RhoBisquare::SecondDerivativeStd(double x, const double scale) const noexcept {
  return BisquareSecondDerivativeValue(x, cc_ * scale) / UpperBound();
}

vec RhoBisquare::SecondDerivativeStd(const vec& x, const double scale) const noexcept {
  vec out;
  this->SecondDerivativeStd(x, scale, &out);
  return out;
}

void RhoBisquare::SecondDerivativeStd(const vec& x, const double scale, vec* out) const noexcept {
  const double cc_scaled = cc_ * scale;
  const double rho_inf = UpperBound();
  auto read_it = x.cbegin();
  out->copy_size(x);
  for (auto write_it = out->begin(), end = out->end(); write_it != end; ++write_it, ++read_it) {
    *write_it = BisquareSecondDerivativeValue(*read_it, cc_scaled) / rho_inf;
  }
}

double RhoBisquare::Weight(double x, const double scale) const noexcept {
  return BisquareWeightValue(x, scale * cc_);
}

vec RhoBisquare::Weight(const vec& x, const double scale) const noexcept {
  vec out;
  this->Weight(x, scale, &out);
  return out;
}

void RhoBisquare::Weight(const vec& x, const double scale, vec* out) const noexcept {
  const double cc_scaled = cc_ * scale;
  auto read_it = x.cbegin();
  out->copy_size(x);
  for (auto write_it = out->begin(), end = out->end(); write_it != end; ++write_it, ++read_it) {
    *write_it = BisquareWeightValue(*read_it, cc_scaled);
  }
}

double RhoBisquare::WeightStd(double x, const double scale) const noexcept {
  return BisquareWeightValue(x, scale * cc_) / UpperBound();
}

vec RhoBisquare::WeightStd(const vec& x, const double scale) const noexcept {
  vec out;
  this->WeightStd(x, scale, &out);
  return out;
}

void RhoBisquare::WeightStd(const vec& x, const double scale, vec* out) const noexcept {
  const double cc_scaled = cc_ * scale;
  const double rho_inf = UpperBound();
  auto read_it = x.cbegin();
  out->copy_size(x);
  for (auto write_it = out->begin(), end = out->end(); write_it != end; ++write_it, ++read_it) {
    *write_it = BisquareWeightValue(*read_it, cc_scaled) / rho_inf;
  }
}
}  // namespace pense

namespace {
inline double BisquareFunctionValueStd(double x, const double cc_scaled) noexcept {
  if (std::abs(x) > cc_scaled) {
    return 1.;
  }
  x /= cc_scaled;
  x *= x;
  return x * (3. + x * (-3. + x));
}

inline double BisquareDerivativeValue(double x, const double cc_scaled) noexcept {
  if (std::abs(x) > cc_scaled) {
    return 0.;
  }
  const double a = x / cc_scaled;
  const double u = 1. - a * a;
  return x * u * u;
}

inline double BisquareSecondDerivativeValue(double x, const double cc_scaled) noexcept {
  if (std::abs(x) > cc_scaled) {
    return 0.;
  }
  x /= cc_scaled;
  x *= x;
  return (1. - x) * (1. - 5. * x);
}

inline double BisquareWeightValue(double x, const double cc_scaled) noexcept {
  if (std::abs(x) > cc_scaled) {
    return 0.;
  }
  x /= cc_scaled;
  x = (1 - x) * (1 + x);
  return x * x;
}

inline double HuberFunctionValue(double x, const double scale, const double cc) noexcept {
  x = std::abs(x) / scale;
  if (x > cc) {
    return cc * (x - 0.5 * cc);
  }
  return 0.5 * x * x;
}

inline double HuberDerivativeValue(double x, const double scale_sq, const double cc_inv_scaled) noexcept {
  x /= scale_sq;
  if (x > cc_inv_scaled) {
    return cc_inv_scaled;
  } else if (-x > cc_inv_scaled) {
    return -cc_inv_scaled;
  }
  return x;
}

inline double HuberSecondDerivativeValue(double x, const double cc_scaled) noexcept {
  return std::abs(x) < cc_scaled ? 1 : 0;
}

inline double HuberWeightValue(double x, const double cc_scaled) noexcept {
  x = std::abs(x);
  if (x > cc_scaled) {
    return cc_scaled / x;
  }
  return 1;
}
}  // namespace
