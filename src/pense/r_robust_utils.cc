//
//  r_robust_utils.cc
//  pense
//
//  Created by David Kepplinger on 2019-04-03.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#include "r_robust_utils.hpp"

#include "constants.hpp"
#include "rcpp_integration.hpp"
#include "r_interface_utils.hpp"
#include "alias.hpp"
#include "robust_scale_location.hpp"

using Rcpp::as;
using pense::Mscale;
using pense::RhoBisquare;
using pense::RhoHuber;
using pense::LocationScaleEstimate;
using pense::GetFallback;
using pense::MLocationScale;
using pense::RhoFunctionType;
using pense::kDefaultHuberLocationCc;
using pense::kDefaultBisquareLocationCc;

namespace {
template<typename T>
LocationScaleEstimate GenericMLocationScale(const arma::vec& x, const Mscale<T>& mscale,
                                             const Rcpp::List& location_opts) {
  switch (static_cast<RhoFunctionType>(GetFallback(location_opts, "rho",
                                                   static_cast<int>(RhoFunctionType::kRhoBisquare)))) {
    case RhoFunctionType::kRhoHuber:
      return MLocationScale(x, mscale, RhoHuber(GetFallback(location_opts, "cc", kDefaultHuberLocationCc)));
    case RhoFunctionType::kRhoBisquare:
    default:
      return MLocationScale(x, mscale, RhoBisquare(GetFallback(location_opts, "cc", kDefaultBisquareLocationCc)));
  }
}

constexpr int kDefaultMLocationMaxIt = 100;
}  // namespace

namespace pense {
namespace r_interface {
//! Compute the tau-Scale of Centered Values
//!
//! @param x numeric values.
//! @return the tau-scale of the centered values.
SEXP TauSize(SEXP r_x) noexcept {
  BEGIN_RCPP
  auto x = MakeVectorView(r_x);
  return Rcpp::wrap(pense::TauSize(*x));
  END_RCPP;
}

//! Compute the M-scale of Centered Values
//!
//! @param x numeric values.
//! @param mscale_opts a list of options for the M-scale equation.
//! @return the M-scale of the centered values.
SEXP MScale(SEXP r_x, SEXP r_mscale_opts) noexcept {
  BEGIN_RCPP
  auto x = MakeVectorView(r_x);
  auto mscale_opts = as<Rcpp::List>(r_mscale_opts);
  switch (static_cast<RhoFunctionType>(GetFallback(mscale_opts, "rho",
                                                   static_cast<int>(RhoFunctionType::kRhoBisquare)))) {
    case RhoFunctionType::kRhoBisquare:
    default:
      return Rcpp::wrap(Mscale<RhoBisquare>(mscale_opts)(*x));
  }
  END_RCPP;
}

//! Compute the derivative M-scale function with respect to each coordinate.
//!
//! @param x numeric values.
//! @param mscale_opts a list of options for the M-scale equation.
//! @param order the order of the derivative to compute (1 or 2)
//! @return the derivative of the M-scale function.
SEXP MScaleDerivative(SEXP r_x, SEXP r_mscale_opts, SEXP r_order) noexcept {
  BEGIN_RCPP
  auto x = MakeVectorView(r_x);
  auto mscale_opts = as<Rcpp::List>(r_mscale_opts);
  auto order = as<int>(r_order);
  switch (static_cast<RhoFunctionType>(GetFallback(mscale_opts, "rho",
                                                   static_cast<int>(RhoFunctionType::kRhoBisquare)))) {
    case RhoFunctionType::kRhoBisquare:
    default:
      switch (order) {
      case 2:
        return Rcpp::wrap(Mscale<RhoBisquare>(mscale_opts).GradientHessian(*x));
      case 1:
      default:
        return Rcpp::wrap(Mscale<RhoBisquare>(mscale_opts).Derivative(*x));
      }
  }
  END_RCPP;
}

//! Compute the maximum derivative of M-scale function over a grid of values
//!
//! @param x original numeric values.
//! @param grid grid of values to look for maximal derivative.
//! @param change number of elements in `x` to change.
//! @param mscale_opts a list of options for the M-scale equation.
//! @return the derivative of the M-scale function.
SEXP MaxMScaleDerivative(SEXP r_x, SEXP r_grid, SEXP r_change, SEXP r_mscale_opts) noexcept {
  BEGIN_RCPP
  auto x = as<arma::vec>(r_x);
  auto grid = MakeVectorView(r_grid);
  auto change = as<int>(r_change);
  auto mscale_opts = as<Rcpp::List>(r_mscale_opts);
  switch (static_cast<RhoFunctionType>(GetFallback(mscale_opts, "rho",
                                                   static_cast<int>(RhoFunctionType::kRhoBisquare)))) {
    case RhoFunctionType::kRhoBisquare:
    default:
      auto mscale = Mscale<RhoBisquare>(mscale_opts);
      const auto derivatives = mscale.Derivative(x);
      double max_md = 0;
      if (derivatives.n_elem > 0) {
        max_md = arma::max(arma::abs(derivatives));
      }

      arma::uvec counters(change, arma::fill::zeros);
      int p = 0;
      do {
        for (int i = 0; i < change; ++i) {
          x[i] = grid->at(counters[i]);
        }
        const auto derivatives = mscale.Derivative(x);
        if (derivatives.n_elem > 0) {
          const double md = arma::max(arma::abs(derivatives));
          if (md > max_md) {
            max_md = md;
          }
        }

        p = change - 1;
        while (p >= 0) {
          ++counters[p];
          if (counters[p] >= grid->n_elem) {
            counters[p] = 0;
            --p;
          } else {
            p = -2;
          }
        }
      } while (p == -2);

      return Rcpp::wrap(max_md);
  }
  END_RCPP;
}

//! Compute the maximum entry in the gradient and Hessian of the M-scale
//! function over a grid of values
//!
//! @param x original numeric values.
//! @param grid grid of values to look for maximal derivative.
//! @param change number of elements in `x` to change.
//! @param mscale_opts a list of options for the M-scale equation.
//! @return a vector with 4 elements:
//!   [0] the maximum element of the gradient,
//!   [1] the maximum element of the Hessian,
//!   [2] the M-scale associated with the maximum gradient,
//!   [3] the M-scale associated with the maximum Hessian.
SEXP MaxMScaleGradientHessian(SEXP r_x, SEXP r_grid, SEXP r_change,
                              SEXP r_mscale_opts) noexcept {
  BEGIN_RCPP
  auto x = as<arma::vec>(r_x);
  auto grid = MakeVectorView(r_grid);
  auto change = as<int>(r_change);
  auto mscale_opts = as<Rcpp::List>(r_mscale_opts);
  const auto rho_fun = static_cast<RhoFunctionType>(
    GetFallback(mscale_opts, "rho",
                static_cast<int>(RhoFunctionType::kRhoBisquare)));

  switch (rho_fun) {
  case RhoFunctionType::kRhoBisquare:
  default:
    auto mscale = Mscale<RhoBisquare>(mscale_opts);
    const auto tmp_maxima = mscale.MaxGradientHessian(x);

    if (tmp_maxima.n_elem < 1) {
      return Rcpp::wrap(tmp_maxima);
    }

    arma::vec::fixed<4> maxima = { tmp_maxima[1], tmp_maxima[2],
                                   tmp_maxima[0], tmp_maxima[0] };

    if (change < 1) {
      return Rcpp::wrap(maxima);
    }

    arma::uvec counters(change, arma::fill::zeros);
    int p = 0;
    do {
      for (int i = 0; i < change; ++i) {
        x[i] = grid->at(counters[i]);
      }
      const auto tmp_maxima = mscale.MaxGradientHessian(x);
      if (tmp_maxima.n_elem == 3) {
        if (tmp_maxima[1] > maxima[0]) {
          maxima[0] = tmp_maxima[1];
          maxima[2] = tmp_maxima[0];
        }
        if (tmp_maxima[2] > maxima[1]) {
          maxima[1] = tmp_maxima[2];
          maxima[3] = tmp_maxima[0];
        }
      }

      p = change - 1;
      while (p >= 0) {
        ++counters[p];
        if (counters[p] >= grid->n_elem) {
          counters[p] = 0;
          --p;
        } else {
          p = -2;
        }
      }
    } while (p == -2);

    return Rcpp::wrap(maxima);
  }
  END_RCPP;
}


//! Compute the M-location
//!
//! @param x numeric values.
//! @param scale the scale of the values.
//! @param opts a list of options for the M-estimating equation.
//! @return the M-estimate of location.
SEXP MLocation(SEXP r_x, SEXP r_scale, SEXP r_opts) noexcept {
  BEGIN_RCPP
  auto x = MakeVectorView(r_x);
  auto opts = as<Rcpp::List>(r_opts);
  double const * const scale = REAL(r_scale);
  const int max_it = GetFallback(opts, "max_it", kDefaultMLocationMaxIt);
  const double convergence_tol = GetFallback(opts, "eps", kDefaultConvergenceTolerance);

  switch (static_cast<RhoFunctionType>(GetFallback(opts, "rho", static_cast<int>(RhoFunctionType::kRhoBisquare)))) {
    case RhoFunctionType::kRhoHuber:
      return Rcpp::wrap(MLocation(*x, RhoHuber(GetFallback(opts, "cc", kDefaultHuberLocationCc)), *scale,
                                  convergence_tol, max_it));
    case RhoFunctionType::kRhoBisquare:
    default:
      return Rcpp::wrap(MLocation(*x, RhoBisquare(GetFallback(opts, "cc", kDefaultBisquareLocationCc)), *scale,
                                  convergence_tol, max_it));
  }
  END_RCPP;
}

//! Compute the M-estimate of the Location and Scale
//!
//! @param x numeric values.
//! @param mscale_opts a list of options for the scale rho-function and the M-estimating equation.
//! @param location_opts a list of options for the location rho-function
//! @return a vector with 2 elements: the location and the scale estimate.
SEXP MLocationScale(SEXP r_x, SEXP r_mscale_opts, SEXP r_location_opts) noexcept {
  BEGIN_RCPP
  auto x = MakeVectorView(r_x);
  auto mscale_opts = as<Rcpp::List>(r_mscale_opts);
  auto location_opts = as<Rcpp::List>(r_location_opts);
  LocationScaleEstimate m_loc_scale;

  switch (static_cast<RhoFunctionType>(GetFallback(mscale_opts, "rho",
                                                   static_cast<int>(RhoFunctionType::kRhoBisquare)))) {
    case RhoFunctionType::kRhoBisquare:
    default:
      m_loc_scale = GenericMLocationScale(*x, Mscale<RhoBisquare>(mscale_opts), location_opts);
  }

  Rcpp::NumericVector ret_vec;
  ret_vec["location"] = m_loc_scale.location;
  ret_vec["scale"] = m_loc_scale.scale;
  return Rcpp::wrap(ret_vec);
  END_RCPP;
}
}  // namespace r_interface
}  // namespace pense
