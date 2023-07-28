//
//  r_en_regression.cc
//  pense
//
//  Created by David Kepplinger on 2019-04-03.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#include "r_en_regression.hpp"

#include <type_traits>

#include "constants.hpp"
#include "rcpp_integration.hpp"
#include "r_interface_utils.hpp"
#include "alias.hpp"
#include "regularization_path.hpp"

using Rcpp::as;
using Rcpp::wrap;

using nsoptim::Metrics;
using nsoptim::LsRegressionLoss;
using nsoptim::WeightedLsRegressionLoss;
using nsoptim::EnPenalty;
using nsoptim::AdaptiveEnPenalty;
using nsoptim::RidgePenalty;
using nsoptim::LinearizedAdmmOptimizer;
using nsoptim::AugmentedLarsOptimizer;
using nsoptim::DalEnOptimizer;
using nsoptim::CoordinateDescentOptimizer;
using SparseCoefs = nsoptim::RegressionCoefficients<arma::sp_vec>;
using DenseCoefs = nsoptim::RegressionCoefficients<arma::vec>;

template<class LossFunction>
using AugmentedRidgeOptimizer = nsoptim::AugmentedLarsOptimizer<LossFunction, nsoptim::RidgePenalty,
                                                                nsoptim::RegressionCoefficients<arma::vec>>;

using pense::alias::ConstRegressionDataPtr;
using pense::alias::FwdList;
using pense::r_interface::MakePredictorResponseData;
using pense::r_interface::MakeVectorView;
using pense::r_interface::MakePenalties;
using pense::r_interface::MakeOptimizer;

namespace {
//! Create a loss with the observation weights taken from the list of optional arguments.
template<typename Optimizer>
typename Optimizer::LossFunction GetLoss(SEXP r_x, SEXP r_y, SEXP r_include_intercept, const Rcpp::List& optional_args,
                                         std::true_type /* is_weighted */) {
  using LossFunction = typename Optimizer::LossFunction;
  ConstRegressionDataPtr data(MakePredictorResponseData(r_x, r_y));
  return LossFunction(data, MakeVectorView(optional_args["obs_weights"]), as<bool>(r_include_intercept));
}

//! Create a loss function without observation weights.
template<typename Optimizer>
typename Optimizer::LossFunction GetLoss(SEXP r_x, SEXP r_y, SEXP r_include_intercept, const Rcpp::List& optional_args,
                                         std::false_type /* is_weighted */) {
  using LossFunction = typename Optimizer::LossFunction;
  ConstRegressionDataPtr data(MakePredictorResponseData(r_x, r_y));
  return LossFunction(data, as<bool>(r_include_intercept));
}

//! Compute the EN regression regularization path.
template<typename Optimizer>
SEXP LsEnRegressionImpl(SEXP r_x, SEXP r_y, SEXP r_penalties, SEXP r_include_intercept,
                        const Rcpp::List& optional_args) {
  using is_weighted_tag = typename nsoptim::traits::is_weighted<typename Optimizer::LossFunction>::type;
  const auto en_options = pense::GetFallback(optional_args, "en_options", Rcpp::List());

  auto loss = GetLoss<Optimizer>(r_x, r_y, r_include_intercept, optional_args, is_weighted_tag{});
  auto penalties = MakePenalties<Optimizer>(r_penalties, optional_args);
  auto optimizer = MakeOptimizer<Optimizer>(en_options);

  Rcpp::List output_reg_path;
  Metrics all_metrics("reg_path");
  pense::RegPath0<Optimizer> reg_path(optimizer, loss, penalties);
  while (!reg_path.End()) {
    auto next_optimum = reg_path.Next();
    if (next_optimum.metrics) {
      all_metrics.AddSubMetrics(std::move(*next_optimum.metrics));
      next_optimum.metrics.reset();
    }
    output_reg_path.push_back(pense::WrapOptimum(next_optimum));
  }
  Rcpp::List en_results = Rcpp::List::create(Rcpp::Named("estimates") = output_reg_path,
                                             Rcpp::Named("metrics") = wrap(all_metrics));

  return wrap(en_results);
}

template<typename LossFunction, typename PenaltyFunction>
SEXP LsEnRegressionDispatch(SEXP x, SEXP y, SEXP penalties, SEXP include_intercept,
                            const Rcpp::List& optional_args) {
  const auto en_options = pense::GetFallback(optional_args, "en_options", Rcpp::List());
  const bool use_sparse_coefs = pense::GetFallback(en_options, "sparse", pense::kDefaultUseSparse);
  switch (pense::GetFallback(en_options, "algorithm", pense::kDefaultEnAlgorithm)) {
    case pense::EnAlgorithm::kDal:
      // If using the DAL optimizer, always use sparse coefficients.
      return LsEnRegressionImpl<DalEnOptimizer<LossFunction, PenaltyFunction>>(
        x, y, penalties, include_intercept, optional_args);
    case pense::EnAlgorithm::kRidge:
      // If using the ridge optimizer, always use the EN penalty (there's not adaptiveness for Ridge) and dens
      // coefficients.
        return LsEnRegressionImpl<AugmentedRidgeOptimizer<LossFunction>>(
          x, y, penalties, include_intercept, optional_args);
    case pense::EnAlgorithm::kLars:
      if (use_sparse_coefs) {
        return LsEnRegressionImpl<AugmentedLarsOptimizer<LossFunction, PenaltyFunction, SparseCoefs>>(
          x, y, penalties, include_intercept, optional_args);
      } else {
        return LsEnRegressionImpl<AugmentedLarsOptimizer<LossFunction, PenaltyFunction, DenseCoefs>>(
          x, y, penalties, include_intercept, optional_args);
      }
    case pense::EnAlgorithm::kCoordinateDescent:
      if (use_sparse_coefs) {
        return LsEnRegressionImpl<CoordinateDescentOptimizer<LossFunction, PenaltyFunction, SparseCoefs>>(
          x, y, penalties, include_intercept, optional_args);
      } else {
        return LsEnRegressionImpl<CoordinateDescentOptimizer<LossFunction, PenaltyFunction, DenseCoefs>>(
          x, y, penalties, include_intercept, optional_args);
      }
    case pense::EnAlgorithm::kLinearizedAdmm:
    default:
      if (use_sparse_coefs) {
        return LsEnRegressionImpl<LinearizedAdmmOptimizer<LossFunction, PenaltyFunction, SparseCoefs>>(
          x, y, penalties, include_intercept, optional_args);
      } else {
        return LsEnRegressionImpl<LinearizedAdmmOptimizer<LossFunction, PenaltyFunction, DenseCoefs>>(
          x, y, penalties, include_intercept, optional_args);
      }
  }
}
}  // namespace

namespace pense {
namespace r_interface {
//! Compute the EN Regularization Path.
//!
//! @param x numeric predictor matrix with `n` rows and `p` columns.
//! @param y numeric response vector with `n` elements.
//! @param penalties a list of EN penalties with decreasing values of the lambda hyper-parameter.
//! @param include_intercept include an intercept in the loss function?
//! @param optional_args a list containing the following named items:
//!                       `en_options` ... control options for the EN algorithm
//!                       `obs_weights` ... optional vector of length `n` with non-negative observation weights.
//!                       `pen_loadings` ... optional vector of length `p` with non-negative penalty loadings
SEXP LsEnRegression(SEXP x, SEXP y, SEXP penalties, SEXP include_intercept, SEXP r_optional_args) noexcept {
  BEGIN_RCPP
  Rcpp::List optional_args = as<Rcpp::List>(r_optional_args);
  if (optional_args.containsElementNamed("obs_weights")) {
    if (optional_args.containsElementNamed("pen_loadings")) {
      return LsEnRegressionDispatch<WeightedLsRegressionLoss, AdaptiveEnPenalty>(x, y, penalties, include_intercept,
                                                                       optional_args);
    }

    return LsEnRegressionDispatch<WeightedLsRegressionLoss, EnPenalty>(x, y, penalties, include_intercept, optional_args);
  } else if (optional_args.containsElementNamed("pen_loadings")) {
    return LsEnRegressionDispatch<LsRegressionLoss, AdaptiveEnPenalty>(x, y, penalties, include_intercept, optional_args);
  }
  
  return LsEnRegressionDispatch<LsRegressionLoss, EnPenalty>(x, y, penalties, include_intercept, optional_args);

  END_RCPP
}

}  // namespace r_interface
}  // namespace pense
