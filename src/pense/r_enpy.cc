//
//  r_enpy.cc
//  pense
//
//  Created by David Kepplinger on 2019-04-03.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#include "r_enpy.hpp"

#include "rcpp_integration.hpp"
#include "r_interface_utils.hpp"
#include "alias.hpp"
#include "robust_scale_location.hpp"
#include "s_loss.hpp"
#include "enpy_initest.hpp"

using Rcpp::as;
using Rcpp::wrap;
using nsoptim::LsRegressionLoss;
using nsoptim::EnPenalty;
using nsoptim::AdaptiveEnPenalty;
using nsoptim::RidgePenalty;
using nsoptim::MMOptimizer;
using nsoptim::RegressionCoefficients;
using nsoptim::AugmentedLarsOptimizer;
using nsoptim::DalEnOptimizer;
using nsoptim::LinearizedAdmmOptimizer;
template<class LossFunction>
using AugmentedRidgeOptimizer = nsoptim::AugmentedLarsOptimizer<LossFunction, nsoptim::RidgePenalty,
                                                                nsoptim::RegressionCoefficients<arma::vec>>;

using pense::alias::ConstRegressionDataPtr;
using pense::alias::FwdList;
using pense::r_interface::MakePredictorResponseData;
using pense::r_interface::MakePenalties;
using pense::r_interface::MakeOptimizer;
using pense::PenaYohaiInitialEstimators;
using pense::Mscale;
using pense::RhoBisquare;
using pense::SLoss;

using SparseCoefs = nsoptim::RegressionCoefficients<arma::sp_vec>;
using DenseCoefs = nsoptim::RegressionCoefficients<arma::vec>;

namespace {
//! Generic ENPY Initial Estimator
//!
//! For parameter documentation, see the specializations below.
template<typename Optimizer>
SEXP PenPyInitialEstimatorImpl(SEXP r_x, SEXP r_y, SEXP r_penalties, SEXP r_sloss_params, const Rcpp::List& enpy_opts,
                               const Rcpp::List& en_options, const Rcpp::List& optional_args) {
  ConstRegressionDataPtr data(MakePredictorResponseData(r_x, r_y));
  const auto sloss_params = as<Rcpp::List>(r_sloss_params);
  const auto penalties = MakePenalties<Optimizer>(r_penalties, optional_args);

  Mscale<RhoBisquare> mscale(as<Rcpp::List>(sloss_params["mscale"]));
  SLoss loss(data, mscale, as<bool>(sloss_params["intercept"]));

  Optimizer optimizer = MakeOptimizer<Optimizer>(en_options);
  auto py_res = PenaYohaiInitialEstimators(loss, penalties, optimizer, enpy_opts);
  return wrap(py_res);
}

template<typename Optimizer>
SEXP WrapPscResults(const FwdList<pense::PscResult<Optimizer>>& psc_results) {
  Rcpp::List out;

  for (auto&& psc_res : psc_results) {
    out.push_back(psc_res.pscs);
  }

  return wrap(out);
}

//! Generic PSCs
//!
//! For parameter documentation, see the specializations below.
template<typename Optimizer>
SEXP PscImpl(SEXP r_x, SEXP r_y, SEXP r_penalties, const Rcpp::List& en_options, const Rcpp::List& optional_args) {
  ConstRegressionDataPtr data(MakePredictorResponseData(r_x, r_y));
  const auto penalties = MakePenalties<Optimizer>(r_penalties, optional_args);
  const bool include_intercept = pense::GetFallback(optional_args, "intercept", true);
  const int num_threads = pense::GetFallback(optional_args, "num_threads", 1);

  LsRegressionLoss loss(data, include_intercept);
  Optimizer optimizer = MakeOptimizer<Optimizer>(en_options);
  auto psc_res = pense::PrincipalSensitiviyComponents(loss, penalties, optimizer, num_threads);
  return WrapPscResults(psc_res);
}

template<typename PenaltyFunction>
SEXP PenPyInitialEstimatorDispatch(SEXP x, SEXP y, SEXP penalties, SEXP sloss_params, SEXP r_enpy_opts,
                                   const Rcpp::List& optional_args) {
  const auto enpy_opts = Rcpp::as<Rcpp::List>(r_enpy_opts);
  const auto en_options = Rcpp::as<Rcpp::List>(enpy_opts["en_options"]);
  const bool use_sparse_coefs = pense::GetFallback(en_options, "sparse", pense::kDefaultUseSparse);

  switch (pense::GetFallback(en_options, "algorithm", pense::kDefaultEnAlgorithm)) {
    case pense::EnAlgorithm::kDal:
      // If using the DAL optimizer, always use sparse coefficients.
      return PenPyInitialEstimatorImpl<DalEnOptimizer<LsRegressionLoss, PenaltyFunction>>(
        x, y, penalties, sloss_params, enpy_opts, en_options, optional_args);
    case pense::EnAlgorithm::kRidge:
      // If using the ridge optimizer, always use the Ridge penalty (there's not adaptiveness for Ridge) and dense
      // coefficients.
      return PenPyInitialEstimatorImpl<AugmentedRidgeOptimizer<LsRegressionLoss>>(
        x, y, penalties, sloss_params, enpy_opts, en_options, optional_args);
    case pense::EnAlgorithm::kLars:
      if (use_sparse_coefs) {
        return PenPyInitialEstimatorImpl<AugmentedLarsOptimizer<LsRegressionLoss, PenaltyFunction, SparseCoefs>>(
          x, y, penalties, sloss_params, enpy_opts, en_options, optional_args);
      } else {
        return PenPyInitialEstimatorImpl<AugmentedLarsOptimizer<LsRegressionLoss, PenaltyFunction, DenseCoefs>>(
          x, y, penalties, sloss_params, enpy_opts, en_options, optional_args);
      }
    case pense::EnAlgorithm::kLinearizedAdmm:
    default:
      if (use_sparse_coefs) {
        return PenPyInitialEstimatorImpl<LinearizedAdmmOptimizer<LsRegressionLoss, PenaltyFunction, SparseCoefs>>(
          x, y, penalties, sloss_params, enpy_opts, en_options, optional_args);
      } else {
        return PenPyInitialEstimatorImpl<LinearizedAdmmOptimizer<LsRegressionLoss, PenaltyFunction, DenseCoefs>>(
          x, y, penalties, sloss_params, enpy_opts, en_options, optional_args);
      }
  }
}

template<typename PenaltyFunction>
SEXP PscDispatch(SEXP x, SEXP y, SEXP penalties, SEXP r_en_options, const Rcpp::List& optional_args) {
  const auto en_options = Rcpp::as<Rcpp::List>(r_en_options);
  const bool use_sparse_coefs = pense::GetFallback(en_options, "sparse", pense::kDefaultUseSparse);

  switch (pense::GetFallback(en_options, "algorithm", pense::kDefaultEnAlgorithm)) {
    case pense::EnAlgorithm::kDal:
      // If using the DAL optimizer, always use sparse coefficients.
      return PscImpl<DalEnOptimizer<LsRegressionLoss, PenaltyFunction>>(x, y, penalties,  en_options, optional_args);
    case pense::EnAlgorithm::kRidge:
      // If using the ridge optimizer, always use the Ridge penalty (there's not adaptiveness for Ridge) and dense
      // coefficients.
      return PscImpl<AugmentedRidgeOptimizer<LsRegressionLoss>>(x, y, penalties, en_options, optional_args);
    case pense::EnAlgorithm::kLars:
      if (use_sparse_coefs) {
        return PscImpl<AugmentedLarsOptimizer<LsRegressionLoss, PenaltyFunction, SparseCoefs>>(
          x, y, penalties,  en_options, optional_args);
      } else {
        return PscImpl<AugmentedLarsOptimizer<LsRegressionLoss, PenaltyFunction, DenseCoefs>>(
          x, y, penalties, en_options, optional_args);
      }
    case pense::EnAlgorithm::kLinearizedAdmm:
    default:
      if (use_sparse_coefs) {
        return PscImpl<LinearizedAdmmOptimizer<LsRegressionLoss, PenaltyFunction, SparseCoefs>>(
          x, y, penalties, en_options, optional_args);
      } else {
        return PscImpl<LinearizedAdmmOptimizer<LsRegressionLoss, PenaltyFunction, DenseCoefs>>(
          x, y, penalties, en_options, optional_args);
      }
  }
}

}  // namespace

namespace pense {
namespace r_interface {
//! Compute the penalized PY Initial Estimators.
//!
//! @param x numeric predictor matrix with `n` rows and `p` columns.
//! @param y numeric response vector with `n` elements.
//! @param penalties a list of EN penalties with decreasing values of the lambda hyper-parameter.
//! @param sloss_params parameters for the M-scale.
//! @param enpy_opts a list of options for the ENPY algorithm.
//! @param optional_args a list containing the following named items:
//!                       `pen_loadings` ... optional vector of length `p` with non-negative penalty loadings.
SEXP PenPyInitialEstimator(SEXP x, SEXP y, SEXP penalties, SEXP sloss_params, SEXP enpy_opts,
                           SEXP r_optional_args) noexcept {
  BEGIN_RCPP
  auto optional_args = as<Rcpp::List>(r_optional_args);

  if (optional_args.containsElementNamed("pen_loadings")) {
    return PenPyInitialEstimatorDispatch<AdaptiveEnPenalty>(x, y, penalties, sloss_params, enpy_opts, optional_args);
  }
  return PenPyInitialEstimatorDispatch<EnPenalty>(x, y, penalties, sloss_params, enpy_opts, optional_args);
  END_RCPP
}

//! Compute the Principal Sensitivity Components.
//!
//! @param x numeric predictor matrix with `n` rows and `p` columns.
//! @param y numeric response vector with `n` elements.
//! @param penalties a list of EN penalties with decreasing values of the lambda hyper-parameter.
//! @param en_options options for the EN algorithm.
//! @param optional_args a list containing the following named items:
//!                       `intercept` ... boolean determining if an intercept should be included.
//!                       `num_threads` ... number of threads.
//!                       `pen_loadings` ... optional vector of length `p` with non-negative penalty loadings.
SEXP PrincipalSensitivityComponents(SEXP x, SEXP y, SEXP penalties, SEXP en_options, SEXP r_optional_args) noexcept {
  BEGIN_RCPP
  auto optional_args = as<Rcpp::List>(r_optional_args);

  if (optional_args.containsElementNamed("pen_loadings")) {
    return PscDispatch<AdaptiveEnPenalty>(x, y, penalties, en_options, optional_args);
  }
  return PscDispatch<EnPenalty>(x, y, penalties, en_options, optional_args);
  END_RCPP
}

}  // namespace r_interface
}  // namespace pense
