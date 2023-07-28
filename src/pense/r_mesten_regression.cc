//
//  r_mesten_regression.cc
//  pense
//
//  Created by David Kepplinger on 2020-06-08
//  Copyright © 2020 David Kepplinger. All rights reserved.
//

#include "r_mesten_regression.hpp"

#include "rcpp_integration.hpp"
#include "r_interface_utils.hpp"
#include "alias.hpp"
#include "rho.hpp"
#include "robust_scale_location.hpp"
#include "regularization_path.hpp"
#include "constants.hpp"
#include "m_loss.hpp"

using Rcpp::as;

namespace {
constexpr double kDefaultRhoCc = 4.5;
constexpr bool kDefaultIncludeIntercept = true;
constexpr double kDefaultExploreTol = 1e-3;
constexpr double kDefaultExploreIt = 20;
constexpr double kDefaultComparisonTol = 1e-3;
constexpr int kDefaultMaxOptima = 1;
constexpr int kDefaultTracks = 1;
constexpr bool kDefaultStrategy0 = true;
constexpr int kDefaultNumberOfThreads = 1;

//! Compute the maximum lambda without penalty loadings.
double MestEnMaxGradient(const arma::mat& x, const arma::vec& weights);
//! Compute the maximum lambda with penalty loadings.
double MestEnMaxGradient(const arma::mat& x, const arma::vec& weights, std::unique_ptr<const arma::vec> loadings);

//! Compute the regularization path for the provided penalty function.
//! This function determines the appropriate LS-EN algorithm for the MM algorithm.
//! @return the regularization path.
template<typename PenaltyFunction>
SEXP MMAlgorithmDispatch(SEXP x, SEXP y, SEXP scale, SEXP penalties, SEXP mest_opts, const Rcpp::List& optional_args);

//! Compute the regularization path for the provided penalty function.
//! This function determines the appropriate LS-EN algorithm for the MM algorithm.
//! @return the regularization path.
template<typename MMOptimizer>
SEXP MestEnRegressionImpl(MMOptimizer optimizer, SEXP x, SEXP y, SEXP scale, SEXP penalties,
                          const Rcpp::List& mest_opts, const Rcpp::List& optional_args);

//! Compute the maximum lambda without penalty loadings.
double MestEnMaxGradient(const arma::mat& x, const arma::vec& weights) {
  double max_gradient = 0.;
  for (arma::uword j = 0; j < x.n_cols; ++j) {
    const double gradient_j = std::abs(arma::mean(x.col(j) % weights));
    if (gradient_j > max_gradient) {
      max_gradient = gradient_j;
    }
  }
  return max_gradient;
}

//! Compute the maximum lambda with penalty loadings.
double MestEnMaxGradient(const arma::mat& x, const arma::vec& weights, std::unique_ptr<const arma::vec> loadings) {
  double max_gradient = 0.;
  for (arma::uword j = 0; j < x.n_cols; ++j) {
    const double gradient_j = std::abs(arma::mean(x.col(j) % weights)) / loadings->at(j);
    if (gradient_j > max_gradient) {
      max_gradient = gradient_j;
    }
  }
  return max_gradient;
}

//! Compute the PENSE Regularization Path for the provided penalty function.
//! This dispatcher inspects `pense_opts` to determine the PENSE algorithm to use. Can be one of the following:
//!  * ADMM ... use ADMM directly for the S-loss and the penalty function.
//!  * MM   ... use an MM algorithm on the convex surrogate of the S-loss.
//!
//! See `PenseEnRegression` for parameter documentation.
//! @return the regularization path.
template<typename PenaltyFunction>
SEXP MMAlgorithmDispatch(SEXP x, SEXP y, SEXP scale, SEXP penalties, SEXP r_mest_opts,
                         const Rcpp::List& optional_args) {
  using pense::GetFallback;
  using Rcpp::List;
  using SparseCoefs = nsoptim::RegressionCoefficients<arma::sp_vec>;
  using DenseCoefs = nsoptim::RegressionCoefficients<arma::vec>;
  using pense::r_interface::MakeOptimizer;
  using MLoss = pense::MLoss<pense::RhoBisquare>;
  using SurrogateLoss = typename MLoss::ConvexSurrogateType;

  const auto mest_opts = as<List>(r_mest_opts);
  const auto mm_options = GetFallback(mest_opts, "algo_opts", List());
  const auto en_options = GetFallback(mm_options, "en_options", List());
  const bool use_sparse_coefs = GetFallback(mest_opts, "sparse", pense::kDefaultUseSparse);

  switch (GetFallback(en_options, "algorithm", pense::kDefaultEnAlgorithm)) {
    case pense::EnAlgorithm::kDal: {
      using InnerOptimizer = nsoptim::DalEnOptimizer<SurrogateLoss, PenaltyFunction>;
      using MMOptimizer = nsoptim::MMOptimizer<MLoss, typename InnerOptimizer::PenaltyFunction, InnerOptimizer,
                                               typename InnerOptimizer::Coefficients>;
      return MestEnRegressionImpl(MakeOptimizer<MMOptimizer>(mm_options, en_options), x, y, scale, penalties,
                                  mest_opts, optional_args);
    }
    case pense::EnAlgorithm::kRidge: {
      using InnerOptimizer = nsoptim::AugmentedLarsOptimizer<SurrogateLoss, nsoptim::RidgePenalty, DenseCoefs>;
      using MMOptimizer = nsoptim::MMOptimizer<MLoss, typename InnerOptimizer::PenaltyFunction, InnerOptimizer,
                                               typename InnerOptimizer::Coefficients>;
      return MestEnRegressionImpl(MakeOptimizer<MMOptimizer>(mm_options, en_options), x, y, scale, penalties,
                                  mest_opts, optional_args);
    }
    case pense::EnAlgorithm::kLinearizedAdmm:
      if (use_sparse_coefs) {
        using InnerOptimizer = nsoptim::LinearizedAdmmOptimizer<SurrogateLoss, PenaltyFunction, SparseCoefs>;
        using MMOptimizer = nsoptim::MMOptimizer<MLoss, typename InnerOptimizer::PenaltyFunction, InnerOptimizer,
                                                 typename InnerOptimizer::Coefficients>;
        return MestEnRegressionImpl(MakeOptimizer<MMOptimizer>(mm_options, en_options), x, y, scale, penalties,
                                    mest_opts, optional_args);
      } else {
        using InnerOptimizer = nsoptim::LinearizedAdmmOptimizer<SurrogateLoss, PenaltyFunction, DenseCoefs>;
        using MMOptimizer = nsoptim::MMOptimizer<MLoss, typename InnerOptimizer::PenaltyFunction, InnerOptimizer,
                                                 typename InnerOptimizer::Coefficients>;
        return MestEnRegressionImpl(MakeOptimizer<MMOptimizer>(mm_options, en_options), x, y, scale, penalties,
                                    mest_opts, optional_args);
      }
    case pense::EnAlgorithm::kLars:
    default:
      if (use_sparse_coefs) {
        using InnerOptimizer = nsoptim::AugmentedLarsOptimizer<SurrogateLoss, PenaltyFunction, SparseCoefs>;
        using MMOptimizer = nsoptim::MMOptimizer<MLoss, typename InnerOptimizer::PenaltyFunction, InnerOptimizer,
                                                 typename InnerOptimizer::Coefficients>;
        return MestEnRegressionImpl(MakeOptimizer<MMOptimizer>(mm_options, en_options), x, y, scale, penalties,
                                    mest_opts, optional_args);
      } else {
        using InnerOptimizer = nsoptim::AugmentedLarsOptimizer<SurrogateLoss, PenaltyFunction, DenseCoefs>;
        using MMOptimizer = nsoptim::MMOptimizer<MLoss, typename InnerOptimizer::PenaltyFunction, InnerOptimizer,
                                                 typename InnerOptimizer::Coefficients>;
        return MestEnRegressionImpl(MakeOptimizer<MMOptimizer>(mm_options, en_options), x, y, scale, penalties,
                                    mest_opts, optional_args);
      }
  }
}


//! Compute the M-EN regularization path
template<typename MMOptimizer>
SEXP MestEnRegressionImpl(MMOptimizer optimizer, SEXP r_x, SEXP r_y, SEXP r_scale, SEXP r_penalties,
                          const Rcpp::List& mest_opts, const Rcpp::List& optional_args) {
  using Rcpp::List;
  using Rcpp::wrap;
  using nsoptim::Metrics;
  using pense::r_interface::MakePredictorResponseData;
  using pense::r_interface::MakePenalties;
  using pense::r_interface::MakeOptimizer;
  using pense::GetFallback;
  using pense::alias::ConstRegressionDataPtr;
  using MLoss = pense::MLoss<pense::RhoBisquare>;
  using CoefficientsList = pense::alias::FwdList<typename MMOptimizer::Coefficients>;

  ConstRegressionDataPtr data(MakePredictorResponseData(r_x, r_y));
  const pense::RhoBisquare rho(GetFallback(mest_opts, "cc", kDefaultRhoCc));
  const double scale = as<double>(r_scale);
  MLoss loss(data, rho, scale, GetFallback(mest_opts, "intercept", kDefaultIncludeIntercept));

  auto penalties = MakePenalties<MMOptimizer>(r_penalties, optional_args);

  const double eps = GetFallback(mest_opts, "eps", pense::kDefaultConvergenceTolerance);
  optimizer.convergence_tolerance(eps);

  Metrics metrics("mstep");
  pense::RegPathCombined<MMOptimizer> reg_paths(optimizer, loss, penalties,
                                                GetFallback(mest_opts, "max_optima", kDefaultMaxOptima),
                                                GetFallback(mest_opts, "nr_tracks", kDefaultTracks),
                                                GetFallback(mest_opts, "explore_tol", kDefaultExploreTol),
                                                GetFallback(mest_opts, "explore_it", kDefaultExploreIt),
                                                GetFallback(mest_opts, "comparison_tol", kDefaultComparisonTol),
                                                GetFallback(mest_opts, "num_threads", kDefaultNumberOfThreads));

  // Enable computation of the 0-based solutions, if requested.
  if (GetFallback(mest_opts, "strategy_0", kDefaultStrategy0)) {
    reg_paths.Add();
  }

  // Enable computation of the regularization paths using shared starting points (i.e., the same starting
  // point at every penalty).
  if (optional_args.containsElementNamed("shared_starts")) {
    auto shared_starts = Rcpp::as<CoefficientsList>(optional_args["shared_starts"]);
    for (auto&& start : shared_starts) {
      reg_paths.Add(start);
    }
  }

  List combined_reg_path;
  while (!reg_paths.End()) {
    List solutions;
    Metrics& sub_metrics = metrics.CreateSubMetrics("lambda");

    // Compute the optima at the next penalty level.
    for (auto&& optimum : reg_paths.Next()) {
      if (optimum.metrics) {
        sub_metrics.AddSubMetrics(*optimum.metrics);
      }
      solutions.push_back(WrapOptimum(optimum));
    }

    combined_reg_path.push_back(solutions);

    Rcpp::checkUserInterrupt();
  }

  List mstep_results = List::create(Rcpp::Named("estimates") = combined_reg_path,
                                    Rcpp::Named("metrics") = wrap(metrics));

  return wrap(mstep_results);
}
}  // namespace

namespace pense {
namespace r_interface {
//! Compute the (Adaptive) M-EN Regularization Path.
//!
//! @param x numeric predictor matrix with `n` rows and `p` columns.
//! @param y numeric response vector with `n` elements.
//! @param scale scalar, numeric, auxiliary scale.
//! @param penalties a list of EN penalties with decreasing values of the lambda hyper-parameter.
//! @param mest_opts a list of options for the M-estimation algorithm.
//! @param optional_args a list containing the following named items:
//!                      `shared_starts` ... optional list of coefficients to start at every penalty.
//!                      `individual_starts` ... optional list the same length as `penalties` with a list coefficients
//!                                               to start at the corresponding penalty.
//!                      `pen_loadings` ... optional vector of length `p` with non-negative penalty loadings.
SEXP MestEnRegression(SEXP x, SEXP y, SEXP scale, SEXP penalties, SEXP mest_opts, SEXP r_optional_args) noexcept {
  BEGIN_RCPP
  const auto optional_args = as<Rcpp::List>(r_optional_args);

  if (optional_args.containsElementNamed("pen_loadings")) {
    return MMAlgorithmDispatch<nsoptim::AdaptiveEnPenalty>(x, y, scale, penalties, mest_opts, optional_args);
  }
  return MMAlgorithmDispatch<nsoptim::EnPenalty>(x, y, scale, penalties, mest_opts, optional_args);
  END_RCPP
}

//! Get the smallest lambda such that the (Adaptive) M-EN-estimate gives the empty model.
//!
//! @param x numeric predictor matrix with `n` rows and `p` columns.
//! @param y numeric response vector with `n` elements.
//! @param scale scalar, numeric, auxiliary scale.
//! @param mest_opts a list of options for the M-estimation algorithm.
//! @param optional_args a list containing the following named items:
//!                       `pen_loadings` ... optional vector of length `p` with non-negative penalty loadings.
SEXP MestEnMaxLambda(SEXP r_x, SEXP r_y, SEXP r_scale, SEXP r_mest_opts, SEXP r_optional_args) noexcept {
  BEGIN_RCPP
  auto data = MakePredictorResponseData(r_x, r_y);
  const auto mest_opts = as<Rcpp::List>(r_mest_opts);
  const auto optional_args = as<Rcpp::List>(r_optional_args);
  const double scale = as<double>(r_scale);
  RhoBisquare rhofun(as<double>(mest_opts["cc"]));
  const double location = MLocation(data->cy(), rhofun, scale, as<double>(mest_opts["eps"]),
                                    as<int>(mest_opts["max_it"]));
  const arma::vec residuals = data->cy() - location;
  arma::vec weights = residuals % rhofun.Weight(residuals, scale);  // dividing by the square scale can be done after!

  if (optional_args.containsElementNamed("pen_loadings")) {
    return Rcpp::wrap(MestEnMaxGradient(data->cx(), weights, MakeVectorView(optional_args["pen_loadings"])) /
                      (scale * scale));
  } else {
    return Rcpp::wrap(MestEnMaxGradient(data->cx(), weights) / (scale * scale));
  }
  END_RCPP
}
}  // namespace r_interface
}  // namespace pense
