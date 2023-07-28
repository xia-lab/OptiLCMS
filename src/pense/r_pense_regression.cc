//
//  r_pense_regression.cc
//  pense
//
//  Created by David Kepplinger on 2019-04-03.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#include "r_pense_regression.hpp"

#include <unordered_set>

#include "rcpp_integration.hpp"
#include "r_interface_utils.hpp"
#include "alias.hpp"
#include "enpy_initest.hpp"
#include "robust_scale_location.hpp"
#include "s_loss.hpp"
#include "cd_pense.hpp"
#include "regularization_path_new.hpp"
#include "constants.hpp"

using Rcpp::as;
using nsoptim::LsRegressionLoss;
using nsoptim::Metrics;
template<class LossFunction>
using AugmentedRidgeOptimizer = nsoptim::AugmentedLarsOptimizer<
  LossFunction, nsoptim::RidgePenalty, nsoptim::RegressionCoefficients<arma::vec>>;

using pense::GetFallback;
using pense::PyResult;
using pense::alias::Optima;
using pense::alias::FwdList;
using pense::alias::ConstRegressionDataPtr;
using pense::r_interface::MakePredictorResponseData;
using pense::r_interface::MakePenalties;
using pense::r_interface::MakeOptimizer;
using pense::SLoss;

using SparseCoefs = nsoptim::RegressionCoefficients<arma::sp_vec>;
using DenseCoefs = nsoptim::RegressionCoefficients<arma::vec>;

template<typename Optimizer>
using CoefficientsList = FwdList<typename Optimizer::Coefficients>;

template<typename Optimizer>
using StartCoefficientsList = FwdList<CoefficientsList<Optimizer>>;

template<typename Optimizer>
using PenaltyList = FwdList<typename Optimizer::PenaltyFunction>;

namespace {
constexpr double kDefaultExploreTol = 0.1;
constexpr double kDefaultComparisonTol = 1e-3;
constexpr double kDefaultExploreIt = 20;
constexpr int kDefaultMaxOptima = 10;
constexpr int kDefaultExploreSolutions = 10;
constexpr bool kDefaultUseWarmStarts = true;
constexpr bool kDefaultStrategy0 = true;
constexpr bool kDefaultStrategyEnpyShared = true;
constexpr bool kDefaultStrategyEnpyIndividual = false;
constexpr bool kDefaultStrategyOtherShared = true;
constexpr bool kDefaultStrategyOtherIndividual = false;
constexpr int kDefaultNumberOfThreads = 1;

//! Expand a list of PY Results to a list of start coefficients.
//!
//! @param py_res the results of the ENPY algorithm.
//! @param penalties a list of the penalties.
//! @param r_indices vector of 1-based indices in the penalties list at which the PY estimates
//!     were computed.
//! @return a list the same length and order as `penalties`.
template<typename Optimizer>
StartCoefficientsList<Optimizer> PyResultToStartCoefficients(
     const FwdList<PyResult<Optimizer>>& py_res, const PenaltyList<Optimizer>& penalties,
     SEXP r_indices) {
  const Rcpp::IntegerVector indices(r_indices);
  StartCoefficientsList<Optimizer> start_coefs;

  auto start_coefs_it = start_coefs.before_begin();

  auto py_res_it = py_res.begin();
  auto index_it = indices.cbegin();
  auto index_end = indices.cend();
  int list_index = 1;
  for (auto it = penalties.cbegin(), end = penalties.cend(); it != end; ++it, ++list_index) {
    // Add empty list of coefficients at the end.
    start_coefs_it = start_coefs.emplace_after(start_coefs_it);

    if ((index_it != index_end) && (list_index == *index_it)) {
      // Copy coefficients
      auto coefs_list_it = start_coefs_it->before_begin();
      for (auto&& optimum : py_res_it->initial_estimates) {
        coefs_list_it = start_coefs_it->insert_after(coefs_list_it, optimum.coefs);
      }
      ++index_it;
      ++py_res_it;
    }
  }
  return start_coefs;
}

//! Add an estimate to the list of solutions.
//!
//! The coefficients at `iterator` are added to the list of solutions and the metrics at
//! `iterator` are added to the given Metrics object. The iterator will not be incremented
//! past the end.
//!
//! @param iterator iterator pointing to the estimate(s) to add.
//! @param end end marker for the given iterator.
//! @param metrics metrics object to which to add the estimate's metrics.
//! @param solutions a list to add the estimate to.
//! @return the updated iterator (if not yet pointing to the end).
template<typename InputIterator>
InputIterator AddEstimate(InputIterator iterator, InputIterator end, Metrics* metrics,
                          Rcpp::List* solutions, const Rcpp::String& origin) {
  if (iterator != end) {
    if (iterator->metrics) {
      metrics->AddSubMetrics(std::move(*(iterator->metrics)));
      iterator->metrics.reset();
    }
    solutions->push_back(WrapOptimum(*iterator));
    return ++iterator;
  }
  return iterator;
}

//! Add a list of estimates to the list of solutions.
//!
//! The optima at `iterator` are added to the list of solutions and the metrics aof all optima
//! `iterator` are added to the given Metrics object. The iterator will not be incremented
//! past the end.
//!
//! @param iterator iterator pointing to the estimates to add.
//! @param end end marker for the given iterator.
//! @param metrics metrics object to which to add the optimas' metrics.
//! @param solutions a list to add the estimate to.
//! @return the updated iterator (if not yet pointing to the end).
template<typename InputIterator>
InputIterator AddEstimates(InputIterator iterator, InputIterator end, Metrics* metrics, Rcpp::List* solutions,
                           const Rcpp::String& origin) {
  if (iterator != end) {
    for (auto list_it = iterator->begin(), list_end = iterator->end(); list_it != list_end;) {
      list_it = AddEstimate(list_it, list_end, metrics, solutions, origin);
    }
    return ++iterator;
  }
  return iterator;
}

//! Compute the ENPY initial estimates using the specified `LsOptimizer` class.
//! This implementation of the function is used if the `LsOptimizer` class can not handle the desired penalty function
//! or the desired coefficients type.
//!
//! @return an empty list of start coefficients.
template<typename LsOptimizer, typename SOptimizer>
StartCoefficientsList<SOptimizer> EnpyInitialEstimatesImpl(const SLoss&,
                                                           const FwdList<typename SOptimizer::PenaltyFunction>&,
                                                           SEXP, SEXP, const Rcpp::List&, const Rcpp::List&,
                                                           const Rcpp::List&, Metrics * const, double) {
  return StartCoefficientsList<SOptimizer>();
}

//! Compute the ENPY initial estimates using the specified `LsOptimizer` class.
//! This implementation of the function is used if the `LsOptimizer` can handle the desired penalty function and the
//! desired coefficients type.
//!
//! @return a list of start coefficients. The returned list is either empty or contains lists of start coefficients for
//!         *each* penalty parameter.
template<typename LsOptimizer, typename SOptimizer, typename = typename
         std::enable_if<std::is_same<typename LsOptimizer::Coefficients, typename SOptimizer::Coefficients>::value &&
                        std::is_same<typename LsOptimizer::PenaltyFunction,
                                     typename SOptimizer::PenaltyFunction>::value>::type>
StartCoefficientsList<SOptimizer> EnpyInitialEstimatesImpl(
    const SLoss& loss, const FwdList<typename SOptimizer::PenaltyFunction>& penalties, SEXP r_penalties,
    SEXP r_enpy_inds, const Rcpp::List& enpy_opts, const Rcpp::List& en_options, const Rcpp::List& optional_args,
    Metrics * const metrics, int) {
  const auto enpy_penalties = MakePenalties<LsOptimizer>(r_penalties, r_enpy_inds, optional_args);

  if (enpy_penalties.empty()) {
    return StartCoefficientsList<LsOptimizer>();
  } else {
    auto py_res = PenaYohaiInitialEstimators(loss, enpy_penalties, MakeOptimizer<LsOptimizer>(en_options), enpy_opts);

    // Move metrics from the PY results.
    auto&& enpy_metrics = metrics->CreateSubMetrics("enpy_initest");
    for (auto&& single_py_res : py_res) {
      enpy_metrics.AddSubMetrics(std::move(single_py_res.metrics));
    }

    return PyResultToStartCoefficients(py_res, penalties, r_enpy_inds);
  }
}

//! Compute the ENPY initial estimates using the options provided in `enpy_opts`.
//! This function inspects `enpy_opts["en_options"]` to determine which LS-EN algorithm to use and only returns a
//! non-empty list of start coefficients if LS-EN algorithm is compatible with both the coefficients and the penalty
//! function used by the specified `SOptimizer` class.
//!
//! @return a list of start coefficients. The returned list is either empty or contains lists of start coefficients for
//!         *each* penalty parameter.
template<typename SOptimizer>
StartCoefficientsList<SOptimizer> EnpyInitialEstimates(const SLoss& loss,
                                                       const FwdList<typename SOptimizer::PenaltyFunction>& penalties,
                                                       SEXP r_penalties, SEXP r_enpy_inds, SEXP r_enpy_opts,
                                                       const Rcpp::List& optional_args, Metrics * metrics) {
  using PenaltyFunction = typename SOptimizer::PenaltyFunction;
  using Coefficients = typename SOptimizer::Coefficients;
  using LsEnDal = nsoptim::DalEnOptimizer<LsRegressionLoss, PenaltyFunction>;
  using LsRidge = AugmentedRidgeOptimizer<LsRegressionLoss>;
  using LsEnAdmm = nsoptim::LinearizedAdmmOptimizer<LsRegressionLoss, PenaltyFunction, Coefficients>;
  using LsEnLars = nsoptim::AugmentedLarsOptimizer<LsRegressionLoss, PenaltyFunction, Coefficients>;

  const auto enpy_opts = as<Rcpp::List>(r_enpy_opts);
  const auto en_options = as<Rcpp::List>(enpy_opts["en_options"]);

  switch (GetFallback(en_options, "algorithm", pense::kDefaultEnAlgorithm)) {
    case pense::EnAlgorithm::kDal:
      return EnpyInitialEstimatesImpl<LsEnDal, SOptimizer>(loss, penalties, r_penalties, r_enpy_inds, enpy_opts,
                                                           en_options, optional_args, metrics, 1);
    case pense::EnAlgorithm::kRidge:
      return EnpyInitialEstimatesImpl<LsRidge, SOptimizer>(loss, penalties, r_penalties, r_enpy_inds, enpy_opts,
                                                           en_options, optional_args, metrics, 1);
    case pense::EnAlgorithm::kLinearizedAdmm:
      return EnpyInitialEstimatesImpl<LsEnAdmm, SOptimizer>(loss, penalties, r_penalties, r_enpy_inds, enpy_opts,
                                                            en_options, optional_args, metrics, 1);
    case pense::EnAlgorithm::kLars:
    default:
      return EnpyInitialEstimatesImpl<LsEnLars, SOptimizer>(loss, penalties, r_penalties, r_enpy_inds, enpy_opts,
                                                            en_options, optional_args, metrics, 1);
  }
}

//! Compute the PENSE Regularization Path using the provided optimizer for the PENSE objective function.
//! The optimizer determines:
//!   * the type of penalty function (EN vs. adaptive EN)
//!   * the coefficients type (sparse vs. dense)
//!
//! See `PenseEnRegression` for parameter documentation.
//! @return the regularization path.
template<class SOptimizer>
SEXP PenseRegressionImpl(SOptimizer optimizer, SEXP r_x, SEXP r_y, SEXP r_penalties,
                         SEXP r_enpy_inds, const Rcpp::List& pense_opts, SEXP r_enpy_opts,
                         const Rcpp::List& optional_args) {
  using RList = Rcpp::List;

  ConstRegressionDataPtr data(MakePredictorResponseData(r_x, r_y));

  pense::Mscale<pense::RhoBisquare> mscale(as<RList>(pense_opts["mscale"]));
  SLoss loss(data, mscale, as<bool>(pense_opts["intercept"]));
  auto penalties = MakePenalties<SOptimizer>(r_penalties, optional_args);

  const double eps = GetFallback(pense_opts, "eps", pense::kDefaultConvergenceTolerance);
  optimizer.convergence_tolerance(eps);
  optimizer.loss(loss);

  Metrics metrics("pense");
  pense::RegularizationPath<SOptimizer> reg_path(
    optimizer, penalties,
    GetFallback(pense_opts, "max_optima", kDefaultMaxOptima),
    GetFallback(pense_opts, "comparison_tol", kDefaultComparisonTol),
    GetFallback(pense_opts, "num_threads", kDefaultNumberOfThreads));

  reg_path.ExplorationOptions(GetFallback(pense_opts, "explore_it", kDefaultExploreIt),
                              GetFallback(pense_opts, "explore_tol", kDefaultExploreTol),
                              GetFallback(pense_opts, "nr_tracks", kDefaultExploreSolutions));

  reg_path.EnableWarmStarts(GetFallback(pense_opts, "warm_starts", kDefaultUseWarmStarts));

  // Compute the initial estimators
  auto&& cold_starts = EnpyInitialEstimates<SOptimizer>(loss, penalties, r_penalties, r_enpy_inds, r_enpy_opts,
                                                        optional_args, &metrics);

  Rcpp::checkUserInterrupt();

  // Enable computation of EN-PY-based solutions
  if (!cold_starts.empty()) {
    if (GetFallback(pense_opts, "strategy_enpy_individual", kDefaultStrategyEnpyIndividual)) {
      // Use the EN-PY solutions only for the penalty they were computed for.
      reg_path.EmplaceIndividualStartingPoints(std::move(cold_starts));
    }
    if (GetFallback(pense_opts, "strategy_enpy_shared", kDefaultStrategyEnpyShared)) {
      // Use every EN-PY solution for all penalties.
      for (auto&& starts_at_lambda : cold_starts) {
        for (auto&& start : starts_at_lambda) {
          reg_path.EmplaceSharedStartingPoint(std::move(start));
        }
      }
    }
  }

  // Enable computation of the 0-based solutions, if requested.
  if (GetFallback(pense_opts, "strategy_0", kDefaultStrategy0)) {
    StartCoefficientsList<SOptimizer> zeros;
    zeros.emplace_front(1, loss.ZeroCoefficients<typename SOptimizer::Coefficients>());
    reg_path.EmplaceIndividualStartingPoints(std::move(zeros));
  }

  // Enable computation of the regularization paths using shared starting points (i.e., the same starting
  // point at every penalty).
  if (GetFallback(pense_opts, "strategy_other_shared", kDefaultStrategyOtherShared) &&
      optional_args.containsElementNamed("shared_starts")) {
    auto shared_starts = as<CoefficientsList<SOptimizer>>(optional_args["shared_starts"]);
    for (auto&& start : shared_starts) {
      reg_path.EmplaceSharedStartingPoint(std::move(start));
    }
  }

  // Enable computation of the regularization paths using individual starting points (i.e., a list of starting
  // points, different for every penalty).
  if (GetFallback(pense_opts, "strategy_other_individual", kDefaultStrategyOtherIndividual) &&
      optional_args.containsElementNamed("individual_starts")) {
    reg_path.EmplaceIndividualStartingPoints(
      as<StartCoefficientsList<SOptimizer>>(optional_args["individual_starts"]));
  }

  RList combined_reg_path;
  while (!reg_path.End()) {
    RList solutions;
    Metrics& sub_metrics = metrics.CreateSubMetrics("lambda");

    // Compute the optima at the next penalty level.
    auto next = reg_path.Next();

    sub_metrics.AddMetric("alpha", next.penalty.alpha());
    sub_metrics.AddMetric("lambda", next.penalty.lambda());

    for (auto&& optimum : next.optima) {
      if (optimum.metrics) {
        optimum.metrics->AddDetail("objf_value", optimum.objf_value);
        sub_metrics.AddSubMetrics(*optimum.metrics);
      }
      solutions.push_back(WrapOptimum(optimum));
    }

    combined_reg_path.push_back(solutions);

    Rcpp::checkUserInterrupt();
  }

  RList pense_results = RList::create(Rcpp::Named("estimates") = combined_reg_path,
                                      Rcpp::Named("metrics") = Rcpp::wrap(metrics));

  return Rcpp::wrap(pense_results);
}

//! Compute the PENSE Regularization Path for the provided penalty function using the MM algorithm with the provided
//! inner optimizer.
//! By default, i.e., unless the inner optimizer can handle the given PenaltyFunction class, returns `R_NilValue`.
//!
//! See `PenseEnRegression` for parameter documentation.
//! @return `R_NilValue`.
template<typename InnerOptimizer, typename PenaltyFunction>
SEXP PenseMMPenaltyImpl(SEXP, SEXP, SEXP, SEXP, const Rcpp::List&, SEXP, const Rcpp::List&,
                        const Rcpp::List&, const Rcpp::List&, double) {
  return R_NilValue;
}

//! Compute the PENSE Regularization Path for the provided penalty function using the MM algorithm with the provided
//! inner optimizer.
//! This function is only enabled if the inner optimizer can handle the given PenaltyFunction class.
//!
//! See `PenseEnRegression` for parameter documentation.
//! @return the regularization path.
template<typename InnerOptimizer, typename PenaltyFunction, typename = typename
         std::enable_if<std::is_same<PenaltyFunction, typename InnerOptimizer::PenaltyFunction>::value>::type >
SEXP PenseMMPenaltyImpl(SEXP x, SEXP y, SEXP penalties, SEXP enpy_inds,
                        const Rcpp::List& pense_opts, SEXP enpy_opts, const Rcpp::List& optional_args,
                        const Rcpp::List& mm_options, const Rcpp::List& en_options, int) {
  using MMSOptimizer = nsoptim::MMOptimizer<SLoss, PenaltyFunction, InnerOptimizer,
                                            typename InnerOptimizer::Coefficients>;
  return PenseRegressionImpl(MakeOptimizer<MMSOptimizer>(mm_options, en_options),
                             x, y, penalties, enpy_inds, pense_opts, enpy_opts, optional_args);
}

//! Compute the PENSE Regularization Path for the provided penalty function using the MM algorithm
//! This dispatcher inspects `pense_opts["mm_options"]` to determine what EN algorithm to use as the inner optimizer.
//! If the penalty function is not supported by the inner optimizer (e.g., the Ridge optimizer can only handle the
//! Ridge penalty and thus no penalty loadings), returns `R_NilValue`!
//!
//! See `PenseEnRegression` for parameter documentation.
//! @return the regularization path.
template<typename PenaltyFunction>
SEXP PenseMMDispatch(SEXP x, SEXP y, SEXP penalties, SEXP enpy_inds, const Rcpp::List& pense_opts,
                     SEXP enpy_opts, const Rcpp::List& optional_args) {
  using SurrogateLoss = typename SLoss::ConvexSurrogateType;
  const auto mm_options = GetFallback(pense_opts, "algo_opts", Rcpp::List());
  const auto en_options = GetFallback(mm_options, "en_options", Rcpp::List());
  const bool use_sparse_coefs = GetFallback(pense_opts, "sparse", pense::kDefaultUseSparse);
  const auto algorithm = GetFallback(en_options, "algorithm", pense::kDefaultEnAlgorithm);
  switch (algorithm) {
    case pense::EnAlgorithm::kDal: {
      using Optimizer = nsoptim::DalEnOptimizer<SurrogateLoss, PenaltyFunction>;
      return PenseMMPenaltyImpl<Optimizer, PenaltyFunction>(x, y, penalties, enpy_inds, pense_opts,
                                                            enpy_opts, optional_args, mm_options, en_options, 1);
    }
    case pense::EnAlgorithm::kRidge: {
      using Optimizer = AugmentedRidgeOptimizer<SurrogateLoss>;
      return PenseMMPenaltyImpl<Optimizer, nsoptim::RidgePenalty>(x, y, penalties, enpy_inds, pense_opts,
                                                                  enpy_opts, optional_args, mm_options, en_options, 1);
    }
    case pense::EnAlgorithm::kLinearizedAdmm:
      if (use_sparse_coefs) {
        using Optimizer = nsoptim::LinearizedAdmmOptimizer<SurrogateLoss, PenaltyFunction, SparseCoefs>;
        return PenseMMPenaltyImpl<Optimizer, PenaltyFunction>(x, y, penalties, enpy_inds, pense_opts,
                                                              enpy_opts, optional_args, mm_options, en_options, 1);
      } else {
        using Optimizer = nsoptim::LinearizedAdmmOptimizer<SurrogateLoss, PenaltyFunction, DenseCoefs>;
        return PenseMMPenaltyImpl<Optimizer, PenaltyFunction>(x, y, penalties, enpy_inds, pense_opts,
                                                              enpy_opts, optional_args, mm_options, en_options, 1);
      }
    case pense::EnAlgorithm::kCoordinateDescent:
      if (use_sparse_coefs) {
        using Optimizer = nsoptim::CoordinateDescentOptimizer<SurrogateLoss, PenaltyFunction, SparseCoefs>;
        return PenseMMPenaltyImpl<Optimizer, PenaltyFunction>(x, y, penalties, enpy_inds, pense_opts,
                                                              enpy_opts, optional_args, mm_options, en_options, 1);
      } else {
        using Optimizer = nsoptim::CoordinateDescentOptimizer<SurrogateLoss, PenaltyFunction, DenseCoefs>;
        return PenseMMPenaltyImpl<Optimizer, PenaltyFunction>(x, y, penalties, enpy_inds, pense_opts,
                                                              enpy_opts, optional_args, mm_options, en_options, 1);
      }
    case pense::EnAlgorithm::kLars:
    default:
      if (use_sparse_coefs) {
        using Optimizer = nsoptim::AugmentedLarsOptimizer<SurrogateLoss, PenaltyFunction, SparseCoefs>;
        return PenseMMPenaltyImpl<Optimizer, PenaltyFunction>(x, y, penalties, enpy_inds, pense_opts,
                                                              enpy_opts, optional_args, mm_options, en_options, 1);
      } else {
        using Optimizer = nsoptim::AugmentedLarsOptimizer<SurrogateLoss, PenaltyFunction, DenseCoefs>;
        return PenseMMPenaltyImpl<Optimizer, PenaltyFunction>(x, y, penalties, enpy_inds, pense_opts,
                                                              enpy_opts, optional_args, mm_options, en_options, 1);
      }
  }
}

//! Compute the PENSE Regularization Path for the provided penalty function using the CD algorithm.
//! Unless the CD optimizer can handle the given PenaltyFunction class, returns `R_NilValue`.
//!
//! See `PenseEnRegression` for parameter documentation.
//! @return `R_NilValue`.
template<typename PenaltyFunction, typename Coefficients>
SEXP PenseCDPenaltyImpl(SEXP, SEXP, SEXP, SEXP, const Rcpp::List&, SEXP, const Rcpp::List&,
                        const Rcpp::List&, double) {
  return R_NilValue;
}

//! Compute the PENSE Regularization Path for the provided penalty function using the CD algorithm.
//! This function is only enabled if the CD optimizer can handle the requested PenaltyFunction class.
//!
//! See `PenseEnRegression` for parameter documentation.
//! @return the regularization path.
template<typename PenaltyFunction, typename Coefficients, typename = typename
         std::enable_if<!std::is_same<PenaltyFunction, nsoptim::RidgePenalty>::value>::type >
SEXP PenseCDPenaltyImpl(SEXP x, SEXP y, SEXP penalties, SEXP enpy_inds,
                        const Rcpp::List& pense_opts, SEXP enpy_opts, const Rcpp::List& optional_args,
                        const Rcpp::List& cd_options, int) {
  using Optimizer = pense::CDPense<PenaltyFunction, Coefficients>;
  return PenseRegressionImpl(MakeOptimizer<Optimizer>(cd_options),
                             x, y, penalties, enpy_inds, pense_opts, enpy_opts, optional_args);
}


//! Compute the PENSE Regularization Path for the provided penalty function using the CD algorithm direcly on
//! the S-loss.
//! If the penalty function is not supported by CD (e.g., the Ridge penalty), returns `R_NilValue`!
//!
//! See `PenseEnRegression` for parameter documentation.
//! @return the regularization path.
template<typename PenaltyFunction>
SEXP PenseCDDispatch(SEXP x, SEXP y, SEXP penalties, SEXP enpy_inds, const Rcpp::List& pense_opts,
                     SEXP enpy_opts, const Rcpp::List& optional_args) {
  const auto cd_options = GetFallback(pense_opts, "algo_opts", Rcpp::List());
  const bool use_sparse_coefs = GetFallback(pense_opts, "sparse", pense::kDefaultUseSparse);

  if (use_sparse_coefs) {
    return PenseCDPenaltyImpl<PenaltyFunction, SparseCoefs>(
      x, y, penalties, enpy_inds, pense_opts, enpy_opts, optional_args, cd_options, 1);
  } else {
    return PenseCDPenaltyImpl<PenaltyFunction, DenseCoefs>(
      x, y, penalties, enpy_inds, pense_opts, enpy_opts, optional_args, cd_options, 1);
  }
}

//! Compute the PENSE Regularization Path for the provided penalty function.
//! This dispatcher inspects `pense_opts` to determine the PENSE algorithm to use. Can be one of the following:
//!
//!  * ADMM ... use ADMM directly for the S-loss and the penalty function.
//!  * MM   ... use an MM algorithm on the convex surrogate of the S-loss.
//!  * CD   ... use coordinate descent directly on the S-loss.
//!
//! See `PenseEnRegression` for parameter documentation.
//! @return the regularization path.
template<typename PenaltyFunction>
SEXP PensePenaltyDispatch(SEXP x, SEXP y, SEXP penalties, SEXP enpy_inds,
                          SEXP r_pense_opts, SEXP enpy_opts, const Rcpp::List& optional_args) {
  const auto pense_options = as<Rcpp::List>(r_pense_opts);
  switch (GetFallback(pense_options, "algorithm", pense::kDefaultPenseAlgorithm)) {
    case pense::PenseAlgorithm::kCoordinateDescent:
      return PenseCDDispatch<PenaltyFunction>(x, y, penalties, enpy_inds, pense_options,
                                              enpy_opts, optional_args);
    case pense::PenseAlgorithm::kAdmm:
     // Currently not implemented! Fall through to default MM algorithm.
    case pense::PenseAlgorithm::kMm:
    default:
      return PenseMMDispatch<PenaltyFunction>(x, y, penalties, enpy_inds, pense_options,
                                              enpy_opts, optional_args);
  }
}

//! Compute the maximum lambda without penalty loadings.
SEXP PenseMaxGradient(const arma::mat& x, const arma::vec& weights) {
  double max_gradient = 0.;
  for (arma::uword j = 0; j < x.n_cols; ++j) {
    const double gradient_j = std::abs(arma::mean(x.col(j) % weights));
    if (gradient_j > max_gradient) {
      max_gradient = gradient_j;
    }
  }
  return Rcpp::wrap(max_gradient);
}

//! Compute the maximum lambda with penalty loadings.
SEXP PenseMaxGradient(const arma::mat& x, const arma::vec& weights, std::unique_ptr<const arma::vec> loadings) {
  double max_gradient = 0.;
  for (arma::uword j = 0; j < x.n_cols; ++j) {
    const double gradient_j = std::abs(arma::mean(x.col(j) % weights)) / loadings->at(j);
    if (gradient_j > max_gradient) {
      max_gradient = gradient_j;
    }
  }
  return Rcpp::wrap(max_gradient);
}

}  // namespace

namespace pense {
namespace r_interface {
//! Compute the (Adaptive) PENSE Regularization Path.
//!
//! @param x numeric predictor matrix with `n` rows and `p` columns.
//! @param y numeric response vector with `n` elements.
//! @param penalties a list of EN penalties with decreasing values of the lambda hyper-parameter.
//! @param enpy_inds a vector of 1-based indices for the `penalties` list, at which initial ENPY estimates should be
//!                  computed.
//! @param pense_opts a list of options for the PENSE algorithm.
//! @param enpy_opts a list of options for the ENPY algorithm.
//! @param optional_args a list containing the following named items:
//!                       `shared_starts` ... optional list of coefficients to start at every penalty.
//!                       `individual_starts` ... optional list the same length as `penalties` with a list coefficients
//!                                               to start at the corresponding penalty.
//!                       `pen_loadings` ... optional vector of length `p` with non-negative penalty loadings.
SEXP PenseEnRegression(SEXP x, SEXP y, SEXP penalties, SEXP enpy_inds,
                       SEXP pense_opts, SEXP enpy_opts, SEXP r_optional_args) noexcept {
  BEGIN_RCPP
  const auto optional_args = as<Rcpp::List>(r_optional_args);

  if (optional_args.containsElementNamed("pen_loadings")) {
    return PensePenaltyDispatch<nsoptim::AdaptiveEnPenalty>(x, y, penalties, enpy_inds, pense_opts,
                                                            enpy_opts, optional_args);
  }

  return PensePenaltyDispatch<nsoptim::EnPenalty>(x, y, penalties, enpy_inds, pense_opts, enpy_opts, optional_args);

  END_RCPP
}

//! Get the smallest lambda such that the (Adaptive) PENSE estimate gives the empty model.
//!
//! @param x numeric predictor matrix with `n` rows and `p` columns.
//! @param y numeric response vector with `n` elements.
//! @param pense_opts a list of options for the PENSE algorithm.
//! @param optional_args a list containing the following named items:
//!                       `pen_loadings` ... optional vector of length `p` with non-negative penalty loadings.
SEXP PenseMaxLambda(SEXP r_x, SEXP r_y, SEXP r_pense_opts, SEXP r_optional_args) noexcept {
  BEGIN_RCPP
  auto data = MakePredictorResponseData(r_x, r_y);
  const auto pense_opts = as<Rcpp::List>(r_pense_opts);
  const auto optional_args = as<Rcpp::List>(r_optional_args);
  pense::Mscale<pense::RhoBisquare> mscale(as<Rcpp::List>(pense_opts["mscale"]));
  const auto locscale = MLocationScale(data->cy(), mscale, mscale.rho());
  const arma::vec residuals = data->cy() - locscale.location;
  arma::vec weights = residuals % mscale.rho().Weight(residuals, locscale.scale);
  const double denom = arma::mean(weights % residuals);
  weights *= locscale.scale * locscale.scale / denom;

  if (optional_args.containsElementNamed("pen_loadings")) {
    return PenseMaxGradient(data->cx(), weights, MakeVectorView(optional_args["pen_loadings"]));
  } else {
    return PenseMaxGradient(data->cx(), weights);
  }
  END_RCPP
}

}  // namespace r_interface
}  // namespace pense
