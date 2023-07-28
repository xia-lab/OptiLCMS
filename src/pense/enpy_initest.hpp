//
//  enpy_initest.hpp
//  pense
//
//  Created by David Kepplinger on 2019-01-30.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifndef ENPY_INITEST_HPP_
#define ENPY_INITEST_HPP_

#include <exception>
#include <string>
#include <utility>
#include <type_traits>
#include <algorithm>
#include <unordered_set>
#include <functional>

#include "nsoptim.hpp"
#include "alias.hpp"
#include "robust_scale_location.hpp"
#include "s_loss.hpp"
#include "enpy_psc.hpp"
#include "omp_utils.hpp"
#include "container_utility.hpp"
#include "enpy_types.hpp"

namespace pense {
//! Compute the Pena-Yohai Initial Estimators for the S-loss.
//!
//! @param loss the S loss object for which to obtain initial estimates.
//! @param penalties a list of penalties to compute the initial estimators for.
//! @param optim the optimizer to compute the estimates.
//! @param config list with configuration options.
//! @return a list of ENPY results, one for each given penalty, in the same order as `penalties`.
template<typename Optimizer>
alias::FwdList<PyResult<Optimizer>> PenaYohaiInitialEstimators(
    const SLoss& loss, const alias::FwdList<typename Optimizer::PenaltyFunction>& penalties,
    const Optimizer& optim, const Rcpp::List& config);



namespace enpy_initest_internal {
constexpr arma::uword kMinObs = 3;  //!< Mininum number of observations in the PSC-filtered data.

template<class Optimizer>
class CandidateComparator {
  using Optimum = typename Optimizer::Optimum;
 public:
  inline bool operator()(const Optimum& a, const Optimum& b) const noexcept {
    return a.objf_value < b.objf_value;
  }
};

//! Get a list of subsets from the given PSCs
//!
//! @param pscs a matrix of principal sensitivity components.
//! @param subset_size the desired size of the subsets.
//! @return a list of unique subsets of length `subset_size`.
alias::FwdList<arma::uvec> GetSubsetList(const arma::mat& pscs, const arma::uword subset_size) noexcept;

//! Get a list of subsets from the given PSCs
//!
//! @param pscs a matrix of principal sensitivity components.
//! @param indices the "true" indices of the rows in `pscs`. Must be sorted.
//! @param subset_size the desired size of the subsets.
//! @return a list of unique subsets of length `subset_size`.
alias::FwdList<arma::uvec> GetSubsetList(const arma::mat& pscs, const arma::uvec& indices,
                                         const arma::uword subset_size) noexcept;


//! Merge the metrics and data from *psc_result* with *metrics*.
template<typename T>
void AppendPscMetrics(PscResult<T>&& psc_result, nsoptim::Metrics* metrics) noexcept {
  metrics->AddSubMetrics(std::move(psc_result.metrics));
  metrics->AddMetric("num_pscs", static_cast<int>(psc_result.pscs.n_cols));
  metrics->AddMetric("psc_status", static_cast<int>(psc_result.status));
  metrics->AddDetail("psc_warnings", psc_result.warnings);
  if (psc_result.message.size() > 0) {
    metrics->AddMetric("psc_message", psc_result.message);
  }
}

//! Construct a Metrics object with the given name, using data from *psc_result*.
template<typename T>
inline nsoptim::Metrics CreatePscMetrics(std::string&&name, PscResult<T>&& psc_result) noexcept {
  nsoptim::Metrics metrics("enpy_initest");
  AppendPscMetrics(std::move(psc_result), &(metrics.CreateSubMetrics(std::move(name))));
  return metrics;
}

struct PyConfiguration {
  int max_it;  //!< Maximum number of iterations.
  double eps;  //!< Numerical tolerance level to determine convergence.
  double keep_psc_proportion;  //!< Proportion of observations to keep based on PSCs.
  bool use_residual_threshold;  //!< Use a fixed threshold instead of a proportion to screen observations based on
                                //!< their residuals.
  double keep_residuals_proportion;  //!< Proportion of observations to keep based on the residuals.
  double keep_residuals_threshold;  //!< Fixed threshold to keep observations based on the residuals.
  double retain_best_factor;  //!< Retain not only the candidates from the last iteration, but also those that are
                              //!< within this factor of the best candidate.
  int retain_max;  //!< Retain at most this number of candidates. Candidates are ordered by the objective function.
                   //!< If negative, all candidates are retained.
  int num_threads;  //!< Number of concurrent threads to use.
};

//! Parse an Rcpp::List into the PyConfiguration structure.
PyConfiguration ParseConfiguration(const Rcpp::List& config) noexcept;

//! Compute the hash of a sorted index vector.
//!
//! @param vector a sorted vector of indices.
//! @return the hash of `vector`.
arma::uword HashIndexVector(const arma::uvec& vector) noexcept;

//! Compute the hash of of the sequence 0 ... `to`.
//!
//! This is identical to calling `HashIndexVector(arma::regspace<arma::uvec>(0, to))`;
//!
//! @param to upper end of the sequence
//! @return the hash of the sequence.
arma::uword HashSequence(const arma::uword to) noexcept;

//! Get an index vector of observations to keep based on the given residuals and configuration.
//! The returned vector of indices is always sorted (smallest to largest)!
//!
//! @param residuals vector of residuals.
//! @param mscale_est the M-scale of the given residuals.
//! @param config configuration.
//! @param all_indices a vector of all indices in unspecified order.
//! @return a sorted vector of indices (smallest to largest).
arma::uvec GetResidualKeepIndices(const arma::vec& residuals, const double mscale_est,
                                  const PyConfiguration& config, arma::uvec* all_indices);

//! Enumeration of different results that can occur when computing the minimizer on a subset of the data.
enum class SubsetEstimateResult {
  kOk, kDuplicate, kOptimizerWarning, kOptimizerError
};

//! Compute the Pena-Yohai initial estimator
//!
//! @param loss the S-loss for which to obtain initial estimates.
//! @param penalty the penalty to compute the initial estimators for.
//! @param psc_result the PSC result on the full data.
//! @param optim the optimizer to compute the estimates with.
//! @param pyconfig configuration object.
//! @param num_threads the number of threads to use. Overrides the configuration from *pyconfig*.
template<typename Optimizer>
PyResult<Optimizer> PYIterations(SLoss loss, const typename Optimizer::PenaltyFunction& penalty,
                                 PscResult<Optimizer>&& full_psc_result, Optimizer optim,
                                 const PyConfiguration& pyconfig, const int num_threads);

//! Compute the Pena-Yohai Initial Estimators for the S-loss if OpenMP is enabled and needed.
//!
//! @param loss the S loss object for which to obtain initial estimates.
//! @param penalties a list of penalties to compute the initial estimators for.
//! @param optim the optimizer to compute the estimates.
//! @param pyconfig Configuration options.
//! @return a list of ENPY results, one for each given penalty, in the same order as `penalties`.
template<typename Optimizer>
alias::FwdList<PyResult<Optimizer>> ComputeENPY(const SLoss& loss,
                                                const alias::FwdList<typename Optimizer::PenaltyFunction>& penalties,
                                                const Optimizer& optim, const PyConfiguration& pyconfig,
                                                int num_threads) {
  using arma::uword;
  // For each penalty, compute the optimizer and PSCs on the full data.
  nsoptim::LsRegressionLoss full_ls_loss(loss.SharedData(), loss.IncludeIntercept());
  pense::utility::OrderedList<double, PyResult<Optimizer>, std::greater<double>> py_initest_results;
  auto psc_results = PrincipalSensitiviyComponents(full_ls_loss, penalties, optim, num_threads);

  // The PY iterations are done separately for each penalty in parallel.
  #pragma omp parallel num_threads(num_threads) default(none) \
    shared(py_initest_results, psc_results, penalties, optim, loss, pyconfig)
  {
    #pragma omp single nowait
    {
      auto penalty_it = penalties.begin();
      for (auto psc_res_it = psc_results.begin(), psc_res_end = psc_results.end(); psc_res_it != psc_res_end;
           ++psc_res_it, ++penalty_it) {
        if (psc_res_it->status != PscStatusCode::kError) {
          #pragma omp task default(none) firstprivate(psc_res_it, penalty_it) \
            shared(py_initest_results, pyconfig, loss, optim)
          {
            // The following line copies the *optim* optimzer, so it is okay that *optim* is shared among threads!
            auto pyit_res = PYIterations(loss, *penalty_it, std::move(*psc_res_it), optim, pyconfig, 1);
            #pragma omp critical(emplace_pyit_res)
            py_initest_results.emplace(penalty_it->lambda(), std::move(pyit_res));
          }
        } else {
          py_initest_results.emplace(penalty_it->lambda(), CreatePscMetrics("full_data", std::move(*psc_res_it)));
        }
      }
    }
  }

  return py_initest_results.items();
}


//! Compute the Pena-Yohai Initial Estimators for the S-loss if OpenMP is disabled or not needed.
//!
//! @param loss the S loss object for which to obtain initial estimates.
//! @param penalties a list of penalties to compute the initial estimators for.
//! @param optim the optimizer to compute the estimates.
//! @param pyconfig Configuration options.
//! @return a list of ENPY results, one for each given penalty, in the same order as `penalties`.
template<typename Optimizer>
alias::FwdList<PyResult<Optimizer>> ComputeENPY(const SLoss& loss,
                                                const alias::FwdList<typename Optimizer::PenaltyFunction>& penalties,
                                                const Optimizer& optim, const PyConfiguration& pyconfig) {
  using arma::uword;

  // For each penalty, compute the optimizer and PSCs on the full data.
  nsoptim::LsRegressionLoss full_ls_loss(loss.SharedData(), loss.IncludeIntercept());
  alias::FwdList<PyResult<Optimizer>> py_initest_results;
  auto psc_results = PrincipalSensitiviyComponents(full_ls_loss, penalties, optim, 1);

  // The PY iterations are done separately.
  auto penalty_it = penalties.begin();
  auto py_initest_res_it = py_initest_results.before_begin();
  for (auto&& psc_result : psc_results) {
    if (psc_result.status != PscStatusCode::kError) {
      py_initest_res_it = py_initest_results.insert_after(
        py_initest_res_it, PYIterations(loss, *penalty_it, std::move(psc_result), optim, pyconfig, 1));
    } else {
      py_initest_res_it = py_initest_results.emplace_after(py_initest_res_it,
                                                          CreatePscMetrics("full_data", std::move(psc_result)));
    }
    ++penalty_it;
  }
  return py_initest_results;
}

template<typename Optimizer>
PyResult<Optimizer> PYIterations(SLoss loss, const typename Optimizer::PenaltyFunction& penalty,
                                 PscResult<Optimizer>&& full_psc_result, Optimizer pyinit_optim,
                                 const PyConfiguration& pyconfig, const int num_threads) {
  using arma::uvec;
  using arma::vec;
  using arma::uword;
  using arma::regspace;
  using nsoptim::PredictorResponseData;
  using nsoptim::Metrics;
  using HashSet = std::unordered_set<uword>;
  using SubsetList = alias::FwdList<arma::uvec>;

  const PredictorResponseData& data = loss.data();
  PyResult<Optimizer> py_result(CreatePscMetrics("full_data", std::move(full_psc_result)));

  // Add the optimum on the full data set. This will be the first "best candidate".
  py_result.initial_estimates.emplace_front(full_psc_result.optimum);

  // Add details about the PSCs from the full data set.
  Metrics* iter_metrics = &py_result.metrics.CreateSubMetrics("full_data_iter");

  uvec all_indices = regspace<uvec>(0, data.n_obs() - 1);
  uvec residuals_keep_ind = all_indices;
  HashSet residual_filtered_data_hashes { HashSequence(all_indices.n_elem - 1) };
  auto best_candidate_it = py_result.initial_estimates.begin();
  nsoptim::LsRegressionLoss ls_loss(loss.SharedData(), loss.IncludeIntercept());
  SubsetList psc_subsets = GetSubsetList(full_psc_result.pscs,
                                         std::max<uword>(pyconfig.keep_psc_proportion * data.n_obs(), kMinObs));
  SubsetList* current_psc_subsets = &psc_subsets;

  // The initial "best candidate" comes from the PSC and thus has the wrong objective function value.
  best_candidate_it->objf_value = loss.Evaluate(best_candidate_it->coefs) + penalty.Evaluate(best_candidate_it->coefs);

  // Set the correct penalty for the iterations.
  pyinit_optim.penalty(penalty);

  // Start the PY iterations.
  int iter = 0;
  decltype(best_candidate_it) insert_candidate_it;
  while (true) {
    insert_candidate_it = best_candidate_it;

    for (auto&& subset : *current_psc_subsets) {
      auto&& psc_metric = iter_metrics->CreateSubMetrics("psc_subset");
      psc_metric.AddDetail("n_obs", static_cast<int>(subset.n_elem));
      pyinit_optim.loss(nsoptim::LsRegressionLoss(
        std::make_shared<nsoptim::PredictorResponseData>(loss.data().Observations(subset)), loss.IncludeIntercept()));
      auto subset_optimum = pyinit_optim.Optimize();
      // Remove the reference to the subset loss.
      subset_optimum.loss = ls_loss;
      psc_metric.AddDetail("estimate_status", static_cast<int>(subset_optimum.status));
      if (subset_optimum.status != nsoptim::OptimumStatus::kOk) {
        psc_metric.AddDetail("estimate_status_message", subset_optimum.message);
        if (subset_optimum.status == nsoptim::OptimumStatus::kError) {
          // Don't add subsets which result in an error!
          continue;
        }
      }
      if (subset_optimum.metrics) {
        psc_metric.AddSubMetrics(std::move(*subset_optimum.metrics));
        subset_optimum.metrics.reset();
      }
      insert_candidate_it = py_result.initial_estimates.insert_after(insert_candidate_it, std::move(subset_optimum));
    }

    // Reset the loss to the
    pyinit_optim.loss(ls_loss);

    // Check if we are at the end of our iterations
    if (++iter >= pyconfig.max_it) {
      break;
    }

    // Find the estimator with the lowest objective function
    auto new_best_candidate_it = py_result.initial_estimates.begin();
    // No need to initialize the following two variables. They will be assigned if a new best candidate is found.
    // If there is no new best candidate, the PY iterations stop right after and the variables are never used.
    arma::vec best_candidate_residuals;
    double best_candidate_mscale = 1.;

    // Skip the first element since this is the `best_candidate_it` and already evaluated and iterate just after the
    // last inserted element.
    auto end_check_candidate_it = insert_candidate_it;
    ++end_check_candidate_it;
    for (auto cand_it = ++py_result.initial_estimates.begin(); cand_it != end_check_candidate_it; ++cand_it) {
      // Compute the S-loss of the candidate
      const arma::vec candidate_residuals = data.cy() - data.cx() * cand_it->coefs.beta - cand_it->coefs.intercept;
      const auto candidate_eval = loss.EvaluateResiduals(candidate_residuals);
      cand_it->objf_value = candidate_eval.loss + penalty.Evaluate(cand_it->coefs);

      if (cand_it->objf_value < new_best_candidate_it->objf_value) {
        new_best_candidate_it = cand_it;
        best_candidate_residuals = candidate_residuals;
        best_candidate_mscale = candidate_eval.scale;
      }
    }

    // Check if the best candidate is still the first element, i.e., didn't change.
    if (new_best_candidate_it == best_candidate_it) {
      iter_metrics->AddMetric("stop", "no change in best candidate");
      break;
    }

    // Check if the new best candidate differs from the previous best candidate.
    // The obj. function value is already computed using the original S-loss!
    const double rel_diff = std::abs(best_candidate_it->objf_value - new_best_candidate_it->objf_value);

    iter_metrics->AddDetail("rel_diff", rel_diff);
    if (rel_diff < pyconfig.eps) {
      // The previously best candidate is (almost) identical to the newly best candidate. We are done!
      iter_metrics->AddMetric("stop", "best candidate is almost identical");
      break;
    }

    iter_metrics = &py_result.metrics.CreateSubMetrics("py_iteration");

    // Move the new best candidate to the front of the list.
    if (insert_candidate_it != new_best_candidate_it) {
      std::iter_swap(new_best_candidate_it, py_result.initial_estimates.begin());
      best_candidate_it = py_result.initial_estimates.begin();
    } else {
      // If the last inserted element is the new best candidate, we need to ensure that `insert_candidate_it`
      // still points to the last element for this iteration after the swap!
      insert_candidate_it = py_result.initial_estimates.begin();
      std::iter_swap(new_best_candidate_it, insert_candidate_it);
      best_candidate_it = py_result.initial_estimates.begin();
    }

    // Use the new best solution to thin-out the data based on the residuals and compute new PSCs.
    residuals_keep_ind = GetResidualKeepIndices(best_candidate_residuals, best_candidate_mscale, pyconfig,
                                                &all_indices);

    // Store the hashed index vector to avoid unnecessary iterations.
    const auto hash_insert = residual_filtered_data_hashes.insert(HashIndexVector(residuals_keep_ind));
    if (!hash_insert.second) {
      // The resulting subset of observations was already used before.
      // Therefore, the next iteration would result in duplicate solutions and we can stop.
      iter_metrics->AddMetric("stop", "filtered residuals are duplicates");
      break;
    }

    // This subset of observations was not yet considered. Continue.
    nsoptim::LsRegressionLoss filtered_ls_loss(
      std::make_shared<PredictorResponseData>(data.Observations(residuals_keep_ind)), loss.IncludeIntercept());

    pyinit_optim.loss(filtered_ls_loss);
    iter_metrics->AddDetail("n_obs", static_cast<int>(residuals_keep_ind.n_elem));

    // Compute the optimzer and the PSCs on the residual-filtered data
    const arma::uword new_subsets_size = std::max<uword>(pyconfig.keep_psc_proportion * residuals_keep_ind.n_elem,
                                                         kMinObs);

    PscResult<Optimizer> psc_result = PrincipalSensitiviyComponents(filtered_ls_loss, pyinit_optim, num_threads);
    psc_subsets = GetSubsetList(psc_result.pscs, residuals_keep_ind, new_subsets_size);

    AppendPscMetrics(std::move(psc_result), iter_metrics);
    if (psc_result.status == PscStatusCode::kError) {
      iter_metrics->AddMetric("stop", "Could not compute PSCs");
      break;
    }

    // Check if the LS estimate, computed on the filtered data is better than the previously best candidate.
    psc_result.optimum.objf_value = loss(psc_result.optimum.coefs) + penalty(psc_result.optimum.coefs);
    if (psc_result.optimum.objf_value < best_candidate_it->objf_value) {
      // If so, replace the best candidate by this estimate.
      py_result.initial_estimates.push_front(std::move(psc_result.optimum));
      best_candidate_it = py_result.initial_estimates.begin();
    }

    current_psc_subsets = &psc_subsets;
  }

  // Only retain the estimates that are within a certain factor of the best solution.
  // `insert_candidate_it` points to the last element inserted in the last iteration.
  if (pyconfig.retain_best_factor > 1) {
    const double objf_value_cutoff = pyconfig.retain_best_factor * best_candidate_it->objf_value;
    insert_candidate_it = best_candidate_it;
    while (true) {
      const auto current_end = py_result.initial_estimates.end();
      if (insert_candidate_it == current_end) {
        break;
      }

      auto tmp_it = insert_candidate_it;
      // Search for either the last element or the next element that is "good enough".
      int removing = 0;
      while (++tmp_it != current_end && tmp_it->objf_value > objf_value_cutoff) { ++removing; }
      if (removing > 0) {
        py_result.initial_estimates.erase_after(insert_candidate_it, tmp_it);
      }
      insert_candidate_it = tmp_it;
    }
  } else {
    // Only retain the estimates from the last iteration.
    py_result.initial_estimates.erase_after(insert_candidate_it, py_result.initial_estimates.end());
  }

  if (pyconfig.retain_max > 0) {
    // Retain only a certain number of the best candidates.
    if (pyconfig.retain_max < std::distance(py_result.initial_estimates.begin(), py_result.initial_estimates.end())) {
      // There are more than `retain_max` candidates.
      py_result.initial_estimates.sort(CandidateComparator<Optimizer>());
      auto delete_rest_it = py_result.initial_estimates.before_begin();
      std::advance(delete_rest_it, pyconfig.retain_max);
      py_result.initial_estimates.erase_after(delete_rest_it, py_result.initial_estimates.end());
    }
  }

  return py_result;
}
}  // namespace enpy_initest_internal


//! Compute the Pena-Yohai Initial Estimators for the S-loss.
//!
//! @param loss the S loss object for which to obtain initial estimates.
//! @param penalties a list of penalties to compute the initial estimators for.
//! @param optim the optimizer to compute the estimates.
//! @param config list with configuration options.
//! @return a list of ENPY results, one for each given penalty, in the same order as `penalties`.
template<typename Optimizer>
alias::FwdList<PyResult<Optimizer>> PenaYohaiInitialEstimators(
    const SLoss& loss, const alias::FwdList<typename Optimizer::PenaltyFunction>& penalties,
    const Optimizer& optim, const Rcpp::List& config) {
  const auto pyconfig = enpy_initest_internal::ParseConfiguration(config);
  if (omp::Enabled(pyconfig.num_threads)) {
    return enpy_initest_internal::ComputeENPY(loss, penalties, optim, pyconfig, pyconfig.num_threads);
  } else {
    return enpy_initest_internal::ComputeENPY(loss, penalties, optim, pyconfig);
  }
}

}  // namespace pense

#endif  // ENPY_INITEST_HPP_
