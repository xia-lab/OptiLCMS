//
//  enpy_psc.hpp
//  pense
//
//  Created by David Kepplinger on 2019-01-30.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifndef ENPY_PSC_HPP_
#define ENPY_PSC_HPP_

#include <algorithm>
#include <exception>
#include <functional>
#include <string>
#include <utility>
#include <type_traits>
#include "nsoptim.hpp"

#include "alias.hpp"
#include "omp_utils.hpp"
#include "container_utility.hpp"

namespace pense {
//! PSC status code.
enum class PscStatusCode {
  kOk,
  kWarning,
  kError
};

namespace enpy_psc_internal {

using DirectRidgeOptimizer = nsoptim::AugmentedLarsOptimizer<nsoptim::LsRegressionLoss, nsoptim::RidgePenalty,
                                                             nsoptim::RegressionCoefficients<arma::vec>>;
//! Type trait to determine if PSCs should be computed "manually" using the linear algebra identities for
//! the Ridge penalty.
template<typename Optimizer>
using EnableDirectRidge = std::is_same<Optimizer, DirectRidgeOptimizer>;

//! Determine the worst of the two given status codes.
inline PscStatusCode WorstStatusCode(const PscStatusCode a, PscStatusCode b) noexcept {
  if (a == PscStatusCode::kError || b == PscStatusCode::kError) {
    return PscStatusCode::kError;
  } else if (a == PscStatusCode::kWarning || b == PscStatusCode::kWarning) {
    return PscStatusCode::kWarning;
  } else {
    return PscStatusCode::kOk;
  }
}

//! Determine the worst of the two given status codes.
inline PscStatusCode WorstStatusCode(const PscStatusCode a, nsoptim::OptimumStatus b) noexcept {
  if (a == PscStatusCode::kError || b == nsoptim::OptimumStatus::kError) {
    return PscStatusCode::kError;
  } else if (a == PscStatusCode::kWarning || b == nsoptim::OptimumStatus::kWarning) {
    return PscStatusCode::kWarning;
  } else {
    return PscStatusCode::kOk;
  }
}

//! Status from several LOO fits.
struct LooStatus {
  explicit LooStatus(const PscStatusCode _status) noexcept : status(_status), warnings(0) {}

  alias::FwdList<nsoptim::Metrics> metrics;
  PscStatusCode status;
  int warnings;
};

//! PSC result class.
class PscResult {
 public:
  PscResult() noexcept : metrics("psc"), status(PscStatusCode::kOk), warnings(0) {}

  explicit PscResult(LooStatus&& loo_status) noexcept : metrics("psc"), status(loo_status.status),
                                                        warnings(loo_status.warnings) {
    for (auto&& loo_metric : loo_status.metrics) {
      metrics.AddSubMetrics(std::move(loo_metric));
    }
    loo_status.metrics.clear();
  }

  void SetLooStatus(LooStatus&& loo_status) {
    status = WorstStatusCode(status, loo_status.status);
    warnings += loo_status.warnings;
    for (auto&& loo_metric : loo_status.metrics) {
      metrics.AddSubMetrics(std::move(loo_metric));
    }
    loo_status.metrics.clear();
  }

  nsoptim::Metrics metrics;
  PscStatusCode status;
  int warnings;
  std::string message;
  arma::mat pscs;
};
}  // namespace enpy_psc_internal

//! Container for PSC results, holding the matrix with principal sensitivity components and the LS-EN estimate
//! computed on the full data.
template<typename T>
class PscResult : public enpy_psc_internal::PscResult {
 public:
  explicit PscResult(const typename T::Optimum& optim) noexcept : optimum(optim) {}
  typename T::Optimum optimum;
};

namespace enpy_psc_internal {
//! Finalize the PSCs, i.e., compute the eigen decomposition of the sensitivity matrix and add the status messages.
void FinalizePSC(const arma::mat& sensitivity_matrix, enpy_psc_internal::PscResult* psc_result);

//! Concatenate two lists containing LooStatus objects, one for each penalty.
//! The metrics and status from *single* are concatenated to the metrics and the status from *combined*.
//!
//! @param single List to concatenate to the other list.
//! @param combined Output of the concatenated lists.
void ConcatenateLooStatus(alias::FwdList<LooStatus>* single, alias::FwdList<LooStatus>* combined) noexcept;

//! Compute the LOO residuals for rows ``[loo_start_index; loo_end_index)``.
//!
//! @param loss Loss object to compute the LOO residuals for.
//! @param penalties List of penalties for which the LOO residuals should be computed at once.
//! @param loo_start_index Lower bound for row indices to leave out. The lower bound is inclusive.
//! @param loo_end_index Upper bound for the row indices to leave out. The upper bound is exclusive.
//! @param optimizer In/Out. Optimizer to use to compute the leave-one-out residuals.
//! @param sensitivity_matrices Out. A list, the same length as *penalties*, with matrices from which columns the
//!                             LOO residuals are subtracted.
//! @return a list of LooStatus objects, one for each penalty.
template<typename T>
alias::FwdList<LooStatus> ComputeLoo(const nsoptim::LsRegressionLoss& loss,
                                     const alias::FwdList<typename T::PenaltyFunction>& penalties,
                                     arma::uword loo_start_index, const arma::uword loo_end_index, T* optimizer,
                                     alias::FwdList<arma::mat>* sensitivity_matrices) {
  const nsoptim::PredictorResponseData& data = loss.data();
  // Create the LOO data set by removing the observation at `loo_start_index`.
  auto loo_data = std::make_shared<nsoptim::PredictorResponseData>(data.RemoveObservation(loo_start_index));

  alias::FwdList<LooStatus> loo_statuses;
  bool fill_loo_statuses = true;

  // This `loo_loss` holds a shared pointer to the LOO data, therefore, the data used by the loss changes
  // in accordance to `loo_data`!
  // It is nevertheless important to change the loss for the optimizer whenever the data changes to inform the
  // optimizer about the changes!
  nsoptim::LsRegressionLoss loo_loss(loo_data, loss.IncludeIntercept());

  while (loo_start_index < loo_end_index) {
    // Set the loss to the loss with the LOO data.
    optimizer->loss(loo_loss);

    // Compute the LOO optima for all the penalties.
    auto sens_mat_it = sensitivity_matrices->begin();
    auto loo_status_it = fill_loo_statuses ? loo_statuses.before_begin() : loo_statuses.before_begin();
    for (auto&& penalty : penalties) {
      if (fill_loo_statuses) {
        loo_status_it = loo_statuses.emplace_after(loo_status_it, PscStatusCode::kOk);
      } else {
        ++loo_status_it;
      }

      // Only compute LOO residuals if the sensitivity matrix is initialized (i.e., the LS-EN estimate on the full
      // data was computed).
      if (sens_mat_it->n_elem > 0) {
        optimizer->penalty(penalty);
        auto loo_optimum = optimizer->Optimize();

        // This write does not need any protection because this thread is guarantueed to be the only one writing
        // to this column!
        sens_mat_it->col(loo_start_index) -= data.cx() * loo_optimum.coefs.beta + loo_optimum.coefs.intercept;

        loo_status_it->metrics.emplace_front("loo_fit");
        auto&& loo_fit_metric = loo_status_it->metrics.front();
        loo_fit_metric.AddMetric("loo_index", static_cast<int>(loo_start_index));

        if (loo_optimum.metrics) {
          loo_fit_metric.AddSubMetrics(std::move(*loo_optimum.metrics));
          loo_optimum.metrics.reset();
        }

        if (loo_optimum.status != nsoptim::OptimumStatus::kOk) {
          loo_fit_metric.AddMetric("lsen_status", static_cast<int>(loo_optimum.status));
          loo_fit_metric.AddMetric("lsen_message", loo_optimum.message);

          loo_status_it->status = WorstStatusCode(loo_status_it->status, loo_optimum.status);
        }
      }
      ++sens_mat_it;
    }

    // "Hide" next row if there are any rows left.
    if (loo_start_index < loo_end_index - 1) {
      loo_data->x().row(loo_start_index) = data.cx().row(loo_start_index);
      loo_data->y()[loo_start_index] = data.cy()[loo_start_index];
    }
    ++loo_start_index;
    fill_loo_statuses = false;
  }

  return loo_statuses;
}

//! Compute the Principal Sensitivity Components if OpenMP is enabled and needed by computing the LOO
//! residuals one by one.
//!
//! @param loss the LS-loss object to compute the PSCs for.
//! @param penalties a lisf of penalties for which the PSCs should be computed at once.
//! @param optimizer Optimizer to use to compute leave-one-out residuals.
//! @param num_threads number of threads to use.
//! @return A list of PSC structures, one for each given penalty, in the same order as `penalties`.
template<typename Optimizer, typename = typename std::enable_if<!EnableDirectRidge<Optimizer>::value>::type >
alias::FwdList<pense::PscResult<Optimizer>> ComputePscs(
    const nsoptim::LsRegressionLoss& loss, const alias::FwdList<typename Optimizer::PenaltyFunction>& penalties,
    Optimizer optimizer, int num_threads) {
  // using PenaltyFunction = typename Optimizer::PenaltyFunction;
  using arma::uword;
  using alias::FwdList;
  using LooStatusList = FwdList<LooStatus>;

  // Note: the initializer for the following reduction uses the default constructor!
  #pragma omp declare reduction(c:FwdList<LooStatus>:ConcatenateLooStatus(&omp_in, &omp_out))

  const nsoptim::PredictorResponseData& data = loss.data();
  // A list of PscResult objects and the sensitivity matrices, i.e., matrices `R` in the paper, one tuple for each
  // penalty.
  pense::utility::OrderedList<double, pense::PscResult<Optimizer>, std::greater<double>> psc_results;
  pense::utility::OrderedList<double, arma::mat, std::greater<double>> sensitivity_matrices;

  // First optimize with respect to the full data set and compute the predictions.
  // This should not be parallelized because it is more efficient to use the same optimizer for all penalties!
  optimizer.loss(loss);
  for (auto&& penalty : penalties) {
    optimizer.penalty(penalty);
    auto psc_result_it = psc_results.emplace(penalty.lambda(), optimizer.Optimize());
    // Retain metrics from all LOO fits.
    auto&& full_fit_metrics = psc_result_it->metrics.CreateSubMetrics("full_fit");
    if (psc_result_it->optimum.metrics) {
      full_fit_metrics.AddSubMetrics(std::move(*(psc_result_it->optimum.metrics)));
      psc_result_it->optimum.metrics.reset();
    }

    switch (psc_result_it->optimum.status) {
      case nsoptim::OptimumStatus::kError:
        psc_result_it->status = PscStatusCode::kError;
        psc_result_it->message = std::string("Can not compute LS-EN residuals: ") + psc_result_it->optimum.message;
        sensitivity_matrices.emplace(penalty.lambda(), arma::mat());
        break;
      case nsoptim::OptimumStatus::kWarning:
        ++(psc_result_it->warnings);
        psc_result_it->status = PscStatusCode::kWarning;
        psc_result_it->message = std::string("LS-EN residuals are not reliable: ") + psc_result_it->optimum.message +
                                 "; ";
        // Fall through the case of no error or warning!
      default:
        // If status != kError, fill the sensitivity matrix with the LS-EN residuals.
        sensitivity_matrices.emplace(penalty.lambda(), arma::repmat(
          data.cx() * psc_result_it->optimum.coefs.beta + psc_result_it->optimum.coefs.intercept, 1, data.n_obs()));
        break;
    }
  }

  const uword block_size = data.n_obs() / num_threads + static_cast<uword>(data.n_obs() % num_threads > 0);
  LooStatusList loo_statuses;
  #pragma omp parallel num_threads(num_threads) default(none) firstprivate(block_size) \
    shared(data, loss, penalties, loo_statuses, sensitivity_matrices, psc_results, optimizer)
  {
    #pragma omp for reduction(c:loo_statuses)
    for (uword start_index = 0; start_index < data.n_obs(); start_index += block_size) {
      const uword upper_index = std::min(start_index + block_size, data.n_obs());
      Optimizer thread_private_optimizer(optimizer);
      loo_statuses = ComputeLoo(loss, penalties, start_index, upper_index, &thread_private_optimizer,
                                &sensitivity_matrices.items());
    }

    #pragma omp single nowait
    {
      auto loo_status_it = loo_statuses.begin();
      auto sens_mat_it = sensitivity_matrices.begin();
      for (auto psc_result_it = psc_results.begin(), end = psc_results.end(); psc_result_it != end;
          ++psc_result_it, ++sens_mat_it, ++loo_status_it) {
        // Skip results with errors.
        if (psc_result_it->status == PscStatusCode::kError || loo_status_it->status == PscStatusCode::kError) {
          psc_result_it->SetLooStatus(std::move(*loo_status_it));
          continue;
        }
        #pragma omp task firstprivate(sens_mat_it, psc_result_it, loo_status_it) default(none)
        {
          psc_result_it->SetLooStatus(std::move(*loo_status_it));
          enpy_psc_internal::FinalizePSC(*sens_mat_it, &(*psc_result_it));
        }
      }
    }
  }

  return psc_results.items();
}

//! Compute the Principal Sensitivity Components if OpenMP is disabled by computing the LOO
//! residuals one by one.
//!
//! @param loss the LS-loss object to compute the PSCs for.
//! @param penalties a lisf of penalties for which the PSCs should be computed at once.
//! @param optimizer Optimizer to use to compute leave-one-out residuals.
//! @return A list of PSC structures, one for each given penalty, in the same order as `penalties`.
template<typename Optimizer, typename = typename std::enable_if<!EnableDirectRidge<Optimizer>::value>::type>
alias::FwdList<pense::PscResult<Optimizer>> ComputePscs(
    const nsoptim::LsRegressionLoss& loss, const alias::FwdList<typename Optimizer::PenaltyFunction>& penalties,
    Optimizer optimizer) {
  using arma::uword;
  using enpy_psc_internal::ComputeLoo;
  using LooStatusList = alias::FwdList<enpy_psc_internal::LooStatus>;

  const nsoptim::PredictorResponseData& data = loss.data();
  // A list of PscResult objects and the sensitivity matrices, i.e., matrices `R` in the paper, one per penalty.
  alias::FwdList<pense::PscResult<Optimizer>> psc_results;
  alias::FwdList<arma::mat> sensitivity_matrices;

  // auto psc_state_it = psc_states.before_begin();
  auto psc_result_it = psc_results.before_begin();
  auto sens_mat_it = sensitivity_matrices.before_begin();

  // First optimize with respect to the full data set and compute the predictions.
  optimizer.loss(loss);
  for (auto&& penalty : penalties) {
    optimizer.penalty(penalty);
    psc_result_it = psc_results.emplace_after(psc_result_it, optimizer.Optimize());

    // Retain metrics from all LOO fits.
    auto&& full_fit_metrics = psc_result_it->metrics.CreateSubMetrics("full_fit");
    if (psc_result_it->optimum.metrics) {
      full_fit_metrics.AddSubMetrics(std::move(*(psc_result_it->optimum.metrics)));
      psc_result_it->optimum.metrics.reset();
    }

    switch (psc_result_it->optimum.status) {
      case nsoptim::OptimumStatus::kError:
        psc_result_it->status = PscStatusCode::kError;
        psc_result_it->message = std::string("Can not compute LS-EN residuals: ") + psc_result_it->optimum.message;
        sens_mat_it = sensitivity_matrices.emplace_after(sens_mat_it);
        break;
      case nsoptim::OptimumStatus::kWarning:
        ++(psc_result_it->warnings);
        psc_result_it->status = PscStatusCode::kWarning;
        psc_result_it->message = std::string("LS-EN residuals are not reliable: ") + psc_result_it->optimum.message +
                                 "; ";
        // Fall through the case of no error or warning!
      default:
        // If status != kError, fill the sensitivity matrix with the LS-EN residuals.
        sens_mat_it = sensitivity_matrices.emplace_after(sens_mat_it, arma::repmat(
          data.cx() * psc_result_it->optimum.coefs.beta + psc_result_it->optimum.coefs.intercept, 1, data.n_obs()));
        break;
    }
  }

  LooStatusList loo_statuses = ComputeLoo(loss, penalties, 0, data.n_obs(), &optimizer, &sensitivity_matrices);
  auto loo_status_it = loo_statuses.begin();
  sens_mat_it = sensitivity_matrices.begin();
  for (auto psc_result_it = psc_results.begin(), end = psc_results.end(); psc_result_it != end;
      ++psc_result_it, ++sens_mat_it, ++loo_status_it) {
    psc_result_it->SetLooStatus(std::move(*loo_status_it));

    // Skip results with errors.
    if (psc_result_it->status == PscStatusCode::kError) {
      continue;
    }
    enpy_psc_internal::FinalizePSC(*sens_mat_it, &(*psc_result_it));
  }
  return psc_results;
}

//! Compute the Principal Sensitivity Components using identities for the Ridge penalty
//! if OpenMP is disabled.
//!
//! @param loss the LS regression loss object to compute the PSCs for.
//! @param penalties a lisf of Ridge penalties for which the PSCs should be computed at once.
//! @param optimizer Optimizer to use to compute leave-one-out residuals.
//! @return A list of PSC structures, one for each given penalty, in the same order as `penalties`.
alias::FwdList<pense::PscResult<DirectRidgeOptimizer>> ComputeRidgePscs(const nsoptim::LsRegressionLoss& loss,
    const alias::FwdList<nsoptim::RidgePenalty>& penalties, DirectRidgeOptimizer optimizer);

//! Compute the Principal Sensitivity Components using identities for the Ridge penalty
//! if OpenMP is enabled and necessary.
//!
//! @param loss the LS regression loss object to compute the PSCs for.
//! @param penalties a lisf of Ridge penalties for which the PSCs should be computed at once.
//! @param optimizer Optimizer to use to compute leave-one-out residuals.
//! @param num_threads number of OpenMP threads.
//! @return A list of PSC structures, one for each given penalty, in the same order as `penalties`.
alias::FwdList<pense::PscResult<DirectRidgeOptimizer>> ComputeRidgePscs(const nsoptim::LsRegressionLoss& loss,
    const alias::FwdList<nsoptim::RidgePenalty>& penalties, const DirectRidgeOptimizer& optimizer, int num_threads);

//! Compute the Principal Sensitivity Components using identities for the Ridge penalty
//! if OpenMP is disabled.
//!
//! @param loss the LS regression loss object to compute the PSCs for.
//! @param penalties a lisf of Ridge penalties for which the PSCs should be computed at once.
//! @param optimizer Optimizer to use to compute leave-one-out residuals.
//! @return A list of PSC structures, one for each given penalty, in the same order as `penalties`.
template<typename Optimizer, typename = typename std::enable_if<EnableDirectRidge<Optimizer>::value>::type>
alias::FwdList<pense::PscResult<Optimizer>> ComputePscs(const nsoptim::LsRegressionLoss& loss,
    const alias::FwdList<nsoptim::RidgePenalty>& penalties, const Optimizer& optimizer) {
  return ComputeRidgePscs(loss, penalties, optimizer);
}

//! Compute the Principal Sensitivity Components using identities for the Ridge penalty
//! if OpenMP is enabled and needed.
//!
//! @param loss the LS regression loss object to compute the PSCs for.
//! @param penalties a lisf of Ridge penalties for which the PSCs should be computed at once.
//! @param optimizer Optimizer to use to compute leave-one-out residuals.
//! @param num_threads number of OpenMP threads.
//! @return A list of PSC structures, one for each given penalty, in the same order as `penalties`.
template<typename Optimizer, typename = typename std::enable_if<EnableDirectRidge<Optimizer>::value>::type>
alias::FwdList<pense::PscResult<Optimizer>> ComputePscs(const nsoptim::LsRegressionLoss& loss,
    const alias::FwdList<nsoptim::RidgePenalty>& penalties, const Optimizer& optimizer, int num_threads) {
  return ComputeRidgePscs(loss, penalties, optimizer, num_threads);
}

}  // namespace enpy_psc_internal

//! Compute the Principal Sensitivity Components as defined in
//! Pena, D. & Yohai, V. (1999). A Fast Procedure for Outlier Diagnostics in Large Regression
//! Problems. JASA.
//!
//! @param loss the S-loss object to compute the PSCs for.
//! @param penalties a lisf of penalties for which the PSCs should be computed at once.
//! @param optimizer Optimizer to use to compute leave-one-out residuals.
//! @param num_threads number of threads to use.
//! @return A list of PSC structures, one for each given penalty, in the same order as `penalties`.
template<typename Optimizer>
alias::FwdList<PscResult<Optimizer>> PrincipalSensitiviyComponents(
    const nsoptim::LsRegressionLoss& loss, const alias::FwdList<typename Optimizer::PenaltyFunction>& penalties,
    const Optimizer& optimizer, const int num_threads) {
  if (omp::Enabled(num_threads)) {
    return enpy_psc_internal::ComputePscs(loss, penalties, optimizer, num_threads);
  } else {
    return enpy_psc_internal::ComputePscs(loss, penalties, optimizer);
  }
}

//! Compute the Principal Sensitivity Components as defined in
//! Pena, D. & Yohai, V. (1999). A Fast Procedure for Outlier Diagnostics in Large Regression
//! Problems. JASA.
//!
//! @param loss the S-loss object to compute the PSCs for.
//! @param optim the optimizer to use to compute leave-one-out residuals.
//! @return a matrix of PSCs.
template<typename Optimizer>
PscResult<Optimizer> PrincipalSensitiviyComponents(const nsoptim::LsRegressionLoss& loss, const Optimizer& optim,
                                                   const int num_threads) {
  const alias::FwdList<typename Optimizer::PenaltyFunction> penalties { optim.penalty() };

  if (omp::Enabled(num_threads)) {
    return enpy_psc_internal::ComputePscs(loss, penalties, optim, num_threads).front();
  } else {
    return enpy_psc_internal::ComputePscs(loss, penalties, optim).front();
  }
}
}  // namespace pense

#endif  // ENPY_PSC_HPP_
