//
//  enpy_psc.cc
//  pense
//
//  Created by David Kepplinger on 2019-01-30.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#include "constants.hpp"
#include "enpy_psc.hpp"

using arma::mat;
using arma::uword;
using arma::vec;
using arma::uvec;
using arma::eig_sym;
using arma::find;

namespace pense {
namespace enpy_psc_internal {
//! Compute the Principal Sensitivity Components using identities for the Ridge penalty
//! if OpenMP is enabled and necessary.
//!
//! @param loss the LS regression loss object to compute the PSCs for.
//! @param penalties a lisf of Ridge penalties for which the PSCs should be computed at once.
//! @param optim Optimizer to use to compute estimate on full data.
//! @param num_threads number of OpenMP threads.
//! @return A list of PSC structures, one for each given penalty, in the same order as `penalties`.
alias::FwdList<pense::PscResult<DirectRidgeOptimizer>> ComputeRidgePscs(const nsoptim::LsRegressionLoss& loss,
    const alias::FwdList<nsoptim::RidgePenalty>& penalties, const DirectRidgeOptimizer& optim, int num_threads) {

  const nsoptim::PredictorResponseData& data = loss.data();
  // A list of PscResult objects, one per penalty.
  pense::utility::OrderedList<double, pense::PscResult<DirectRidgeOptimizer>, std::greater<double>> psc_results;
  arma::mat x_int;
  const arma::mat& x = loss.IncludeIntercept() ? x_int : data.cx();
  arma::mat ridge_gram;
  double gram_diag_int = 0;

  if (loss.IncludeIntercept()) {
    x_int = arma::join_rows(arma::ones(data.n_obs()), data.cx());
    ridge_gram = x_int.t() * x_int;
    gram_diag_int = ridge_gram.at(0, 0);
  } else {
    ridge_gram = data.cx().t() * data.cx();
  }

  // Computing PSCs can be done in parallel for each penalty. (default(none) does not work in gcc 9 and up)
  #pragma omp parallel num_threads(num_threads) \
    shared(psc_results, penalties, loss, data, x, optim) firstprivate(gram_diag_int, ridge_gram)
  {
    #pragma omp single nowait
    {
      for (auto pen_it = penalties.cbegin(), pen_end = penalties.cend(); pen_it != pen_end; ++pen_it) {
        // Compute optimum on full data.
        #pragma omp task firstprivate(pen_it, ridge_gram, gram_diag_int) \
          shared(psc_results, loss, data, x, optim)
        {
          auto optimizer = optim;
          optimizer.loss(loss);
          optimizer.penalty(*pen_it);
          auto psc_result_it = psc_results.begin();

          #pragma omp critical(emplace_psc_result)
          psc_result_it = psc_results.emplace(pen_it->lambda(), optimizer.Optimize());

          // Compute hat matrix for LOO residuals manually
          ridge_gram.diag() += (data.n_obs() - 1) * pen_it->lambda();
          if (loss.IncludeIntercept()) {
            ridge_gram.at(0, 0) = gram_diag_int;
          }

          arma::mat hat = x * arma::solve(ridge_gram, x.t());
          // Fitted y from all data:
          const arma::vec y_hat = data.cx() * psc_result_it->optimum.coefs.beta +
            psc_result_it->optimum.coefs.intercept;
          // Fitted y from LOO
          const arma::vec y_hat_loo = hat * data.cy();

          // Compute sensitivity matrix
          hat.each_row() %= arma::trans((data.cy() - y_hat_loo) / (1 - hat.diag()));
          hat.each_col() += y_hat - y_hat_loo;
          enpy_psc_internal::FinalizePSC(hat, &(*psc_result_it));
        }
      }
    }
  }

  return psc_results.items();
}

//! Compute the Principal Sensitivity Components using identities for the Ridge penalty
//! if OpenMP is disabled.
//!
//! @param loss the LS regression loss object to compute the PSCs for.
//! @param penalties a lisf of Ridge penalties for which the PSCs should be computed at once.
//! @param optimizer Optimizer to use to compute leave-one-out residuals.
//! @return A list of PSC structures, one for each given penalty, in the same order as `penalties`.
alias::FwdList<pense::PscResult<DirectRidgeOptimizer>> ComputeRidgePscs(const nsoptim::LsRegressionLoss& loss,
    const alias::FwdList<nsoptim::RidgePenalty>& penalties, DirectRidgeOptimizer optimizer) {
  const nsoptim::PredictorResponseData& data = loss.data();
  // A list of PscResult objects, one per penalty.
  alias::FwdList<pense::PscResult<DirectRidgeOptimizer>> psc_results;

  auto psc_result_it = psc_results.before_begin();
  arma::mat x_int;
  const arma::mat& x = loss.IncludeIntercept() ? x_int : data.cx();
  arma::mat ridge_gram;
  double gram_diag_int = 0;

  if (loss.IncludeIntercept()) {
    x_int = arma::join_rows(arma::ones(data.n_obs()), data.cx());
    ridge_gram = x_int.t() * x_int;
    gram_diag_int = ridge_gram.at(0, 0);
  } else {
    ridge_gram = data.cx().t() * data.cx();
  }
  double diagonal_offset = 0;

  // First optimize with respect to the full data set and compute the predictions.
  optimizer.loss(loss);
  for (auto&& penalty : penalties) {
    // Compute optimum on full data.
    optimizer.penalty(penalty);
    psc_result_it = psc_results.emplace_after(psc_result_it, optimizer.Optimize());

    // Compute hat matrix for LOO residuals manually
    const double prev_diagonal_offset = diagonal_offset;
    diagonal_offset = (data.n_obs() - 1) * penalty.lambda();
    ridge_gram.diag() += (diagonal_offset - prev_diagonal_offset);
    if (loss.IncludeIntercept()) {
      ridge_gram.at(0, 0) = gram_diag_int;
    }

    arma::mat hat = x * arma::solve(ridge_gram, x.t());
    // Fitted y from all data:
    const arma::vec y_hat = data.cx() * psc_result_it->optimum.coefs.beta + psc_result_it->optimum.coefs.intercept;
    // Fitted y from LOO
    const arma::vec y_hat_loo = hat * data.cy();

    // Compute sensitivity matrix
    hat.each_row() %= arma::trans((data.cy() - y_hat_loo) / (1 - hat.diag()));
    hat.each_col() += y_hat - y_hat_loo;
    enpy_psc_internal::FinalizePSC(hat, &(*psc_result_it));
  }

  return psc_results;
}

void FinalizePSC(const mat& sensitivity_matrix, PscResult* psc_result) {
  if (psc_result->warnings > 0) {
    psc_result->status = PscStatusCode::kWarning;
    psc_result->message.append("Some LOO residuals are unreliable; ");
  }

  // Get the Eigen decomposition of the outer product of the sensitivity vectors, i.e., the eigenvectors
  // of matrix `P` (eq. 9 in the paper)
  vec eigenvalues;
  const bool success = eig_sym(eigenvalues, psc_result->pscs, sensitivity_matrix * sensitivity_matrix.t());
  if (!success) {
    psc_result->pscs.reset();
    psc_result->status = PscStatusCode::kError;
    psc_result->message.append("Eigendecomposition failed");
    return;
  }

  // Only use the Eigenvectors with "non-zero" Eigenvalue.
  // "Non-zero" means Eigenvalues larger than the numerical tolerance times the largest Eigenvalue.
  uword cutoff_index = eigenvalues.n_elem - 1;

  if (eigenvalues[cutoff_index] < kNumericZero) {
    psc_result->pscs.reset();
    psc_result->status = PscStatusCode::kError;
    psc_result->message.append("All Eigenvalues are zero");
    return;
  }

  const double cutoff_threshold = kNumericZero * eigenvalues[cutoff_index];
  // Determine the index at which the Eigenvalues are too small (the largest one is apparently fine, c.f. line 42).
  while (cutoff_index > 0 && eigenvalues[--cutoff_index] > cutoff_threshold) {}

  // Hide the Eigenvectors for zero Eigenvalues
  if (cutoff_index > 0) {
    // Eigenvalues are ordered ascending, so we can simply strip all the eigenvectors for the first `k` eigenvalues.
    psc_result->pscs.shed_cols(0, cutoff_index);
  }
}

//! Concatenate two lists containing LooStatus objects, one for each penalty.
//! The metrics and status from *single* are concatenated to the metrics and the status from *combined*.
//!
//! @param single List to concatenate to the other list.
//! @param combined Output of the concatenated lists.
void ConcatenateLooStatus(alias::FwdList<LooStatus>* single, alias::FwdList<LooStatus>* combined) noexcept {
  auto in_it = single->begin();
  auto out_it = combined->before_begin();
  const auto in_end = single->end();
  const auto out_end = combined->end();

  while (in_it != in_end) {
    auto out_it_insert = out_it++;
    if (out_it == out_end) {
      out_it = combined->emplace_after(out_it_insert, std::move(*in_it));
      out_it->warnings += in_it->warnings;
    } else {
      out_it->warnings += in_it->warnings;
      out_it->metrics.splice_after(out_it->metrics.before_begin(), in_it->metrics);
      if (in_it->status == PscStatusCode::kError) {
        out_it->status = PscStatusCode::kError;
      } else if (in_it->status == PscStatusCode::kWarning && out_it->status != PscStatusCode::kError) {
        out_it->status = PscStatusCode::kWarning;
      }
    }
    ++in_it;
  }
}

}  // namespace enpy_psc_internal
}  // namespace pense
