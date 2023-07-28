//
//  r_interface_utils.cc
//  pense
//
//  Created by David Kepplinger on 2019-01-30.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#include "r_interface_utils.hpp"

#include <memory>
#include <exception>

#include "rcpp_integration.hpp"

using arma::sp_vec;
using arma::vec;

using Rcpp::as;

using nsoptim::AdaptiveLassoPenalty;
using nsoptim::AdaptiveEnPenalty;
using nsoptim::PredictorResponseData;
using SpRegressionCoefficients = nsoptim::RegressionCoefficients<sp_vec>;

namespace {
struct MatrixDimensions {
  int n_rows;
  int n_cols;
};

//! Get dimensions of an R matrix
//! @param matrix the R matrix object
//! @param nrows the variable to store the number of rows in
//! @param ncols the variable to store the number of columns in
MatrixDimensions GetMatrixDimensions(SEXP matrix) noexcept;
}  // namespace

namespace pense {
namespace r_interface {
std::unique_ptr<const PredictorResponseData> MakePredictorResponseData(SEXP x, SEXP y) {
  const int nobs_y = Rf_length(y);
  const MatrixDimensions x_dims = GetMatrixDimensions(x);

  // Check that x and y have the same nr. of observations
  if (nobs_y != x_dims.n_rows) {
    throw std::invalid_argument("y and x have differing number of observations");
  }

  // Check that x and y have the required data type
  if (TYPEOF(x) != REALSXP || TYPEOF(y) != REALSXP) {
    throw std::invalid_argument("y and x must be numeric");
  }

  return std::make_unique<const PredictorResponseData>(arma::mat(REAL(x), x_dims.n_rows, x_dims.n_cols, false, true),
                                                       vec(REAL(y), nobs_y, false, true));
}

AdaptiveEnPenalty MakeAdaptiveEnPenalty(SEXP r_penalty, std::shared_ptr<const vec> loadings) {
  const Rcpp::List penalty(r_penalty);
  return AdaptiveEnPenalty(loadings, as<double>(penalty["alpha"]), as<double>(penalty["lambda"]));
}

//! Create a list of adaptive EN penalties.
//! Note: we are not using the as path here as the loadings are shared among all penalties.
//!
//! @param penalty an R-list object containing information about the hyper-parameters.
//! @param loadings vector with penalty loadings.
//! @return a list of adaptive EN penalty objects, in the same order as `indices`.
std::forward_list<AdaptiveEnPenalty> MakeAdaptiveEnPenaltyList(SEXP r_penalties, SEXP r_loadings) {
  std::shared_ptr<const vec> loadings = MakeVectorView(r_loadings);
  std::forward_list<AdaptiveEnPenalty> ret_list;
  auto ret_list_it = ret_list.before_begin();

  for (auto&& penalty_sexp : Rcpp::List(r_penalties)) {
    const Rcpp::List penalty(penalty_sexp);
    ret_list_it = ret_list.emplace_after(ret_list_it, loadings, as<double>(penalty["alpha"]),
                                         as<double>(penalty["lambda"]));
  }
  return ret_list;
}

//! Create a list of adaptive EN penalties, extracting only the penalties at the given indices.
//! Note: we are not using the as path here as the loadings are shared among all penalties.
//!
//! @param penalty an R-list object containing information about the hyper-parameters.
//! @param loadings vector with penalty loadings.
//! @param indices vector of 1-based indices to extract.
//! @return a list of adaptive EN penalty objects, in the same order as `indices`.
std::forward_list<AdaptiveEnPenalty> MakeAdaptiveEnPenaltyList(SEXP r_penalties, SEXP r_indices, SEXP r_loadings) {
  const Rcpp::List penalties(r_penalties);
  std::shared_ptr<const vec> loadings = MakeVectorView(r_loadings);
  std::forward_list<AdaptiveEnPenalty> ret_list;
  auto ret_list_it = ret_list.before_begin();

  for (auto&& index : Rcpp::IntegerVector(r_indices)) {
    const Rcpp::List penalty(penalties[index - 1]);
    ret_list_it = ret_list.emplace_after(ret_list_it, loadings, as<double>(penalty["alpha"]),
                                         as<double>(penalty["lambda"]));
  }

  return ret_list;
}

AdaptiveLassoPenalty MakeAdaptiveLassoPenalty(SEXP r_penalty, std::shared_ptr<const vec> loadings) {
  const Rcpp::List penalty(r_penalty);
  return AdaptiveLassoPenalty(loadings, as<double>(penalty["lambda"]));
}

//! Create a list of adaptive LASSO penalties.
//! Note: we are not using the as path here as the loadings are shared among all penalties.
//!
//! @param penalty an R-list object containing information about the hyper-parameters.
//! @param loadings vector with penalty loadings.
//! @return a list of adaptive EN penalty objects, in the same order as `indices`.
std::forward_list<AdaptiveLassoPenalty> MakeAdaptiveLassoPenaltyList(SEXP r_penalties, SEXP r_loadings) {
  std::shared_ptr<const vec> loadings = MakeVectorView(r_loadings);
  std::forward_list<AdaptiveLassoPenalty> ret_list;
  auto ret_list_it = ret_list.before_begin();

  for (auto&& penalty_sexp : Rcpp::List(r_penalties)) {
    const Rcpp::List penalty(penalty_sexp);
    ret_list_it = ret_list.emplace_after(ret_list_it, loadings, as<double>(penalty["lambda"]));
  }
  return ret_list;
}

//! Create a list of adaptive LASSO penalties, extracting only the penalties at the given indices.
//! Note: we are not using the as path here as the loadings are shared among all penalties.
//!
//! @param penalty an R-list object containing information about the hyper-parameters.
//! @param loadings vector with penalty loadings.
//! @param indices vector of 1-based indices to extract.
//! @return a list of adaptive EN penalty objects, in the same order as `indices`.
std::forward_list<AdaptiveLassoPenalty> MakeAdaptiveLassoPenaltyList(SEXP r_penalties, SEXP r_loadings,
                                                                     SEXP r_indices) {
  const Rcpp::List penalties(r_penalties);
  std::shared_ptr<const vec> loadings = MakeVectorView(r_loadings);
  std::forward_list<AdaptiveLassoPenalty> ret_list;
  auto ret_list_it = ret_list.before_begin();

  for (auto&& index : Rcpp::IntegerVector(r_indices)) {
    const Rcpp::List penalty(penalties[index - 1]);
    ret_list_it = ret_list.emplace_after(ret_list_it, loadings, as<double>(penalty["lambda"]));
  }

  return ret_list;
}

// pense::EstimateWrapper<sp_vec, AdaptiveLassoPenalty> MakeAdaptiveLassoEstimate(SEXP r_obj,
//                                                                                std::shared_ptr<const vec> loadings) {
//   return pense::EstimateWrapper<sp_vec, AdaptiveLassoPenalty>(
//     as<SpRegressionCoefficients>(r_obj), MakeAdaptiveLassoPenalty(r_obj, loadings));
// }

// alias::FwdList<pense::EstimateWrapper<sp_vec, AdaptiveLassoPenalty>> MakeAdaptiveLassoEstimateList(
//     SEXP r_list, std::shared_ptr<const vec> loadings) {
//   alias::FwdList<pense::EstimateWrapper<sp_vec, AdaptiveLassoPenalty>> ret_list;
//   auto ret_list_it = ret_list.before_begin();

//   for (auto&& list_item : Rcpp::List(r_list)) {
//     ret_list_it = ret_list.emplace_after(ret_list_it, as<SpRegressionCoefficients>(list_item),
//                                          MakeAdaptiveLassoPenalty(list_item, loadings));
//   }
//   return ret_list;
// }

//! Get an unsafe view to the given R vector without copying any data.
//!
//! @param numeric_vector a numeric R vector
std::unique_ptr<const vec> MakeVectorView(SEXP numeric_vector) noexcept {
  if (TYPEOF(numeric_vector) != REALSXP) {
    return std::make_unique<const vec>();
  }
  return std::make_unique<const vec>(REAL(numeric_vector), Rf_length(numeric_vector), false, true);
}
}  // namespace r_interface
}  // namespace pense

namespace {
//! Get dimensions of an R matrix
MatrixDimensions GetMatrixDimensions(SEXP matrix) noexcept {
  SEXP r_dims;
  PROTECT(r_dims = Rf_getAttrib(matrix, R_DimSymbol));
  int const * const dims = INTEGER(r_dims);
  MatrixDimensions mat_dims {dims[0], dims[1]};
  UNPROTECT(1);
  return mat_dims;
}
}  // namespace
