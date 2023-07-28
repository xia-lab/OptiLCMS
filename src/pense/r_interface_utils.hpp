//
//  r_interface_utils.hpp
//  pense
//
//  Created by David Kepplinger on 2019-04-03.
//  Copyright © 2019 David Kepplinger. All rights reserved.
//

#ifndef R_INTERFACE_UTILS_HPP_
#define R_INTERFACE_UTILS_HPP_

#include "rcpp_integration.hpp"
#include "constants.hpp"
#include "alias.hpp"
#include "cd_pense.hpp"

namespace pense {
namespace r_interface {
//! Map an R-list into a forward list, extracting only the elements at the given indices.
//!
//! @param r_list an R-list with elements of type `T`.
//! @param indices vector of 1-based indices.
//! @return a `forward_list` with elements of type `T`, in the same order the given indices.
template <typename T>
alias::FwdList<T> ExtractListSubset(SEXP r_list, SEXP r_indices) {
  const Rcpp::List list(r_list);
  alias::FwdList<T> ret_list;
  auto ret_list_it = ret_list.before_begin();

  for (auto&& index : Rcpp::IntegerVector(r_indices)) {
    ret_list_it = ret_list.insert_after(ret_list_it, Rcpp::as<T>(list[index - 1]));
  }

  return ret_list;
}

//! Create a PredictorResponseData object from the R predictor matrix and the R response vector.
//!
//! This creates a read-only view to the memory managed by R.
//!
//! @param x numeric predictor matrix with `n` observation (rows) and `p` predictors (columns)
//! @param y numeric response vector with `n` elements
//! @return a pointer to the read-only predictor-response data. If the data is invalid, an exception is thrown.
std::unique_ptr<const nsoptim::PredictorResponseData> MakePredictorResponseData(SEXP x, SEXP y);

//! Create an adaptive EN penalty.
//! Note: we are not using the Rcpp::as path here as the loadings are shared among all penalties.
//!
//! @param penalty an R-list object containing information about the hyper-parameters.
//! @param loadings vector with penalty loadings.
//! @return the adaptive EN penalty object.
nsoptim::AdaptiveEnPenalty MakeAdaptiveEnPenalty(SEXP penalty, std::shared_ptr<const arma::vec> loadings);

//! Create a list of adaptive EN penalties.
//! Note: we are not using the Rcpp::as path here as the loadings are shared among all penalties.
//!
//! @param penalty an R-list object containing information about the hyper-parameters.
//! @param loadings vector with penalty loadings.
//! @return a list of adaptive EN penalty objects.
std::forward_list<nsoptim::AdaptiveEnPenalty> MakeAdaptiveEnPenaltyList(SEXP penalties, SEXP loadings);

//! Create a list of adaptive EN penalties, extracting only the penalties at the given indices.
//! Note: we are not using the as path here as the loadings are shared among all penalties.
//!
//! @param penalty an R-list object containing information about the hyper-parameters.
//! @param indices vector of 1-based indices to extract.
//! @param loadings vector with penalty loadings.
//! @return a list of adaptive EN penalty objects.
std::forward_list<nsoptim::AdaptiveEnPenalty> MakeAdaptiveEnPenaltyList(SEXP penalties, SEXP indices, SEXP loadings);

//! Create an adaptive LASSO penalty.
//! Note: we are not using the Rcpp::as path here as the loadings are shared among all penalties.
//!
//! @param penalty an R-list object containing information about the hyper-parameter.
//! @param loadings vector with penalty loadings.
//! @return the adaptive LASSO penalty object.
nsoptim::AdaptiveLassoPenalty MakeAdaptiveLassoPenalty(SEXP penalty, std::shared_ptr<const arma::vec> loadings);

//! Create a list of adaptive LASSO penalties.
//! Note: we are not using the as path here as the loadings are shared among all penalties.
//!
//! @param penalty an R-list object containing information about the hyper-parameters.
//! @param loadings vector with penalty loadings.
//! @return a list of adaptive LASSO penalty objects.
std::forward_list<nsoptim::AdaptiveLassoPenalty> MakeAdaptiveLassoPenaltyList(SEXP penalties, SEXP loadings);

//! Create a list of adaptive LASSO penalties, extracting only the penalties at the given indices.
//! Note: we are not using the as path here as the loadings are shared among all penalties.
//!
//! @param penalty an R-list object containing information about the hyper-parameters.
//! @param loadings vector with penalty loadings.
//! @param indices vector of 1-based indices to extract.
//! @return a list of adaptive LASSO penalty objects.
std::forward_list<nsoptim::AdaptiveLassoPenalty> MakeAdaptiveLassoPenaltyList(SEXP penalties, SEXP loadings,
                                                                              SEXP indices);

//! Get an unsafe view to the given R vector without copying any data.
//!
//! @param numeric_vector a numeric R vector
std::unique_ptr<const arma::vec> MakeVectorView(SEXP numeric_vector) noexcept;

namespace utils_internal {
//! Create a list of penalties with the penalty loadings taken from the list of optional arguments.
//!
//! @param r_penalties a list of adaptive EN-penalties.
//! @param optional_args a list with element `pen_loadings`, a numeric vector with *p* positive penalty loadings.
//! @return a list of penalty functions.
template<typename T>
alias::FwdList<T> MakePenalties(SEXP r_penalties, const Rcpp::List& optional_args, std::true_type) {
  return MakeAdaptiveEnPenaltyList(r_penalties, optional_args["pen_loadings"]);
}

//! Create a list of penalties without penalty loadings.
//!
//! @param r_penalties a list of adaptive EN-penalties.
//! @return a list of penalty functions.
template<typename T>
alias::FwdList<T> MakePenalties(SEXP r_penalties, const Rcpp::List&, std::false_type) {
  return Rcpp::as<alias::FwdList<T>>(r_penalties);
}

//! Create a list of penalties with the penalty loadings taken from the list of optional arguments.
//! Only penalties with indices in `r_indices` are extracted.
//!
//! @param r_penalties a list of adaptive EN-penalties.
//! @param indices vector of 1-based indices.
//! @param optional_args a list with element `pen_loadings`, a numeric vector with *p* positive penalty loadings.
//! @return a list of penalty functions.
template<typename T>
alias::FwdList<T> MakePenalties(SEXP r_penalties, SEXP r_indices, const Rcpp::List& optional_args, std::true_type) {
  return MakeAdaptiveEnPenaltyList(r_penalties, r_indices, optional_args["pen_loadings"]);
}

//! Create a list of penalties without penalty loadings.
//! Only penalties with indices in `r_indices` are extracted.
//!
//! @param r_penalties a list of adaptive EN-penalties.
//! @param indices vector of 1-based indices.
//! @return a list of penalty functions.
template<typename T>
alias::FwdList<T> MakePenalties(SEXP r_penalties, SEXP r_indices, const Rcpp::List&, std::false_type) {
  return ExtractListSubset<T>(r_penalties, r_indices);
}

//! Conditional alias for the `CDPense` class.
//! It is only defined if `T` actually is a specialization of the `CDPense` templated class.
template<typename T>
using CDPenseOptimizer = std::enable_if<
  std::is_same<T, pense::CDPense<typename T::PenaltyFunction, typename T::Coefficients>>::value, T>;

//! Conditional alias for the `nsoptim::CoordinateDescentOptimizer`.
//! It is only defined if `T` actually is a specialization of the `nsoptim::CoordinateDescentOptimizer` templated class.
template<typename T>
using CoordinateDescentOptimizer = std::enable_if<
  std::is_same<T, nsoptim::CoordinateDescentOptimizer<typename T::LossFunction, typename T::PenaltyFunction,
                                                      typename T::Coefficients>>::value, T>;

//! Conditional alias for the `nsoptim::AugmentedLarsOptimizer`.
//! It is only defined if `T` actually is a specialization of the `nsoptim::AugmentedLarsOptimizer` templated class.
template<typename T>
using AugmentedLarsOptimizer = std::enable_if<
  std::is_same<T, nsoptim::AugmentedLarsOptimizer<typename T::LossFunction, typename T::PenaltyFunction,
                                                  typename T::Coefficients>>::value, T>;

//! Conditional alias for the `nsoptim::DalEnOptimizer`.
//! It is only defined if `T` actually is a specialization of the `nsoptim::DalEnOptimizer` templated class.
template<typename T>
using DalOptimizer = std::enable_if<
  std::is_same<T, nsoptim::DalEnOptimizer<typename T::LossFunction, typename T::PenaltyFunction>>::value, T>;

//! Conditional alias for the `nsoptim::MMOptimizer`.
//! It is only defined if `T` actually is a specialization of the `nsoptim::MMOptimizer` templated class.
template<typename T>
using MMOptimizer = std::enable_if<
  std::is_same<T, nsoptim::MMOptimizer<typename T::LossFunction, typename T::PenaltyFunction,
                                       typename T::InnerOptimizer, typename T::Coefficients>>::value, T>;

//! Conditional alias for the `nsoptim::LinearizedAdmmOptimizer`.
//! It is only defined if `T` actually is a specialization of the `nsoptim::LinearizedAdmmOptimizer` templated class.
template<typename T>
using LinearizedAdmmOptimizer = std::enable_if<
  std::is_same<T, nsoptim::LinearizedAdmmOptimizer<typename T::ProximalOperator, typename T::PenaltyFunction,
                                                   typename T::Coefficients>>::value, T>;

//! Templated wrapper around the configuration object for the proximal operator class `T`.
template<typename T>
struct ProximalOperatorConfig {
  using type = Rcpp::List;
};

//! Templated wrapper around the configuration object for the proximal operator, specialized for `LsProximalOperator`.
template<>
struct ProximalOperatorConfig<nsoptim::LsProximalOperator> {
  using type = double;
};

//! Templated wrapper around the configuration object for the proximal operator,
//! specialized for `WeightedLsProximalOperator`.
template<>
struct ProximalOperatorConfig<nsoptim::WeightedLsProximalOperator> {
  using type = double;
};

//! Create an object of the given arbitrary Optimizer class.
//!
//! @param options... Ignored for optimizers without specialization.
template<typename Optimizer, typename... Ts>
Optimizer MakeOptimizer(double, Ts&&... /* options */) {
  return Optimizer();
}

//! Create an object of the CD optimizer for the S-loss class.
//!
//! @param options a list of options for the CDPense optimizer.
template<typename Optimizer>
typename CDPenseOptimizer<Optimizer>::type MakeOptimizer(int, const Rcpp::List& options);

//! Create an object of the CD-LS optimizer class.
//!
//! @param options a list of options for the CD-LS optimizer.
template<typename Optimizer>
typename CoordinateDescentOptimizer<Optimizer>::type MakeOptimizer(int, const Rcpp::List& options);

//! Create an object of the augmented LARS optimizer class.
//!
//! @param options a list of options for the variable augmented LARS optimizer.
template<typename Optimizer>
typename AugmentedLarsOptimizer<Optimizer>::type MakeOptimizer(int, const Rcpp::List& options);

//! Create an object of the linearized ADMM optimizer class.
//!
//! @param options a list of options for the linearized ADMM optimizer.
//! @param prox_opts options for the proximal operator.
template<typename Optimizer>
typename LinearizedAdmmOptimizer<Optimizer>::type MakeOptimizer(int, const Rcpp::List& options,
                                                                const Rcpp::List& prox_opts);

//! Create an object of the linearized ADMM optimizer class.
//!
//! @param options a list of options for the linearized ADMM optimizer.
template<typename Optimizer>
typename LinearizedAdmmOptimizer<Optimizer>::type MakeOptimizer(int, const Rcpp::List& options);

//! Create an object of the given DAL optimizer class.
//!
//! @param options a list of options for the DAL optimizer.
template<typename Optimizer>
typename DalOptimizer<Optimizer>::type MakeOptimizer(int, const Rcpp::List& options);

//! Create an object of the given MM optimizer class.
//!
//! @param options a list of options for the MM optimizer and it's inner optimizer type in "inner_opts".
template<typename Optimizer>
typename MMOptimizer<Optimizer>::type MakeOptimizer(int, const Rcpp::List& options);

//! Create an object of the given MM optimizer class.
//!
//! @param options a list of options for the MM optimizer and it's inner optimizer type.
//! @param other_options... further options passed on to the inner optimizer.
template<typename Optimizer, typename... Ts>
typename MMOptimizer<Optimizer>::type MakeOptimizer(int, const Rcpp::List& options, Ts&&... other_options);

//! Create an object of the CD optimizer for the S-loss class.
//!
//! @param options a list of options for the CDPense optimizer.
template<typename Optimizer>
typename CDPenseOptimizer<Optimizer>::type MakeOptimizer(int, const Rcpp::List& options) {
  Optimizer optim(Rcpp::as<CDPenseConfiguration>(options));
  optim.convergence_tolerance(pense::GetFallback(options, "eps", pense::kDefaultConvergenceTolerance));
  return optim;
}

//! Create an object of the CD-LS optimizer class.
//!
//! @param options a list of options for the CD-LS optimizer.
template<typename Optimizer>
typename CoordinateDescentOptimizer<Optimizer>::type MakeOptimizer(int, const Rcpp::List& options) {
  Optimizer optim(Rcpp::as<nsoptim::CDConfiguration>(options));
  optim.convergence_tolerance(pense::GetFallback(options, "eps", pense::kDefaultConvergenceTolerance));
  return optim;
}

//! Create an object of the augmented LARS optimizer class.
//!
//! @param options a list of options for the variable augmented LARS optimizer.
template<typename Optimizer>
typename AugmentedLarsOptimizer<Optimizer>::type MakeOptimizer(int, const Rcpp::List& options) {
  return Optimizer();
}

//! Make the configuration object for the given proximal operator class `T`.
//! The default implementation doesn't do anything but return its argument.
//!
//! @param prox_opts user-supplied options for the proximal operator.
//! @return options for the proximal operator as the correct type.
template<typename T>
inline typename ProximalOperatorConfig<T>::type MakeProximalOperatorConfig(const Rcpp::List& prox_opts) noexcept {
  return prox_opts;
}

//! Make the configuration object for the weighted LS proximal operator.
//! @param prox_opts options for the weighted LS proximal operator.
template<>
inline ProximalOperatorConfig<nsoptim::LsProximalOperator>::type
MakeProximalOperatorConfig<nsoptim::LsProximalOperator>(const Rcpp::List& prox_opts) noexcept {
  return pense::GetFallback(prox_opts, "tau", -1.);
}

//! Make the configuration object for the weighted LS proximal operator.
//! @param prox_opts options for the weighted LS proximal operator.
template<>
inline ProximalOperatorConfig<nsoptim::WeightedLsProximalOperator>::type
MakeProximalOperatorConfig<nsoptim::WeightedLsProximalOperator>(const Rcpp::List& prox_opts) noexcept {
  return pense::GetFallback(prox_opts, "tau", -1.);
}

//! Create an object of the linearized ADMM optimizer class.
//!
//! @param options a list of options for the linearized ADMM optimizer.
//! @param prox_opts options for the proximal operator.
template<typename Optimizer>
typename LinearizedAdmmOptimizer<Optimizer>::type MakeOptimizer(int, const Rcpp::List& options,
                                                                const Rcpp::List& prox_opts) {
  const auto prox_opts_obj = MakeProximalOperatorConfig<typename Optimizer::ProximalOperator>(prox_opts);
  Optimizer optim(Rcpp::as<nsoptim::AdmmLinearConfiguration>(options), prox_opts_obj);
  optim.convergence_tolerance(pense::GetFallback(options, "eps", pense::kDefaultConvergenceTolerance));
  return optim;
}

//! Create an object of the linearized ADMM optimizer class.
//!
//! @param options a list of options for the linearized ADMM optimizer.
template<typename Optimizer>
typename LinearizedAdmmOptimizer<Optimizer>::type MakeOptimizer(int, const Rcpp::List& options) {
  if (options.containsElementNamed("prox_opts")) {
    const auto prox_opts = Rcpp::as<Rcpp::List>(options["prox_opts"]);
    return MakeOptimizer<Optimizer>(1, options, prox_opts);
  } else {
    Optimizer optim(Rcpp::as<nsoptim::AdmmLinearConfiguration>(options));
    optim.convergence_tolerance(pense::GetFallback(options, "eps", pense::kDefaultConvergenceTolerance));
    return optim;
  }
}

//! Create an object of the given DAL optimizer class.
//!
//! @param options a list of options for the DAL optimizer.
template<typename Optimizer>
typename DalOptimizer<Optimizer>::type MakeOptimizer(int, const Rcpp::List& options) {
  Optimizer optim(Rcpp::as<nsoptim::DalEnConfiguration>(options));
  optim.convergence_tolerance(pense::GetFallback(options, "eps", pense::kDefaultConvergenceTolerance));
  return optim;
}

//! Create an object of the given MM optimizer class.
//!
//! @param options a list of options for the MM optimizer and it's inner optimizer type.
template<typename Optimizer>
typename MMOptimizer<Optimizer>::type MakeOptimizer(int, const Rcpp::List& options) {
  const Rcpp::List inner_optimizer_options = pense::GetFallback(options, "inner_opts", Rcpp::List());
  Optimizer optim(MakeOptimizer<typename Optimizer::InnerOptimizer>(inner_optimizer_options),
                  Rcpp::as<nsoptim::MMConfiguration>(options));
  optim.convergence_tolerance(pense::GetFallback(options, "eps", pense::kDefaultConvergenceTolerance));
  return optim;
}

//! Create an object of the given MM optimizer class.
//!
//! @param options a list of options for the MM optimizer and it's inner optimizer type.
template<typename Optimizer, typename... Ts>
typename MMOptimizer<Optimizer>::type MakeOptimizer(int, const Rcpp::List& options, Ts&&... other_options) {
  const auto mm_config = Rcpp::as<nsoptim::MMConfiguration>(options);
  Optimizer optim(MakeOptimizer<typename Optimizer::InnerOptimizer>(1, std::forward<Ts>(other_options)...),
                  mm_config);
  optim.convergence_tolerance(pense::GetFallback(options, "eps", pense::kDefaultConvergenceTolerance));
  return optim;
}
}  // namespace utils_internal

//! Create a list of penalties. If the penalties require penalty loadings, they are taken from the list of optional
//! arguments.
//! If the optimizer-type requires penalty loadings, they are taken from the list of optional arguments.
//!
//! @param penalties a list of (adaptive) EN-penalties.
//! @param optional_args a list with element `pen_loadings`, a numeric vector with *p* positive penalty loadings.
//! @return a list of penalty functions.
template<typename Optimizer>
alias::FwdList<typename Optimizer::PenaltyFunction> MakePenalties(SEXP penalties, const Rcpp::List& optional_args) {
  using is_adaptive_tag = typename nsoptim::traits::is_adaptive<typename Optimizer::PenaltyFunction>::type;
  using PenaltyFunction = typename Optimizer::PenaltyFunction;
  return utils_internal::MakePenalties<PenaltyFunction>(penalties, optional_args, is_adaptive_tag{});
}

//! Create a list of penalties with the penalty loadings taken from the list of optional arguments.
//! Only penalties with indices in `r_indices` are extracted.
//! If the optimizer-type requires penalty loadings, they are taken from the list of optional arguments.
//!
//! @param penalties a list of adaptive EN-penalties.
//! @param indices vector of 1-based indices.
//! @param optional_args a list with element `pen_loadings`, a numeric vector with *p* positive penalty loadings.
//! @return a list of penalty functions.
template<typename Optimizer>
alias::FwdList<typename Optimizer::PenaltyFunction> MakePenalties(SEXP penalties, SEXP indices,
                                                                  const Rcpp::List& optional_args) {
  using is_adaptive_tag = typename nsoptim::traits::is_adaptive<typename Optimizer::PenaltyFunction>::type;
  using PenaltyFunction = typename Optimizer::PenaltyFunction;
  return utils_internal::MakePenalties<PenaltyFunction>(penalties, indices, optional_args, is_adaptive_tag{});
}

//! Create an object of the given Optimizer.
//!
//! @param options 0 or more lists of options for the optimizer.
template<typename Optimizer, typename... Ts>
Optimizer MakeOptimizer(Ts&&... options) {
  return utils_internal::MakeOptimizer<Optimizer>(1, std::forward<Ts>(options)...);
}

}  // namespace r_interface
}  // namespace pense

#endif  // R_INTERFACE_UTILS_HPP_
